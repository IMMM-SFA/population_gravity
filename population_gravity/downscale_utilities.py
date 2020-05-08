"""
Helper functions for population downscaling.

@author: Hamidreza Zoraghein
@date: Created on Mon Dec 11 11:24:13 2017

@history:
SVN: $Id: pop_downscaling_module.py 224 2018-04-06 22:57:02Z kauff $
SVN: $URL: https://svn-iam-thesis.cgd.ucar.edu/population_spatial/trunk/src/pop_downscaling_module.py $

"""

import multiprocessing
import rasterio

import numpy as np
import pandas as pd

from rasterio.merge import merge
from collections import deque
from pathos.multiprocessing import ProcessingPool as Pool
from scipy.spatial import cKDTree


def mosaic(raster_list, out_raster, source_metadata):
    """Create a raster mosiac from multiple rasters and save to file.

    :param raster_list:             List of full path to rasters with file name and extensions
    :param out_raster:              Full path with file name and extension to write the raster to
    :param source_metadata:         Metadata rasterio object from the target states init raster

    :return:                        Mosaicked raster file name

    """
    # build list of raster objects
    raster_objects = [rasterio.open(i) for i in raster_list]

    # create mosaic
    mosaic, out_transform = merge(raster_objects)

    # update metadata with mosiac values
    source_metadata.update({"height": mosaic.shape[1], "width": mosaic.shape[2], "transform": out_transform})

    # write output
    with rasterio.open(out_raster, 'w', **source_metadata) as dest:
        dest.write(out_raster)

    return out_raster


def raster_to_array(raster):
    """Convert a raster to a NumPy array.

    :param raster:                  Full path to file with name and extension to the input raster

    :return:                        NumPy array

    """
    # read the input raster before converting it to an array
    with rasterio.open(raster) as src_raster:
        return src_raster.read(1)
    

def array_to_raster(input_raster, input_array, within_indices, output_raster):
    """Save NumPy array to a raster

    :param input_raster:            ?
    :param input_array:             ?
    :param within_indices:          ?
    :param output_raster:           ?

    """
    # read the template raster to be filled by output array values later
    with rasterio.open(input_raster) as src_raster:
        band = src_raster.read(1)
        src_profile = src_raster.profile
        row_count = band.shape[0]
        col_count = band.shape[1]
        flat_array = band.flatten()
    
    # replace initial array values with those from the input array
    flat_array[within_indices] = input_array
    
    array = flat_array.reshape(row_count, col_count)
    
    with rasterio.open(output_raster, "w", **src_profile) as dst:
        dst.write_band(1, array)
    

def all_index_retriever(array, columns, row_col='row', column_col='column', all_index_col='all_index'):
    """Build data frame in the shape of the input array.

    :param array:                   Input 2D array from input raster
    :param columns:                 list; Target columns
    :param row_col:                 Name of row column
    :param column_col:              Name of column column
    :param all_index_col:           Name of `all_index` column

    :return:                        Typed data frame of indicies

    """
    
    # put the row, column and linear indices of all elements in a dataframe
    shape = array.shape
    index = pd.MultiIndex.from_product([range(s)for s in shape], names=columns)
    df = pd.DataFrame({all_index_col: array.flatten()}, index=index).reset_index()
    df[all_index_col] = df.index
    
    return df.astype({row_col: np.int32, column_col: np.int32, all_index_col: np.int32})


def suitability_estimator(pop_dist_params):
    """Estimate suitability based on neighbor characteristics

    :param pop_dist_params:                 [0] target grid cell
                                            [1] ind_diffs
                                            [2] total_population_1st
                                            [3] alpha_parameter
                                            [4] exp_xx_inv_beta_dist

    """
    # the id of the current focal point
    element = pop_dist_params[0]
    
    # construct nearby indices
    neigh_indices = element + pop_dist_params[1]

    # population of close points
    pop = pop_dist_params[2][neigh_indices]

    # if 0.0 then inf will be generated which will raise a divide by zero runtime warning for np.power()
    pop = np.where(pop == 0.0, np.nan, pop)

    # calculation of the other elements of the suitability equation
    pop_xx_alpha = np.power(pop, pop_dist_params[3])

    # convert back to 0.0
    pop_xx_alpha = np.where(np.isnan(pop_xx_alpha), 0.0, pop_xx_alpha)
    pop_xx_alpha = np.where(np.isinf(pop_xx_alpha), 0.0, pop_xx_alpha)

    pop_dist_factor = pop_xx_alpha * pop_dist_params[4]
    
    return np.sum(pop_dist_factor)


def dist_matrix_calculator(first_index, cut_off_meters, all_indices, coordinate_array, row_col='row',
                           column_col='column', all_index_col='all_index', dis_col='dis',
                           ind_diff_col='ind_diff', row_diff_col='row_diff', col_diff_col='col_diff'):
    """TODO: Fill in description


    :param first_index:             ?
    :param cut_off_meters:          ?
    :param all_indices:             ?
    :param coordinate_array:        ?
    :param row_col:                 ?
    :param column_col:              ?
    :param all_index_col:           ?
    :param dis_col:                 ?
    :param ind_diff_col:            ?
    :param row_diff_col:            ?
    :param col_diff_col:            ?

    :return:                        ?

    """
    # calculate distances between the first point and all other points within a 100km neighborhood
    cut_off_metres = cut_off_meters + 1
    tree_1 = cKDTree(coordinate_array[first_index:first_index + 1, [0, 1]])
    tree_2 = cKDTree(coordinate_array[:, [0, 1]])
    tree_dist = cKDTree.sparse_distance_matrix(tree_1, tree_2, cut_off_metres, output_type='dict', p=2)
    
    # put distances and indices of neighboring in a data frame
    dist_df = pd.DataFrame(columns=["near_id", dis_col])

    # add list for zip to accommodate Python 3
    dist_df["near_id"] = coordinate_array[list(zip(*tree_dist))[1], 2].astype(np.int32)
    dist_df[dis_col] = tree_dist.values()
    dist_df = dist_df.loc[dist_df.loc[:, dis_col] != 0, :] # Remove the distance to itself
    
    # bring row and column indices of neighboring points by a join
    dist_df = dist_df.join(all_indices, on="near_id")
    
    # add to columns holding the relative difference in rows and colums beween focal point and its neighbors
    foc_indices = all_indices.loc[first_index, [row_col, column_col]].values
    dist_df[ind_diff_col] = dist_df["near_id"] - first_index
    dist_df[row_diff_col] = dist_df[row_col] - foc_indices[0]
    dist_df[col_diff_col] = dist_df[column_col] - foc_indices[1]
    
    # drop unwanted columns
    dist_df = dist_df.drop([row_col, column_col, "near_id", all_index_col], axis=1)
        
    return dist_df.astype({ind_diff_col: np.int32, row_diff_col: np.int32, col_diff_col: np.int32})


def pop_min_function(z, arr_pop_1st, arr_pop_2nd, arr_tot_pop_1st, points_mask,
                     dist_matrix, within_indices, ind_diff_col='ind_diff', dist_col='dis'):
    """TODO: fill in description

    :param z:                       ?
    :param arr_pop_1st:
    :param arr_pop_2nd:
    :param arr_tot_pop_1st:
    :param points_mask:             Mask values of points
    :param dist_matrix:             Template distance matrix
    :param within_indices:         Indices of points within the state boundary (subset of the above)
    :param ind_diff_col:            Column name for index difference between focal and nearby points
    :param dist_col:                Column containing distance

    :return:                        ?

    """
    # unpack z tuple
    a, b = z

    # calculate aggregate urban/rural population at times 1 and 2
    pop_t1 = np.sum(arr_pop_1st[within_indices])
    pop_t2 = np.sum(arr_pop_2nd[within_indices])
    
    # population change between the two reference years
    pop_change = pop_t2 - pop_t1    
    if pop_change < 0:
        negative_mod = 1
    else:
        negative_mod = 0
    
    # differences in index between focal and nearby points as a template
    ind_diffs = dist_matrix[ind_diff_col].values
    
    # distances between current point and its close points
    ini_dist = dist_matrix[dist_col].values/1000.0
    dist = -b * ini_dist

    exp_xx_inv_beta_dist = np.exp(dist)

    # initialize the parallelization
    pool = Pool(processes=multiprocessing.cpu_count())

    # provide the inputs for the parallelized function
    parallel_elements = deque([(i, ind_diffs, arr_tot_pop_1st, a, exp_xx_inv_beta_dist) for i in within_indices])

    # derive suitability estimates
    suitability_estimates = np.array(pool.map(suitability_estimator, parallel_elements))

    # extract only the necessary mask values that fall within the state boundary
    points_mask = points_mask[within_indices]
    
    # population in the first year
    pop_first_year = arr_pop_1st[within_indices]
    
    # in case of population decline, suitability estimates are reciprocated for non-zero values
    if negative_mod:
        
        # find those whose mask is 0 but have population, they should decline anyway 
        mask_zero = np.where(points_mask == 0)[0]
        pop_non_zero = np.where(pop_first_year != 0)[0]
        
        # those cells with mask value of zero and population are the intersection of the two above arrays
        pop_mask = np.intersect1d(mask_zero, pop_non_zero, assume_unique=True)
        
        # change the mask value of the above cells to the mean so that they also lose population
        points_mask[pop_mask] = points_mask.mean()
        
        # adjust suitability values by applying mask values
        suitability_estimates = points_mask * suitability_estimates
        
        # inverse current mask values for a better reflection of population decline
        suitability_estimates[suitability_estimates != 0] = 1.0/suitability_estimates[suitability_estimates != 0]
    
    else:
        # adjust suitability values by applying mask values
        suitability_estimates = points_mask * suitability_estimates
    
    # total suitability for the whole area, which is the summation of all individual suitability values
    tot_suitability = suitability_estimates.sum()
    
    # final population estimate for each point if negative mode is off
    pop_estimates = suitability_estimates/tot_suitability * pop_change + pop_first_year
    
    if negative_mod:
        while any(pop < 0 for pop in pop_estimates): # To ensure there is no negative population
            
            # treating negative population values
            extra_pop_mod = abs(pop_estimates[pop_estimates < 0].sum())
            pop_estimates[pop_estimates < 0] = 0

            # calculate the new total suitability value based on those points whose projected population is positive
            new_tot_suitability = suitability_estimates[pop_estimates > 0].sum()       
            
            # adjust non-negative population values to maintain the total aggregated population
            pop_estimates[pop_estimates > 0] = pop_estimates[pop_estimates > 0] - (suitability_estimates[pop_estimates > 0]/new_tot_suitability) * extra_pop_mod

    # produce the total error compared to observed values
    tot_error = abs(arr_pop_2nd[within_indices] - pop_estimates).sum()
    
    return tot_error
