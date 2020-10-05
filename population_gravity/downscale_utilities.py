"""
Helper functions for population downscaling.

@author: Hamidreza Zoraghein, Chris R. Vernon
@date: Created on Mon Dec 11 11:24:13 2017

"""

import gzip

import multiprocessing
import rasterio

import numpy as np
import pandas as pd

from rasterio.merge import merge
from rasterio.io import MemoryFile
from collections import deque
from pathos.multiprocessing import ProcessingPool as Pool
from scipy.spatial import cKDTree


def get_raster_with_metadata(source_raster):
    """Get a copy of the metadata from a source raster.

    :param source_raster:           Full path with file name and extension of an input raster

    :return:                        Rasterio object, Rasterio metadata object

    """
    raster_obj = rasterio.open(source_raster)

    return raster_obj, raster_obj.meta.copy()


def mosaic(raster_list, out_raster, source_metadata):
    """Create a raster mosiac from multiple rasters and save to file.

    :param raster_list:             List of full path to rasters with file name and extensions
    :param out_raster:              Full path with file name and extension to write the raster to
    :param source_metadata:         Metadata rasterio object from the target states init raster

    :return:                        Mosaicked rasterio object

    """

    # build list of raster objects
    raster_objects = [rasterio.open(i) for i in raster_list]

    # create mosaic
    mosaic, out_transform = merge(raster_objects)

    # update metadata with mosiac values
    source_metadata.update({"height": mosaic.shape[1], "width": mosaic.shape[2], "transform": out_transform})

    # write output
    with rasterio.open(out_raster, 'w', **source_metadata) as dest:
        dest.write(mosaic)

    # open output object with updated metadata
    return rasterio.open(out_raster)


def mosaic_memory(raster_objects, source_metadata):
    """Create a raster mosiac from multiple rasters and save to file.

    :param raster_list:             List of raster objects
    :param source_metadata:         Metadata rasterio object from the target states init raster

    :return:                        Mosaicked rasterio object

    """

    # create mosaic
    mosaic, out_transform = merge(raster_objects)

    # update metadata with mosiac values
    source_metadata.update({"height": mosaic.shape[1], "width": mosaic.shape[2], "transform": out_transform})

    # write output to memory
    with MemoryFile() as memfile:

        dataset = memfile.open(**source_metadata)
        dataset.write(mosaic)

        return dataset


def create_bbox(raster_object):
    """Create a bounding box array from a rasterio object.

    :param raster_object:           Rasterio raster object

    :return:                        bounding box

    """

    return np.array([raster_object.bounds.left,
                     raster_object.bounds.bottom,
                     raster_object.bounds.right,
                     raster_object.bounds.top])


def raster_to_array(raster):
    """Convert a raster to a NumPy array.

    :param raster:                  Full path to file with name and extension to the input raster

    :return:                        NumPy array

    """
    # read the input raster before converting it to an array
    with rasterio.open(raster) as src_raster:
        return src_raster.read(1)
    

def array_to_raster(template_raster_object, input_array, within_indices, output_raster):
    """Save NumPy array to a raster.

    :param template_raster_object:          Object containing [2D array, 1D array, row count, column count, profile]
    :param input_array:                     1D numpy array from run
    :param within_indices:                  Full path with file name and extension to the text file
                                            containing a file structured as a Python list (e.g. [0, 1]) that
                                            contains the index of each grid cell when flattened from a 2D array to
                                            a 1D array for the target state.
    :param output_raster:                   Full path with file name and extension to the output raster

    """

    flat_array = template_raster_object[1].copy()

    # replace initial array values with those from the input array
    flat_array[within_indices] = input_array

    # reshape array to raster size
    array = flat_array.reshape(template_raster_object[2], template_raster_object[3])
    
    with rasterio.open(output_raster, "w", **template_raster_object[4]) as dst:
        dst.write(array, 1)


def reshape_array_to_raster(template_raster_object, input_array, within_indices, output_array):
    """Reshape and save a flattened array to the shape of an input raster.

    :param template_raster_object:          Object containing [2D array, 1D array, row count, column count, profile]
    :param input_array:                     1D numpy array from run
    :param within_indices:                  Full path with file name and extension to the text file
                                            containing a file structured as a Python list (e.g. [0, 1]) that
                                            contains the index of each grid cell when flattened from a 2D array to
                                            a 1D array for the target state.

    :param output_array:                    Full path with file name and extension to the output array

    """

    flat_array = template_raster_object[1].copy()

    # replace initial array values with those from the input array
    flat_array[within_indices] = input_array

    # reshape array to raster size
    array = flat_array.reshape(template_raster_object[2], template_raster_object[3])

    # set in raster grid space
    with MemoryFile() as memfile:

        with memfile.open(**template_raster_object[4]) as dataset:
            dataset.write_band(1, array)

            np.save(output_array, dataset.read(1))


def array_to_raster_memory(template_raster_object, input_array, within_indices):
    """Save NumPy array to a raster in memory.

    :param template_raster_object:          Object containing [2D array, 1D array, row count, column count, profile]
    :param input_array:                     1D numpy array from run
    :param within_indices:                  Full path with file name and extension to the text file
                                            containing a file structured as a Python list (e.g. [0, 1]) that
                                            contains the index of each grid cell when flattened from a 2D array to
                                            a 1D array for the target state.

    """

    flat_array = template_raster_object[1].copy()

    # replace initial array values with those from the input array
    flat_array[within_indices] = input_array

    # reshape array to raster size
    array = flat_array.reshape(template_raster_object[2], template_raster_object[3])

    # write output to memory
    with MemoryFile() as memfile:

        dataset = memfile.open(**template_raster_object[4])
        dataset.write_band(1, array)

        return dataset


def array2d_to_raster_memory(array, raster_profile):
    """Write 2D array to raster in memory.

    :param array:                   2D numpy array
    :param raster_profile:          Profile from template raster.

    :return:                        rasterio dataset in memory

    """

    with MemoryFile() as memfile:

        dataset = memfile.open(**raster_profile)
        dataset.write_band(1, array)

        return dataset


def join_coords_to_value(vaild_coordinates_csv, valid_raster_values_csv, out_csv=None):
    """Join non-NODATA CSV raster value outputs to their corresponding X, Y coordinates.

    :param vaild_coordinates_csv:               Full path with file name and extension to the input CSV file containing
                                                non-NODATA grid cell coordinates where fields are [`x_coord`, `y_coord`]
                                                and represent the X, Y coordinates of the projected coordinate system.
                                                The key of this file for joining purposes is the order position (index).

    :param valid_raster_values_csv:             Full path with file name and extension to the input CSV file containing
                                                non-NODATA grid cell values from the output raster where fields are
                                                [`value`].  The key of this file for joining purposes is the order
                                                position (index).

    :param out_csv:                             Full path with file name and extension to write the output CSV file to.
                                                Output fields are [`x_coord`, `y_coord`, `value`].

    :return:                                    Joined Pandas DataFrame


    Example:

        >>> import population_gravity as pgr
        >>>
        >>> vaild_coordinates_csv = "<path to file>"
        >>> valid_raster_values_csv = "<path to file>"
        >>> out_csv = "<path to file>"
        >>>
        >>> df = pgr.join_coords_to_value(vaild_coordinates_csv, valid_raster_values_csv, out_csv)

    """

    # read in data
    df_coords = pd.read_csv(vaild_coordinates_csv)
    df_values = pd.read_csv(valid_raster_values_csv)

    # create key from index
    df_coords['key'] = df_coords.index
    df_values['key'] = df_values.index

    # conduct left join
    df_join = pd.merge(left=df_coords, right=df_values, how='left', on='key')

    # drop key
    df_join.drop(columns=['key'], inplace=True)

    # save to CSV if desired
    if out_csv is not None:
        df_join.to_csv(out_csv, index=False)

    return df_join


def raster_to_csv(input_array, grid_coordinates_array, out_csv, compress=True, export_value_only=True):
    """Create a CSV file with ['x_coord', 'y_coord', 'value'] from an input raster omitting cells with NODATA.

    :param input_array:                         Input 1D array of valid grid cells

    :param grid_coordinates_array:              Numpy Array containing the coordinates for each 1 km grid cell
                                                within the target state. File includes a header with the fields
                                                XCoord, YCoord, FID.
                                                Where data types and field descriptions are as follows:
                                                (XCoord, float, X coordinate in meters),
                                                (YCoord, float, Y coordinate in meters),
                                                (FID, int, Unique feature id)

    :param out_csv:                             str.  Full path with file name and extension to the output CSV file.

    :param compress:                            Boolean.  Compress CSV using GNU zip (gzip) compression; Default True

    :param export_value_only:                   Boolean.  Export only the value column; Default True

    """

    # just get x, y from array; match dtype of input raster
    coordinate_array = grid_coordinates_array[:, :2]

    # build data frame
    if export_value_only:
        df = pd.DataFrame({'value': input_array})
    else:
        df = pd.DataFrame({'x_coord': coordinate_array[:, 0], 'y_coord': coordinate_array[:, 1], 'value': input_array})

    # generate output file name
    if compress:
        out_csv = f'{out_csv}.gz'

    # write output
    df.to_csv(out_csv, index=False, compression='infer')


def gzip_csv_to_array(gzip_csv_file, header_rows_to_skip=1):
    """Read GZIP CSV file to a NumPy array.

    :param gzip_csv_file:                       str.  Full path with file name and extension to an GZIP CSV file.
    :param header_rows_to_skip:                 int.  The number of header rows to skip.

    :return:                                    NumPy array

    """

    with gzip.open(gzip_csv_file, mode='rt') as f:
        return np.genfromtxt(f, delimiter=',', skip_header=header_rows_to_skip)


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
    # try:
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
    # first_index = 502287

    # calculate distances between the first point and all other points within the user-defined neighborhood
    cut_off_metres = cut_off_meters + 1
    tree_1 = cKDTree(coordinate_array[first_index:first_index + 1, [0, 1]])

    tree_2 = cKDTree(coordinate_array[:, [0, 1]])
    tree_dist = cKDTree.sparse_distance_matrix(tree_1, tree_2, cut_off_metres, output_type='dict', p=2)

    tree_dist = {k: tree_dist[k] for k in sorted(tree_dist)}
    
    # put distances and indices of neighboring in a data frame
    dist_df = pd.DataFrame(columns=["near_id", dis_col])

    # add list for zip to accommodate Python 3
    dist_df["near_id"] = coordinate_array[list(zip(*tree_dist))[1], 2].astype(np.int32)
    dist_df[dis_col] = tree_dist.values()
    dist_df = dist_df.loc[dist_df.loc[:, dis_col] != 0, :]  # Remove the distance to itself
    
    # bring row and column indices of neighboring points by a join
    dist_df = dist_df.join(all_indices, on="near_id")

    # add to columns holding the relative difference in rows and columns between focal point and its neighbors
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
