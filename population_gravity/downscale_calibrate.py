"""
Population calibration function for the population_gravity model

@author   Hamid Zoraghein, Chris R. Vernon
@email:   hzoraghein@popcouncil.org, chris.vernon@pnnl.gov

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

"""

import os
import logging
import simplejson
import pickle

import numpy as np
import pandas as pd
import scipy.optimize
from scipy.spatial import cKDTree

from population_gravity.sampler import equal_interval
import population_gravity.downscale_utilities as pdm


def all_index_retriever(array, columns):
    """Build data frame in the shape of the input array.

    :param array:                   Input 2D array from input raster
    :param columns:                 list; Target columns
    :param row_col:                 Name of row column
    :param column_col:              Name of column column
    :param all_index_col:           Name of `all_index` column

    :return:                        Typed data frame of indicies

    """

    # Put the row, column and linear indices of all elements in a dataframe
    shape = array.shape
    index = pd.MultiIndex.from_product([range(s) for s in shape], names=columns)
    df = pd.DataFrame({'all_index': array.flatten()}, index=index).reset_index()
    df["all_index"] = df.index

    df = df.astype({"row": np.int32, "column": np.int32, "all_index": np.int32})
    return df


def dist_matrix_calculator(first_index, cut_off_meters, all_indices, coors_csv_file):

    # Read all points with their coordinates
    points = np.genfromtxt(coors_csv_file, delimiter=',', skip_header=1, usecols=(0, 1, 2), dtype=float)

    # Calculate distances between the first point and all other points within a 100km neighborhood
    cut_off_metres = cut_off_meters + 1
    tree_1 = cKDTree(points[first_index:first_index + 1, [0, 1]])
    tree_2 = cKDTree(points[:, [0, 1]])
    tree_dist = cKDTree.sparse_distance_matrix(tree_1, tree_2, cut_off_metres, output_type='dict', p=2)

    # Put distances and indices of neighboring in a dataframe
    dist_df = pd.DataFrame(columns=["near_id", "dis"])
    dist_df["near_id"] = points[list(zip(*tree_dist))[1], 2].astype(np.int32)
    dist_df["dis"] = tree_dist.values()
    dist_df = dist_df.loc[dist_df.loc[:, "dis"] != 0, :]  # Remove the distance to itself

    # Bring row and column indices of neighboring points by a join
    dist_df = dist_df.join(all_indices, on="near_id")

    # Add to columns holding the relative difference in rows and colums beween focal point and its neighbors
    foc_indices = all_indices.loc[first_index, ["row", "column"]].values
    dist_df["ind_diff"] = dist_df["near_id"] - first_index
    dist_df["row_diff"] = dist_df["row"] - foc_indices[0]
    dist_df["col_diff"] = dist_df["column"] - foc_indices[1]

    # Drop unwanted columns
    dist_df = dist_df.drop(["row", "column", "near_id", "all_index"], axis=1)

    dist_df = dist_df.astype({"ind_diff": np.int32, "row_diff": np.int32, "col_diff": np.int32})

    return dist_df


def calibration(cfg):
    """Calibrate alpha and beta parameters."""

    # Define local variables
    all_rasters = {}  # Dictionary storing initial urban and rural population grids
    rur_pop_files = []  # List containing rural population grids
    urb_pop_files = []  # List containing urban population grids
    population_1st = {}  # Dictionary containing population of each point in year 1
    population_2nd = {}  # Dictionary containing population of each point in year 2
    parameters_dict = {}  # Dictionary storing urban/rural calibration parameters

    # Rural
    rur_pop_files.append(cfg.calibration_rural_year_one_raster)
    rur_pop_files.append(cfg.calibration_rural_year_two_raster)

    # Urban
    urb_pop_files.append(cfg.calibration_urban_year_one_raster)
    urb_pop_files.append(cfg.calibration_urban_year_two_raster)

    # Urban and rural
    all_rasters["Rural"] = rur_pop_files
    all_rasters["Urban"] = urb_pop_files

    # Populate the array containing mask values
    points_mask = pdm.raster_to_array(cfg.historical_suitability_raster)

    # Read historical urban and rural population grids into arrays
    for setting in all_rasters:

        # Create the dictionary containing population of each point in year 1
        population_1st[setting] = pdm.raster_to_array(all_rasters[setting][0])

        # Create the dictionary containing population of each point in year 2
        population_2nd[setting] = pdm.raster_to_array(all_rasters[setting][1])

    # Create an array containing total population values in the first historical year
    total_population_1st = population_1st["Rural"] + population_1st["Urban"]

    # A csv file to write the parameters into
    out_cal = cfg.output_directory + cfg.state_name + "_" + cfg.scenario + "_Params.csv"

    # All indices
    arr_upop_2d = pdm.raster_to_array(cfg.calibration_urban_year_two_raster)
    all_indices = all_index_retriever(arr_upop_2d, ["row", "column"])

    # Read indices of points that fall within the state boundary
    with open(cfg.one_dimension_indices_file, 'r') as r:
        within_indices = simplejson.load(r)

    # Calculate a distance matrix that serves as a template
    dist_matrix = dist_matrix_calculator(within_indices[0], cfg.kernel_distance_meters, all_indices, cfg.grid_coordinates_file)

    # Parameter calculation for rural and urban
    for setting in all_rasters:

        if setting.lower() == 'urban':
            a_lower = cfg.urban_alpha_lower_bound
            a_upper = cfg.urban_alpha_upper_bound
            b_lower = cfg.urban_beta_lower_bound
            b_upper = cfg.urban_beta_upper_bound

        else:
            a_lower = cfg.rural_alpha_lower_bound
            a_upper = cfg.rural_alpha_upper_bound
            b_lower = cfg.rural_beta_lower_bound
            b_upper = cfg.rural_beta_upper_bound

        # build sample
        a_list, b_list, fst_results = equal_interval(alpha_lower=a_lower,
                                                     alpha_upper=a_upper,
                                                     beta_lower=b_lower,
                                                     beta_upper=b_upper,
                                                     alpha_intervals=8,
                                                     beta_intervals=10)

        # extract arrays
        arr_1st = population_1st[setting]
        arr_2nd = population_2nd[setting]

        params = (arr_1st.flatten(), arr_2nd.flatten(), total_population_1st.flatten(), points_mask.flatten(), dist_matrix, within_indices)

        # Run brute force to calculate optimization per grid point
        i = 0
        for a in a_list:
            for b in b_list:
                local_result = pdm.pop_min_function((a, b), *params)
                fst_results.loc[(fst_results["alpha_param"] == a) & (fst_results["beta_param"] == b), "estimate"] = local_result
                i += 1

        # Pickle the current optimization file
        pickle_file = os.path.join(cfg.output_directory, f"local_optimization_{cfg.state_name.lower()}_{setting.lower()}.pkl")
        fst_results.to_pickle(pickle_file)

        # Use the point with the minimum value as an initial guess for the second optimizer
        (a0, b0) = fst_results.loc[fst_results["estimate"].idxmin(), ["alpha_param", "beta_param"]]

        # Final optimization
        rranges = ((a_lower, a_upper), (b_lower, b_upper))
        parameters = scipy.optimize.minimize(pdm.pop_min_function, x0=(a0, b0), args=params, method="SLSQP",
                                             tol=0.001, options={"disp": True}, bounds=rranges)

        parameters_dict[setting] = parameters["x"]

    out_string = f"{cfg.state_name},{cfg.scenario},{parameters_dict['Rural'][0]},{parameters_dict['Rural'][1]},{parameters_dict['Urban'][0]},{parameters_dict['Urban'][1]}"
    logging.info(f"Final optimization:  {out_string}")

    # Write the parameters to the designated csv file
    with open(out_cal, 'w') as out:
        out.write("Region,SSP,Alpha_Rural,Beta_Rural,Alpha_Urban,Beta_Urban,Comments\n")
        out.write(f"{out_string}\n")

    return parameters_dict['Urban'][0], parameters_dict['Rural'][0], parameters_dict['Urban'][1], parameters_dict['Rural'][1]
