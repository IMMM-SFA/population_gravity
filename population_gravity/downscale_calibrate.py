"""
Population calibration function for the population_gravity model

@author   Hamid Zoraghein, Chris R. Vernon
@email:   hzoraghein@popcouncil.org, chris.vernon@pnnl.gov

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

"""

import os
import logging

import numpy as np
import pandas as pd
import scipy.optimize

import population_gravity.downscale_utilities as utils


def build_iterator(obj_x, obj_y):
    """Build a cross-tabulated iterator"""
    l = []
    for x in obj_x:
        for y in obj_y:
            l.append((x, y))

    return l


def final_optimization(df, params, parameters_dict, bounds, setting):
    """
    TODO: Fill in docs

    """
    # use the point with the minimum value as an initial guess for the second optimizer
    (a0, b0) = df.loc[df["estimate"].idxmin(), ["a", "b"]]

    # final optimization
    parameters = scipy.optimize.minimize(utils.pop_min_function, x0=(a0, b0), args=params, method="SLSQP",
                                         tol=0.001, options={"disp": True, 'eps': 0.01, 'ftol': 0.01},
                                         bounds=bounds)

    logging.info("\tOptimization Outcomes for {}: {}, {}".format(setting, parameters["x"], parameters["fun"]))

    parameters_dict[setting] = parameters["x"]

    return parameters_dict


def calibration(cfg):
    """Calibrate alpha and beta parameters for the gravity model for a target state.

    :param cfg:                             Configuration object

    :returns:                               [0] alpha_urban
                                            [1] alpha_rural
                                            [2] beta_urban
                                            [3] beta_rural

    """

    # A csv file to write the parameters into
    out_cal = os.path.join(cfg.output_directory, '{}_{}_cablibration_parameters.csv'.format(cfg.state_name, cfg.scenario))

    # Dictionary storing initial urban and rural population grids
    all_rasters = {'Rural': [cfg.calibration_rural_year_one_raster, cfg.calibration_rural_year_two_raster],
                   'Urban': [cfg.calibration_urban_year_one_raster, cfg.calibration_urban_year_two_raster]}

    # dictionary storing urban and rural calibration parameters
    parameters_dict = {}

    # Create an array containing total population values in the first historical year
    arr_pop_rur_1st = utils.raster_to_array(cfg.calibration_rural_year_one_raster).flatten()
    arr_pop_rur_2nd = utils.raster_to_array(cfg.calibration_rural_year_two_raster).flatten()
    arr_pop_urb_1st = utils.raster_to_array(cfg.calibration_urban_year_one_raster).flatten()

    arr_pop_urb_2nd_2D = utils.raster_to_array(cfg.calibration_urban_year_two_raster)
    arr_pop_urb_2nd = arr_pop_urb_2nd_2D.flatten()

    arr_pop_tot_1st = arr_pop_rur_1st + arr_pop_urb_1st

    # All indices
    df_all_indices = utils.all_index_retriever(arr_pop_urb_2nd_2D, ["row", "column"])

    dist_matrix = utils.dist_matrix_calculator(cfg.one_dimension_indices[0], cfg.kernel_distance_meters, df_all_indices,
                                               cfg.grid_coordinates_array)

    # initial alpha values
    a_lower = cfg.pass_one_alpha_lower
    a_upper = cfg.pass_one_alpha_upper

    # initial beta values
    b_lower = cfg.pass_one_beta_lower
    b_upper = cfg.pass_one_beta_upper

    # parameters to be used in optimization evenly distributed from lower to upper bound
    a_list = np.linspace(a_lower, a_upper, cfg.brute_n_alphas)
    b_list = np.linspace(b_lower, b_upper, cfg.brute_n_betas)
    ab_iter = build_iterator(a_list, b_list)

    # Parameter calculation for rural and urban
    for setting in all_rasters.keys():  # this is a dictionary of lists

        logging.info("Processing: {}".format(setting))

        if setting == 'Urban':
            arr_1st = arr_pop_urb_1st
            arr_2nd = arr_pop_urb_2nd
        else:
            arr_1st = arr_pop_rur_1st
            arr_2nd = arr_pop_rur_2nd

        params = (setting, arr_1st, arr_2nd, arr_pop_tot_1st, cfg.historical_suitability_array, dist_matrix,
                  cfg.one_dimension_indices)

        # initialize the data frame that will hold values of the brute force
        fst_results = pd.DataFrame(data={"a": np.repeat(a_list, cfg.brute_n_betas).astype(np.float32),
                                         "b": np.tile(b_list, cfg.brute_n_alphas).astype(np.float32),
                                         "estimate": np.full((cfg.brute_n_alphas * cfg.brute_n_betas), np.nan, dtype=np.float64)
                                         })

        # run brute force to calculate optimization per grid point
        for index, a_b in enumerate(ab_iter):

            estimate = utils.pop_min_function(a_b, arr_1st, arr_2nd, arr_pop_tot_1st, cfg.historical_suitability_array,
                                              dist_matrix, cfg.one_dimension_indices)

            fst_results.loc[(fst_results["a"] == a_b[0]) & (fst_results["b"] == a_b[1]), "estimate"] = estimate

        bounds = ((cfg.pass_two_alpha_lower, cfg.pass_two_alpha_upper), (cfg.pass_two_beta_lower, cfg.pass_two_beta_upper))
        parameters_dict = final_optimization(fst_results, params[1:], parameters_dict, bounds, setting)

    # write the parameters to the designated csv file
    logging.info("\tWriting parameterization file:  {}".format(out_cal))

    alpha_urban = parameters_dict["Urban"][0]
    alpha_rural = parameters_dict["Rural"][0]
    beta_urban = parameters_dict["Urban"][1]
    beta_rural = parameters_dict["Rural"][1]

    with open(out_cal, 'w') as out_csv:
        out_csv.write("Region,SSP,Alpha_Rural,Beta_Rural,Alpha_Urban,Beta_Urban\n")
        out_csv.write('{},{},{},{},{},{}\n'.format(cfg.state_name, cfg.scenario, alpha_rural, beta_rural, alpha_urban, beta_urban))

    return alpha_urban, alpha_rural, beta_urban, beta_rural
