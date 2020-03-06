"""
Calibration module for population downscaling.

@author: Hamidreza Zoraghein

"""


import logging
import os
import pickle

import simplejson
import numpy as np
import pandas as pd
import scipy.optimize

import population_gravity.downscale_utilities as pdm


def build_iterator(obj_x, obj_y):
    """Build a cross-tabulated iterator"""
    l = []
    for x in obj_x:
        for y in obj_y:
            l.append((x, y))

    return l


def final_optimization(df, params, parameters_dict, bounds, setting):

    # use the point with the minimum value as an initial guess for the second optimizer
    (a0, b0) = df.loc[df["estimate"].idxmin(), ["a", "b"]]

    # final optimization
    parameters = scipy.optimize.minimize(pdm.pop_min_function, x0=(a0, b0), args=params, method="SLSQP",
                                         tol=0.001, options={"disp": True}, bounds=bounds)

    logging.info("\tOptimization Outcomes for {}: {}, {}".format(setting, parameters["x"], parameters["fun"]))

    parameters_dict[setting] = parameters["x"]

    return parameters_dict


def calibration(cfg):
    """
    Calibration for population ?

    :param cfg:                         Configuration object containing all settings
    :param urb_pop_fst_year:            ?
    :param urb_pop_snd_year:            ?
    :param rur_pop_fst_year:            ?
    :param rur_pop_snd_year:            ?
    :param mask_raster:                 ?
    :param region_code:                 ?
    :param ssp_code:                    ?
    :param point_indices:               ?
    :param point_coors:                 ?
    :param datadir_output:              ?

    """

    # set variables
    region_code = cfg.region_code
    ssp_code = cfg.ssp_code

    urb_outfile = os.path.join(cfg.datadir_output, '{}_{}_urban_calib_output.pkl'.format(cfg.region_code.replace(' ', '_'), cfg.ssp_code))
    rur_outfile = os.path.join(cfg.datadir_output, '{}_{}_rural_calib_output.pkl'.format(cfg.region_code.replace(' ', '_'), cfg.ssp_code))

    # Read indices of points that fall within the state boundary
    with open(cfg.point_indices, 'r') as r:
        within_indices = simplejson.load(r)

    arr_coords = np.genfromtxt(cfg.point_coors, delimiter=',', skip_header=1, usecols=(0, 1, 2), dtype=float)

    # A csv file to write the parameters into
    out_cal = os.path.join(cfg.datadir_output, '{}_{}_Parmas.csv'.format(region_code, ssp_code))

    # Dictionary storing initial urban and rural population grids
    all_rasters = {'Rural': [cfg.rur_pop_fst_year, cfg.rur_pop_snd_year],
                   'Urban': [cfg.urb_pop_fst_year, cfg.urb_pop_fst_year]}

    # define local variables
    parameters_dict = {}  # Dictionary storing urban/rural calibration parameters

    # populate the array containing mask values
    points_mask = pdm.raster_to_array(cfg.mask_raster).flatten()

    # Create an array containing total population values in the first historical year
    arr_pop_rur_1st = pdm.raster_to_array(cfg.rur_pop_fst_year).flatten()
    arr_pop_rur_2nd = pdm.raster_to_array(cfg.rur_pop_snd_year).flatten()
    arr_pop_urb_1st = pdm.raster_to_array(cfg.urb_pop_fst_year).flatten()
    arr_pop_urb_2nd_2D = pdm.raster_to_array(cfg.urb_pop_snd_year)
    arr_pop_urb_2nd = arr_pop_urb_2nd_2D.flatten()
    arr_pop_tot_1st = arr_pop_rur_1st + arr_pop_urb_1st

    # All indices
    df_all_indices = pdm.all_index_retriever(arr_pop_urb_2nd_2D, ["row", "column"])

    # Calculate a distance matrix that serves as a template
    DIST_THRESHOLD_METERS = 100000

    dist_matrix = pdm.dist_matrix_calculator(within_indices[0], DIST_THRESHOLD_METERS, df_all_indices, arr_coords)

    # initial alpha values
    a_lower = -2.0
    a_upper = 2.0

    # initial beta values
    b_lower = -2.0
    b_upper = 4.0

    bounds = ((a_lower, a_upper), (b_lower, b_upper))

    # parameters to be used in optimization evenly distributed from lower to upper bound
    a_list = np.linspace(a_lower, a_upper, 8)
    b_list = np.linspace(b_lower, b_upper, 10)
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

        params = (setting, arr_1st, arr_2nd, arr_pop_tot_1st, points_mask, dist_matrix, within_indices)

        # initialize the dataframe that will hold values of the brute force
        fst_results = pd.DataFrame(data={"a": np.repeat(a_list, 10).astype(np.float32),
                                         "b": np.tile(b_list, 8).astype(np.float32),
                                         "estimate": np.full(80, np.nan, dtype=np.float64)
                                         })

        # run brute force to calculate optimization per grid point
        for index, a_b in enumerate(ab_iter):

            estimate = pdm.pop_min_function(a_b, arr_1st, arr_2nd, arr_pop_tot_1st, points_mask, dist_matrix, within_indices)

            fst_results.loc[(fst_results["a"] == a_b[0]) & (fst_results["b"] == a_b[1]), "estimate"] = estimate

            # logging.info("\t{}th iteration done".format(index))

        # pickle the current optimization files
        # pickle_file = '/Users/d3y010/Desktop/fst_results_urb.pik'
        # with open(pickle_file, "wb") as pickle_params:
        #     pickle.dump(fst_results, pickle_params)

        parameters_dict = final_optimization(fst_results, params[1:], parameters_dict, bounds, setting)

    # write the parameters to the designated csv file
    logging.info("\tWriting parameterization file:  {}".format(out_cal))

    with open(out_cal, 'w') as out_csv:
        out_csv.write("Region,SSP,Alpha_Rural,Beta_Rural,Alpha_Urban,Beta_Urban,Comments\n")
        out_csv.write('{},{},{},{},{},{},{}\n'.format(region_code, ssp_code, parameters_dict["Rural"][0],
                                                      parameters_dict["Rural"][1], parameters_dict["Urban"][0],
                                                      parameters_dict["Urban"][1], ""))
