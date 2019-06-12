import csv
import logging
import os
import pickle

import simplejson
import numpy as np
import pandas as pd
import scipy.optimize

import spatial_downscale.downscale_utilities as pdm


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

    :return:

    """

    # set variables
    urb_pop_fst_year = cfg.urb_pop_fst_year
    urb_pop_snd_year = cfg.urb_pop_snd_year
    rur_pop_fst_year = cfg.rur_pop_fst_year
    rur_pop_snd_year = cfg.rur_pop_snd_year
    mask_raster = cfg.mask_raster
    region_code = cfg.region_code
    ssp_code = cfg.ssp_code
    point_indices = cfg.point_indices
    point_coors = cfg.point_coors
    datadir_output = cfg.datadir_output

    # define local variables
    all_rasters = {}  # Dictionary storing initial urban and rural population grids
    rur_pop_files = []  # List containing rural population grids
    urb_pop_files = []  # List containing urban population grids
    population_1st = {}  # Dictionary containing population of each point in year 1
    population_2nd = {}  # Dictionary containing population of each point in year 2
    parameters_dict = {}  # Dictionary storing urban/rural calibration parameters

    # rural
    rur_pop_files.append(rur_pop_fst_year)
    rur_pop_files.append(rur_pop_snd_year)

    # urban
    urb_pop_files.append(urb_pop_fst_year)
    urb_pop_files.append(urb_pop_snd_year)

    # urban and rural
    all_rasters["Rural"] = rur_pop_files
    all_rasters["Urban"] = urb_pop_files

    # populate the array containing mask values
    points_mask = pdm.raster_to_array(mask_raster)

    # read historical urban and rural population grids into arrrays
    for setting in all_rasters:

        # create the dictionary containing population of each point in year 1
        population_1st[setting] = pdm.raster_to_array(all_rasters[setting][0])

        # create the dictionary containing population of each point in year 2
        population_2nd[setting] = pdm.raster_to_array(all_rasters[setting][1])

    # Create an array containing total population values in the first historical year
    total_population_1st = population_1st["Rural"] + population_1st["Urban"]

    # A csv file to write the parameters into
    out_cal = datadir_output + region_code + "_" + ssp_code + "_Params.csv"

    # All indices
    all_indices = pdm.all_index_retriever(urb_pop_snd_year, ["row", "column"])

    # Read indices of points that fall within the state boundary
    with open(point_indices, 'r') as r:
        within_indices = simplejson.load(r)

    # Calculate a distance matrix that serves as a template
    cut_off_meters = 100000
    dist_matrix = pdm.dist_matrix_calculator(within_indices[0], cut_off_meters, all_indices, point_coors)

    # Parameter calculation for rural and urban
    for setting in all_rasters:
        # initial alpha values
        a_lower = -2.0
        a_upper = 2.0
        # initial beta values
        b_lower = -2.0
        b_upper = 4.0

        rranges = ((a_lower, a_upper), (b_lower, b_upper))

        # parameters to be used in optimization
        a_list = np.linspace(a_lower, a_upper, 8)
        b_list = np.linspace(b_lower, b_upper, 10)
        params = (setting, population_1st, population_2nd, total_population_1st,
                  points_mask, dist_matrix, within_indices)

        # initialize the dataframe that will hold values of the brute force
        fst_results = pd.DataFrame(data={"a": np.repeat(a_list, 10).astype(np.float32),
                                         "b": np.tile(b_list, 8).astype(np.float32),
                                         "estimate": np.full(80, np.nan, dtype=np.float64)
                                         }
                                   )

        # run brute force to calculate optimization per grid point
        # TODO:  we should be able to vectorize this
        i = 0
        for a in a_list:
            for b in b_list:
                fst_results.loc[(fst_results["a"] == a) & (fst_results["b"] == b),
                                "estimate"] = pdm.pop_min_function((a, b), *params)

                logging.info("\t{}th iteration done".format(i))

                i += 1

        # pickle the current optimization file
        pickle_file = os.path.join(datadir_output, setting)
        with open(pickle_file, "wb") as pickle_params:
            pickle.dump(fst_results, pickle_params)

        # use the point with the minimum value as an initial guess for the second optimizer
        (a0, b0) = fst_results.loc[fst_results["estimate"].idxmin(), ["a", "b"]]

        # final optimization
        parameters = scipy.optimize.minimize(pdm.pop_min_function, x0=(a0, b0), args=params, method="SLSQP",
                                             tol=0.001, options={"disp": True}, bounds=rranges)

        logging.info("\tOptimization Outcomes for {}: {}, {}".format(setting, parameters["x"], parameters["fun"]))

        parameters_dict[setting] = parameters["x"]

    # write the parameters to the designated csv file
    logging.info("\tWriting parameterization file:  {}".format(out_cal))

    with open(out_cal, 'wb') as out_csv:
        writer = csv.writer(out_csv)
        writer.writerow(["Region", "SSP", "Alpha_Rural", "Beta_Rural", "Alpha_Urban", "Beta_Urban", "Comments"])
        writer.writerow([region_code, ssp_code, parameters_dict["Rural"][0], parameters_dict["Rural"][1],
                         parameters_dict["Urban"][0], parameters_dict["Urban"][1], ""])