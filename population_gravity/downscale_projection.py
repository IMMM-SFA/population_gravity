"""
Population downscaling function for the population_gravity model

@author   Hamid Zoraghein, Chris R. Vernon
@email:   hzoraghein@popcouncil.org, chris.vernon@pnnl.gov

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

"""

import logging
import os

import numpy as np
import multiprocessing
from pathos.multiprocessing import ProcessingPool as Pool
from collections import deque

import population_gravity.downscale_utilities as utils


def pop_projection(cfg, urban_raster, rural_raster, alpha_urban, beta_urban, alpha_rural, beta_rural, rural_pop_proj_n,
                   urban_pop_proj_n, yr, cut_off_meters=100000):
    """Downscale population from state-level projections for urban and rural to 1 km gridded data.

    :param cfg:                             Configuration object
    :param urban_raster:                    string, Urban raster from historic data for base year
    :param rural_raster:                    string, Rural raster from historic data for base year
    :param alpha_urban:                     Alpha parameter for urban
    :param beta_urban:                      Beta parameter for urban
    :param alpha_rural:                     Alpha parameter for rural
    :param beta_rural:                      Beta parameter for rural
    :param rural_pop_proj_n:                Population number for rural for the projection year
    :param urban_pop_proj_n:                Population number for urban for the projection year
    :param cut_off_meters:                  Distance kernel in meters; default 100,000 meters

    """

    # Define local variables
    time_one_data = {}  # Dictionary storing the base year population grids
    population_1st = {}  # Dictionary storing the base year population arrays
    time_one_data['Rural'] = rural_raster  # Rural
    time_one_data['Urban'] = urban_raster  # Urban
    final_arrays = {}  # Dictionary containing final projected arrays
    final_raster = os.path.join(cfg.output_directory, "{}_1km_{}_Total_{}.tif".format(cfg.state_name, cfg.scenario, yr))

    # Calculate a distance matrix that serves as a template
    dist_matrix = utils.dist_matrix_calculator(cfg.one_dimension_indices[0], cut_off_meters, cfg.df_indicies, cfg.grid_coordinates_array)

    # Read historical urban and rural population grids into arrays
    for setting in time_one_data:

        rast_array = utils.raster_to_array(time_one_data[setting]).flatten()

        # this ensures that the no data value is filtered out otherwise the array will overflow on add
        rast_array = np.where(rast_array < 0.0, 0.0, rast_array)

        # create the dictionary containing population of each point in year 0
        population_1st[setting] = rast_array

    # create an array containing total population values in the first historical year
    total_population_1st = population_1st["Rural"] + population_1st["Urban"]

    # number of columns of population grid required to derive linear indices
    ind_diffs = dist_matrix["ind_diff"].values

    # distances between current point and its close points
    ini_dist = dist_matrix["dis"].values / 1000.0

    # downscale population projection
    for setting in time_one_data:

        # output raster path and file name with extension
        output = os.path.join(cfg.output_directory, "{}_1km_{}_{}_{}.tif".format(cfg.state_name, cfg.scenario, setting, yr))

        # calculate aggregate urban/rural population at time 1
        pop_first_year = population_1st[setting][cfg.one_dimension_indices]
        pop_t1 = pop_first_year.sum()

        # load the SSP file to retrieve the aggregated projected population at time 2 for downscaling
        if cfg.df_projected is None:
            if setting == "Urban":
                pop_t2 = urban_pop_proj_n
            else:
                pop_t2 = rural_pop_proj_n

        else:

            if setting == "Urban":
                pop_t2 = cfg.df_projected.loc[(cfg.df_projected["Year"] == yr) & (cfg.df_projected["Scenario"] == cfg.scenario),
                                      "UrbanPop"].values[0]
            else:
                pop_t2 = cfg.df_projected.loc[(cfg.df_projected["Year"] == yr) & (cfg.df_projected["Scenario"] == cfg.scenario),
                                      "RuralPop"].values[0]

        # population change between years 1 and 2
        pop_change = pop_t2 - pop_t1

        if pop_change < 0:
            negative_mod = 1
        else:
            negative_mod = 0

        if setting == "Urban":
            alpha_parameter = alpha_urban
            beta_parameter = -beta_urban  # note negative multiplier
        else:
            alpha_parameter = alpha_rural
            beta_parameter = -beta_rural  # note negative multiplier

        dist = beta_parameter * ini_dist
        exp_xx_inv_beta_dist = np.exp(dist)

        # initialize the parallelization
        pool = Pool(processes=multiprocessing.cpu_count())

        # provide the inputs for the parallelized function
        parallel_elements = deque([(i, ind_diffs, total_population_1st, alpha_parameter, exp_xx_inv_beta_dist)
                                   for i in cfg.one_dimension_indices])

        # derive suitability estimates
        suitability_estimates = pool.map(utils.suitability_estimator, parallel_elements)

        # change suitability estimates to a numpy array
        suitability_estimates = np.array(suitability_estimates)

        # extract only the necessary mask values that fall within the state boundary
        cur_points_mask = cfg.historical_suitability_array[cfg.one_dimension_indices]

        # in case of population decline, suitability estimates are reciprocated
        if negative_mod:

            # find those whose mask is 0 but have population, they should decline anyway
            mask_zero = np.where(cur_points_mask == 0)[0]
            pop_non_zero = np.where(pop_first_year != 0)[0]

            # those cells with mask value of zero and population are the intersection of the two above arrays
            pop_mask = np.intersect1d(mask_zero, pop_non_zero, assume_unique=True)

            # change the mask value of the above cells to 1 so that they also lose population
            cur_points_mask[pop_mask] = cur_points_mask.mean()

            # adjust suitability values by applying mask values
            suitability_estimates = cur_points_mask * suitability_estimates

            # inverse current mask values for a better reflection of population decline
            suitability_estimates[suitability_estimates != 0] = 1.0 / suitability_estimates[suitability_estimates != 0]

        else:
            # adjust suitability values by applying its mask values
            suitability_estimates = cur_points_mask * suitability_estimates

        # total suitability for the whole area, which is the summation of all individual suitability values
        tot_suitability = suitability_estimates.sum()

        # tinal population estimates if negative mode is off
        pop_estimates = suitability_estimates / tot_suitability * pop_change + pop_first_year

        # adjust the projection so that no cell can have less than 0 individuals.
        if negative_mod:

            # to make sure that there is no negative population
            while any(pop < 0 for pop in pop_estimates):

                # treating negative population values
                extra_pop_mod = abs(pop_estimates[pop_estimates < 0].sum())
                pop_estimates[pop_estimates < 0] = 0

                # calculate the new total suitability value based on points with positive projected population
                new_tot_suitability = suitability_estimates[pop_estimates > 0].sum()

                # adjust non-negative population values to maintain the total aggregated population
                pop_estimates[pop_estimates > 0] = pop_estimates[pop_estimates > 0] - (
                            suitability_estimates[pop_estimates > 0] / new_tot_suitability) * extra_pop_mod

        # save the projection array of the current setting
        final_arrays[setting] = pop_estimates

        # save the final urban, rural raster
        logging.info("Saving {} raster to:  {}".format(setting, output))

        utils.array_to_raster(cfg.historical_suitability_raster, pop_estimates, cfg.one_dimension_indices, output)

    # calculate the total population array
    total_array = final_arrays["Rural"] + final_arrays["Urban"]

    # save the final total raster
    logging.info("Saving total population raster to:  {}".format(final_raster))

    utils.array_to_raster(cfg.historical_suitability_raster, total_array, cfg.one_dimension_indices, final_raster)
