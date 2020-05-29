"""
Process population projection by step for the popultation_gravity model

@author   Chris R. Vernon
@email:   chris.vernon@pnnl.gov

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

"""

import logging
import os
import pkg_resources
import time
import rasterio

import pandas as pd

import population_gravity.downscale_utilities as utils
from population_gravity.downscale_projection import pop_projection


class ProcessStep:
    """Process population downscaling for a single time step

    :param cfg:                 Configuration file object
    :param yr:                  Target year (YYYY) for a time step

    """

    def __init__(self, cfg, yr):

        # start time
        td = time.time()

        logging.info("Downscaling year:  {}".format(yr))

        self.cfg = cfg
        self.yr = yr

        # find which states to run based on the target state
        self.states_list = self.get_neighboring_states_list()

        if self.yr == self.cfg.projection_start_year:

            # run downscaling
            pop_projection(self.cfg, self.cfg.historical_urban_pop_raster, self.cfg.historical_rural_pop_raster,
                           self.cfg.alpha_urban, self.cfg.beta_urban, self.cfg.alpha_rural, self.cfg.beta_rural,
                           self.cfg.rural_pop_proj_n, self.cfg.urban_pop_proj_n, self.yr,
                           self.cfg.kernel_distance_meters)

        else:

            prev_step = self.yr - self.cfg.time_step

            # build raster list for all neighboring raster outputs from the previous time step
            urban_raster_list = [self.construct_file_name(i, prev_step, 'Urban') for i in self.cfg.neighbors]
            rural_raster_list = [self.construct_file_name(i, prev_step, 'Rural') for i in self.cfg.neighbors]

            # construct output mosaic file names
            urban_mosaic_file = self.construct_file_name(self.cfg.state_name, prev_step, 'Urban', suffix='_mosaic')
            rural_mosaic_file = self.construct_file_name(self.cfg.state_name, prev_step, 'Rural', suffix='_mosaic')

            # construct output mask file names
            urban_mask_file = self.construct_file_name(self.cfg.state_name, prev_step, 'Urban', suffix='_mask')
            rural_mask_file = self.construct_file_name(self.cfg.state_name, prev_step, 'Rural', suffix='_mask')

            # build mosaics
            urban_mosaic = utils.mosaic(urban_raster_list, urban_mosaic_file, self.cfg.metadata.copy())
            rural_mosaic = utils.mosaic(rural_raster_list, rural_mosaic_file, self.cfg.metadata.copy())

            # create a rasterio.windows.Window object from the bounding box of the target for the full mosaic
            target_window = urban_mosaic.window(*self.cfg.bbox)

            # write the urban masked mosiac raster to file
            self.mask_raster(urban_mask_file, urban_mosaic, target_window)

            # write the urban masked mosiac raster to file
            self.mask_raster(rural_mask_file, rural_mosaic, target_window)

            pop_projection(self.cfg, urban_mask_file, rural_mask_file, self.cfg.alpha_urban, self.cfg.beta_urban,
                           self.cfg.alpha_rural, self.cfg.beta_rural, self.cfg.rural_pop_proj_n,
                           self.cfg.urban_pop_proj_n, self.yr, self.cfg.kernel_distance_meters)

        logging.info("Downscaling for year {} completed in {} minutes.".format(yr, (time.time() - td) / 60))

    def mask_raster(self, file_name, raster_object, source_window):
        """Subset the mosiac using a window of the source historical raster and write to file.

        :param file_name:                   Full path with file name and extension to the output file
        :type file_name:                    str

        :param raster_object:               Rasterio object

        :param source_window:               A rasterio.windows.Window object from the bounding box of the target for
                                            the full mosaic

        """
        # read only the window of the mosaic to an array
        masked = raster_object.read(1, window=source_window)

        # write the mosiac raster to file
        with rasterio.open(file_name, 'w', **self.cfg.metadata) as dest:
            dest.write(masked, indexes=1)

        return masked

    def get_neighboring_states_list(self):
        """Get a list of all states within a 100 km distance from the target state border, including the target state,
        as a list.

        :return:                            List of state names that are lower case and underscore separated

        """
        # load neighboring state reference
        near_states_file = pkg_resources.resource_filename('population_gravity', 'data/neighboring_states_100km.csv')

        # find which states to run based on the target state
        states_df = pd.read_csv(near_states_file).apply(lambda x: x.str.lower().str.replace(' ', '_'))

        return states_df.groupby('target_state')['near_state'].apply(list).to_dict()[self.cfg.state_name]

    def construct_file_name(self, state_name, prev_step, designation, suffix='', extension='tif'):
        """Construct output file name.

        :param state_name:                  Name of state
        :type state_name:                   str

        :param prev_step:                   Previous time step.
        :type prev_step:                    int, str

        :param designation:                 Either 'urban', 'rural', or 'total'
        :type designation:                  str

        :param suffix:                      String to append to end of file name before extension
        :type suffix:                       str

        :param extension:                   Raster extension without the dot; default 'tif'
        :type extension:                    str

        :returns:                           Full file path for an output raster

        """
        fle = f"{state_name}_1km_{self.cfg.scenario}_{designation}_{prev_step}{suffix}.{extension}"

        return os.path.join(self.cfg.output_directory, fle)
