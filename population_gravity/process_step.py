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

import numpy as np
import pandas as pd

from rasterio.io import MemoryFile

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
            urban_output, rural_output = pop_projection(self.cfg, self.cfg.historical_urban_pop_raster,
                                                        self.cfg.historical_rural_pop_raster,
                                                        self.cfg.alpha_urban, self.cfg.beta_urban, self.cfg.alpha_rural,
                                                        self.cfg.beta_rural,
                                                        self.cfg.rural_pop_proj_n, self.cfg.urban_pop_proj_n, self.yr,
                                                        self.cfg.kernel_distance_meters)

        else:

            prev_step = self.yr - self.cfg.time_step

            # build raster list for all neighboring raster outputs from the previous time step
            urban_file_list = self.construct_file_list(prev_step, 'Urban')
            rural_file_list = self.construct_file_list(prev_step, 'Rural')

            # build mosaics
            urban_mosaic = utils.mosaic_memory(urban_file_list, self.cfg.metadata.copy())
            rural_mosaic = utils.mosaic_memory(rural_file_list, self.cfg.metadata.copy())

            # create a rasterio.windows.Window object from the bounding box of the target for the full mosaic
            target_window = urban_mosaic.window(*self.cfg.bbox)

            # write the urban masked mosiac raster to file
            urban_mask_file = self.mask_raster_memory(urban_mosaic, target_window)
            urban_mosaic.close()

            # write the urban masked mosiac raster to file
            rural_mask_file = self.mask_raster_memory(rural_mosaic, target_window)
            rural_mosaic.close()

            urban_output, rural_output = pop_projection(self.cfg, urban_mask_file, rural_mask_file,
                                                        self.cfg.alpha_urban, self.cfg.beta_urban,
                                                        self.cfg.alpha_rural, self.cfg.beta_rural,
                                                        self.cfg.rural_pop_proj_n,
                                                        self.cfg.urban_pop_proj_n, self.yr,
                                                        self.cfg.kernel_distance_meters)

            # close in memory mask objects
            urban_mask_file.close()
            rural_mask_file.close()

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

    def mask_raster_memory(self, raster_object, source_window):
        """Subset the mosiac using a window of the source historical raster and write to memory.

        :param raster_object:               Rasterio object

        :param source_window:               A rasterio.windows.Window object from the bounding box of the target for
                                            the full mosaic

        """
        # read only the window of the mosaic to an array
        masked = raster_object.read(1, window=source_window)

        # write output to memory
        with MemoryFile() as memfile:
            dataset = memfile.open(**self.cfg.metadata)
            dataset.write(masked, indexes=1)

            return dataset

    def get_neighboring_states_list(self):
        """Get a list of all states within a 100 km distance from the target state border, including the target state,
        as a list.

        :return:                            List of state names that are lower case and underscore separated

        """
        # load neighboring state reference
        near_states_file = os.path.join(os.path.dirname(self.cfg.historical_rural_pop_raster), 'neighboring_states_100km.csv')

        # find which states to run based on the target state
        states_df = pd.read_csv(near_states_file).apply(lambda x: x.str.lower().str.replace(' ', '_'))

        return states_df.groupby('target_state')['near_state'].apply(list).to_dict()[self.cfg.state_name]

    def construct_file_name(self, state_name, prev_step, designation, suffix='', extension='.tif'):
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

        if (len(suffix) > 0) or (self.cfg.state_name != state_name):

            fle = f"{state_name}_1km_{self.cfg.scenario}_{designation}_{prev_step}{suffix}{extension}"

        else:

            if type(self.cfg.run_number) == int:
                delim = '_'
            else:
                delim = ''

            fle = f"{state_name}_1km_{self.cfg.scenario}_{designation}_{prev_step}{delim}{self.cfg.run_number}{suffix}{extension}"

        return os.path.join(self.cfg.output_directory, fle)

    def construct_file_list(self, prev_step, setting):
        """Construct a list of arrays from rasters or arrays.

        :param prev_step:                       int.  Previous time step; e.g., year
        :param setting:                         str.  Either 'Urban' or 'Rural'

        #TODO:  load prev years CSV files

        """

        out_list = []

        for i in self.cfg.neighbors:

            # check for either the 'tif', 'npy', or 'csv' files in the output directory, use 'tif' first
            tif = self.construct_file_name(i, prev_step, setting, extension='.tif')
            npy1d = self.construct_file_name(i, prev_step, setting, extension='_1d.npy')
            npy2d = self.construct_file_name(i, prev_step, setting, extension='_2d.npy')
            csv_gz = self.construct_file_name(i, prev_step, setting, extension='.csv.gz')

            if os.path.isfile(tif):
                logging.info(f"Using file '{tif}' for previous time step data.")
                out_list.append(rasterio.open(tif))

            elif os.path.isfile(npy1d):
                logging.info(f"Using file '{npy1d}' for previous time step data.")

                # load npy file to array
                array1d = np.load(npy1d)
                out_list.append(utils.array_to_raster_memory(self.cfg.template_raster, array1d, self.cfg.one_dimension_indices))

            elif os.path.isfile(npy2d):
                logging.info(f"Using file '{npy2d}' for previous time step data.")

                # load npy file to array
                array2d = np.load(npy2d)
                out_list.append(utils.array2d_to_raster_memory(array2d, raster_profile=self.cfg.template_raster[4]))

            elif os.path.isfile(csv_gz):
                logging.info(f"Using file '{csv_gz}' for previous time step data.")

                array1d = pd.read_csv(csv_gz, compression='gzip', sep=',')['value'].values
                out_list.append(utils.array_to_raster_memory(self.cfg.template_raster, array1d, self.cfg.one_dimension_indices))

        return out_list
