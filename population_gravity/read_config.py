import os
import yaml

import numpy as np
import pandas as pd

import population_gravity.downscale_utilities as utils


class ReadConfig:

    def __init__(self, config_file=None, datadir_histdata=None, ssp_data_directory=None, ssp_code=None,
                 region_code=None, output_directory=None, calibration_parameters_file=None, future_start_year=None,
                 future_end_year=None, time_step=None, alpha_urban=None, beta_urban=None, alpha_rural=None, beta_rural=None,
                 rural_pop_proj_n=None, urban_pop_proj_n=None, historic_base_year=None):

        if config_file is None:

            self.datadir_histdata = datadir_histdata
            self.datadir_future = ssp_data_directory
            self.ssp_code = ssp_code
            self.region_code = region_code
            self.datadir_output = output_directory
            self.calibration_parameters_file = calibration_parameters_file
            self.future_end_year = future_end_year
            self.historic_base_year = historic_base_year
            self.future_start_year = future_start_year

            self.time_step = time_step

            # urban population projection number
            self.urban_pop_proj_n = urban_pop_proj_n

            # urban population projection number
            self.rural_pop_proj_n = rural_pop_proj_n

            if self.calibration_parameters_file is None:
                # calibration parameters
                self.alpha_urban = alpha_urban
                self.beta_urban = beta_urban
                self.alpha_rural = alpha_rural
                self.beta_rural = beta_rural

            else:
                # unpack calibration parameters
                self.alpha_urban, self.beta_urban, self.alpha_rural, self.beta_rural = self.unpack_calibration_parameters()
        else:

            # extract config file to YAML object
            cfg = self.get_yaml(config_file)

            # directory of main inputs
            self.datadir_histdata = cfg['datadir_histdata']

            # directory of projections
            self.datadir_future = cfg['ssp_data_directory']

            # SSP
            self.ssp_code = cfg['ssp_code']

            # region code
            self.region_code = cfg['region_code']

            # output directory
            self.datadir_output = cfg['output_directory']

            self.historic_base_year = cfg['historic_base_year']

            # start year
            self.future_start_year = cfg['future_start_year']

            # end year
            self.future_end_year = cfg['future_end_year']

            # interval
            self.time_step = cfg['time_step']

            # calibration parameters file
            self.calibration_parameters_file = cfg['calibration_parameters_file']

            # unpack calibration parameters
            self.alpha_urban, self.beta_urban, self.alpha_rural, self.beta_rural = self.unpack_calibration_parameters()

        self.steps = range(self.future_start_year, self.future_end_year + self.time_step, self.time_step)

        # calibration inputs -- derived from user inputs
        self.urb_pop_fst_year = os.path.join(self.datadir_histdata, "{}_urban_{}_1km.tif".format(self.region_code, self.historic_base_year))
        self.urb_pop_snd_year = os.path.join(self.datadir_histdata, "{}_urban_{}_1km.tif".format(self.region_code, self.historic_base_year + self.time_step))
        self.rur_pop_fst_year = os.path.join(self.datadir_histdata, "{}_rural_{}_1km.tif".format(self.region_code, self.historic_base_year))
        self.rur_pop_snd_year = os.path.join(self.datadir_histdata, "{}_rural_{}_1km.tif".format(self.region_code, self.historic_base_year + self.time_step))
        self.mask_raster_file = os.path.join(self.datadir_histdata, "{}_mask_short_term.tif".format(self.region_code))
        self.point_indices = os.path.join(self.datadir_histdata, "{}_within_indices.txt".format(self.region_code))
        self.point_coordinates_file = os.path.join(self.datadir_histdata, "{}_coordinates.csv".format(self.region_code))
        self.point_coordinates_array = np.genfromtxt(self.point_coordinates_file, delimiter=',', skip_header=1, usecols=(0, 1, 2), dtype=float)

        # downscaling inputs -- derived from user inputs
        self.mask_raster = utils.raster_to_array(self.mask_raster_file).flatten()

        self.urb_pop_init_year = os.path.join(self.datadir_histdata, "{}_urban_{}_1km.tif".format(self.region_code, self.historic_base_year))
        self.rur_pop_init_year = os.path.join(self.datadir_histdata, "{}_rural_{}_1km.tif".format(self.region_code, self.historic_base_year))

        # if the user wants to pass in the projections by file then use it, if not get params from user
        if (self.urban_pop_proj_n is None) and (self.rural_pop_proj_n is None):
            self.ssp_proj_file = os.path.join(self.datadir_future, "{}_{}_popproj.csv".format(self.region_code, self.ssp_code))
            self.ssp_data = pd.read_csv(self.ssp_proj_file)

        else:
            self.ssp_proj_file = None

        self.params_file = self.calibration_parameters_file

    def unpack_calibration_parameters(self):
        """Extract calibration parameters from file.

        :retrun:            [0]  float, Alpha Urban
                            [1]  float, Beta Urban
                            [2]  float, Alpha Rural
                            [3]  float, Beta Rural

        """

        df = pd.read_csv(self.calibration_parameters_file)

        urban_arr = df.loc[df['SSP'] == self.ssp_code, ["Alpha_Urban", "Beta_Urban"]].values
        rural_arr = df.loc[df['SSP'] == self.ssp_code, ["Alpha_Rural", "Beta_Rural"]].values

        alpha_urban = urban_arr[0][0]
        beta_urban = urban_arr[0][1]

        alpha_rural = rural_arr[0][0]
        beta_rural = rural_arr[0][1]

        return alpha_urban, beta_urban, alpha_rural, beta_rural

    @staticmethod
    def get_yaml(config_file):
        """Read the YAML config file

        :param config_file:         Full path with file name and extension to the input config.yml file

        :return:                    YAML config object

        """
        with open(config_file, 'r') as ymlfile:
            return yaml.load(ymlfile)
