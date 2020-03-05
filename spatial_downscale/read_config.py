import os
import yaml


class ReadConfig:

    def __init__(self, config_file=None, datadir_histdata=None, ssp_data_directory=None, ssp_code=None,
                 region_code=None, output_directory=None, calibration_parameters_file=None,
                 compute_params=False, compute_proj=True, start_year=None, end_year=None,
                 time_step=None):

        if config_file is None:

            self.datadir_histdata = datadir_histdata
            self.datadir_future = ssp_data_directory
            self.ssp_code = ssp_code
            self.region_code = region_code
            self.datadir_output = output_directory
            self.calibration_parameters_file = calibration_parameters_file
            self.compute_params = compute_params
            self.compute_proj = compute_proj
            self.start_year = start_year
            self.end_year = end_year
            self.time_step = time_step

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

            # Boolean compute parameters
            self.compute_params = cfg['compute_params']

            # Boolean compute projections
            self.compute_proj = cfg['compute_proj']

            # output directory
            self.datadir_output = cfg['output_directory']

            # start year
            self.start_year = cfg['start_year']

            # end year
            self.end_year = cfg['end_year']

            # interval
            self.time_step = cfg['time_step']

            # calibration parameters file
            self.calibration_parameters_file = cfg['calibration_parameters_file']

        # calibration inputs -- derived from user inputs
        self.urb_pop_fst_year = os.path.join(self.datadir_histdata, "{}_urban_{}_1km.tif".format(self.region_code, self.start_year))
        self.urb_pop_snd_year = os.path.join(self.datadir_histdata, "{}_urban_{}_1km.tif".format(self.region_code, self.start_year + self.time_step))
        self.rur_pop_fst_year = os.path.join(self.datadir_histdata, "{}_rural_{}_1km.tif".format(self.region_code, self.start_year))
        self.rur_pop_snd_year = os.path.join(self.datadir_histdata, "{}_rural_{}_1km.tif".format(self.region_code, self.start_year + self.time_step))

        self.mask_raster = os.path.join(self.datadir_histdata, "{}_mask_short_term.tif".format(self.region_code))

        # TODO: why are there two point_indicies files with different names - same variable name
        self.point_indices = os.path.join(self.datadir_histdata, "{}_within_indices.txt".format(self.region_code))
        self.point_coors = os.path.join(self.datadir_histdata, "{}_coordinates.csv".format(self.region_code))

        # downscaling inputs -- derived from user inputs
        self.urb_pop_init_year = os.path.join(self.datadir_histdata, "{}_urban_{}_1km.tif".format(self.region_code, self.start_year))
        self.rur_pop_init_year = os.path.join(self.datadir_histdata, "{}_rural_{}_1km.tif".format(self.region_code, self.start_year))

        self.ssp_dataFn = os.path.join(self.datadir_future, "{}_{}_popproj.csv".format(self.region_code, self.ssp_code))

        self.params_file = self.calibration_parameters_file

    @staticmethod
    def get_yaml(config_file):
        """Read the YAML config file

        :param config_file:         Full path with file name and extension to the input config.yml file

        :return:                    YAML config object

        """
        with open(config_file, 'r') as ymlfile:
            return yaml.load(ymlfile)
