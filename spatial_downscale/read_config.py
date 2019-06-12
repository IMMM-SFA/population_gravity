import os
import yaml


class ReadConfig:

    def __init__(self, config_file):

        # extract config file to YAML object
        cfg = self.get_yaml(config_file)

        # root directory of the data
        self.data_rootdir = cfg['data_rootdir']

        # directory of main inputs
        self.datadir_histdata = cfg['datadir_histdata']

        # directory of projections
        self.datadir_future = cfg['datadir_future']

        # point indicies file
        self.point_indices = cfg['point_indices']

        # SSP
        self.ssp_code = cfg['ssp_code']

        # region code
        self.region_code = cfg['region_code']

        # Boolean compute parameters
        self.compute_params = cfg['compute_params']

        # Boolean compute projections
        self.compute_proj = cfg['compute_proj']

        # output directory
        self.datadir_output = cfg['datadir_output']

        # calibration inputs -- derived from user inputs
        self.urb_pop_fst_year = os.path.join(self.datadir_histdata, "{}_Urban_{}_1km.tif".format(self.region_code, cfg['first_yr']))
        self.urb_pop_snd_year = os.path.join(self.datadir_histdata, "{}_Urban_{}_1km.tif".format(self.region_code, cfg['second_yr']))
        self.rur_pop_fst_year = os.path.join(self.datadir_histdata, "{}_Rural_{}_1km.tif".format(self.region_code, cfg['first_yr']))
        self.rur_pop_snd_year = os.path.join(self.datadir_histdata, "{}_Rural_{}_1km.tif".format(self.region_code, cfg['second_yr']))

        self.mask_raster = os.path.join(self.datadir_histdata, "{}_Mask_short_term.tif".format(self.region_code))

        # TODO: why are there two point_indicies files with different names - same variable name
        self.point_indices = os.path.join(self.datadir_histdata, "{}_Within_Indices.txt".format(self.region_code))
        self.point_coors = os.path.join(self.datadir_histdata, "{}_Coors.csv".format(self.region_code))

        # downscaling inputs -- derived from user inputs
        self.urb_pop_init_year = os.path.join(self.datadir_histdata, "{}_Urban_{}_1km.tif".format(self.region_code, cfg['first_yr']))
        self.rur_pop_init_year = os.path.join(self.datadir_histdata, "{}_Rural_{}_1km.tif".format(self.region_code, cfg['first_yr']))

        self.ssp_dataFn = os.path.join(self.datadir_future, "{}_{}_popproj.csv".format(self.region_code, self.ssp_code))

        self.params_file = os.path.join(self.datadir_output, "{}_{}_Params.csv".format(self.region_code, self.ssp_code))

    @staticmethod
    def get_yaml(config_file):
        """Read the YAML config file

        :param config_file:         Full path with file name and extension to the input config.yml file

        :return:                    YAML config object

        """
        with open(config_file, 'r') as ymlfile:
            return yaml.load(ymlfile)
