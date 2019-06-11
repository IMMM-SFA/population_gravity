import yaml


class ReadConfig:

    def __init__(self, config_file):

        # extract config file to YAML object
        cfg = self.get_yaml(config_file)

        # root directory of the data
        self.data_rootdir = cfg['data_rootdir']

        # directory of main inputs
        self.datadir_histdata = cfg['data_rootdir']

        # directory of projections
        self.datadir_future = cfg['data_rootdir']

        # point indicies file
        self.point_indices = cfg['data_rootdir']

        # SSP
        self.ssp_code = cfg['data_rootdir']

        # region code
        self.region_code = cfg['data_rootdir']

        # Boolean compute parameters
        self.compute_params = cfg['data_rootdir']

        # Boolean compute projections
        self.compute_proj = cfg['data_rootdir']

        # output directory
        self.datadir_output = cfg['data_rootdir']

    @staticmethod
    def get_yaml(config_file):
        """Read the YAML config file

        :param config_file:         Full path with file name and extension to the input config.yml file

        :return:                    YAML config object

        """
        with open(config_file, 'r') as ymlfile:

            return yaml.load(ymlfile)
