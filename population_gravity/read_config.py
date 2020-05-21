"""
Configuration reader for the population_gravity model

@author   Chris R. Vernon
@email:   chris.vernon@pnnl.gov

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

"""

import pkg_resources
import simplejson

import yaml

import numpy as np
import pandas as pd

import population_gravity.downscale_utilities as utils


class ReadConfig:
    """Read configuration data either provided in the configuration YAML file or as passed in via arguments.

    :param config_file:                         string. Full path to configuration YAML file with file name and
                                                extension. If not provided by the user, the code will default to the
                                                expectation of alternate arguments.

    :param grid_coordinates_file:               string. Full path with file name and extension to the CSV file
                                                containing the coordinates for each 1 km grid cell within the target
                                                state. File includes a header with the fields XCoord, YCoord, FID.
                                                Where data types and field descriptions are as follows:
                                                (XCoord, float, X coordinate in meters),
                                                (YCoord, float, Y coordinate in meters),
                                                (FID, int, Unique feature id)

    :param historical_suitability_raster:       string. Full path with file name and extension to the suitability
                                                raster containing values from 0.0 to 1.0 for each 1 km grid cell
                                                representing suitability depending on topographic and land use and
                                                land cover characteristics within the target state.

    :param historical_rural_pop_raster:         string. Full path with file name and extension to a raster containing
                                                rural population counts for each 1 km grid cell for the historical
                                                base time step.

    :param historical_urban_pop_raster:         string. Full path with file name and extension to a raster containing
                                                urban population counts for each 1 km grid cell for the historical
                                                base time step.

    :param projected_population_file:           string. Full path with file name and extension to a CSV file containing
                                                population projections per year separated into urban and rural
                                                categories.  Field descriptions for require fields as follows:
                                                (Year, integer, four digit year),
                                                (UrbanPop, float, population count for urban),
                                                (RuralPop, float, population count for rural),
                                                (Scenario, string, scenario as set in the `scenario` variable)

    :param one_dimension_indices_file:          string. Full path with file name and extension to the text file
                                                containing a file structured as a Python list (e.g. [0, 1]) that
                                                contains the index of each grid cell when flattened from a 2D array to
                                                a 1D array for the target state.

    :param output_directory:                    string. Full path with file name and extension to the output directory
                                                where outputs and the log file will be written.

    :param alpha_urban:                         float. Alpha parameter for urban. Represents the degree to which the
                                                population size of surrounding cells translates into the suitability
                                                of a focal cell.  A positive value indicates that the larger the
                                                population that is located within the 100 km neighborhood, the more
                                                suitable the focal cell is.  More negative value implies less suitable.
                                                Acceptable range:  -2.0 to 2.0

    :param beta_urban:                          float. Beta parameter for urban. Reflects the significance of distance
                                                to surrounding cells on the suitability of a focal cell.  Within 100 km,
                                                beta determines how distance modifies the effect on suitability.
                                                Acceptable range:  -0.5 to 2.0

    :param alpha_rural:                         float. Alpha parameter for rural. Represents the degree to which the
                                                population size of surrounding cells translates into the suitability
                                                of a focal cell.  A positive value indicates that the larger the
                                                population that is located within the 100 km neighborhood, the more
                                                suitable the focal cell is.  More negative value implies less suitable.
                                                Acceptable range:  -2.0 to 2.0

    :param beta_rural:                          float. Beta parameter for rural. Reflects the significance of distance
                                                to surrounding cells on the suitability of a focal cell.  Within 100 km,
                                                beta determines how distance modifies the effect on suitability.
                                                Acceptable range:  -0.5 to 2.0

    :param scenario:                            string. String representing the scenario with no spaces. Must match
                                                what is in the `projected_population_file` if passing population
                                                projections in using a file.

    :param state_name:                          string. Target state name with no spaces separated by an underscore.

    :param historic_base_year:                  int. Four digit historic base year.

    :param projection_start_year:               int. Four digit first year to process for the projection.

    :param projection_end_year:                 int. Four digit last year to process for the projection.

    :param time_step:                           int. Number of steps (e.g. number of years between projections)

    :param rural_pop_proj_n:                    float.  Rural population projection count for the projected year being
                                                calculated.  These can be read from the `projected_population_file`
                                                instead.

    :param urban_pop_proj_n:                    float.  Urban population projection count for the projected year being
                                                calculated.  These can be read from the `projected_population_file`
                                                instead.

    :param calibration_urban_year_one_raster:   string. Only used for running calibration.  Full path with file name and
                                                extension to a raster containing urban population counts for each 1 km
                                                grid cell for year one of the calibration.

    :param calibration_urban_year_two_raster:   string. Only used for running calibration.  Full path with file name and
                                                extension to a raster containing urban population counts for each 1 km
                                                grid cell for year two of the calibration.

    :param calibration_rural_year_one_raster:   string. Only used for running calibration.  Full path with file name and
                                                extension to a raster containing rural population counts for each 1 km
                                                grid cell for year one of the calibration.

    :param calibration_rural_year_two_raster:   string. Only used for running calibration.  Full path with file name and
                                                extension to a raster containing rural population counts for each 1 km
                                                grid cell for year two of the  calibration.

    :param kernel_distance_meters:              float. Distance kernel in meters; default 100,000 meters.

    :param raster_to_csv:                       boolean. Optionally export raster as a CSV file without nodata values

    """

    def __init__(self, config_file=None, grid_coordinates_file=None, historical_suitability_raster=None,
                 historical_rural_pop_raster=None, historical_urban_pop_raster=None, projected_population_file=None,
                 one_dimension_indices_file=None, output_directory=None, alpha_urban=None, beta_urban=None,
                 alpha_rural=None, beta_rural=None, scenario=None, state_name=None, historic_base_year=None,
                 projection_start_year=None,  projection_end_year=None, time_step=None, rural_pop_proj_n=None,
                 urban_pop_proj_n=None, calibration_urban_year_one_raster=None, calibration_urban_year_two_raster=None,
                 calibration_rural_year_one_raster=None, calibration_rural_year_two_raster=None,
                 kernel_distance_meters=None, raster_to_csv=False):

        if config_file is None:

            self.grid_coordinates_file = grid_coordinates_file
            self.historical_suitability_raster = historical_suitability_raster
            self.historical_rural_pop_raster = historical_rural_pop_raster
            self.historical_urban_pop_raster = historical_urban_pop_raster
            self.projected_population_file = projected_population_file
            self.one_dimension_indices_file = one_dimension_indices_file
            self.output_directory = output_directory
            self.alpha_urban = alpha_urban
            self.beta_urban = beta_urban
            self.alpha_rural = alpha_rural
            self.beta_rural = beta_rural
            self.scenario = scenario
            self.state_name = state_name.lower()
            self.historic_base_year = historic_base_year
            self.projection_start_year = projection_start_year
            self.projection_end_year = projection_end_year
            self.time_step = time_step
            self.rural_pop_proj_n = rural_pop_proj_n
            self.urban_pop_proj_n = urban_pop_proj_n
            self.kernel_distance_meters = kernel_distance_meters
            self.raster_to_csv = raster_to_csv

            # specific to calibration run
            self.calibration_urban_year_one_raster = calibration_urban_year_one_raster
            self.calibration_urban_year_two_raster = calibration_urban_year_two_raster
            self.calibration_rural_year_one_raster = calibration_rural_year_one_raster
            self.calibration_rural_year_two_raster = calibration_rural_year_two_raster

        else:

            # extract config file to YAML object
            cfg = self.get_yaml(config_file)

            self.grid_coordinates_file = self.validate_key(cfg, 'grid_coordinates_file')
            self.historical_suitability_raster = self.validate_key(cfg, 'historical_suitability_raster')
            self.historical_rural_pop_raster = self.validate_key(cfg, 'historical_rural_pop_raster')
            self.historical_urban_pop_raster = self.validate_key(cfg, 'historical_urban_pop_raster')
            self.projected_population_file = self.validate_key(cfg, 'projected_population_file')
            self.one_dimension_indices_file = self.validate_key(cfg, 'one_dimension_indices_file')
            self.output_directory = self.validate_key(cfg, 'output_directory')
            self.alpha_urban = self.validate_key(cfg, 'alpha_urban')
            self.beta_urban = self.validate_key(cfg, 'beta_urban')
            self.alpha_rural = self.validate_key(cfg, 'alpha_rural')
            self.beta_rural = self.validate_key(cfg, 'beta_rural')
            self.scenario = self.validate_key(cfg, 'scenario')
            self.state_name = self.validate_key(cfg, 'state_name').lower()
            self.historic_base_year = self.validate_key(cfg, 'historic_base_year')
            self.projection_start_year = self.validate_key(cfg, 'projection_start_year')
            self.projection_end_year = self.validate_key(cfg, 'projection_end_year')
            self.time_step = self.validate_key(cfg, 'time_step')
            self.rural_pop_proj_n = self.validate_key(cfg, 'rural_pop_proj_n')
            self.urban_pop_proj_n = self.validate_key(cfg, 'urban_pop_proj_n')
            self.kernel_distance_meters = self.validate_key(cfg, 'kernel_distance_meters')
            self.raster_to_csv = self.validate_key(cfg, 'raster_to_csv')

            # specific to calibration run
            self.calibration_urban_year_one_raster = self.validate_key(cfg, 'calibration_urban_year_one_raster')
            self.calibration_urban_year_two_raster = self.validate_key(cfg, 'calibration_urban_year_two_raster')
            self.calibration_rural_year_one_raster = self.validate_key(cfg, 'calibration_rural_year_one_raster')
            self.calibration_rural_year_two_raster = self.validate_key(cfg, 'calibration_rural_year_two_raster')

        # list of time steps in projection
        self.steps = range(self.projection_start_year, self.projection_end_year + self.time_step, self.time_step)

        # read in historical suitability mask as an array
        historical_suitability_2darray = utils.raster_to_array(self.historical_suitability_raster)
        self.historical_suitability_array = historical_suitability_2darray.flatten()

        # build data frame in the shape of the raster array
        self.df_indicies = utils.all_index_retriever(historical_suitability_2darray, ["row", "column"])

        # read grid indices of points that fall within the state boundary
        with open(self.one_dimension_indices_file, 'r') as r:
            self.one_dimension_indices = simplejson.load(r)

        # read in grid coordinates and feature ids
        self.grid_coordinates_array = np.genfromtxt(self.grid_coordinates_file, delimiter=',', skip_header=1, usecols=(0, 1, 2), dtype=float)

        # if the user wants to pass in the projections by file then use it, if not get params from user
        if (self.urban_pop_proj_n is None) and (self.rural_pop_proj_n is None):
            self.df_projected = pd.read_csv(self.projected_population_file)

        else:
            self.projected_population_file = None

        # get all neighboring states including the target state as a list
        self.neighbors = self.get_state_neighbors(state_name)

        # get a copy of the raster metadata from a states input raster
        self.template_raster_object, self.metadata = utils.get_raster_with_metadata(self.historical_suitability_raster)

        # get a bounding box from the historical raster
        self.bbox = utils.create_bbox(self.template_raster_object)

    @staticmethod
    def get_state_neighbors(state_name):
        """Get all neighboring states and the target state from lookup file as a list"""

        df = pd.read_csv(pkg_resources.resource_filename('population_gravity', 'data/neighboring_states_100km.csv'))

        # get the actual state name from the near states because they are not lower case like what is being passed
        state_find = df.loc[(df['target_state'] == state_name) & (df['near_state'].str.lower() == state_name)].values[0][-1]

        # extract a list of all neighboring states including the target state
        state_list = df.loc[df['target_state'] == state_name]['near_state'].to_list()

        # ensure that the target state comes first to prevent any issue with the reverse painter's algorithm for merge
        state_list.insert(0, state_list.pop(state_list.index(state_find)))

        # make all lower case
        return [i.lower() for i in state_list]

    @staticmethod
    def validate_key(yaml_object, key):
        """Check to see if key is in YAML file, if not return None"""
        try:
            return yaml_object[key]
        except KeyError:
            return None

    @staticmethod
    def get_yaml(config_file):
        """Read the YAML config file

        :param config_file:         Full path with file name and extension to the input config.yml file

        :return:                    YAML config object

        """
        with open(config_file, 'r') as yml:
            return yaml.load(yml)
