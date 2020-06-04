"""
Configuration reader for the population_gravity model

@author   Chris R. Vernon
@email:   chris.vernon@pnnl.gov

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

"""

import datetime
import os
import pkg_resources
import simplejson

import rasterio
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

    :param write_raster:                        boolean. Optionally export raster output; Default True

    :param write_csv:                           boolean. Optionally export raster as a CSV file without nodata values

    :param write_array2d:                       boolean. Optionally export a NumPy 2D array for each output in the shape
                                                of the template raster

    :param write_array1d:                       boolean. Optionally export a Numpy 1D flattened array of only grid cells
                                                within the target state

    :param run_number:                          int. Add on for the file name when running sensitivity analysis

    :param write_logfile:                       boolean.  Optionally write log to file.; Default True

    :param compress_csv:                        boolean.  Optionally compress CSV file to GZIP if outputting in CSV

    :param output_total:                        boolean.  Choice to output total (urban + rural) dataset; Defualt True

    """

    # format for datetime string
    DATETIME_FORMAT = '%Y-%m-%d_%Hh%Mm%Ss'

    # key names from YAML config file
    OUT_DIR_KEY = 'output_directory'
    START_STEP_KEY = 'start_step'
    THROUGH_STEP_KEY = 'through_step'
    TIME_STEP_KEY = 'time_step'
    ALPHA_KEY = 'alpha_param'
    BETA_KEY = 'beta_param'

    # definition of acceptable range of values for parameters
    MAX_PARAM_VALUE = 2.0
    MIN_PARAM_VALUE = -2.0

    def __init__(self, config_file=None, grid_coordinates_file=None, historical_suitability_raster=None,
                 historical_rural_pop_raster=None, historical_urban_pop_raster=None, projected_population_file=None,
                 one_dimension_indices_file=None, output_directory=None, alpha_urban=None, beta_urban=None,
                 alpha_rural=None, beta_rural=None, scenario=None, state_name=None, historic_base_year=None,
                 projection_start_year=None,  projection_end_year=None, time_step=None, rural_pop_proj_n=None,
                 urban_pop_proj_n=None, calibration_urban_year_one_raster=None, calibration_urban_year_two_raster=None,
                 calibration_rural_year_one_raster=None, calibration_rural_year_two_raster=None,
                 kernel_distance_meters=None, write_raster=True, write_csv=False, write_array1d=False,
                 write_array2d=False, run_number='', write_logfile=True, compress_csv=True, output_total=True):

        self._config_file = config_file
        self._alpha_urban = alpha_urban
        self._alpha_rural = alpha_rural
        self._beta_urban = beta_urban
        self._beta_rural = beta_rural
        self._output_directory = output_directory
        self._grid_coordinates_file = grid_coordinates_file
        self._historical_suitability_raster = historical_suitability_raster
        self._historical_rural_pop_raster = historical_rural_pop_raster
        self._historical_urban_pop_raster = historical_urban_pop_raster
        self._projected_population_file = projected_population_file
        self._one_dimension_indices_file = one_dimension_indices_file
        self._scenario = scenario
        self._state_name = state_name.lower()
        self._historic_base_year = historic_base_year
        self._projection_start_year = projection_start_year
        self._projection_end_year = projection_end_year
        self._time_step = time_step
        self._rural_pop_proj_n = rural_pop_proj_n
        self._urban_pop_proj_n = urban_pop_proj_n
        self._kernel_distance_meters = kernel_distance_meters
        self._write_raster = write_raster
        self._write_csv = write_csv
        self._write_array1d = write_array1d
        self._write_array2d = write_array2d
        self._run_number = run_number
        self._write_logfile = write_logfile
        self._compress_csv = compress_csv
        self._output_total = output_total

        # specific to calibration run
        self._calibration_urban_year_one_raster = calibration_urban_year_one_raster
        self._calibration_urban_year_two_raster = calibration_urban_year_two_raster
        self._calibration_rural_year_one_raster = calibration_rural_year_one_raster
        self._calibration_rural_year_two_raster = calibration_rural_year_two_raster

        # get a copy of the raster metadata from a states input raster
        self._template_raster_object, self._metadata = utils.get_raster_with_metadata(self.historical_suitability_raster)

    @property
    def output_total(self):
        """Choice to output total dataset (urban + rural)."""

        return self._output_total

    @property
    def run_number(self):
        """An integer add on for the file name when running sensitivity analysis."""

        return self._run_number

    @property
    def compress_csv(self):
        """Compress CSV to GZIP option."""

        return self._compress_csv

    @property
    def write_logfile(self):
        """Optionally write log outputs to a file."""

        return self._write_logfile

    @property
    def write_raster(self):
        """Optionally save outputs to a raster."""

        return self._write_raster

    @property
    def write_array1d(self):
        """Optionally save outputs to a 1D array for cells within the target state."""

        return self._write_array1d

    @property
    def write_array2d(self):
        """Optionally save outputs to a 1D array for cells within the target state."""

        return self._write_array2d

    @property
    def write_csv(self):
        """Optionally export raster as a CSV file without nodata values;  option set to compress CSV using gzip.
        Exports values for non-NODATA grid cells as field name `value`.

        """

        return self._write_csv

    @property
    def historical_urban_pop_raster(self):
        """Full path with file name and extension to a raster containing urban population counts for each 1 km
        grid cell for the historical base time step.

        """

        return self.validate_file(self._historical_urban_pop_raster)

    @property
    def historical_rural_pop_raster(self):
        """Full path with file name and extension to a raster containing rural population counts for each 1 km
        grid cell for the historical base time step.

        """

        return self.validate_file(self._historical_rural_pop_raster)

    @property
    def kernel_distance_meters(self):
        """Distance kernel in meters; default 100,000 meters."""

        return self.validate_float(self._kernel_distance_meters)

    @kernel_distance_meters.setter
    def kernel_distance_meters(self, value):
        """Setter for kernel_distance_meters."""

        self._kernel_distance_meters = value

    @property
    def bbox(self):
        """Get a bounding box from the historical raster."""

        return utils.create_bbox(self.template_raster_object)

    @property
    def neighbors(self):
        """Get all neighboring states including the target state as a list."""

        return self.get_state_neighbors(self.state_name)

    @property
    def metadata(self):
        """Get a copy of the raster metadata from a states input raster."""

        return self._metadata

    @property
    def template_raster_object(self):
        """Get a copy of the raster metadata from a states input raster."""

        return self._template_raster_object

    @property
    def df_projected(self):
        """From population projection file if exists."""
        if (self.urban_pop_proj_n is None) and (self.rural_pop_proj_n is None):
            return pd.read_csv(self.projected_population_file)
        else:
            return None

    @property
    def date_time_string(self):
        """Get a current time in a string matching the specified datetime format."""

        return datetime.datetime.now().strftime(self.DATETIME_FORMAT)

    @property
    def datetime_format(self):
        """Convenience wrapper for the DATETIME_FORMAT class attribute."""

        return self.DATETIME_FORMAT

    @property
    def output_directory(self):
        """Validate output directory."""

        if self.config is None:
            return self.validate_directory(self._output_directory)
        else:
            key = self.validate_key(self.config, self.OUT_DIR_KEY)
            return self.validate_directory(key)

    @property
    def scenario(self):
        """Target scenario name."""

        return self._scenario

    @property
    def state_name(self):
        """Target state name."""

        return self._state_name

    @property
    def logfile(self):
        """Full path with file name and extension to the logfile."""

        return os.path.join(self.output_directory, 'logfile_{}_{}_{}.log'.format(self.scenario,
                                                                                 self.state_name,
                                                                                 self.date_time_string))
    @property
    def time_step(self):
        """Number of time steps."""

        return self.validate_step(self._time_step, self.TIME_STEP_KEY)

    @property
    def projection_start_year(self):
        """Four digit first year to process for the projection."""

        return self.validate_step(self._projection_start_year, 'projection_start_year')

    @property
    def projection_end_year(self):
        """Four digit last year to process for the projection."""

        return self.validate_step(self._projection_end_year, 'projection_end_year')

    @property
    def historic_base_year(self):
        """Four digit historic base year."""

        return self.validate_step(self._historic_base_year, 'historic_base_year')

    @property
    def steps(self):
        """Create a list of time steps from the start and through steps by the step interval."""

        return range(self.projection_start_year, self.projection_end_year + self.time_step, self.time_step)

    @property
    def historical_suitability_raster(self):
        """Full path with file name and extension to the suitability raster containing values from 0.0 to 1.0
        for each 1 km grid cell representing suitability depending on topographic and land use and land cover
        characteristics within the target state.

        """

        return self.validate_file(self._historical_suitability_raster)

    @property
    def historical_suitability_2darray(self):
        """Read in historical suitability mask as an array"""

        return utils.raster_to_array(self.historical_suitability_raster)

    @property
    def historical_suitability_array(self):
        """Flatten historical suitability array."""

        return self.historical_suitability_2darray.flatten()

    @property
    def df_indicies(self):
        """Build data frame in the shape of the raster array."""

        return utils.all_index_retriever(self.historical_suitability_2darray, ["row", "column"])

    @property
    def one_dimension_indices_file(self):
        """File that describe grid indices of points that fall within the state boundary."""

        return self.validate_file(self._one_dimension_indices_file)

    @property
    def one_dimension_indices(self):
        """Grid indices for the state to an array."""

        with open(self.one_dimension_indices_file, 'r') as r:
            return simplejson.load(r)

    @property
    def grid_coordinates_file(self):
        """File with grid coordinates and feature ids."""

        return self.validate_file(self._grid_coordinates_file)

    @property
    def grid_coordinates_array(self):
        """Grid coordinates to array."""

        return np.genfromtxt(self.grid_coordinates_file, delimiter=',', skip_header=1, usecols=(0, 1, 2), dtype=float)

    @property
    def urban_pop_proj_n(self):
        """Urban population projection count for the projected year being calculated.  These can be read from
        the `projected_population_file` instead.

        """

        return self.validate_float(self._urban_pop_proj_n)

    @property
    def rural_pop_proj_n(self):
        """Rural population projection count for the projected year being calculated.  These can be read from
        the `projected_population_file` instead.

        """

        return self.validate_float(self._rural_pop_proj_n)

    @property
    def projected_population_file(self):
        """Full path with file name and extension to a CSV file containing population projections per year
        separated into urban and rural categories.

        """

        return self.validate_file(self._projected_population_file)

    @property
    def config(self):
        """Read the YAML config file object"""

        if self._config_file is None:
            return None

        else:
            with open(self._config_file, 'r') as yml:
                return yaml.load(yml)

    @property
    def template_raster(self):
        """Generate template raster specifications.

        :return:                        [0] 2D array of template raster values
                                        [1] 1D flattened array
                                        [2] row count
                                        [3] column count
                                        [4] profile

        """

        with rasterio.open(self.historical_suitability_raster) as src_raster:
            profile = src_raster.profile
            array2d = src_raster.read(1)
            row_count = array2d.shape[0]
            col_count = array2d.shape[1]
            array1d = array2d.flatten()

        return array2d, array1d, row_count, col_count, profile

    @property
    def alpha_urban(self):
        """Alpha urban parameter for model."""

        return self.validate_parameter(self._alpha_urban, 'alpha_urban')

    @alpha_urban.setter
    def alpha_urban(self, value):
        """Setter for alpha urban parameter."""

        self._alpha_urban = self.validate_parameter(value, 'alpha_urban')

    @property
    def alpha_rural(self):
        """Alpha rural parameter for model."""

        return self.validate_parameter(self._alpha_rural, 'alpha_rural')

    @alpha_rural.setter
    def alpha_rural(self, value):
        """Setter for alpha rural parameter."""

        self._alpha_rural = self.validate_parameter(value, 'alpha_rural')

    @property
    def beta_urban(self):
        """Beta urban parameter for model."""

        return self.validate_parameter(self._beta_urban, 'beta_urban')

    @beta_urban.setter
    def beta_urban(self, value):
        """Setter for beta urban parameter."""

        self._beta_urban = self.validate_parameter(value, 'beta_urban')

    @property
    def beta_rural(self):
        """Beta rural parameter for model."""

        return self.validate_parameter(self._beta_rural, 'beta_rural')

    @beta_rural.setter
    def beta_rural(self, value):
        """Setter for beta rural parameter."""

        self._beta_rural = self.validate_parameter(value, 'beta_rural')

    def validate_parameter(self, param, key):
        """Validate parameter existence and range.

        :param param:               Parameter value
        :type param:                float

        :param key:                 Configuration key from YAML file
        :type key:                  str

        :return:                    int; parameter

        """

        if self.config is None:
            is_float = self.validate_float(param)
            return self.validate_range(is_float)
        else:
            is_key = self.validate_key(self.config, key)
            is_float = self.validate_float(is_key)
            return self.validate_range(is_float)

    def validate_range(self, value):
        """Ensure value falls within an acceptable range."""

        if (value >= self.MIN_PARAM_VALUE) and (value <= self.MAX_PARAM_VALUE):
            return value
        else:
            raise ValueError(f"Parameter value '{value}' is not within the valid range of {self.MIN_PARAM_VALUE} - {self.MAX_PARAM_VALUE}.")

    @staticmethod
    def validate_float(val):
        """Ensure parameter value is type float"""

        if val is not None:
            try:
                return float(val)
            except TypeError:
                raise TypeError(f"Parameter value '{val}' is not a float.")
        else:
            return None

    @staticmethod
    def validate_directory(directory):
        """Validate directory to ensure it exists.

        :param directory:                       Full path to the target directory.
        :type directory:                        str

        :return:                                Full path of a valid directory

        """
        if (directory is not None) and (os.path.isdir(directory)):
            return directory
        elif (directory is not None) and (os.path.isdir(directory is False)):
            raise NotADirectoryError(f"Directory: {directory} does not exist.")
        else:
            return None

    @staticmethod
    def validate_file(file):
        """Validate file to ensure it exists.

        :param file:                            Full path to the target file.
        :type file:                             str

        :return:                                Full path of a valid file

        """
        if (file is not None) and (os.path.isfile(file)):
            return file
        elif (file is not None) and (os.path.isfile(file) is False):
            raise FileNotFoundError(f"File: {file} does not exist.")
        else:
            return None

    def validate_step(self, step, key):
        """Validate step existence and value.

        :param step:                Time step value
        :type step:                 int

        :param key:                 Configuration key from YAML file
        :type key:                  str

        :return:                    int; time step

        """

        if self.config is None:
            return self.validate_int(step)
        else:
            is_key = self.validate_key(self.config, key)
            return self.validate_int(is_key)

    @staticmethod
    def validate_int(n):
        """Ensure time step is type int"""

        if n is not None:
            try:
                return int(n)
            except TypeError:
                raise TypeError(f"Value '{n}' is not an integer.")
        else:
            return None

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
        """Check to see if key is in YAML file, if not return None.

        :param yaml_object:                     YAML object for the configuration file

        :param key:                             Target key name from the configuration file.
        :type key:                              str

        :return:                                Value from configuration file matching the key. If no key present,
                                                return None.
        """
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
