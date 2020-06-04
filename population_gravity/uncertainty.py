import os
import pkg_resources

import numpy as np
from SALib.sample import latin
from SALib.analyze import delta

from population_gravity import Model


STATE_NAME = 'vermont'
SCENARIO = 'SSP2'
GRID_COORD_FILE = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_coordinates.csv'.format(STATE_NAME))
HIST_RURAL_RASTER = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_rural_2010_1km.tif'.format(STATE_NAME))
HIST_URBAN_RASTER = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_urban_2010_1km.tif'.format(STATE_NAME))
PROJ_POP_FILE = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_{}_popproj.csv'.format(STATE_NAME, SCENARIO))
ONE_D_IND_FILE = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_within_indices.txt'.format(STATE_NAME))
HIST_SUITABILITY = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_mask_short_term.tif'.format(STATE_NAME))
OUTPUT_DIRECTORY = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs')




d = {'alpha_urban': [-2.0, 2.0],
     'alpha_rural': [-2.0, 2.0],
     'beta_urban': [-2.0, 2.0],
     'beta_rural': [-2.0, 2.0],
     'kernel_distance_meters': [90000, 100000]}

class LhsMoment:

    # parameter value limits
    LOWER_ALPHA_LIMIT = -2.0
    UPPER_ALPHA_LIMIT = 2.0
    LOWER_BETA_LIMIT = -2.0
    UPPER_BETA_LIMIT = 2.0

    # kernel distance in meters
    LOWER_KD_LIMIT = 100
    UPPER_KD_LIMIT = 100000

    def __init__(self, alpha_urban_bounds=None, alpha_rural_bounds=None, beta_urban_bounds=None, beta_rural_bounds=None,
                 kernel_distance_meters_bounds=None):

        self._alpha_urban_bounds = alpha_urban_bounds
        self._alpha_rural_bounds = alpha_rural_bounds
        self._beta_urban_bounds = beta_urban_bounds
        self._beta_rural_bounds = beta_rural_bounds
        self._kernel_distance_meters_bounds = kernel_distance_meters_bounds


    @property
    def alpha_urban_bounds(self):
        """Ensure values are within acceptable bounds."""

        if self._alpha_urban_bounds is not None:

            # ensure lower is not >= upper
            valid_lower, valid_upper = self.validate_min_max(self._alpha_urban_bounds, 'alpha_urban')

            return self.validate_limits(valid_upper, valid_lower, self.UPPER_ALPHA_LIMIT, self.LOWER_ALPHA_LIMIT, 'alpha_urban')

        else:
            return self._alpha_urban_bounds

    @property
    def alpha_rural_bounds(self):
        """Ensure values are within acceptable bounds."""

        if self._alpha_rural_bounds is not None:

            # ensure lower is not >= upper
            valid_lower, valid_upper = self.validate_min_max(self._alpha_rural_bounds, 'alpha_rural')

            return self.validate_limits(valid_upper, valid_lower, self.UPPER_ALPHA_LIMIT, self.LOWER_ALPHA_LIMIT, 'alpha_rural')

        else:
            return self._alpha_rural_bounds

    @property
    def beta_urban_bounds(self):
        """Ensure values are within acceptable bounds."""

        if self._beta_urban_bounds is not None:

            # ensure lower is not >= upper
            valid_lower, valid_upper = self.validate_min_max(self._beta_urban_bounds, 'beta_urban')

            return self.validate_limits(valid_upper, valid_lower, self.UPPER_BETA_LIMIT, self.LOWER_BETA_LIMIT, 'beta_urban')

        else:
            return self._beta_urban_bounds

    @property
    def beta_rural_bounds(self):
        """Ensure values are within acceptable bounds."""

        if self._beta_rural_bounds is not None:

            # ensure lower is not >= upper
            valid_lower, valid_upper = self.validate_min_max(self._beta_rural_bounds, 'beta_rural')

            return self.validate_limits(valid_upper, valid_lower, self.UPPER_BETA_LIMIT, self.LOWER_BETA_LIMIT, 'beta_rural')

        else:
            return self._beta_rural_bounds

    @property
    def kernel_distance_meters_bounds(self):
        """Ensure values are withing acceptable bounds."""

        if self._kernel_distance_meters_bounds is not None:

            # ensure lower is not >= upper
            valid_lower, valid_upper = self.validate_min_max(self._beta_rural_bounds, 'kernel_distance_meters')

            return self.validate_limits(valid_upper, valid_lower, self.UPPER_KD_LIMIT, self.LOWER_KD_LIMIT, 'kernel_distance_meters')

        else:
            return self._kernel_distance_meters_bounds

    @staticmethod
    def validate_limits(upper_value, lower_value, upper_limit, lower_limit, variable):
        """Ensure the values for bounds are not outside of the acceptable limits.

        :param upper_value:                     upper value for a parameter
        :type:                                  int; float

        :param lower_value:                     lower value for a parameter
        :type:                                  int; float

        :param upper_limit:                     upper acceptable value for a parameter
        :type:                                  int; float

        :param lower_limit:                     lower acceptable value for a parameter
        :type:                                  int; float

        """

        if lower_value < lower_limit:
            raise ValueError(f"Lower limit '{lower_value}'' for '{variable}' cannot be less than '{lower_limit}'")

        elif upper_value > upper_limit:
            raise ValueError(f"Upper limit '{upper_value}'' for '{variable}' cannot be less than '{upper_limit}'")

        return lower_value, upper_value

    @staticmethod
    def validate_min_max(bounds, variable):
        """Check to make sure min is not >= max

        :param bounds:                      [min_value, max_value]
        :type bounds:                       list

        :param variable:                    variable name
        :type variable:                     str

        """

        min_bound, max_bound = bounds

        if min_bound >= max_bound:
            raise ValueError(f"Minimum bound '{min_bound}' for '{variable}' is >= '{max_bound}'")
        else:
            return bounds






# Define the model inputs
problem = {
    'num_vars': 5,
    'names': ['alpha_urban', 'alpha_rural', 'beta_urban', 'beta_rural', 'kernel_distance_meters'],
    'bounds': [[-2.0, 2.0],
               [-2.0, 2.0],
               [-2.0, 2.0],
               [-2.0, 2.0],
               [90000, 100000]]
}

# create a latin hypercube for the problem
param_values = latin.sample(problem, 10)

for index, i in enumerate(param_values):

    run = Model(grid_coordinates_file=GRID_COORD_FILE,
                historical_rural_pop_raster=HIST_RURAL_RASTER,
                historical_urban_pop_raster=HIST_URBAN_RASTER,
                historical_suitability_raster=HIST_SUITABILITY,
                projected_population_file=PROJ_POP_FILE,
                one_dimension_indices_file=ONE_D_IND_FILE,
                output_directory=OUTPUT_DIRECTORY,
                alpha_urban=i[0],
                alpha_rural=i[1],
                beta_urban=i[2],
                beta_rural=i[3],
                kernel_distance_meters=i[4],
                scenario=SCENARIO,
                state_name=STATE_NAME,
                historic_base_year=2010,
                projection_start_year=2020,
                projection_end_year=2020,
                time_step=10,
                write_raster=False,
                write_array1d=True,
                write_csv=False,
                write_logfile=False,
                run_number=index)

    run.downscale()

# get urban
outputs = [os.path.join(OUTPUT_DIRECTORY, i) for i in os.listdir(OUTPUT_DIRECTORY) if (i.split('_')[0] == 'vermont') and (os.path.splitext(i)[-1] == '.npy') and ('Urban' in i)]

# read arrays

# shape (n_runs, grid cells)
arr = np.zeros(shape=(len(outputs), 24905))

for index, i in enumerate(outputs):

    arr[index, :] = np.load(i)


with open('/Users/d3y010/Desktop/delta_result.csv', 'w') as out:

    out.write('param,delta,delta_conf,S1,S1_conf,gridcell\n')

    for i in range(24905):

        print(i)

        Y = arr[:, i]

        # Print the first-order sensitivity indices
        Si = delta.analyze(problem, param_values, Y, print_to_console=False)

        out.write(f"alpha_urban,{Si['delta'][0]},{Si['delta_conf'][0]},{Si['S1'][0]},{Si['S1_conf'][0]},{i}\n")
        out.write(f"alpha_rural,{Si['delta'][1]},{Si['delta_conf'][1]},{Si['S1'][1]},{Si['S1_conf'][1]},{i}\n")
        out.write(f"beta_urban,{Si['delta'][2]},{Si['delta_conf'][2]},{Si['S1'][2]},{Si['S1_conf'][2]},{i}\n")
        out.write(f"beta_rural,{Si['delta'][3]},{Si['delta_conf'][3]},{Si['S1'][3]},{Si['S1_conf'][3]},{i}\n")
        out.write(f"kernel_distance_meters,{Si['delta'][4]},{Si['delta_conf'][4]},{Si['S1'][4]},{Si['S1_conf'][4]},{i}\n")
