import pkg_resources

from SALib.sample import saltelli
from SALib.sample import latin
from SALib.analyze import sobol
from SALib.analyze import delta
from SALib.test_functions import Ishigami
import numpy as np

from population_gravity import Model

from population_gravity.downscale_utilities import raster_to_array

STATE_NAME = 'vermont'
SCENARIO = 'SSP2'
GRID_COORD_FILE = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_coordinates.csv'.format(STATE_NAME))
HIST_RURAL_RASTER = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_rural_2010_1km.tif'.format(STATE_NAME))
HIST_URBAN_RASTER = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_urban_2010_1km.tif'.format(STATE_NAME))
PROJ_POP_FILE = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_{}_popproj.csv'.format(STATE_NAME, SCENARIO))
ONE_D_IND_FILE = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_within_indices.txt'.format(STATE_NAME))
HIST_SUITABILITY = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_mask_short_term.tif'.format(STATE_NAME))
OUTPUT_DIRECTORY = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs')

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

param_values = latin.sample(problem, 4)

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
                projection_end_year=2030,
                time_step=10,
                raster_to_csv=False,
                run_number=index)

    run.downscale()

# aggregate results to numpy arrays
# Y = ...

# Print the first-order sensitivity indices
# Si = delta.analyze(problem, param_values, Y, print_to_console=True)

