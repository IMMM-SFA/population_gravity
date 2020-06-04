import os
import pkg_resources

from SALib.sample import latin
from SALib.analyze import delta
import numpy as np

from population_gravity import Model
import os

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


param_values = latin.sample(problem, 100)

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
