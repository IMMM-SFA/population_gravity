from population_gravity.sensitivity import Lhs
from population_gravity.sensitivity import BatchModelRun
from population_gravity.sensitivity import DeltaMomentIndependent

import os
import pkg_resources

STATE_NAME = 'vermont'
SCENARIO = 'SSP2'

GRID_COORD_FILE = pkg_resources.resource_filename('population_gravity', f'tests/data/inputs/{STATE_NAME}_coordinates.csv')
HIST_RURAL_RASTER = pkg_resources.resource_filename('population_gravity', f'tests/data/inputs/{STATE_NAME}_rural_2010_1km.tif')
HIST_URBAN_RASTER = pkg_resources.resource_filename('population_gravity', f'tests/data/inputs/{STATE_NAME}_urban_2010_1km.tif')
PROJ_POP_FILE = pkg_resources.resource_filename('population_gravity', f'tests/data/inputs/{STATE_NAME}_{SCENARIO}_popproj.csv')
ONE_D_IND_FILE = pkg_resources.resource_filename('population_gravity', f'tests/data/inputs/{STATE_NAME}_within_indices.txt')
HIST_SUITABILITY = pkg_resources.resource_filename('population_gravity', f'tests/data/inputs/{STATE_NAME}_mask_short_term.tif')
OUTPUT_DIRECTORY = pkg_resources.resource_filename('population_gravity', f'tests/data/outputs')


# generate latin hypercube sample
lhs = Lhs(alpha_urban_bounds=[-2.0, 2.0],
          alpha_rural_bounds=[-2.0, 2.0],
          beta_urban_bounds=[-2.0, 2.0],
          beta_rural_bounds=[-2.0, 2.0],
          kernel_distance_meters_bounds=[90000, 100000],
          n_samples=10,
          problem_dict_outfile=os.path.join(OUTPUT_DIRECTORY, 'problem_dict.p'),
          sample_outfile=os.path.join(OUTPUT_DIRECTORY, 'sample.npy'))

# create batch model run
x = BatchModelRun(grid_coordinates_file=GRID_COORD_FILE,
                  historical_rural_pop_raster=HIST_RURAL_RASTER,
                  historical_urban_pop_raster=HIST_URBAN_RASTER,
                  historical_suitability_raster=HIST_SUITABILITY,
                  projected_population_file=PROJ_POP_FILE,
                  one_dimension_indices_file=ONE_D_IND_FILE,
                  output_directory=OUTPUT_DIRECTORY,
                  alpha_urban=1.99999999995073, # these are default values that will be overridden by samples
                  alpha_rural=0.0750326293181678,
                  beta_urban=1.77529986067379,
                  beta_rural=1.42410799449511,
                  kernel_distance_meters=100000,
                  scenario=SCENARIO, # shared socioeconomic pathway abbreviation
                  state_name=STATE_NAME,
                  historic_base_year=2010,
                  projection_start_year=2020,
                  projection_end_year=2020,
                  time_step=1,
                  write_raster=False,
                  write_csv=False,
                  compress_csv=True,
                  write_array1d=True,
                  write_logfile=False,
                  output_total=False,
                  sample=lhs.sample,
                  problem_dict=lhs.problem_dict)

# generate runs using LHS parameter values
x.run_batch()

# instantiate delta model
delta_run = DeltaMomentIndependent(problem_dict=lhs.problem_dict,
                                   file_directory=OUTPUT_DIRECTORY,
                                   state_name=x.state_name,
                                   sample=lhs.sample,
                                   setting='Urban',  # either 'Urban' or 'Rural'
                                   file_extension='.npy',  # file extension matching the output format from run files
                                   output_file=os.path.join(OUTPUT_DIRECTORY, 'output.csv'))

# run analysis
output_list = delta_run.run_analysis()
