from population_gravity import Model


# state_name = 'rhode_island'
# scenario = 'SSP2'
#
# run = Model(grid_coordinates_file='/Users/d3y010/repos/github/spatial_population_downscaling_model/population_gravity/tests/data/inputs/{}_coordinates.csv'.format(state_name),
#             historical_rural_pop_raster='/Users/d3y010/repos/github/spatial_population_downscaling_model/population_gravity/tests/data/inputs/{}_rural_2010_1km.tif'.format(state_name),
#             historical_urban_pop_raster='/Users/d3y010/repos/github/spatial_population_downscaling_model/population_gravity/tests/data/inputs/{}_urban_2010_1km.tif'.format(state_name),
#             historical_suitability_raster='/Users/d3y010/repos/github/spatial_population_downscaling_model/population_gravity/tests/data/inputs/{}_mask_short_term.tif'.format(state_name),
#             projected_population_file='/Users/d3y010/repos/github/spatial_population_downscaling_model/population_gravity/tests/data/inputs/{}_{}_popproj.csv'.format(state_name, scenario),
#             one_dimension_indices_file='/Users/d3y010/repos/github/spatial_population_downscaling_model/population_gravity/tests/data/inputs/{}_within_indices.txt'.format(state_name),
#             output_directory='/Users/d3y010/repos/github/spatial_population_downscaling_model/population_gravity/tests/data/outputs',
#             alpha_urban=-2.0,
#             alpha_rural=-0.34,
#             beta_urban=0.46,
#             beta_rural=1.0,
#             scenario=scenario,
#             state_name=state_name,
#             historic_base_year=2010,
#             projection_start_year=2020,
#             projection_end_year=2030,
#             time_step=10)

run = Model(config_file='/Users/d3y010/repos/github/spatial_population_downscaling_model/example/config.yml')

run.downscale()
