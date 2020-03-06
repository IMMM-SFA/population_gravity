from population_gravity import Model

run = Model(
    datadir_histdata='/Users/d3y010/repos/github/spatial_population_downscaling_model/population_gravity/tests/data/inputs/main_inputs',
    ssp_data_directory='/Users/d3y010/repos/github/spatial_population_downscaling_model/population_gravity/tests/data/inputs/projection',
    ssp_code='SSP2',
    region_code='district_of_columbia',
    output_directory='/Users/d3y010/Desktop/outputs',
    historic_base_year=2000,
    future_start_year=2010,
    future_end_year=2020,
    time_step=10,
    alpha_urban=-0.606381394,
    alpha_rural=0,
    beta_urban=1.999999534,
    beta_rural=0)

# config_file = '/Users/d3y010/repos/github/spatial_population_downscaling_model/example/config.yml'
# run = Model(config_file)

# run.calibrate()

run.downscale()
