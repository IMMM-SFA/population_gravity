from population_gravity import Model

run = Model(
    datadir_histdata='/Users/d3y010/projects/population/rhode_island/inputs',
    ssp_data_directory='/Users/d3y010/projects/population/rhode_island/inputs/projection',
    ssp_code='SSP2',
    region_code='rhode_island',
    output_directory='/Users/d3y010/Desktop/outputs',
    historic_base_year=2000,
    future_start_year=2010,
    future_end_year=2020,
    time_step=10,
    alpha_urban=-2,
    alpha_rural=-0.341883268,
    beta_urban=0.461313537,
    beta_rural=0.996417737)

run.downscale()
