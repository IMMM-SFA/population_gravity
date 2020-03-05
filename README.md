# spatial_population_downscaling_model

### A model to allocate urban and rural populations for a defined region to a grid

## Overview
The Spatial Population Downscaling Model allocates aggregate urban and rural populations for a defined region to grid cells within that region. In the IM3 application, the model is applied to US state-level urban and rural populations, which are allocated to 1 km grid cells within states. Allocation to grid cells is based on the relative suitability of each cell. Suitability is calculated using a gravity-based approach, in which the suitability of a given cell is determined by the population of surrounding cells (100 km radius), their distance away, and two parameters, namely alpha and beta. Alpha and beta are estimated based on historical population data and indicate the importance of returns to scale and distance in determining suitability values of cells, respectively. The model is composed of two components: calibration and projection. The calibration component uses historical urban/rural population grids of each state in 2000 and 2010 and an optimization algorithm to estimate the alpha and beta parameters that minimize the absolute difference between the actual population grid in 2010 and the one derived from the model. The two parameters can be modified to reflect distinctive forms of population development that may be desired in different socio-economic scenarios. Once the parameters are defined, the projection component downscales state-level urban/rural population aggregates of each state from 2020 to 2100 under different scenarios to grid cells within the state.

## Getting Started

Fill in once we get the package built

### Setting up a run

Fill in with DOI link once we get it setup

## Examples

### Run population downscaling for all years by passing argument values
```python
from spatial_downscale import Model

run = Model(
    datadir_histdata='<full path to historical input file directory>',
    ssp_data_directory='<full path to projected input directory>',
    ssp_code='SSP2', # shared socioeconomic pathway abbreviation
    region_code='district_of_columbia',  # US State name
    output_directory='<full path to the output directory>',
    start_year=2000,
    end_year=2020,
    time_step=10,
    alpha_urban=-0.606381394, # calibrated alpha parameter for urban
    alpha_rural=0, # calibrated alpha parameter for rural
    beta_urban=1.999999534, # calibrated beta parameter for urban
    beta_rural=0) # calibrated beta parameter for rural

# run downscaling for all years
run.downscale()
```

### Modify values per time step
```python
from spatial_downscale import Model

run = Model(
    datadir_histdata='<full path to historical input file directory>',
    ssp_data_directory='<full path to projected input directory>',
    ssp_code='SSP2', # shared socioeconomic pathway abbreviation
    region_code='district_of_columbia',  # US State name
    output_directory='<full path to the output directory>',
    start_year=2000,
    end_year=2020,
    time_step=10,
    alpha_urban=-0.606381394, # calibrated alpha parameter for urban
    alpha_rural=0, # calibrated alpha parameter for rural
    beta_urban=1.999999534, # calibrated beta parameter for urban
    beta_rural=0) # calibrated beta parameter for rural

# initialize model
run.initialize()

# downscale year 0
run.advance_step()

# modify the calibrated alpha parameter value for urban
run.cfg.alpha_urban = -0.1

# run next step with modified parameters
run.advance_step()

# close out run
run.close()

```

## Contact
For questions please contact:

Hamidreza Zoraghein @ <Hamidreza.Zoraghein@du.edu>

Brian O'Neill @ <Brian.ONeill@du.edu>
