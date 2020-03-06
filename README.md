# population_gravity

### A model to allocate urban and rural populations for a defined region to a grid

## Overview
The Spatial Population Downscaling Model allocates aggregate urban and rural populations for a defined region to grid cells within that region. In the IM3 application, the model is applied to US state-level urban and rural populations, which are allocated to 1 km grid cells within states. Allocation to grid cells is based on the relative suitability of each cell. Suitability is calculated using a gravity-based approach, in which the suitability of a given cell is determined by the population of surrounding cells (100 km radius), their distance away, and two parameters, namely alpha and beta. Alpha and beta are estimated based on historical population data and indicate the importance of returns to scale and distance in determining suitability values of cells, respectively. The model is composed of two components: calibration and projection. The calibration component uses historical urban/rural population grids of each state in 2000 and 2010 and an optimization algorithm to estimate the alpha and beta parameters that minimize the absolute difference between the actual population grid in 2010 and the one derived from the model. The two parameters can be modified to reflect distinctive forms of population development that may be desired in different socio-economic scenarios. Once the parameters are defined, the projection component downscales state-level urban/rural population aggregates of each state from 2020 to 2100 under different scenarios to grid cells within the state.

## Getting Started
The `population_gravity` package uses only Python 3.3 and up.

### Step 1:
You can install `population_gravity` by running the following from your cloned directory (NOTE: ensure that you are using the desired `python` instance):

`pip install git+https://github.com/IMMM-SFA/spatial_population_downscaling_model.git --user`

### Step 2:
Confirm that the module and its dependencies have been installed by running from your prompt:

```python
from population_gravity import Model
```

If no error is returned then you are ready to go!

### Setting up a run

The following are required input files to start a model run for a target state:
| Input File                   | Description                                                                                                                                                                                                                                                                                                       | Source             |
|------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------|
| <state_name>_coordinates.csv | Coordinates for each 1 km grid cell within the target state.  File includes a header with the fields `XCoord`, `YCoord`, `FID`. Where data types and field descriptions are as follows: (`XCoord`, `float`, X coordinate in meters), (`YCoord`, `float`, Y coordinate in meters), (`FID`, `int`, Unique feature id) | Generated in a GIS |

### Key variables
Users can modify any key variables after model initialization.  This includes updating values between time steps.

| variable         | type  | description                                                 |
|------------------|-------|-------------------------------------------------------------|
| alpha_urban      | float | Alpha parameter for urban                                   |
| beta_urban       | float | Beta parameter for urban                                    |
| alpha_rural      | float | Alpha parameter for rural                                   |
| beta_rural       | float | Beta parameter for rural                                    |
| urban_pop_proj_n | float | Urban population projection number for the second time step |
| rural_pop_proj_n | float | Rural population projection number for the second time step |


## Examples

### Run population downscaling for all years by passing argument values
```python
from population_gravity import Model

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
from population_gravity import Model

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
