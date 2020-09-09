[![Build Status](https://travis-ci.org/IMMM-SFA/population_gravity.svg?branch=master)](https://travis-ci.org/IMMM-SFA/population_gravity) 
[![codecov](https://codecov.io/gh/IMMM-SFA/population_gravity/branch/master/graph/badge.svg)](https://codecov.io/gh/IMMM-SFA/population_gravity)

# population_gravity

### A model to downscale urban and rural population for a defined state to a 1 km grid

## Overview
The `population_gravity` model allocates aggregate urban and rural populations for a defined region to grid cells within that region. In the IM3 application, the model is applied to US state-level urban and rural populations, which are allocated to 1 km grid cells within states. Allocation to grid cells is based on the relative suitability of each cell. Suitability is calculated using a gravity-based approach, in which the suitability of a given cell is determined by the population of surrounding cells (100 km radius), their distance away, and two parameters, namely alpha and beta. Alpha and beta are estimated based on historical population data and indicate the importance of returns to scale and distance in determining suitability values of cells, respectively. The model is composed of two components: calibration and projection. The calibration component uses historical urban/rural population grids of each state in 2000 and 2010 and an optimization algorithm to estimate the alpha and beta parameters that minimize the absolute difference between the actual population grid in 2010 and the one derived from the model. The two parameters can be modified to reflect distinctive forms of population development that may be desired in different socio-economic scenarios. Once the parameters are defined, the projection component downscales state-level urban/rural population aggregates of each state from 2020 to 2100 under different scenarios to grid cells within the state.

## Associated publications
Zoraghein, H., & O’Neill, B. (2020a). The methodological foundation of a gravity-based model to downscale U.S. state-level populations to high-resolution distributions for integrated human-environment analysis. *Privileged*.

Zoraghein, H., & O’Neill, B. (2020b). U.S. state-level projections of the spatial distribution of population consistent with
Shared Socioeconomic Pathways. *Privileged*.


## Getting Started
The `population_gravity` package uses only **Python 3.6** and up.

### Step 1:
You can install `population_gravity` by running the following from your cloned directory (NOTE: ensure that you are using the desired `pip` instance that matches your Python3 distribution):

`pip3 install git+https://github.com/IMMM-SFA/spatial_population_downscaling_model.git --user`

### Step 2:
Confirm that the module and its dependencies have been installed by running from your prompt:

```python
from population_gravity import Model
```

If no error is returned then you are ready to go!

## Setting up a run

### Expected arguments
See examples below for how to pass into the `Model` class

| Argument | Type | Description |
|----|----|----|
| `config_file` | string | Full path to configuration YAML file with file name and extension. If not provided by the user, the code will default to the expectation of alternate arguments. |
| `grid_coordinates_file` | string | Full path with file name and extension to the CSV file containing the coordinates for each 1 km grid cell within the target state. File includes a header with the fields XCoord, YCoord, FID.,Where data types and field descriptions are as follows: (XCoord, float, X coordinate in meters),(YCoord, float, Y coordinate in meters),(FID, int, Unique feature id) |
| `historical_suitability_raster` | string | Full path with file name and extension to the suitability raster containing values from 0.0 to 1.0 for each 1 km grid cell representing suitability depending on topographic and land use and land cover characteristics within the target state. |
| `historical_rural_pop_raster` | string | Full path with file name and extension to a raster containing rural population counts for each 1 km grid cell for the historical base time step. |
| `historical_urban_pop_raster` | string | Full path with file name and extension to a raster containing urban population counts for each 1 km grid cell for the historical base time step. |
| `projected_population_file` | string | Full path with file name and extension to a CSV file containing population projections per year separated into urban and rural categories.,Field descriptions for require fields as follows: (Year, integer, four digit year), (UrbanPop, float, population count for urban), (RuralPop, float, population count for rural), (Scenario, string, scenario as set in the `scenario` variable) |
| `one_dimension_indices_file` | string | Full path with file name and extension to the text file containing a file structured as a Python list (e.g. [0, 1]) that contains the index of each grid cell when flattened from a 2D array to a 1D array for the target state. |
| `output_directory` | string | Full path with file name and extension to the output directory where outputs and the log file will be written. |
| `alpha_urban` | float | Alpha parameter for urban. Represents the degree to which the population size of surrounding cells translates into the suitability of a focal cell.,A positive value indicates that the larger the population that is located within the 100 km neighborhood, the more suitable the focal cell is.,More negative value implies less suitable. Acceptable range:,-2.0 to 2.0 |
| `beta_urban` | float | Beta parameter for urban. Reflects the significance of distance to surrounding cells on the suitability of a focal cell.,Within 100 km, beta determines how distance modifies the effect on suitability. Acceptable range:,-2.0 to 2.0 |
| `alpha_rural` | float | Alpha parameter for rural. Represents the degree to which the population size of surrounding cells translates into the suitability of a focal cell.,A positive value indicates that the larger the population that is located within the 100 km neighborhood, the more suitable the focal cell is.,More negative value implies less suitable. Acceptable range:,-2.0 to 2.0 |
| `beta_rural` | float | Beta parameter for rural. Reflects the significance of distance to surrounding cells on the suitability of a focal cell.,Within 100 km, beta determines how distance modifies the effect on suitability. Acceptable range:,-2.0 to 2.0 |
| `scenario` | string | String representing the scenario with no spaces. Must match what is in the `projected_population_file` if passing population projections in using a file. |
| `state_name` | string | Target state name with no spaces separated by an underscore. |
| `historic_base_year` | integer | Four digit historic base year. |
| `projection_start_year` | integer | Four digit first year to process for the projection. |
| `projection_end_year` | integer | Four digit last year to process for the projection. |
| `time_step` | integer | Number of steps (e.g. number of years between projections) |
| `rural_pop_proj_n` | float | Rural population projection count for the projected year being calculated. These can be read from the `projected_population_file` instead. |
| `urban_pop_proj_n` | float | Urban population projection count for the projected year being calculated. These can be read from the `projected_population_file` instead. |
| `calibration_urban_year_one_raster` | string | Only used for running calibration.  Full path with file name and extension to a raster containing urban population counts for each 1 km grid cell for year one of the calibration. |
| `calibration_urban_year_two_raster` | string | Only used for running calibration.  Full path with file name and extension to a raster containing urban population counts for each 1 km grid cell for year two of the calibration. |
| `calibration_rural_year_one_raster` | string | Only used for running calibration.  Full path with file name and extension to a raster containing rural population counts for each 1 km grid cell for year one of the calibration. |
| `calibration_rural_year_two_raster` | string | Only used for running calibration.  Full path with file name and extension to a raster containing rural population counts for each 1 km grid cell for year two of the calibration. |
| `kernel_distance_meters` | float | Distance kernel in meters; default 100,000 meters. |
| `write_csv` | boolean | Optionally export raster as a CSV file without nodata values; option set to compress CSV using gzip.  Exports values for non-NODATA grid cells as field name `value`. |
| `compress_csv` | boolean | Optionally compress CSV file to GZIP if outputting in CSV; Default True |
| `write_raster` | boolean | Optionally export raster output; Default True |
| `write_array2d` | boolean | Optionally export a NumPy 2D array for each output in the shape of the template raster |
| `write_array1d` | boolean | Optionally export a Numpy 1D flattened array of only grid cells within the target state |
| `run_number` | int | Add on for the file name when running sensitivity analysis |
| `output_total` | boolean | Choice to output total (urban + rural) dataset; Defualt True |

### Variable arguments
Users can update variable argument values after model initialization; this includes updating values between time steps (see **Example 3**).  The following are variable arguments:
- `alpha_urban`
- `beta_urban`
- `alpha_rural`
- `beta_rural`
- `urban_pop_proj_n`
- `rural_pop_proj_n`
- `kernel_distance_meters`

### YAML configuration file option (e.g., config.yml)
Arguments can be passed into the `Model` class using a YAML configuration file as well (see **Example 1**):

```yaml
# Example configuration file setup
grid_coordinates_file: '<Full path with file name and extension to the file>'
historical_rural_pop_raster: '<Full path with file name and extension to the file>'
historical_urban_pop_raster: '<Full path with file name and extension to the file>'
historical_suitability_raster: '<Full path with file name and extension to the file>'
projected_population_file: '<Full path with file name and extension to the file>'
one_dimension_indices_file: '<Full path with file name and extension to the file>'
output_directory: '<Full path with file name and extension to the file>'
alpha_urban: 2.0
alpha_rural: 0.08
beta_urban: 1.78
beta_rural: 1.42
kernel_distance_meters: 100000
scenario: 'SSP2'
state_name: 'vermont'
historic_base_year: 2010
projection_start_year: 2020
projection_end_year: 2050
time_step: 10
write_raster: True
output_total: True

# calibration specific entries
calibration_urban_year_one_raster: '<Full path with file name and extension to the file>'
calibration_urban_year_two_raster: '<Full path with file name and extension to the file>'
calibration_rural_year_one_raster: '<Full path with file name and extension to the file>'
calibration_rural_year_two_raster: '<Full path with file name and extension to the file>'
```

### Generate calibration parameters
If the calibration has not yet been conducted, follow **Example 4** to generate calibration parameters for a target state.

### Expected outputs
Each downscaling run will output a raster for urban, rural, and total population count for each 1 km grid cell for the target state.  These will be written to where the `output_directory` has been assigned.

Each calibration run will output a CSV file containing the calibration parameters for the target state and scenario.  These will be written to where the `output_directory` has been assigned.

## Examples

### Example 1:  Run population downscaling for all years using a configuration file
```python
from population_gravity import Model

run = Model(config_file='<Full path with file name and extension to the YAML configuration file (e.g., config.yml)>')

run.downscale()
```

### Example 2:  Run population downscaling for all years by passing argument values
```python
from population_gravity import Model

run = Model(grid_coordinates_file='<Full path with file name and extension to the file>',
            historical_rural_pop_raster='<Full path with file name and extension to the file>',
            historical_urban_pop_raster='<Full path with file name and extension to the file>',
            historical_suitability_raster='<Full path with file name and extension to the file>',
            projected_population_file='<Full path with file name and extension to the file>',
            one_dimension_indices_file='<Full path with file name and extension to the file>',
            output_directory='<Full path with file name and extension to the file>',
            alpha_urban=1.99999999995073,
            alpha_rural=0.0750326293181678,
            beta_urban=1.77529986067379,
            beta_rural=1.42410799449511,
            kernel_distance_meters=100000,
            scenario='SSP2', # shared socioeconomic pathway abbreviation
            state_name='vermont',
            historic_base_year=2010,
            projection_start_year=2020,
            projection_end_year=2030,
            time_step=10,
            write_raster=True)

run.downscale()
```

### Example 3:  Run population downscaling by year by passing argument values; update value in between time step
```python
from population_gravity import Model

run = Model(grid_coordinates_file='<Full path with file name and extension to the file>',
            historical_rural_pop_raster='<Full path with file name and extension to the file>',
            historical_urban_pop_raster='<Full path with file name and extension to the file>',
            historical_suitability_raster='<Full path with file name and extension to the file>',
            projected_population_file='<Full path with file name and extension to the file>',
            one_dimension_indices_file='<Full path with file name and extension to the file>',
            output_directory='<Full path with file name and extension to the file>',
            alpha_urban=1.99999999995073,
            alpha_rural=0.0750326293181678,
            beta_urban=1.77529986067379,
            beta_rural=1.42410799449511,
            kernel_distance_meters=100000,
            scenario='SSP2', # shared socioeconomic pathway abbreviation
            state_name='vermont',
            historic_base_year=2010,
            projection_start_year=2020,
            projection_end_year=2030,
            time_step=10,
            write_raster=True)

# initialize model
run.initialize()

# downscale year 0
run.advance_step()

# modify the calibrated alpha parameter value for urban
run.alpha_urban = -0.1

# run next step with modified parameters
run.advance_step()

# close out run
run.close()
```

### Example 4:  Calibrate downscaling parameters for a target state using a configuration file
```python
from population_gravity import Model

run = Model(config_file='<Full path with file name and extension to the YAML configuration file (e.g., config.yml)>')

run.calibrate()
```

### Example 5: Join raster values CSV file containing non-NODATA grid cell values to valid X, Y coordinates
```python

import population_gravity as pgr

vaild_coordinates_csv = "<path to file>"
valid_raster_values_csv = "<path to file>"
out_csv = "<path to file>"

# optionally choose not to write CSV and just return the data frame by `out_csv=None`
df = pgr.join_coords_to_value(vaild_coordinates_csv, valid_raster_values_csv, out_csv)
```
