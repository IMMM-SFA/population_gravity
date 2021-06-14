[![build](https://github.com/IMMM-SFA/population_gravity/actions/workflows/build.yml/badge.svg)](https://github.com/IMMM-SFA/population_gravity/actions/workflows/build.yml)
[![codecov](https://codecov.io/gh/IMMM-SFA/population_gravity/branch/master/graph/badge.svg)](https://codecov.io/gh/IMMM-SFA/population_gravity)

# population_gravity

### A model to downscale urban and rural population for a defined state to a 1 km grid

## Overview
The `population_gravity` model allocates aggregate urban and rural populations for a defined region to grid cells within that region. In the IM3 application, the model is applied to US state-level urban and rural populations, which are allocated to 1 km grid cells within states. Allocation to grid cells is based on the relative suitability of each cell. Suitability is calculated using a gravity-based approach, in which the suitability of a given cell is determined by the population of surrounding cells (100 km radius), their distance away, and two parameters, namely alpha and beta. Alpha and beta are estimated based on historical population data and indicate the importance of returns to scale and distance in determining suitability values of cells, respectively. The model is composed of two components: calibration and projection. The calibration component uses historical urban/rural population grids of each state in 2000 and 2010 and an optimization algorithm to estimate the alpha and beta parameters that minimize the absolute difference between the actual population grid in 2010 and the one derived from the model. The two parameters can be modified to reflect distinctive forms of population development that may be desired in different socio-economic scenarios. Once the parameters are defined, the projection component downscales state-level urban/rural population aggregates of each state from 2020 to 2100 under different scenarios to grid cells within the state.

## Associated publications

>Zoraghein, H., & O’Neill, B. C. (2020). US State-level Projections of the Spatial Distribution of Population Consistent with Shared Socioeconomic Pathways. Sustainability, 12(8), 3374. https://doi.org/10.3390/su12083374

The input and output data used in this publication can be found here:

>Zoraghein, H., & O'Neill, B. (2020). Data Supplement: U.S. state-level projections of the spatial distribution of population consistent with Shared Socioeconomic Pathways. (Version v0.1.0) [Data set]. Zenodo. http://doi.org/10.5281/zenodo.3756179

State-level population projections were described in this publication:

>Jiang, L., B.C. O'Neill, H. Zoraghein, and S. Dahlke. 2020. Population scenarios for U.S. states consistent with Shared Socioeconomic Pathways. Environmental Research Letters, https://doi.org/10.1088/1748-9326/aba5b1.

The data produced in Jiang et al. (2020) can be downloaded from here:

>Jiang, L., Dahlke, S., Zoraghein, H., & O'Neill, B.C. (2020). Population scenarios for U.S. states consistent with Shared Socioeconomic Pathways (Version v0.1.0) [Data set]. Zenodo. http://doi.org/10.5281/zenodo.3956412

and the state-level model code used in that publication can be found here:

>Zoraghein, H., R. Nawrotzki, L. Jiang, and S. Dahlke (2020). IMMM-SFA/statepop: v0.1.0 (Version v0.1.0). Zenodo. http://doi.org/10.5281/zenodo.3956703


## Getting Started
The `population_gravity` package uses only **Python 3.6** and up.

### Step 1:
You can install `population_gravity` from GitHub by running the following from your terminal:

`python -m pip install -e git://github.com/IMMM-SFA/population_gravity.git@main#egg=population_gravity`

### Step 2:
Confirm that the module and its dependencies have been installed by running from your prompt:

```python
import population_gravity
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
| `base_rural_pop_raster` | string | Full path with file name and extension to a raster containing rural population counts for each 1 km grid cell for the historical base time step. |
| `base_urban_pop_raster` | string | Full path with file name and extension to a raster containing urban population counts for each 1 km grid cell for the historical base time step. |
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
| `projection_year` | integer | Four digit first year to process for the projection. |
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
Users can update variable argument values after model initialization.  The following are variable arguments:
- `alpha_urban`
- `beta_urban`
- `alpha_rural`
- `beta_rural`
- `urban_pop_proj_n`
- `rural_pop_proj_n`
- `kernel_distance_meters`

### YAML configuration file option (e.g., config.yml)
Arguments can be passed into the `Model` class using a YAML configuration file as well (see **Example 1**):

### Generate calibration parameters
If the calibration has not yet been conducted, follow **Example 2** to generate calibration parameters for a target state.

### Expected outputs
Each downscaling run will output a raster for urban, rural, and total population count for each 1 km grid cell for the target state.  These will be written to where the `output_directory` has been assigned.

Each calibration run will output a CSV file containing the calibration parameters for the target state and scenario.  These will be written to where the `output_directory` has been assigned.

## Examples

### Download the original data
Download and unzip the inputs and outputs as archived in Zoraghein and O'Neill (2020) from the following Zenodo archive:  [zoraghein-oneill_population_gravity_inputs_outputs.zip](https://zenodo.org/record/3756179/files/zoraghein-oneill_population_gravity_inputs_outputs.zip?download=1)

### Example 1:  Run population downscaling for Vermont using year 2010 as the base year to downscale population projections for 2020.  Write outputs as GeoTiff files.
```python
from population_gravity import Model

# instantiate model
run = Model(grid_coordinates_file='<Full path with file name and extension to the file>',
            base_rural_pop_raster='<Full path with file name and extension to the file>',
            base_urban_pop_raster='<Full path with file name and extension to the file>',
            historical_suitability_raster='<Full path with file name and extension to the file>',
            projected_population_file='<Full path with file name and extension to the file>',
            one_dimension_indices_file='<Full path with file name and extension to the file>',
            output_directory='<Full path to the desired directory>',
            alpha_urban=alpha_urban,
            alpha_rural=alpha_rural,
            beta_urban=beta_urban,
            beta_rural=beta_rural,
            kernel_distance_meters=kernel_distance_meters,
            scenario=scenario,
            state_name=target_state,
            historic_base_year=historical_year,
            projection_year=projection_year,
            write_raster=write_raster,
            write_logfile=write_logfile,
            output_total=output_total,
            write_array1d=write_array1d,
            run_number=sample_id)

run.downscale()
```

### Example 2:  Calibrate downscaling parameters for a target state using a configuration file
```python
from population_gravity import Model

run = Model(config_file='<Full path with file name and extension to the YAML configuration file (e.g., config.yml)>')

run.calibrate()
```

### Example 3: Join raster values CSV file containing non-NODATA grid cell values to valid X, Y coordinates
```python

import population_gravity as pgr

vaild_coordinates_csv = "<path to file>"
valid_raster_values_csv = "<path to file>"
out_csv = "<path to file>"

# optionally choose not to write CSV and just return the data frame by `out_csv=None`
df = pgr.join_coords_to_value(vaild_coordinates_csv, valid_raster_values_csv, out_csv)
```
