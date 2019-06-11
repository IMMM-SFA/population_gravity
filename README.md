# spatial_population_downscaling_model

### A model to allocate urban and rural populations for a defined region to a grid

## Overview
The Spatial Population Downscaling Model allocates aggregate urban and rural populations for a defined region to grid cells within that region. In the IM3 application, the model is applied to US state-level urban and rural populations, which are allocated to 1 km grid cells within states. Allocation to grid cells is based on the relative suitability of each cell. Suitability is calculated using a gravity-based approach, in which the suitability of a given cell is determined by the population of surrounding cells (100 km radius), their distance away, and two parameters, namely alpha and beta. Alpha and beta are estimated based on historical population data and indicate the importance of returns to scale and distance in determining suitability values of cells, respectively. The model is composed of two components: calibration and projection. The calibration component uses historical urban/rural population grids of each state in 2000 and 2010 and an optimization algorithm to estimate the alpha and beta parameters that minimize the absolute difference between the actual population grid in 2010 and the one derived from the model. The two parameters can be modified to reflect distinctive forms of population development that may be desired in different socio-economic scenarios. Once the parameters are defined, the projection component downscales state-level urban/rural population aggregates of each state from 2020 to 2100 under different scenarios to grid cells within the state.

## Getting Started

Fill in once we get the package built

### Accompanying input data

Fill in with DOI link once we get it setup

## Example

Fill in once we get the package built

## Contact
For questions please contact:

Hamidreza Zoraghein <Hamidreza.Zoraghein@du.edu>
Brian O'Neill <Brian.ONeill@du.edu>
