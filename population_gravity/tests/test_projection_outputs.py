"""test_projection_outputs.py

Tests to ensure high-level functionality and outputs remain consistent.

@author Chris R. Vernon (chris.vernon@pnnl.gov)
@license BSD 2-Clause

"""

import pkg_resources
import unittest

import numpy.testing as np

from population_gravity.downscale_utilities import raster_to_array
from population_gravity import Model


class TestProjectedOutputs(unittest.TestCase):
    """Test configuration integrity."""

    COMP_RURAL_2010 = pkg_resources.resource_filename('population_gravity', 'tests/data/comp_data/rhode_island_1km_SSP2_Rural_2010.tif')
    COMP_RURAL_2020 = pkg_resources.resource_filename('population_gravity', 'tests/data/comp_data/rhode_island_1km_SSP2_Rural_2020.tif')
    COMP_URBAN_2010 = pkg_resources.resource_filename('population_gravity', 'tests/data/comp_data/rhode_island_1km_SSP2_Urban_2010.tif')
    COMP_URBAN_2020 = pkg_resources.resource_filename('population_gravity', 'tests/data/comp_data/rhode_island_1km_SSP2_Urban_2020.tif')
    COMP_TOTAL_2010 = pkg_resources.resource_filename('population_gravity', 'tests/data/comp_data/rhode_island_1km_SSP2_Total_2010.tif')
    COMP_TOTAL_2020 = pkg_resources.resource_filename('population_gravity', 'tests/data/comp_data/rhode_island_1km_SSP2_Total_2020.tif')

    STATE_NAME = 'rhode_island'
    SCENARIO = 'SSP2'
    GRID_COORD_FILE = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_coordinates.csv'.format(STATE_NAME))
    HIST_RURAL_RASTER = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_rural_2000_1km.tif'.format(STATE_NAME))
    HIST_URBAN_RASTER = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_urban_2000_1km.tif'.format(STATE_NAME))
    PROJ_POP_FILE = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/projection/{}_{}_popproj.csv'.format(STATE_NAME, SCENARIO))
    ONE_D_IND_FILE = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_within_indices.txt'.format(STATE_NAME))
    HIST_SUITABILITY = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_mask_short_term.tif'.format(STATE_NAME))
    OUTPUT_DIRECTORY = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs')

    RUN_RURAL_2010 = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs/rhode_island_1km_SSP2_Rural_2010.tif')
    RUN_RURAL_2020 = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs/rhode_island_1km_SSP2_Rural_2020.tif')
    RUN_URBAN_2010 = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs/rhode_island_1km_SSP2_Urban_2010.tif')
    RUN_URBAN_2020 = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs/rhode_island_1km_SSP2_Urban_2020.tif')
    RUN_TOTAL_2010 = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs/rhode_island_1km_SSP2_Total_2010.tif')
    RUN_TOTAL_2020 = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs/rhode_island_1km_SSP2_Total_2020.tif')

    def test_proj_outputs(self):
        """Test for projection outputs"""

        # read in comp data from rasters to arrays
        comp_rural_2010 = raster_to_array(TestProjectedOutputs.COMP_RURAL_2010)
        comp_rural_2020 = raster_to_array(TestProjectedOutputs.COMP_RURAL_2020)
        comp_urban_2010 = raster_to_array(TestProjectedOutputs.COMP_URBAN_2010)
        comp_urban_2020 = raster_to_array(TestProjectedOutputs.COMP_URBAN_2020)
        comp_total_2010 = raster_to_array(TestProjectedOutputs.COMP_TOTAL_2010)
        comp_total_2020 = raster_to_array(TestProjectedOutputs.COMP_TOTAL_2020)

        run = Model(grid_coordinates_file=TestProjectedOutputs.GRID_COORD_FILE,
                    historical_rural_pop_raster=TestProjectedOutputs.HIST_RURAL_RASTER,
                    historical_urban_pop_raster=TestProjectedOutputs.HIST_URBAN_RASTER,
                    historical_suitability_raster=TestProjectedOutputs.HIST_SUITABILITY,
                    projected_population_file=TestProjectedOutputs.PROJ_POP_FILE,
                    one_dimension_indices_file=TestProjectedOutputs.ONE_D_IND_FILE,
                    output_directory=TestProjectedOutputs.OUTPUT_DIRECTORY,
                    alpha_urban=-2,
                    alpha_rural=-0.34,
                    beta_urban=0.46,
                    beta_rural=1.0,
                    scenario=TestProjectedOutputs.SCENARIO,
                    state_name=TestProjectedOutputs.STATE_NAME,
                    historic_base_year=2000,
                    projection_start_year=2010,
                    projection_end_year=2020,
                    time_step=10)

        run.downscale()

        # read in run data from rasters to arrays
        run_rural_2010 = raster_to_array(TestProjectedOutputs.RUN_RURAL_2010)
        run_rural_2020 = raster_to_array(TestProjectedOutputs.RUN_RURAL_2020)
        run_urban_2010 = raster_to_array(TestProjectedOutputs.RUN_URBAN_2010)
        run_urban_2020 = raster_to_array(TestProjectedOutputs.RUN_URBAN_2020)
        run_total_2010 = raster_to_array(TestProjectedOutputs.RUN_TOTAL_2010)
        run_total_2020 = raster_to_array(TestProjectedOutputs.RUN_TOTAL_2020)

        # test equality
        np.assert_array_equal(comp_rural_2010, run_rural_2010)
        np.assert_array_equal(comp_rural_2020, run_rural_2020)
        np.assert_array_equal(comp_urban_2010, run_urban_2010)
        np.assert_array_equal(comp_urban_2020, run_urban_2020)
        np.assert_array_equal(comp_total_2010, run_total_2010)
        np.assert_array_equal(comp_total_2020, run_total_2020)


if __name__ == '__main__':
    unittest.main()
