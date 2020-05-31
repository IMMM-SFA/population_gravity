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

    STATE_NAME = 'vermont'
    SCENARIO = 'SSP2'

    COMP_RURAL_2030 = pkg_resources.resource_filename('population_gravity', 'tests/data/comp_data/{}_1km_{}_Rural_2030.tif'.format(STATE_NAME, SCENARIO))
    COMP_RURAL_2020 = pkg_resources.resource_filename('population_gravity', 'tests/data/comp_data/{}_1km_{}_Rural_2020.tif'.format(STATE_NAME, SCENARIO))
    COMP_URBAN_2030 = pkg_resources.resource_filename('population_gravity', 'tests/data/comp_data/{}_1km_{}_Urban_2030.tif'.format(STATE_NAME, SCENARIO))
    COMP_URBAN_2020 = pkg_resources.resource_filename('population_gravity', 'tests/data/comp_data/{}_1km_{}_Urban_2020.tif'.format(STATE_NAME, SCENARIO))
    COMP_TOTAL_2030 = pkg_resources.resource_filename('population_gravity', 'tests/data/comp_data/{}_1km_{}_Total_2030.tif'.format(STATE_NAME, SCENARIO))
    COMP_TOTAL_2020 = pkg_resources.resource_filename('population_gravity', 'tests/data/comp_data/{}_1km_{}_Total_2020.tif'.format(STATE_NAME, SCENARIO))

    GRID_COORD_FILE = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_coordinates.csv'.format(STATE_NAME))
    HIST_RURAL_RASTER = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_rural_2010_1km.tif'.format(STATE_NAME))
    HIST_URBAN_RASTER = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_urban_2010_1km.tif'.format(STATE_NAME))
    PROJ_POP_FILE = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_{}_popproj.csv'.format(STATE_NAME, SCENARIO))
    ONE_D_IND_FILE = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_within_indices.txt'.format(STATE_NAME))
    HIST_SUITABILITY = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_mask_short_term.tif'.format(STATE_NAME))
    OUTPUT_DIRECTORY = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs')

    RUN_RURAL_2030 = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs/{}_1km_{}_Rural_2030.tif'.format(STATE_NAME, SCENARIO))
    RUN_RURAL_2020 = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs/{}_1km_{}_Rural_2020.tif'.format(STATE_NAME, SCENARIO))
    RUN_URBAN_2030 = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs/{}_1km_{}_Urban_2030.tif'.format(STATE_NAME, SCENARIO))
    RUN_URBAN_2020 = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs/{}_1km_{}_Urban_2020.tif'.format(STATE_NAME, SCENARIO))
    RUN_TOTAL_2030 = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs/{}_1km_{}_Total_2030.tif'.format(STATE_NAME, SCENARIO))
    RUN_TOTAL_2020 = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs/{}_1km_{}_Total_2020.tif'.format(STATE_NAME, SCENARIO))

    def test_proj_outputs(self):
        """Test for projection outputs"""

        # read in comp data from rasters to arrays
        comp_rural_2030 = raster_to_array(TestProjectedOutputs.COMP_RURAL_2030)
        comp_rural_2020 = raster_to_array(TestProjectedOutputs.COMP_RURAL_2020)
        comp_urban_2030 = raster_to_array(TestProjectedOutputs.COMP_URBAN_2030)
        comp_urban_2020 = raster_to_array(TestProjectedOutputs.COMP_URBAN_2020)
        comp_total_2030 = raster_to_array(TestProjectedOutputs.COMP_TOTAL_2030)
        comp_total_2020 = raster_to_array(TestProjectedOutputs.COMP_TOTAL_2020)

        run = Model(grid_coordinates_file=TestProjectedOutputs.GRID_COORD_FILE,
                    historical_rural_pop_raster=TestProjectedOutputs.HIST_RURAL_RASTER,
                    historical_urban_pop_raster=TestProjectedOutputs.HIST_URBAN_RASTER,
                    historical_suitability_raster=TestProjectedOutputs.HIST_SUITABILITY,
                    projected_population_file=TestProjectedOutputs.PROJ_POP_FILE,
                    one_dimension_indices_file=TestProjectedOutputs.ONE_D_IND_FILE,
                    output_directory=TestProjectedOutputs.OUTPUT_DIRECTORY,
                    alpha_urban=1.99999999995073,
                    alpha_rural=0.0750326293181678,
                    beta_urban=1.77529986067379,
                    beta_rural=1.42410799449511,
                    kernel_distance_meters=100000,
                    scenario=TestProjectedOutputs.SCENARIO,
                    state_name=TestProjectedOutputs.STATE_NAME,
                    historic_base_year=2010,
                    projection_start_year=2020,
                    projection_end_year=2030,
                    time_step=10,
                    write_raster=True,
                    write_array=False,
                    write_logfile=False)

        run.downscale()

        # read in run data from rasters to arrays
        run_urban_2020 = raster_to_array(TestProjectedOutputs.RUN_URBAN_2020)
        run_rural_2020 = raster_to_array(TestProjectedOutputs.RUN_RURAL_2020)
        run_total_2020 = raster_to_array(TestProjectedOutputs.RUN_TOTAL_2020)
        run_rural_2030 = raster_to_array(TestProjectedOutputs.RUN_RURAL_2030)
        run_urban_2030 = raster_to_array(TestProjectedOutputs.RUN_URBAN_2030)
        run_total_2030 = raster_to_array(TestProjectedOutputs.RUN_TOTAL_2030)

        # test equality
        np.assert_array_almost_equal(comp_urban_2020, run_urban_2020, decimal=5)
        np.assert_array_almost_equal(comp_rural_2020, run_rural_2020, decimal=5)
        np.assert_array_almost_equal(comp_total_2020, run_total_2020, decimal=5)
        np.assert_array_almost_equal(comp_rural_2030, run_rural_2030, decimal=5)
        np.assert_array_almost_equal(comp_urban_2030, run_urban_2030, decimal=5)
        np.assert_array_almost_equal(comp_total_2030, run_total_2030, decimal=5)


if __name__ == '__main__':
    unittest.main()
