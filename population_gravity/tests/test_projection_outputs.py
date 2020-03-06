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

    COMP_RURAL_2010 = pkg_resources.resource_filename('population_gravity', 'tests/data/comp_data/district_of_columbia_1km_SSP2_Rural_2010.tif')
    COMP_RURAL_2020 = pkg_resources.resource_filename('population_gravity', 'tests/data/comp_data/district_of_columbia_1km_SSP2_Rural_2020.tif')
    COMP_URBAN_2010 = pkg_resources.resource_filename('population_gravity', 'tests/data/comp_data/district_of_columbia_1km_SSP2_Urban_2010.tif')
    COMP_URBAN_2020 = pkg_resources.resource_filename('population_gravity', 'tests/data/comp_data/district_of_columbia_1km_SSP2_Urban_2020.tif')
    COMP_TOTAL_2010 = pkg_resources.resource_filename('population_gravity', 'tests/data/comp_data/district_of_columbia_1km_SSP2_Total_2010.tif')
    COMP_TOTAL_2020 = pkg_resources.resource_filename('population_gravity', 'tests/data/comp_data/district_of_columbia_1km_SSP2_Total_2020.tif')

    HISTORIC_DATA = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/main_inputs')
    PROJECTION_DATA = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/projection')
    OUTPUT_DIRECTORY = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs')
    PARAM_FILE = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/district_of_columbia_SSP2_calibration_parameters.csv')

    RUN_RURAL_2010 = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs/district_of_columbia_1km_SSP2_Rural_2010.tif')
    RUN_RURAL_2020 = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs/district_of_columbia_1km_SSP2_Rural_2020.tif')
    RUN_URBAN_2010 = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs/district_of_columbia_1km_SSP2_Urban_2010.tif')
    RUN_URBAN_2020 = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs/district_of_columbia_1km_SSP2_Urban_2020.tif')
    RUN_TOTAL_2010 = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs/district_of_columbia_1km_SSP2_Total_2010.tif')
    RUN_TOTAL_2020 = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs/district_of_columbia_1km_SSP2_Total_2020.tif')

    def test_proj_outputs(self):
        """Test for projection outputs"""

        # read in comp data from rasters to arrays
        comp_rural_2010 = raster_to_array(TestProjectedOutputs.COMP_RURAL_2010)
        comp_rural_2020 = raster_to_array(TestProjectedOutputs.COMP_RURAL_2020)
        comp_urban_2010 = raster_to_array(TestProjectedOutputs.COMP_URBAN_2010)
        comp_urban_2020 = raster_to_array(TestProjectedOutputs.COMP_URBAN_2020)
        comp_total_2010 = raster_to_array(TestProjectedOutputs.COMP_TOTAL_2010)
        comp_total_2020 = raster_to_array(TestProjectedOutputs.COMP_TOTAL_2020)

        run = Model(datadir_histdata=TestProjectedOutputs.HISTORIC_DATA,
                        ssp_data_directory=TestProjectedOutputs.PROJECTION_DATA,
                        ssp_code='SSP2',
                        region_code='district_of_columbia',
                        output_directory=TestProjectedOutputs.OUTPUT_DIRECTORY,
                        calibration_parameters_file=TestProjectedOutputs.PARAM_FILE,
                        start_year=2000,
                        end_year=2020,
                        time_step=10).downscale()

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
