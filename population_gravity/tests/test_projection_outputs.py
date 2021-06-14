"""test_projection_outputs.py

Tests to ensure high-level functionality and outputs remain consistent.

@author Chris R. Vernon (chris.vernon@pnnl.gov)
@license BSD 2-Clause

"""

import os
import pkg_resources
import tempfile
import unittest

import numpy as np

import population_gravity.downscale_utilities as utils
from population_gravity import Model


class TestProjectedOutputs(unittest.TestCase):
    """Test configuration integrity."""

    # set temporary output dir
    TEMP_DIR = tempfile.TemporaryDirectory()

    STATE_NAME = 'vermont'
    SCENARIO = 'ssp2'

    COMP_RURAL_2020 = utils.raster_to_array(pkg_resources.resource_filename('population_gravity', f'tests/data/comp_data/{STATE_NAME}_1km_{SCENARIO}_rural_2020.tif'))
    COMP_URBAN_2020 = utils.raster_to_array(pkg_resources.resource_filename('population_gravity', f'tests/data/comp_data/{STATE_NAME}_1km_{SCENARIO}_urban_2020.tif'))
    COMP_TOTAL_2020 = utils.raster_to_array(pkg_resources.resource_filename('population_gravity', f'tests/data/comp_data/{STATE_NAME}_1km_{SCENARIO}_total_2020.tif'))

    GRID_COORD_FILE = pkg_resources.resource_filename('population_gravity', f'tests/data/inputs/{STATE_NAME}_coordinates.csv')
    HIST_RURAL_RASTER = pkg_resources.resource_filename('population_gravity', f'tests/data/inputs/{STATE_NAME}_rural_2010_1km.tif')
    HIST_URBAN_RASTER = pkg_resources.resource_filename('population_gravity', f'tests/data/inputs/{STATE_NAME}_urban_2010_1km.tif')
    PROJ_POP_FILE = pkg_resources.resource_filename('population_gravity', f'tests/data/inputs/{STATE_NAME}_{SCENARIO}_popproj.csv')
    ONE_D_IND_FILE = pkg_resources.resource_filename('population_gravity', f'tests/data/inputs/{STATE_NAME}_within_indices.txt')
    HIST_SUITABILITY = pkg_resources.resource_filename('population_gravity', f'tests/data/inputs/{STATE_NAME}_mask_short_term.tif')
    OUTPUT_DIRECTORY = TEMP_DIR.name

    # output rasters
    RUN_RURAL_2020 = os.path.join(OUTPUT_DIRECTORY, f'{STATE_NAME}_1km_{SCENARIO}_rural_2020.tif')
    RUN_URBAN_2020 = os.path.join(OUTPUT_DIRECTORY, f'{STATE_NAME}_1km_{SCENARIO}_urban_2020.tif')
    RUN_TOTAL_2020 = os.path.join(OUTPUT_DIRECTORY, f'{STATE_NAME}_1km_{SCENARIO}_total_2020.tif')

    # output 2D arrays
    RUN_RURAL_2020_2NPY = os.path.join(OUTPUT_DIRECTORY, f'{STATE_NAME}_1km_{SCENARIO}_rural_2020_2d.npy')
    RUN_URBAN_2020_2NPY = os.path.join(OUTPUT_DIRECTORY, f'{STATE_NAME}_1km_{SCENARIO}_urban_2020_2d.npy')
    RUN_TOTAL_2020_2NPY = os.path.join(OUTPUT_DIRECTORY, f'{STATE_NAME}_1km_{SCENARIO}_total_2020_2d.npy')

    # output 1D arrays
    RUN_RURAL_2020_1NPY = os.path.join(OUTPUT_DIRECTORY, f'{STATE_NAME}_1km_{SCENARIO}_rural_2020_1d.npy')
    RUN_URBAN_2020_1NPY = os.path.join(OUTPUT_DIRECTORY, f'{STATE_NAME}_1km_{SCENARIO}_urban_2020_1d.npy')
    RUN_TOTAL_2020_1NPY = os.path.join(OUTPUT_DIRECTORY, f'{STATE_NAME}_1km_{SCENARIO}_total_2020_1d.npy')

    # output CSV.GZ files
    RUN_RURAL_2020_CSV = os.path.join(OUTPUT_DIRECTORY, f'{STATE_NAME}_1km_{SCENARIO}_rural_2020.csv.gz')
    RUN_URBAN_2020_CSV = os.path.join(OUTPUT_DIRECTORY, f'{STATE_NAME}_1km_{SCENARIO}_urban_2020.csv.gz')
    RUN_TOTAL_2020_CSV = os.path.join(OUTPUT_DIRECTORY, f'{STATE_NAME}_1km_{SCENARIO}_total_2020.csv.gz')

    def test_raster_outputs(self):
        """Test for raster outputs."""

        run = Model(grid_coordinates_file=TestProjectedOutputs.GRID_COORD_FILE,
                    base_rural_pop_raster=TestProjectedOutputs.HIST_RURAL_RASTER,
                    base_urban_pop_raster=TestProjectedOutputs.HIST_URBAN_RASTER,
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
                    projection_year=2020,
                    write_raster=True,
                    write_array1d=False,
                    write_array2d=False,
                    write_csv=False,
                    write_logfile=False,
                    write_suitability=False)

        run.downscale()

        # test equality
        self.equality(file_type='tif')

    def test_array1d_outputs(self):
        """Test for 1D array outputs."""

        run = Model(grid_coordinates_file=TestProjectedOutputs.GRID_COORD_FILE,
                    base_rural_pop_raster=TestProjectedOutputs.HIST_RURAL_RASTER,
                    base_urban_pop_raster=TestProjectedOutputs.HIST_URBAN_RASTER,
                    historical_suitability_raster=TestProjectedOutputs.HIST_SUITABILITY,
                    projected_population_file=TestProjectedOutputs.PROJ_POP_FILE,
                    one_dimension_indices_file=TestProjectedOutputs.ONE_D_IND_FILE,
                    output_directory=TestProjectedOutputs.OUTPUT_DIRECTORY,
                    alpha_urban=1.9999999999507383,
                    alpha_rural=0.07503262931816788,
                    beta_urban=1.7752998606737922,
                    beta_rural=1.4241079944951196,
                    kernel_distance_meters=100000,
                    scenario=TestProjectedOutputs.SCENARIO,
                    state_name=TestProjectedOutputs.STATE_NAME,
                    historic_base_year=2010,
                    projection_year=2020,
                    write_raster=False,
                    write_array1d=True,
                    write_array2d=False,
                    write_csv=False,
                    write_suitability=False,
                    write_logfile=False)

        run.downscale()

        # test equality
        self.equality(file_type='1d.npy', run_object=run)

    def test_array2d_outputs(self):
        """Test for 2D array outputs."""

        run = Model(grid_coordinates_file=TestProjectedOutputs.GRID_COORD_FILE,
                    base_rural_pop_raster=TestProjectedOutputs.HIST_RURAL_RASTER,
                    base_urban_pop_raster=TestProjectedOutputs.HIST_URBAN_RASTER,
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
                    projection_year=2020,
                    write_raster=False,
                    write_array1d=False,
                    write_array2d=True,
                    write_csv=False,
                    write_logfile=False)

        run.downscale()

        # test equality
        self.equality(file_type='2d.npy')

    def test_csv_outputs(self):
        """Test for CSV outputs."""

        run = Model(grid_coordinates_file=TestProjectedOutputs.GRID_COORD_FILE,
                    base_rural_pop_raster=TestProjectedOutputs.HIST_RURAL_RASTER,
                    base_urban_pop_raster=TestProjectedOutputs.HIST_URBAN_RASTER,
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
                    projection_year=2020,
                    write_raster=False,
                    write_array1d=False,
                    write_array2d=False,
                    write_csv=True,
                    write_logfile=False)

        run.downscale()

        # test equality
        self.equality(file_type='csv.gz', run_object=run)

    def equality(self, file_type, run_object=None):
        """Generate equality tests."""

        if file_type == 'tif':
            run_urban_2020 = utils.raster_to_array(TestProjectedOutputs.RUN_URBAN_2020)
            run_rural_2020 = utils.raster_to_array(TestProjectedOutputs.RUN_RURAL_2020)
            run_total_2020 = utils.raster_to_array(TestProjectedOutputs.RUN_TOTAL_2020)

        elif file_type == '1d.npy':

            run_urban_2020 = utils.array_to_raster_memory(run_object.template_raster,
                                                          np.load(TestProjectedOutputs.RUN_URBAN_2020_1NPY),
                                                          run_object.one_dimension_indices).read(1)

            run_rural_2020 = utils.array_to_raster_memory(run_object.template_raster,
                                                          np.load(TestProjectedOutputs.RUN_RURAL_2020_1NPY),
                                                          run_object.one_dimension_indices).read(1)

            run_total_2020 = utils.array_to_raster_memory(run_object.template_raster,
                                                          np.load(TestProjectedOutputs.RUN_TOTAL_2020_1NPY),
                                                          run_object.one_dimension_indices).read(1)

        elif file_type == '2d.npy':
            run_urban_2020 = np.load(TestProjectedOutputs.RUN_URBAN_2020_2NPY)
            run_rural_2020 = np.load(TestProjectedOutputs.RUN_RURAL_2020_2NPY)
            run_total_2020 = np.load(TestProjectedOutputs.RUN_TOTAL_2020_2NPY)

        else:
            arr_urban_2020 = utils.gzip_csv_to_array(TestProjectedOutputs.RUN_URBAN_2020_CSV, header_rows_to_skip=1)
            arr_rural_2020 = utils.gzip_csv_to_array(TestProjectedOutputs.RUN_RURAL_2020_CSV, header_rows_to_skip=1)
            arr_total_2020 = utils.gzip_csv_to_array(TestProjectedOutputs.RUN_TOTAL_2020_CSV, header_rows_to_skip=1)

            run_urban_2020 = utils.array_to_raster_memory(run_object.template_raster, arr_urban_2020, run_object.one_dimension_indices).read(1)
            run_rural_2020 = utils.array_to_raster_memory(run_object.template_raster, arr_rural_2020, run_object.one_dimension_indices).read(1)
            run_total_2020 = utils.array_to_raster_memory(run_object.template_raster, arr_total_2020, run_object.one_dimension_indices).read(1)

        np.testing.assert_array_almost_equal(TestProjectedOutputs.COMP_URBAN_2020, run_urban_2020, decimal=4)
        np.testing.assert_array_almost_equal(TestProjectedOutputs.COMP_RURAL_2020, run_rural_2020, decimal=4)
        np.testing.assert_array_almost_equal(TestProjectedOutputs.COMP_TOTAL_2020, run_total_2020, decimal=4)


if __name__ == '__main__':
    unittest.main()
