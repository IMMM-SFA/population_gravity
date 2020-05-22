import pkg_resources
import unittest

from population_gravity import Model


class TestCalibration(unittest.TestCase):

    STATE_NAME = 'vermont'
    SCENARIO = 'SSP2'

    GRID_COORD_FILE = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_coordinates.csv'.format(STATE_NAME))
    PROJ_POP_FILE = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_{}_popproj.csv'.format(STATE_NAME, SCENARIO))
    ONE_D_IND_FILE = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_within_indices.txt'.format(STATE_NAME))
    HIST_SUITABILITY = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_mask_short_term.tif'.format(STATE_NAME))
    OUTPUT_DIRECTORY = pkg_resources.resource_filename('population_gravity', 'tests/data/outputs')
    URBAN_POP_YR1 = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_urban_2000_1km.tif'.format(STATE_NAME))
    URBAN_POP_YR2 = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_urban_2010_1km.tif'.format(STATE_NAME))
    RURAL_POP_YR1 = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_rural_2000_1km.tif'.format(STATE_NAME))
    RURAL_POP_YR2 = pkg_resources.resource_filename('population_gravity', 'tests/data/inputs/{}_rural_2010_1km.tif'.format(STATE_NAME))

    def test_calibration(self):

        run = Model(grid_coordinates_file=TestCalibration.GRID_COORD_FILE,
                    historical_suitability_raster=TestCalibration.HIST_SUITABILITY,
                    projected_population_file=TestCalibration.PROJ_POP_FILE,
                    one_dimension_indices_file=TestCalibration.ONE_D_IND_FILE,
                    output_directory=TestCalibration.OUTPUT_DIRECTORY,
                    calibration_urban_year_one_raster=TestCalibration.URBAN_POP_YR1,
                    calibration_urban_year_two_raster=TestCalibration.URBAN_POP_YR2,
                    calibration_rural_year_one_raster=TestCalibration.RURAL_POP_YR1,
                    calibration_rural_year_two_raster=TestCalibration.RURAL_POP_YR2,
                    kernel_distance_meters=100000,
                    scenario=TestCalibration.SCENARIO,
                    state_name=TestCalibration.STATE_NAME,
                    urban_alpha_lower_bound=-2.0,
                    urban_alpha_upper_bound=2.0,
                    urban_beta_upper_bound=2.0,
                    urban_beta_lower_bound=-0.5,
                    rural_alpha_lower_bound=-2.0,
                    rural_alpha_upper_bound=2.0,
                    rural_beta_upper_bound=2.0,
                    rural_beta_lower_bound=-0.5)

        run.calibrate()


if __name__ == '__main__':
    unittest.main()
