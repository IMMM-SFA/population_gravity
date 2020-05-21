import pkg_resources
import unittest

import pandas as pd

import population_gravity as pgr


class TestUtilities(unittest.TestCase):

    STATE_NAME = 'vermont'
    SSP = 'SSP2'

    VALID_COORDS_CSV = pkg_resources.resource_filename('population_gravity', f'tests/data/inputs/{STATE_NAME}_valid_coordinates.csv')
    VALID_VALUES_CSV = pkg_resources.resource_filename('population_gravity', f'tests/data/comp_data/{STATE_NAME}_1km_{SSP}_Total_2030.csv')
    COMP_CSV = pkg_resources.resource_filename('population_gravity', f'tests/data/comp_data/{STATE_NAME}_1km_{SSP}_Total_2030_coords.csv')

    def test_join_coords_to_value(self):
        """Tests for `join_coords_to_value` function."""

        # read in comp data
        comp_df = pd.read_csv(TestUtilities.COMP_CSV)

        # generate output data frame
        df = pgr.join_coords_to_value(TestUtilities.VALID_COORDS_CSV,
                                      TestUtilities.VALID_VALUES_CSV,
                                      out_csv=None)

        # test for equality in outputs
        pd.testing.assert_frame_equal(comp_df, df)


if __name__ == '__main__':
    unittest.main()
