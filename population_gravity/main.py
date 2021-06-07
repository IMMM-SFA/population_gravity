"""
Model interface for the popultation_gravity model

@author   Chris R. Vernon
@email:   chris.vernon@pnnl.gov

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

"""

import logging
import os
import time

import population_gravity.downscale_calibrate as calib

from population_gravity.process_step import ProcessStep

# Logger inherits ReadConfig
from population_gravity.logger import Logger


class Model(Logger):
    """Run the population downscaling model.

    :param config_file:                         string. Full path to configuration YAML file with file name and
                                                extension. If not provided by the user, the code will default to the
                                                expectation of alternate arguments.

    :param grid_coordinates_file:               string. Full path with file name and extension to the CSV file
                                                containing the coordinates for each 1 km grid cell within the target
                                                state. File includes a header with the fields XCoord, YCoord, FID.
                                                Where data types and field descriptions are as follows:
                                                (XCoord, float, X coordinate in meters),
                                                (YCoord, float, Y coordinate in meters),
                                                (FID, int, Unique feature id)

    :param historical_suitability_raster:       string. Full path with file name and extension to the suitability
                                                raster containing values from 0.0 to 1.0 for each 1 km grid cell
                                                representing suitability depending on topographic and land use and
                                                land cover characteristics within the target state.

    :param base_rural_pop_raster:         string. Full path with file name and extension to a raster containing
                                                rural population counts for each 1 km grid cell for the historical
                                                base time step.

    :param base_urban_pop_raster:         string. Full path with file name and extension to a raster containing
                                                urban population counts for each 1 km grid cell for the historical
                                                base time step.

    :param projected_population_file:           string. Full path with file name and extension to a CSV file containing
                                                population projections per year separated into urban and rural
                                                categories.  Field descriptions for require fields as follows:
                                                (Year, integer, four digit year),
                                                (UrbanPop, float, population count for urban),
                                                (RuralPop, float, population count for rural),
                                                (Scenario, string, scenario as set in the `scenario` variable)

    :param one_dimension_indices_file:          string. Full path with file name and extension to the text file
                                                containing a file structured as a Python list (e.g. [0, 1]) that
                                                contains the index of each grid cell when flattened from a 2D array to
                                                a 1D array for the target state.

    :param output_directory:                    string. Full path with file name and extension to the output directory
                                                where outputs and the log file will be written.

    :param alpha_urban:                         float. Alpha parameter for urban. Represents the degree to which the
                                                population size of surrounding cells translates into the suitability
                                                of a focal cell.  A positive value indicates that the larger the
                                                population that is located within the 100 km neighborhood, the more
                                                suitable the focal cell is.  More negative value implies less suitable.
                                                Acceptable range:  -2.0 to 2.0

    :param beta_urban:                          float. Beta parameter for urban. Reflects the significance of distance
                                                to surrounding cells on the suitability of a focal cell.  Within 100 km,
                                                beta determines how distance modifies the effect on suitability.
                                                Acceptable range:  -0.5 to 2.0

    :param alpha_rural:                         float. Alpha parameter for rural. Represents the degree to which the
                                                population size of surrounding cells translates into the suitability
                                                of a focal cell.  A positive value indicates that the larger the
                                                population that is located within the 100 km neighborhood, the more
                                                suitable the focal cell is.  More negative value implies less suitable.
                                                Acceptable range:  -2.0 to 2.0

    :param beta_rural:                          float. Beta parameter for rural. Reflects the significance of distance
                                                to surrounding cells on the suitability of a focal cell.  Within 100 km,
                                                beta determines how distance modifies the effect on suitability.
                                                Acceptable range:  -0.5 to 2.0

    :param scenario:                            string. String representing the scenario with no spaces. Must match
                                                what is in the `projected_population_file` if passing population
                                                projections in using a file.

    :param state_name:                          string. Target state name with no spaces separated by an underscore.

    :param historic_base_year:                  int. Four digit historic base year.

    :param projection_year:               int. Four digit first year to process for the projection.

    :param projection_end_year:                 int. Four digit last year to process for the projection.

    :param time_step:                           int. Number of steps (e.g. number of years between projections)

    :param rural_pop_proj_n:                    float.  Rural population projection count for the projected year being
                                                calculated.  These can be read from the `projected_population_file`
                                                instead.

    :param urban_pop_proj_n:                    float.  Urban population projection count for the projected year being
                                                calculated.  These can be read from the `projected_population_file`
                                                instead.

    :param calibration_urban_year_one_raster:   string. Only used for running calibration.  Full path with file name and
                                                extension to a raster containing urban population counts for each 1 km
                                                grid cell for year one of the calibration.

    :param calibration_urban_year_two_raster:   string. Only used for running calibration.  Full path with file name and
                                                extension to a raster containing urban population counts for each 1 km
                                                grid cell for year two of the calibration.

    :param calibration_rural_year_one_raster:   string. Only used for running calibration.  Full path with file name and
                                                extension to a raster containing rural population counts for each 1 km
                                                grid cell for year one of the calibration.

    :param calibration_rural_year_two_raster:   string. Only used for running calibration.  Full path with file name and
                                                extension to a raster containing rural population counts for each 1 km
                                                grid cell for year two of the  calibration.

    :param kernel_distance_meters:              float. Distance kernel in meters; default 100,000 meters.

    :param write_raster:                        boolean. Optionally export raster output; Default True

    :param write_csv:                           boolean. Optionally export raster as a CSV file without nodata values

    :param write_array2d:                       boolean. Optionally export a NumPy 2D array for each output in the shape
                                                of the template raster

    :param write_array1d:                       boolean. Optionally export a Numpy 1D flattened array of only grid cells
                                                within the target state

    :param run_number:                          int. Add on for the file name when running sensitivity analysis

    :param write_logfile:                       boolean.  Optionally write log to file.; Default True

    :param compress_csv:                        boolean.  Optionally compress CSV file to GZIP if outputting in CSV

    :param output_total:                        boolean.  Choice to output total (urban + rural) dataset; Defualt True

    """

    def __init__(self, config_file=None, grid_coordinates_file=None, historical_suitability_raster=None,
                 base_rural_pop_raster=None, base_urban_pop_raster=None, projected_population_file=None,
                 one_dimension_indices_file=None, output_directory=None, alpha_urban=None, beta_urban=None,
                 alpha_rural=None, beta_rural=None, scenario=None, state_name=None, historic_base_year=None,
                 projection_year=None, rural_pop_proj_n=None,
                 urban_pop_proj_n=None, calibration_urban_year_one_raster=None, calibration_urban_year_two_raster=None,
                 calibration_rural_year_one_raster=None, calibration_rural_year_two_raster=None,
                 kernel_distance_meters=None, write_raster=True, write_csv=False, write_array1d=False,
                 write_array2d=False, run_number='', write_logfile=True, compress_csv=True, output_total=True,
                 write_suitability=False):

        super(Logger, self).__init__(config_file, grid_coordinates_file, historical_suitability_raster,
                                     base_rural_pop_raster, base_urban_pop_raster,
                                     projected_population_file, one_dimension_indices_file, output_directory,
                                     alpha_urban, beta_urban, alpha_rural, beta_rural, scenario, state_name,
                                     historic_base_year, projection_year,
                                     rural_pop_proj_n, urban_pop_proj_n, calibration_urban_year_one_raster,
                                     calibration_urban_year_two_raster, calibration_rural_year_one_raster,
                                     calibration_rural_year_two_raster, kernel_distance_meters, write_raster, write_csv,
                                     write_array1d, write_array2d, run_number, write_logfile, compress_csv,
                                     output_total, write_suitability)

    @staticmethod
    def make_dir(pth):
        """Create dir if not exists."""

        if not os.path.exists(pth):
            os.makedirs(pth)

    def initialize(self):
        """Setup model."""

        # build output directory first to store logfile and other outputs
        self.make_dir(self.output_directory)

        # initialize logger
        self.initialize_logger()

        logging.info("Start time:  {}".format(time.strftime(self.datetime_format)))

        # log run parameters
        logging.info("Input parameters:")
        logging.info("\tbase_urban_pop_raster = {}".format(self.base_urban_pop_raster))
        logging.info("\tbase_rural_pop_raster = {}".format(self.base_rural_pop_raster))
        logging.info("\tone_dimension_indices_file = {}".format(self.one_dimension_indices_file))

        # for projection
        logging.info("\tprojection_year = {}".format(self.projection_year))
        logging.info("\tgrid_coordinates_file = {}".format(self.grid_coordinates_file))

        # for either
        logging.info("\thistorical_suitability_raster = {}".format(self.historical_suitability_raster))
        logging.info("\tstate_name = {}".format(self.state_name))
        logging.info("\tscenario = {}".format(self.scenario))
        logging.info("\toutput_directory = {}".format(self.output_directory))
        logging.info("\tkernel_distance_meters = {}".format(self.kernel_distance_meters))
        logging.info("\talpha_urban = {}".format(self.alpha_urban))
        logging.info("\talpha_rural = {}".format(self.alpha_rural))
        logging.info("\tbeta_urban = {}".format(self.beta_urban))
        logging.info("\tbeta_rural = {}".format(self.beta_rural))

        if self.projected_population_file is not None:
            logging.info("\tprojected_population_file = {}".format(self.projected_population_file))

        if self.urban_pop_proj_n is not None:
            logging.info("\turban_pop_proj_n = {}".format(self.urban_pop_proj_n))

        if self.rural_pop_proj_n is not None:
            logging.info("\trural_pop_proj_n = {}".format(self.rural_pop_proj_n))


    def calibrate(self):
        """Run the model."""

        # initialize model
        self.initialize()

        # start time
        tc = time.time()

        logging.info("Starting calibration.")

        # run calibration
        calib.calibration(self.config)

        # run time for calibration
        logging.info("Calibration completed in {} minutes.".format((time.time() - tc) / 60))

        self.close()

    def downscale(self):
        """Downscale rural and urban projection for all input years"""

        # initialize model
        self.initialize()

        # start time
        td = time.time()

        logging.info("Starting downscaling.")

        # process projected year
        ProcessStep(self)

        logging.info("Downscaling completed in {} minutes.".format((time.time() - td) / 60))

        # clean logger
        self.close()

    def close(self):
        """End model run and close log files"""

        logging.info("End time:  {}".format(time.strftime("%Y-%m-%d %H:%M:%S")))

        # Remove logging handlers
        self.close_logger()
