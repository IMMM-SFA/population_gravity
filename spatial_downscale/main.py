"""
Model interface for the Spatial Population Downscaling Model.

@author   Hamidreza Zoraghein (for model), Chris R. Vernon (of this file)
@email:   chris.vernon@pnnl.gov

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

"""

import argparse
import logging
import os
import sys
import time

import spatial_downscale.downscale_calibrate as calib

from spatial_downscale.read_config import ReadConfig
from spatial_downscale.process_step import ProcessStep


class Model:
    """Run the spatial population downscaling model.

    :param config_file:             Full path with file name and extension to the config.yml file

    datadir_histdata              : /Users/d3y010/projects/population/district_of_columbia/main_inputs
    ssp_data_directory            : /Users/d3y010/projects/population/district_of_columbia/ssp_projections
    ssp_code                      : SSP2
    region_code                   : district_of_columbia
    output_directory              : /Users/d3y010/projects/population/district_of_columbia/outputs
    calibration_parameters_file   : /Users/d3y010/projects/population/district_of_columbia/district_of_columbia_SSP2_calibration_parameters.csv

    # run settings
    compute_params      : False
    compute_proj        : True
    start_year          : 2000
    end_year            : 2020
    time_step           : 10

    """

    def __init__(self, config_file=None, datadir_histdata=None, ssp_data_directory=None, ssp_code=None,
                 region_code=None, output_directory=None, calibration_parameters_file=None, start_year=None,
                 end_year=None, time_step=None, alpha_urban=None, beta_urban=None, alpha_rural=None, beta_rural=None):

        # read the YAML configuration file
        self.cfg = ReadConfig(config_file=config_file,
                              datadir_histdata=datadir_histdata,
                              ssp_data_directory=ssp_data_directory,
                              ssp_code=ssp_code,
                              region_code=region_code,
                              output_directory=output_directory,
                              calibration_parameters_file=calibration_parameters_file,
                              start_year=start_year,
                              end_year=end_year,
                              time_step=time_step,
                              alpha_urban=alpha_urban,
                              beta_urban=beta_urban,
                              alpha_rural=alpha_rural,
                              beta_rural=beta_rural)

        # logfile path
        self.logfile = os.path.join(self.cfg.datadir_output, 'logfile_{}_{}.log'.format(self.cfg.ssp_code, self.cfg.region_code))

        self.timestep = None

    @staticmethod
    def make_dir(pth):
        """Create dir if not exists."""

        if not os.path.exists(pth):
            os.makedirs(pth)

    def init_log(self):
        """Initialize project-wide logger. The logger outputs to both stdout and a file."""

        log_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        log_level = logging.DEBUG

        logger = logging.getLogger()
        logger.setLevel(log_level)

        # logger console handler
        c_handler = logging.StreamHandler(sys.stdout)
        c_handler.setLevel(log_level)
        c_handler.setFormatter(log_format)
        logger.addHandler(c_handler)

        # logger file handler
        f_handler = logging.FileHandler(self.logfile)
        c_handler.setLevel(log_level)
        c_handler.setFormatter(log_format)
        logger.addHandler(f_handler)

    def initialize(self):
        """Setup model."""

        # build output directory first to store logfile and other outputs
        self.make_dir(self.cfg.datadir_output)

        # initialize logger
        self.init_log()

        logging.info("Start time:  {}".format(time.strftime("%Y-%m-%d %H:%M:%S")))

        # log run parameters
        logging.info("Input parameters:")
        logging.info("\turb_pop_fst_year = {}".format(self.cfg.urb_pop_fst_year))
        logging.info("\turb_pop_snd_year = {}".format(self.cfg.urb_pop_snd_year))
        logging.info("\trur_pop_fst_year = {}".format(self.cfg.rur_pop_fst_year))
        logging.info("\trur_pop_snd_year = {}".format(self.cfg.rur_pop_snd_year))
        logging.info("\tpoint_indices = {}".format(self.cfg.point_indices))

        # for projection
        logging.info("\turb_pop_init_year = {}".format(self.cfg.urb_pop_init_year))
        logging.info("\trur_pop_init_year = {}".format(self.cfg.rur_pop_init_year))
        logging.info("\tparams_file = {}".format(self.cfg.params_file))
        logging.info("\tpoint_coors = {}".format(self.cfg.point_coors))

        # for either
        logging.info("\tmask_raster = {}".format(self.cfg.mask_raster))
        logging.info("\tregion_code = {}".format(self.cfg.region_code))
        logging.info("\tssp_dataFn = {}".format(self.cfg.ssp_dataFn))
        logging.info("\tssp_code = {}".format(self.cfg.ssp_code))
        logging.info("\tdatadir_output = {}".format(self.cfg.datadir_output))

        # set up time step generator
        self.timestep = self.build_step_generator()

    def build_step_generator(self):
        """Build step generator."""

        for step in self.cfg.steps:

            yield ProcessStep(self.cfg, step)

    def advance_step(self):
        """Advance to next time step.

        Python 3 requires the use of `next()` to wrap the generator.

        """

        next(self.timestep)

    def calibrate(self):
        """Run the model."""

        # initialize model
        self.initialize()

        # start time
        tc = time.time()

        logging.info("Starting calibration.")

        # run calibration
        calib.calibration(self.cfg)

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

        # process all years
        for _ in self.cfg.steps:
            self.advance_step()

        logging.info("Downscaling completed in {} minutes.".format((time.time() - td) / 60))

        # clean logger
        self.close()

    def close(self):
        """End model run and close log files"""

        logging.info("End time:  {}".format(time.strftime("%Y-%m-%d %H:%M:%S")))

        # Remove logging handlers
        logger = logging.getLogger()

        for handler in logger.handlers[:]:
            handler.close()
            logger.removeHandler(handler)

        logging.shutdown()
