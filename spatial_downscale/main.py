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
import spatial_downscale.downscale_projection as prj

from spatial_downscale.read_config import ReadConfig


class Downscale:
    """Run the spatial population downscaling model.

    :param config_file:             Full path with file name and extension to the config.yml file

    """

    def __init__(self, config_file):

        # read the YAML configuration file
        self.cfg = ReadConfig(config_file)

        # logfile path
        self.logfile = os.path.join(self.cfg.datadir_output, 'logfile_{}_{}.log'.format(self.cfg.ssp_code, self.cfg.region_code))

    @staticmethod
    def make_dir(pth):
        """Create dir if not exists."""

        if not os.path.exists(pth):
            os.makedirs(pth)

    def init_log(self):
        """Initialize project-wide logger. The logger outputs to both stdout and a file."""

        log_format = logging.Formatter('%(levelname)s: %(message)s')
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

    def stage(self):
        """Set up run."""

        # build output directory first to store logfile and other outputs
        self.make_dir(self.cfg.datadir_output)

        # initialize logger
        self.init_log()

        logging.info("Start time:  {}".format(time.strftime("%Y-%m-%d %H:%M:%S")))

        # log run parameters
        if self.cfg.compute_params:
            logging.info("Calibration inputs:")
            logging.info("\tcompute_params = {}".format(self.cfg.compute_params))
            logging.info("\turb_pop_fst_year = {}".format(self.cfg.urb_pop_fst_year))
            logging.info("\turb_pop_snd_year = {}".format(self.cfg.urb_pop_snd_year))
            logging.info("\trur_pop_fst_year = {}".format(self.cfg.rur_pop_fst_year))
            logging.info("\trur_pop_snd_year = {}".format(self.cfg.rur_pop_snd_year))
            logging.info("\tmask_raster = {}".format(self.cfg.mask_raster))
            logging.info("\tregion_code = {}".format(self.cfg.region_code))
            logging.info("\tssp_code = {}".format(self.cfg.ssp_code))
            logging.info("\tpoint_indices = {}".format(self.cfg.point_indices))
            logging.info("\tpoint_coors = {}".format(self.cfg.point_coors))
            logging.info("\tdatadir_output = {}".format(self.cfg.datadir_output))

        else:
            logging.info("\tCalibration not selected in config file to execute.")

        if self.cfg.compute_proj:
            logging.info("Downscaling inputs:")
            logging.info("\tcompute_proj = {}".format(self.cfg.compute_proj))
            logging.info("\turb_pop_init_year = {}".format(self.cfg.urb_pop_init_year))
            logging.info("\trur_pop_init_year = {}".format(self.cfg.rur_pop_init_year))
            logging.info("\tmask_raster = {}".format(self.cfg.mask_raster))
            logging.info("\tregion_code = {}".format(self.cfg.region_code))
            logging.info("\tssp_dataFn = {}".format(self.cfg.ssp_dataFn))

            # TODO: check for params_file?
            logging.info("\tparams_file = {}".format(self.cfg.params_file))
            logging.info("\tssp_code = {}".format(self.cfg.ssp_code))
            logging.info("\tpoint_coors = {}".format(self.cfg.point_coors))
            logging.info("\tdatadir_output = {}".format(self.cfg.datadir_output))

        else:
            logging.info("\tDownscaling not selecting in config file to execute.")

    def execute(self):
        """Run the model."""

        self.stage()

        if self.cfg.compute_params:

            # start time
            tc = time.time()

            logging.info("Starting calibration.")

            # run calibration
            calib.calibration(self.cfg)

            # run time for calibration
            logging.info("Calibration completed in {} minutes.".format((time.time() - tc) / 60))

        else:
            logging.info("Downscaling will use existing parameter file instead of newly calibrated parameters.")

        if self.cfg.compute_proj:

            # start time
            td = time.time()

            logging.info("Starting downscaling.")

            # iterate through all years
            for yr in range(self.cfg.start_year, self.cfg.end_year + self.cfg.time_step, self.cfg.time_step):

                logging.info("\tDownscaling year:  {}".format(yr))

                if yr == self.cfg.start_year:
                    urban_raster = self.cfg.urb_pop_init_year
                    rural_raster = self.cfg.rur_pop_init_year

                else:
                    urban_raster = os.path.join(self.cfg.datadir_output, "{}_1km_{}_Urban_{}.tif".format(self.cfg.region_code, self.cfg.ssp_code, self.yr))
                    rural_raster = os.path.join(self.cfg.datadir_output, "{}_1km_{}_Rural_{}.tif".format(self.cfg.region_code, self.cfg.ssp_code, self.yr))

                # run downscaling
                prj.pop_projection(self.cfg, urban_raster, rural_raster)

            logging.info("Downscaling completed in {} minutes.".format((time.time() - td) / 60))

        logging.info("End time:  {}".format(time.strftime("%Y-%m-%d %H:%M:%S")))

    @staticmethod
    def cleanup():
        """Close log files."""

        # Remove logging handlers
        logger = logging.getLogger()

        for handler in logger.handlers[:]:
            logger.removeHandler(handler)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', type=str, help='Full path with file name to YAML configuration file.')
    args = parser.parse_args()

    run = Downscale(args.config_file)
    run.execute()
    del run
