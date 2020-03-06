import logging
import os
import time

from population_gravity.downscale_projection import pop_projection


class ProcessStep:
    """Process population downscaling for a single time step

    :param cfg:                 Configuration file object
    :param yr:                  Target year (YYYY) for a time step

    """

    def __init__(self, cfg, yr):

        # start time
        td = time.time()

        logging.info("\tDownscaling year:  {}".format(yr))

        self.cfg = cfg
        self.yr = yr

        if self.yr == self.cfg.start_year:

            # run downscaling
            pop_projection(self.cfg, self.cfg.urb_pop_init_year, self.cfg.rur_pop_init_year)

        else:

            urban_raster = os.path.join(self.cfg.datadir_output,
                                        "{}_1km_{}_Urban_{}.tif".format(self.cfg.region_code, self.cfg.ssp_code, self.yr))
            rural_raster = os.path.join(self.cfg.datadir_output,
                                        "{}_1km_{}_Rural_{}.tif".format(self.cfg.region_code, self.cfg.ssp_code, self.yr))

            pop_projection(self.cfg, urban_raster, rural_raster)

        logging.info("Downscaling for year {} completed in {} minutes.".format(yr, (time.time() - td) / 60))
