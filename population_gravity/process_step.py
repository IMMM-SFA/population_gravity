"""
Process population projection by step for the popultation_gravity model

@author   Chris R. Vernon
@email:   chris.vernon@pnnl.gov

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

"""

import logging
import time

from population_gravity.downscale_projection import pop_projection


class ProcessStep:
    """Process population downscaling for a single time step

    :param cfg:                 Configuration file object
    :param yr:                  Target year (YYYY) for a time step

    """

    def __init__(self, cfg):

        # start time
        td = time.time()

        self.cfg = cfg

        logging.info("Downscaling year:  {}".format(self.cfg.projection_year))

        # run downscaling
        urban_output, rural_output = pop_projection(cfg=self.cfg,
                                                    urban_raster=self.cfg.base_urban_pop_raster,
                                                    rural_raster=self.cfg.base_rural_pop_raster,
                                                    alpha_urban=self.cfg.alpha_urban,
                                                    beta_urban=self.cfg.beta_urban,
                                                    alpha_rural=self.cfg.alpha_rural,
                                                    beta_rural=self.cfg.beta_rural,
                                                    rural_pop_proj_n=self.cfg.rural_pop_proj_n,
                                                    urban_pop_proj_n=self.cfg.urban_pop_proj_n,
                                                    yr=self.cfg.projection_year,
                                                    cut_off_meters=self.cfg.kernel_distance_meters)

        logging.info(f"Downscaling for year {self.cfg.projection_year} completed in {(time.time() - td) / 60} minutes.")
