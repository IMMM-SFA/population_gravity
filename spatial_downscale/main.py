import os
import time
import string

import spatial_downscale.downscale_utilities as pdm
import spatial_downscale.downscale_calibrate as calib
import spatial_downscale.downscale_projection as prj

from spatial_downscale.read_config import ReadConfig

def main(config_file):
    """Run the spatial population downscaling model.

    :param config_file:             Full path with file name and extension to the config.yml file

    """
    # read the YAML configuration file
    cfg = ReadConfig(config_file)

    # target SSP
    ssp_code = cfg.ssp_code

    # target state name
    region_code = cfg.region_code

    point_indices = cfg.point_indices
    data_rootdir = cfg.data_rootdir
    datadir_histdata = cfg.datadir_histdata
    datadir_future = cfg.datadir_future
    datadir_output = cfg.datadir_output
    compute_params = cfg.compute_params
    compute_proj = cfg.compute_proj
    
    # calibration inputs -- derived from user inputs
    urb_pop_fst_year = cfg.urb_pop_fst_year
    urb_pop_snd_year = cfg.urb_pop_snd_year
    rur_pop_fst_year = cfg.rur_pop_fst_year
    rur_pop_snd_year = cfg.rur_pop_snd_year

    # downscaling inputs -- derived from user inputs
    urb_pop_init_year = cfg.urb_pop_init_year
    rur_pop_init_year = cfg.rur_pop_init_year
    mask_raster = cfg.mask_raster
    point_indices = cfg.point_indices
    point_coors = cfg.point_coors
    ssp_dataFn = cfg.ssp_dataFn
    params_file = cfg.params_file
    
    #=========================================================================================
    sub_name = "<main> "

    print sub_name + format(time.strftime("%Y-%m-%d %H:%M:%S"))

    print sub_name + "Calibration inputs:"
    print sub_name + "   compute_params   = " + compute_params
    print sub_name + "   urb_pop_fst_year = " + urb_pop_fst_year
    print sub_name + "   urb_pop_snd_year = " + urb_pop_snd_year
    print sub_name + "   rur_pop_fst_year = " + rur_pop_fst_year
    print sub_name + "   rur_pop_snd_year = " + rur_pop_snd_year
    print sub_name + "   mask_raster      = " + mask_raster
    print sub_name + "   region_code      = " + region_code
    print sub_name + "   ssp_code         = " + ssp_code
    print sub_name + "   point_indices    = " + point_indices
    print sub_name + "   point_coors      = " + point_coors
    print sub_name + "   datadir_output   = " + datadir_output
    print sub_name
    print sub_name + "Downscaling inputs:"
    print sub_name + "   compute_proj      = " + compute_proj
    print sub_name + "   urb_pop_init_year = " + urb_pop_init_year
    print sub_name + "   rur_pop_init_year = " + rur_pop_init_year
    print sub_name + "   mask_raster       = " + mask_raster
    print sub_name + "   ssp_dataFn        = " + ssp_dataFn
    print sub_name + "   params_file       = " + params_file
    print sub_name + "   region_code       = " + region_code
    print sub_name + "   ssp_code          = " + ssp_code
    print sub_name + "   point_coors       = " + point_coors
    print sub_name + "   datadir_output    = " + datadir_output

    print sub_name


    print sub_name + format(time.strftime("%Y-%m-%d %H:%M:%S")) + " call calibration -------"

    if compute_params == "true" :
       print sub_name + "compute historical parameters"
       calib.calibration(urb_pop_fst_year, urb_pop_snd_year, rur_pop_fst_year, rur_pop_snd_year, mask_raster,
                   region_code, ssp_code, point_indices, point_coors, datadir_output)

    else:
       print sub_name + "Note: input downscaling parameters from a file "
       print sub_name + "    : input file = " + params_file
       print sub_name + "    : used to downscale to region code = " + region_code

    print sub_name + format(time.strftime("%Y-%m-%d %H:%M:%S")) + " call pop_projection -------"

    if compute_proj == "true" :
       print sub_name + "compute downscaled future projections"
       prj.pop_projection(urb_pop_init_year, rur_pop_init_year, mask_raster, ssp_dataFn,
                      params_file, region_code, ssp_code, point_coors, datadir_output)

       year = 2010
       while year <= 2050:
          print sub_name + "Projecting starting from year = " + str(year)
          urb_pop_init_year = datadir_output + region_code + "_1km_" + ssp_code + "_Urban_" + str(year) + ".tif"
          rur_pop_init_year = datadir_output + region_code + "_1km_" + ssp_code + "_Rural_" + str(year) + ".tif"
          prj.pop_projection(urb_pop_init_year, rur_pop_init_year, mask_raster, ssp_dataFn,
                         params_file, region_code, ssp_code, point_coors, datadir_output)
          year = year + 10
    else:
       print sub_name + "do not compute downscaled future projections"

    print sub_name + format(time.strftime("%Y-%m-%d %H:%M:%S")) + " finished ------------------"

    
