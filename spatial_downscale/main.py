import os
import time
import string

os.chdir(r"C:\Users\Hamidreza.Zoraghein\Google Drive\Population_Downscaling_1KM\Code")
#import pop_downscaling_module_190201 as pdm
import spatial_downscale.downscale_utilities as pdm
import spatial_downscale.downscale_calibrate as calib

#User inputs -- from stdin "namelist" file 
ssp_code            = "unset"
region_code         = "unset"
point_indices       = "unset"
data_rootdir        = " " # optional functionality
datadir_histdata    = "unset/"
datadir_future      = "unset/"
datadir_output      = "unset/"
compute_params      = "true"
compute_proj        = "false"

#What is the current state?
state = "Pennsylvania"

input_fn = r"C:\Users\Hamidreza.Zoraghein\Google Drive\Population_Downscaling_1KM\Template_Input.txt"
with open(input_fn) as input_list:
   for line in input_list:
      keyword, value = line.partition("=")[::2]
      if (keyword.strip() == "EOF"):
         break
      if keyword.strip() == "ssp_code" :
         ssp_code = value.strip()
      if keyword.strip() == "region_code" :
         region_code = value.strip()
      if keyword.strip() == "data_rootdir" :
         data_rootdir = value.strip()
      if keyword.strip() == "datadir_histdata" :
         datadir_histdata = value.strip()
      if keyword.strip() == "datadir_future" :
         datadir_future = value.strip()
      if keyword.strip() == "datadir_output" :
         datadir_output = value.strip()
      if keyword.strip() == "compute_params" :
         compute_params = value.strip()
      if keyword.strip() == "compute_proj" :
         compute_proj  = value.strip()

#Replace the template with the current state
datadir_histdata = string.replace(datadir_histdata, "Template", state)
datadir_future   = string.replace(datadir_future,    "Template", state)
datadir_output   = string.replace(datadir_output,    "Template", state)
region_code      = state

# replace "ROOTDIR" with data_rootdir
if data_rootdir != "unset" : 
   datadir_histdata = string.replace(datadir_histdata,   "ROOTDIR", data_rootdir)
   datadir_future   = string.replace(datadir_future,     "ROOTDIR", data_rootdir)
   datadir_output   = string.replace(datadir_output,     "ROOTDIR", data_rootdir)

try:
    os.makedirs(datadir_output)
except:
    print "The folder already exists!"
    
#Calibration inputs -- derived from user inputs
urb_pop_fst_year  = datadir_histdata + region_code + "_" + "Urban_2000_1km.tif"
urb_pop_snd_year  = datadir_histdata + region_code + "_" + "Urban_2010_1km.tif"
rur_pop_fst_year  = datadir_histdata + region_code + "_" + "Rural_2000_1km.tif"
rur_pop_snd_year  = datadir_histdata + region_code + "_" + "Rural_2010_1km.tif"
mask_raster       = datadir_histdata + region_code + "_" + "Mask_short_term.tif"
point_indices     = datadir_histdata + region_code + "_" + "Within_Indices.txt"
point_coors       = datadir_histdata + region_code + "_" + "Coors.csv"


#Downscaling inputs -- derived from user inputs
urb_pop_init_year = datadir_histdata + region_code + "_" + "Urban_2000_1km.tif"
rur_pop_init_year = datadir_histdata + region_code + "_" + "Rural_2000_1km.tif"
mask_raster       = datadir_histdata + region_code + "_" + "Mask_short_term.tif"
point_indices     = datadir_histdata + region_code + "_" + "Within_Indices.txt"
point_coors       = datadir_histdata + region_code + "_" + "Coors.csv"
ssp_dataFn        = datadir_future   + region_code + "_" + ssp_code + "_popproj.csv"

if compute_params == "true" : # compute params & put them in this file
    params_file   = datadir_output   + region_code + "_" + ssp_code + "_Params.csv"
    



    
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
   pop_projection(urb_pop_init_year, rur_pop_init_year, mask_raster, ssp_dataFn, 
                  params_file, region_code, ssp_code, point_coors, datadir_output)

   year = 2010
   while year <= 2050:
      print sub_name + "Projecting starting from year = " + str(year)
      urb_pop_init_year = datadir_output + region_code + "_1km_" + ssp_code + "_Urban_" + str(year) + ".tif"
      rur_pop_init_year = datadir_output + region_code + "_1km_" + ssp_code + "_Rural_" + str(year) + ".tif"
      pop_projection(urb_pop_init_year, rur_pop_init_year, mask_raster, ssp_dataFn, 
                     params_file, region_code, ssp_code, point_coors, datadir_output)
      year = year + 10
else:
   print sub_name + "do not compute downscaled future projections"

print sub_name + format(time.strftime("%Y-%m-%d %H:%M:%S")) + " finished ------------------"       

    
