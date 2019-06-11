import os
import time
import string

os.chdir(r"C:\Users\Hamidreza.Zoraghein\Google Drive\Population_Downscaling_1KM\Code")
#import pop_downscaling_module_190201 as pdm
import pop_downscaling_module_parallel as pdm

#User inputs -- from stdin "namelist" file 
ssp_code            = "unset"
region_code         = "unset"
point_indices       = "unset"
data_rootdir        = " "      # optional functionality
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

def calibration(urb_pop_fst_year, urb_pop_snd_year, rur_pop_fst_year, rur_pop_snd_year, mask_raster,
                region_code, ssp_code, point_indices, point_coors, datadir_output):

    #Import the required modules
    import pandas as pd
    import numpy as np
    import scipy.optimize
    import simplejson
    import pickle
    import csv
    
    sub_name = "<calibration> "
    print sub_name + "enter"
    
    #Define local variables    
    all_rasters     = {} ##Dictionary storing initial urban and rural population grids
    rur_pop_files   = [] ##List containing rural population grids
    urb_pop_files   = [] ##List containing urban population grids
    population_1st  = {} ##Dictionary containing population of each point in year 1
    population_2nd  = {} ##Dictionary containing population of each point in year 2
    parameters_dict = {} ##Dictionary storing urban/rural calibration parameters
    
    ##Rural
    rur_pop_files.append(rur_pop_fst_year)
    rur_pop_files.append(rur_pop_snd_year)
    
    ##Urban
    urb_pop_files.append(urb_pop_fst_year)
    urb_pop_files.append(urb_pop_snd_year)
    
    ##Urban and rural
    all_rasters["Rural"] = rur_pop_files
    all_rasters["Urban"] = urb_pop_files
   
    #Populate the array containing mask values 
    points_mask = pdm.raster_to_array(mask_raster)
    
    #Read historical urban and rural population grids into arrrays
    for setting in all_rasters:
        ##Create the dictionary containing population of each point in year 1
        population_1st[setting] = pdm.raster_to_array(all_rasters[setting][0])        
        ##Create the dictionary containing population of each point in year 2
        population_2nd[setting] = pdm.raster_to_array(all_rasters[setting][1])
    
    #Create an array containing total population values in the first historical year    
    total_population_1st = population_1st["Rural"] + population_1st["Urban"] 
    
    #A csv file to write the parameters into
    out_cal  = datadir_output + region_code + "_" + ssp_code + "_Params.csv"
    
    #All indices
    all_indices = pdm.all_index_retriever(urb_pop_snd_year, ["row", "column"])
    
    #Read indices of points that fall within the state boundary
    with open(point_indices, 'r') as r:
        within_indices = simplejson.load(r)
        
    #Calculate a distance matrix that serves as a template
    cut_off_meters = 100000
    dist_matrix = pdm.dist_matrix_calculator(within_indices[0], cut_off_meters, all_indices, point_coors)
        
    #Parameter calculation for rural and urban    
    for setting in all_rasters:                   
        ###Initial alpha values
        a_lower = -2.0
        a_upper = 2.0 
        ###Initial beta values
        b_lower = -2.0
        b_upper = 4.0
                
        rranges = ((a_lower, a_upper), (b_lower, b_upper))
        
        ###Parameters to be used in optimization        
        a_list = np.linspace(a_lower, a_upper, 8)
        b_list = np.linspace(b_lower, b_upper, 10)
        params = (setting, population_1st, population_2nd, total_population_1st,
                  points_mask, dist_matrix, within_indices)
        
        ###Initialize the dataframe that will hold values of the brute force
        fst_results = pd.DataFrame(data={"a"        : np.repeat(a_list, 10).astype(np.float32),
                                         "b"        : np.tile(b_list, 8).astype(np.float32),
                                         "estimate" : np.full(80, np.nan, dtype=np.float64)
                       }
                 )
        
        ###Run brute force to calculate optimization per grid point
        i = 0
        for a in a_list:
            for b in b_list:
                fst_results.loc[(fst_results["a"] == a) & (fst_results["b"] == b),
                                "estimate"] = pdm.pop_min_function((a, b), *params)
                print "{}th iteratin done".format(i)
                i+=1
        
        ###Pickle the current optimization file
        pickle_file = os.path.join(datadir_output, setting)
        with open(pickle_file, "wb") as pickle_params:
            pickle.dump(fst_results, pickle_params)
        
        ###Use the point with the minimum value as an initial guess for the second optimizer
        (a0, b0) = fst_results.loc[fst_results["estimate"].idxmin(), ["a", "b"]] 
        
        ###Final optimization
        parameters = scipy.optimize.minimize(pdm.pop_min_function, x0 = (a0, b0), args = params, method = "SLSQP",
                                             tol = 0.001, options = {"disp" : True}, bounds = rranges)
        
        print sub_name + "Optimization Outcomes for " + setting
        print parameters["x"]
        print parameters["fun"]
        parameters_dict[setting] = parameters["x"]

        
    ###Write the parameters to the designated csv file
    print sub_name + "writing to file: " + out_cal
    with open(out_cal, 'wb') as out_csv:
        writer = csv.writer(out_csv)
        writer.writerow(["Region", "SSP", "Alpha_Rural", "Beta_Rural", "Alpha_Urban", "Beta_Urban", "Comments"])
        writer.writerow([region_code, ssp_code, parameters_dict["Rural"][0], parameters_dict["Rural"][1],
                         parameters_dict["Urban"][0], parameters_dict["Urban"][1], ""])
        
 #=========================================================================================       

def pop_projection(urb_pop_init_year, rur_pop_init_year, mask_raster, ssp_dataFn, 
                  params_file, region_code, ssp_code, point_coors, datadir_output):

    #Import the required modules
    import numpy as np
    import pandas as pd
    import multiprocessing
    from pathos.multiprocessing import ProcessingPool as Pool
    from collections import deque
    import simplejson
    
    sub_name = "<pop_projection> "
    print sub_name + "enter"
    
    #Define local variables
    current_timestep = [int(s) for s in urb_pop_init_year.split("/")[-1][:-4].split("_") if s.isdigit()][0] + 10
    time_one_data    = {} ##Dictionary storing the base year population grids  
    population_1st   = {} ##Dictionary storing the base year population arrays 
    time_one_data['Rural'] = rur_pop_init_year ##Rural
    time_one_data['Urban'] = urb_pop_init_year ##Urban
    final_arrays = {} ##Dictionary containing final projected arrays
    final_raster =  datadir_output + region_code + "_1km_" + ssp_code + "_Total_" + str(current_timestep)  + ".tif"
    
    ##Populate the array containing mask values for all points
    points_mask = pdm.raster_to_array(mask_raster)
    
    #All indices
    all_indices = pdm.all_index_retriever(urb_pop_snd_year, ["row", "column"])
    
    #Read indices of points that fall within the state boundary
    with open(point_indices, 'r') as r:
        within_indices = simplejson.load(r)
        
    #Calculate a distance matrix that serves as a template
    cut_off_meters = 100000
    dist_matrix    = pdm.dist_matrix_calculator(within_indices[0], cut_off_meters, all_indices, point_coors)
    
    #Read historical urban and rural population grids into arrrays
    for setting in time_one_data:
        ##Create the dictionary containing population of each point in year 1
        population_1st[setting] = pdm.raster_to_array(time_one_data[setting]) 
    
    #Create an array containing total population values in the first historical year    
    total_population_1st = population_1st["Rural"] + population_1st["Urban"] 
      
    #Number of columns of population grid required to derive linear indices  
    ind_diffs = dist_matrix["ind_diff"].values
    
    #Distances between current point and its close points
    ini_dist = dist_matrix["dis"].values/1000.0
    
    #Downscale population projection    
    for setting in time_one_data:
        
        ##Define required parameters within the loop
        pop_estimates         = np.zeros(len(within_indices)) ##Population estimates for the second year 
        suitability_estimates = deque() #Suitability estimates in the second year
        
        ##Output raster
        output = datadir_output + region_code + "_1km_" + ssp_code + "_" + setting + "_" + str(current_timestep)  + ".tif"
                        
        ##Calculate aggregate urban/rural population at time 1
        pop_first_year = population_1st[setting][within_indices]
        pop_t1         = pop_first_year.sum()
    
        ##Load the SSP file to retrieve the aggregated projected population at time 2 for downscaling
        ssp_data = pd.read_csv(ssp_dataFn)
        if(setting == "Urban"):
            pop_t2 = ssp_data.loc[(ssp_data["Year"] == current_timestep) & (ssp_data["Scenario"] == ssp_code),
                                      "UrbanPop"].values[0]
        else:
            pop_t2 = ssp_data.loc[(ssp_data["Year"] == current_timestep) & (ssp_data["Scenario"] == ssp_code),
                                      "RuralPop"].values[0] 
        
        ##Population change between years 1 and 2
        pop_change = pop_t2 - pop_t1
        
        if pop_change < 0:
            negative_mod = 1
        else:
            negative_mod = 0
        
        #Extract the alpha and beta values from the calibration files
        calib_params = pd.read_csv(params_file)
        if (setting == "Urban"):
            calib_params = calib_params.loc[calib_params.SSP == ssp_code, ["Alpha_Urban", "Beta_Urban"]].values
        else:
            calib_params = calib_params.loc[calib_params.SSP == ssp_code, ["Alpha_Rural", "Beta_Rural"]].values
        
        a = calib_params[0][0]
        b = calib_params[0][1]
        
        dist                 = -b * ini_dist
        exp_xx_inv_beta_dist = np.exp(dist)
       
        #Initialize the parallelization
        pool = Pool(processes = multiprocessing.cpu_count())
        
        #Provide the inputs for the parallelized function
        parallel_elements = deque([(i, ind_diffs, total_population_1st, a, exp_xx_inv_beta_dist)
                                    for i in within_indices])
        
        #Derive suitability estimates
        suitability_estimates = pool.map(pdm.suitability_estimator, parallel_elements)
        
        #Change suitability estimates to a numpy array
        suitability_estimates = np.array(suitability_estimates)
        
        # Exract only the necessary mask values that fall within the state boundary
        cur_points_mask = points_mask[within_indices]
    
        ##In case of population decline, suitability estimates are reciprocated 
        if negative_mod:
            
            # find those whose mask is 0 but have population, they should decline anyway 
            mask_zero    = np.where(cur_points_mask == 0)[0]
            pop_non_zero = np.where(pop_first_year != 0)[0]
            
            # Those cells with mask value of zero and population are the intersection of the two above arrays
            pop_mask = np.intersect1d(mask_zero, pop_non_zero, assume_unique=True)
            
            # Change the mask value of the above cells to 1 so that they also lose population
            cur_points_mask[pop_mask] = cur_points_mask.mean()
            
            #Adjust suitability values by applying mask values
            suitability_estimates = cur_points_mask * suitability_estimates
            
            # Inverse current mask values for a better reflection of population decline
            suitability_estimates[suitability_estimates != 0] = 1.0/suitability_estimates[suitability_estimates != 0]
        
        else:
            ##Adjust suitability values by applying its mask values
            suitability_estimates = cur_points_mask * suitability_estimates
    
        ##Total suitability for the whole area, which is the summation of all individual suitability values
        tot_suitability = suitability_estimates.sum()

        ##Final population estimates if nagative mode is off
        pop_estimates = suitability_estimates/tot_suitability * pop_change + pop_first_year
       
        #Adjust the projection so that no cell can have less than 0 individuals.    
        if negative_mod:
            ##To make sure that there is no negative population
            while any(pop < 0 for pop in pop_estimates):
                new_tot_suitability = 0 ##Total suitability calculated over points with positive population
                extra_pop_mod = 0 ##For adjusting negative population values
                
                ##Treating negative population values
                extra_pop_mod = abs(pop_estimates[pop_estimates < 0].sum())
                pop_estimates[pop_estimates < 0] = 0
    
                ##Calculate the new total suitability value based on points with positive projected population                
                new_tot_suitability = suitability_estimates[pop_estimates > 0].sum()       
                
                ##Adjust non-negative population values to maintain the total aggregated population 
                pop_estimates[pop_estimates > 0] = pop_estimates[pop_estimates > 0] - (suitability_estimates[pop_estimates > 0]/new_tot_suitability) * extra_pop_mod
        
        #Save the projection array of the current setting
        final_arrays[setting] = pop_estimates
        
        #Save the final urban/rural raster
        print sub_name + "output downscaled data: " + output
        pdm.array_to_raster(mask_raster, pop_estimates, within_indices, output)
    
    #Calculate the total population array    
    total_array = final_arrays["Rural"] + final_arrays["Urban"]
    #Save the final total raster
    print sub_name + "output downscaled data: " + final_raster
    pdm.array_to_raster(mask_raster, total_array, within_indices, final_raster)
    
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
   calibration(urb_pop_fst_year, urb_pop_snd_year, rur_pop_fst_year, rur_pop_snd_year, mask_raster,
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

    
