import os
import logging
import pkg_resources

import rasterio

import numpy as np
import pandas as pd

from rasterio.io import MemoryFile

import population_gravity.downscale_utilities as utils


def mask_raster(metadata, file_name, raster_object, source_window):
    """Subset the mosiac using a window of the source historical raster and write to file.

    :param file_name:                   Full path with file name and extension to the output file
    :type file_name:                    str

    :param raster_object:               Rasterio object

    :param source_window:               A rasterio.windows.Window object from the bounding box of the target for
                                        the full mosaic

    """
    # read only the window of the mosaic to an array
    masked = raster_object.read(1, window=source_window)

    # write the mosiac raster to file
    with rasterio.open(file_name, 'w', **metadata) as dest:
        dest.write(masked, indexes=1)

    return masked


def mask_raster_memory(metadata, raster_object, source_window):
    """Subset the mosiac using a window of the source historical raster and write to memory.

    :param raster_object:               Rasterio object

    :param source_window:               A rasterio.windows.Window object from the bounding box of the target for
                                        the full mosaic

    """
    # read only the window of the mosaic to an array
    masked = raster_object.read(1, window=source_window)

    # write output to memory
    with MemoryFile() as memfile:
        dataset = memfile.open(**metadata)
        dataset.write(masked, indexes=1)

        return dataset


def construct_file_name(state_name, prev_step, designation, scenario, output_directory, run_number='', suffix='',
                        extension='.tif'):
    """Construct output file name.

    :param state_name:                  Name of state
    :type state_name:                   str

    :param prev_step:                   Previous time step.
    :type prev_step:                    int, str

    :param designation:                 Either 'urban', 'rural', or 'total'
    :type designation:                  str

    :param suffix:                      String to append to end of file name before extension
    :type suffix:                       str

    :param extension:                   Raster extension without the dot; default 'tif'
    :type extension:                    str

    :returns:                           Full file path for an output raster

    """

    if (len(suffix) > 0):

        fle = f"{state_name}_1km_{scenario}_{designation}_{prev_step}{suffix}{extension}"

    else:

        if type(run_number) == int:
            delim = '_'
        else:
            delim = ''

        fle = f"{state_name}_1km_{scenario}_{designation}_{prev_step}{delim}{run_number}{suffix}{extension}"

    return os.path.join(output_directory, fle)


def get_state_neighbors(state_name):
    """Get all neighboring states and the target state from lookup file as a list"""

    df = pd.read_csv(pkg_resources.resource_filename('population_gravity', 'data/neighboring_states_150km.csv'))

    # get the actual state name from the near states because they are not lower case like what is being passed
    state_find = df.loc[(df['target_state'] == state_name) & (df['near_state'].str.lower() == state_name)].values[0][-1]

    # extract a list of all neighboring states including the target state
    state_list = df.loc[df['target_state'] == state_name]['near_state'].to_list()

    # ensure that the target state comes first to prevent any issue with the reverse painter's algorithm for merge
    state_list.insert(0, state_list.pop(state_list.index(state_find)))

    # make all lower case
    return [i.lower() for i in state_list]


def construct_file_list(prev_step, setting, state_name, template_raster, one_dimension_indices):
    """Construct a list of arrays from rasters or arrays.

    :param prev_step:                       int.  Previous time step; e.g., year
    :param setting:                         str.  Either 'urban' or 'rural'

    #TODO:  load prev years CSV files

    """

    out_list = []

    # Get all neighboring states including the target state as a list
    neighbors = get_state_neighbors(state_name)

    for i in neighbors:

        # check for either the 'tif', 'npy', or 'csv' files in the output directory, use 'tif' first
        tif = construct_file_name(i, prev_step, setting, extension='.tif')
        npy1d = construct_file_name(i, prev_step, setting, extension='_1d.npy')
        npy2d = construct_file_name(i, prev_step, setting, extension='_2d.npy')
        csv_gz = construct_file_name(i, prev_step, setting, extension='.csv.gz')

        if os.path.isfile(tif):
            logging.info(f"Using file '{tif}' for previous time step mosaic of neighboring states.")
            out_list.append(rasterio.open(tif))

        elif os.path.isfile(npy1d):
            logging.info(f"Using file '{npy1d}' for previous time step mosaic of neighboring states.")

            # load npy file to array
            array1d = np.load(npy1d)
            out_list.append(
                utils.array_to_raster_memory(template_raster, array1d, one_dimension_indices))

        elif os.path.isfile(npy2d):
            logging.info(f"Using file '{npy2d}' for previous time step mosaic of neighboring states.")

            # load npy file to array
            array2d = np.load(npy2d)
            out_list.append(utils.array2d_to_raster_memory(array2d, raster_profile=template_raster[4]))

        elif os.path.isfile(csv_gz):
            logging.info(f"Using file '{csv_gz}' for previous time step mosaic of neighboring states.")

            array1d = pd.read_csv(csv_gz, compression='gzip', sep=',')['value'].values
            out_list.append(
                utils.array_to_raster_memory(template_raster, array1d, one_dimension_indices))

        else:
            raise FileNotFoundError(f"No spatial file found for '{i}' for setting '{setting}' and year '{prev_step}'")

    return out_list


def mosaic_neighbors(yr, metadata, bbox):
    """asdf"""

    # build raster list for all neighboring raster outputs from the previous time step
    urban_file_list = construct_file_list(yr, 'urban')
    rural_file_list = construct_file_list(yr, 'rural')

    # build mosaics
    urban_mosaic = utils.mosaic_memory(urban_file_list, metadata.copy())
    rural_mosaic = utils.mosaic_memory(rural_file_list, metadata.copy())

    # create a rasterio.windows.Window object from the bounding box of the target for the full mosaic
    target_window = urban_mosaic.window(*bbox)

    # write the urban masked mosiac raster to memory
    urban_mask_file = mask_raster_memory(urban_mosaic, target_window)
    urban_mosaic.close()

    # write the urban masked mosiac raster to memory
    rural_mask_file = mask_raster_memory(rural_mosaic, target_window)
    rural_mosaic.close()

    return urban_mask_file, rural_mask_file
