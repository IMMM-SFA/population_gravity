"""Sensitivity protocol for the popultation_gravity model

@author   Chris R. Vernon
@email:   chris.vernon@pnnl.gov

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

"""


import multiprocessing
import os
import pickle

import numpy as np
import pandas as pd
from SALib.sample import latin
from SALib.analyze import delta
from pathos.multiprocessing import ProcessingPool as Pool

from population_gravity.main import Model
from population_gravity.read_config import ReadConfig


class Problem:
    """Generate a validated problem set for SALib."""

    # parameter value limits
    LOWER_ALPHA_LIMIT = -2.0
    UPPER_ALPHA_LIMIT = 2.0
    LOWER_BETA_LIMIT = -2.0
    UPPER_BETA_LIMIT = 2.0

    # kernel distance in meters
    LOWER_KD_LIMIT = 100
    UPPER_KD_LIMIT = 100000

    # variable names
    ALPHA_URBAN_NAME = 'alpha_urban'
    ALPHA_RURAL_NAME = 'alpha_rural'
    BETA_URBAN_NAME = 'beta_urban'
    BETA_RURAL_NAME = 'beta_rural'
    KD_NAME = 'kernel_distance_meters'

    def __init__(self, alpha_urban_bounds=None, alpha_rural_bounds=None, beta_urban_bounds=None, beta_rural_bounds=None,
                 kernel_distance_meters_bounds=None, problem_dict_outfile=None):

        self._alpha_urban_bounds = alpha_urban_bounds
        self._alpha_rural_bounds = alpha_rural_bounds
        self._beta_urban_bounds = beta_urban_bounds
        self._beta_rural_bounds = beta_rural_bounds
        self._kernel_distance_meters_bounds = kernel_distance_meters_bounds
        self._problem_dict_outfile = problem_dict_outfile

    @property
    def alpha_urban_bounds(self):
        """Ensure values are within acceptable bounds."""

        if self._alpha_urban_bounds is not None:

            # ensure lower is not >= upper
            valid_lower, valid_upper = self.validate_min_max(self._alpha_urban_bounds, self.ALPHA_URBAN_NAME)

            return self.validate_limits(valid_upper, valid_lower, self.UPPER_ALPHA_LIMIT, self.LOWER_ALPHA_LIMIT, self.ALPHA_URBAN_NAME)

        else:
            return self._alpha_urban_bounds

    @property
    def alpha_rural_bounds(self):
        """Ensure values are within acceptable bounds."""

        if self._alpha_rural_bounds is not None:

            # ensure lower is not >= upper
            valid_lower, valid_upper = self.validate_min_max(self._alpha_rural_bounds, self.ALPHA_RURAL_NAME)

            return self.validate_limits(valid_upper, valid_lower, self.UPPER_ALPHA_LIMIT, self.LOWER_ALPHA_LIMIT, self.ALPHA_RURAL_NAME)

        else:
            return self._alpha_rural_bounds

    @property
    def beta_urban_bounds(self):
        """Ensure values are within acceptable bounds."""

        if self._beta_urban_bounds is not None:

            # ensure lower is not >= upper
            valid_lower, valid_upper = self.validate_min_max(self._beta_urban_bounds, self.BETA_URBAN_NAME)

            return self.validate_limits(valid_upper, valid_lower, self.UPPER_BETA_LIMIT, self.LOWER_BETA_LIMIT, self.BETA_URBAN_NAME)

        else:
            return self._beta_urban_bounds

    @property
    def beta_rural_bounds(self):
        """Ensure values are within acceptable bounds."""

        if self._beta_rural_bounds is not None:

            # ensure lower is not >= upper
            valid_lower, valid_upper = self.validate_min_max(self._beta_rural_bounds, self.BETA_RURAL_NAME)

            return self.validate_limits(valid_upper, valid_lower, self.UPPER_BETA_LIMIT, self.LOWER_BETA_LIMIT, self.BETA_RURAL_NAME)

        else:
            return self._beta_rural_bounds

    @property
    def kernel_distance_meters_bounds(self):
        """Ensure values are withing acceptable bounds."""

        if self._kernel_distance_meters_bounds is not None:

            # ensure lower is not >= upper
            valid_lower, valid_upper = self.validate_min_max(self._kernel_distance_meters_bounds, self.KD_NAME)

            return self.validate_limits(valid_upper, valid_lower, self.UPPER_KD_LIMIT, self.LOWER_KD_LIMIT, self.KD_NAME)

        else:
            return self._kernel_distance_meters_bounds

    def build_problem_dict(self):
        """Construct a problem dictionary for SALib.

        :return:                    dictionary of problem set.

        """

        names = []
        bounds = []

        if self.alpha_urban_bounds is not None:
            names.append(self.ALPHA_URBAN_NAME)
            bounds.append(self.alpha_urban_bounds)

        if self.alpha_rural_bounds is not None:
            names.append(self.ALPHA_RURAL_NAME)
            bounds.append(self.alpha_rural_bounds)

        if self.beta_urban_bounds is not None:
            names.append(self.BETA_URBAN_NAME)
            bounds.append(self.beta_urban_bounds)

        if self.beta_rural_bounds is not None:
            names.append(self.BETA_RURAL_NAME)
            bounds.append(self.beta_rural_bounds)

        if self.kernel_distance_meters_bounds is not None:
            names.append(self.KD_NAME)
            bounds.append(self.kernel_distance_meters_bounds)

        d = {'num_vars': len(names), 'names': names, 'bounds': bounds}

        # write dict to file
        self.problem_dict_outfile(d)

        return d

    def problem_dict_outfile(self, problem_dict):
        """Optionally write problem_dict to pickled file."""

        if self._problem_dict_outfile is not None:
            with open(self._problem_dict_outfile, 'wb') as out:
                pickle.dump(problem_dict, out)

    @staticmethod
    def validate_limits(upper_value, lower_value, upper_limit, lower_limit, variable):
        """Ensure the values for bounds are not outside of the acceptable limits.

        :param upper_value:                     upper value for a parameter
        :type:                                  int; float

        :param lower_value:                     lower value for a parameter
        :type:                                  int; float

        :param upper_limit:                     upper acceptable value for a parameter
        :type:                                  int; float

        :param lower_limit:                     lower acceptable value for a parameter
        :type:                                  int; float

        """

        if lower_value < lower_limit:
            raise ValueError(f"Lower limit '{lower_value}'' for '{variable}' cannot be less than '{lower_limit}'")

        elif upper_value > upper_limit:
            raise ValueError(f"Upper limit '{upper_value}'' for '{variable}' cannot be less than '{upper_limit}'")

        return lower_value, upper_value

    @staticmethod
    def validate_min_max(bounds, variable):
        """Check to make sure min is not >= max

        :param bounds:                      [min_value, max_value]
        :type bounds:                       list

        :param variable:                    variable name
        :type variable:                     str

        """

        min_bound, max_bound = bounds

        if min_bound >= max_bound:
            raise ValueError(f"Minimum bound '{min_bound}' for '{variable}' is >= '{max_bound}'")
        else:
            return bounds


class Lhs(Problem):
    """Generate model run outputs using Latin Hypercube Sampling."""

    def __init__(self, alpha_urban_bounds=None, alpha_rural_bounds=None, beta_urban_bounds=None, beta_rural_bounds=None,
                 kernel_distance_meters_bounds=None, n_samples=None, sample_outfile=None, problem_dict_outfile=None):

        # initialize problem set
        super(Lhs, self).__init__(alpha_urban_bounds, alpha_rural_bounds, beta_urban_bounds, beta_rural_bounds,
                                  kernel_distance_meters_bounds, problem_dict_outfile)

        self._n_samples = n_samples
        self._sample_outfile = sample_outfile

        self.problem_dict = self.build_problem_dict()
        self.sample = self.generate_sample()

    @property
    def n_samples(self):
        """Validate the number of samples."""

        try:
            return int(self._n_samples)

        except TypeError:
            raise TypeError(f"'n_samples' value '{self._n_samples} must be an integer.")

    def generate_sample(self):
        """Generate a Latin Hypercube sample."""

        sample = latin.sample(self.problem_dict, self.n_samples)

        # write sample to file
        self.sample_outfile(sample)

        return sample

    def sample_outfile(self, sample):
        """Optionally save sample as numpy array."""

        if self._sample_outfile is not None:
            np.save(self._sample_outfile, sample)

        else:
            return self._sample_outfile


class BatchModelRun(ReadConfig):
    """Conduct a batch run of the model in parallel """

    # variable names
    ALPHA_URBAN_NAME = 'alpha_urban'
    ALPHA_RURAL_NAME = 'alpha_rural'
    BETA_URBAN_NAME = 'beta_urban'
    BETA_RURAL_NAME = 'beta_rural'
    KD_NAME = 'kernel_distance_meters'

    def __init__(self, grid_coordinates_file=None, historical_suitability_raster=None, historical_rural_pop_raster=None,
                 historical_urban_pop_raster=None, projected_population_file=None, one_dimension_indices_file=None,
                 output_directory=None, alpha_urban=None, beta_urban=None, alpha_rural=None, beta_rural=None,
                 scenario=None, state_name=None, historic_base_year=None, projection_start_year=None,
                 projection_end_year=None, time_step=None, kernel_distance_meters=None, write_raster=True,
                 write_csv=False, write_array1d=False, write_array2d=False, run_number='', write_logfile=True,
                 compress_csv=True, output_total=True, sample=None, problem_dict=None):

        super(BatchModelRun, self).__init__(grid_coordinates_file=grid_coordinates_file,
                                         historical_rural_pop_raster=historical_rural_pop_raster,
                                         historical_urban_pop_raster=historical_urban_pop_raster,
                                         historical_suitability_raster=historical_suitability_raster,
                                         projected_population_file=projected_population_file,
                                         one_dimension_indices_file=one_dimension_indices_file,
                                         output_directory=output_directory,
                                         alpha_urban=alpha_urban,
                                         alpha_rural=alpha_rural,
                                         beta_urban=beta_urban,
                                         beta_rural=beta_rural,
                                         kernel_distance_meters=kernel_distance_meters,
                                         scenario=scenario,
                                         state_name=state_name,
                                         historic_base_year=historic_base_year,
                                         projection_start_year=projection_start_year,
                                         projection_end_year=projection_end_year,
                                         time_step=time_step,
                                         write_raster=write_raster,
                                         write_array1d=write_array1d,
                                         write_array2d=write_array2d,
                                         write_csv=write_csv,
                                         compress_csv=compress_csv,
                                         write_logfile=write_logfile,
                                         run_number=run_number,
                                         output_total=output_total)

        self.sample = sample
        self.problem_dict = problem_dict

    def run_batch(self):
        """Run all model runs for the parameter set in parallel."""

        for index, i in enumerate(self.sample):

            print(f"Processing sample:  {index}")

            # instantiate model
            run = Model(grid_coordinates_file=self.grid_coordinates_file,
                        historical_rural_pop_raster=self.historical_rural_pop_raster,
                        historical_urban_pop_raster=self.historical_urban_pop_raster,
                        historical_suitability_raster=self.historical_suitability_raster,
                        projected_population_file=self.projected_population_file,
                        one_dimension_indices_file=self.one_dimension_indices_file,
                        output_directory=self.output_directory,
                        alpha_urban=self.set_param_value(self.ALPHA_URBAN_NAME, self.alpha_urban, i),
                        alpha_rural=self.set_param_value(self.ALPHA_RURAL_NAME, self.alpha_rural, i),
                        beta_urban=self.set_param_value(self.BETA_URBAN_NAME, self.beta_urban, i),
                        beta_rural=self.set_param_value(self.BETA_RURAL_NAME, self.beta_rural, i),
                        kernel_distance_meters=self.set_param_value(self.KD_NAME, self.kernel_distance_meters, i),
                        scenario=self.scenario,
                        state_name=self.state_name,
                        historic_base_year=self.historic_base_year,
                        projection_start_year=self.projection_start_year,
                        projection_end_year=self.projection_end_year,
                        time_step=self.time_step,
                        write_raster=self.write_raster,
                        write_array1d=self.write_array1d,
                        write_array2d=self.write_array2d,
                        write_csv=self.write_csv,
                        compress_csv=self.compress_csv,
                        write_logfile=self.write_logfile,
                        run_number=index,
                        output_total=self.output_total)

            run.downscale()

    def set_param_value(self, param_name, default_value, sample_list):
        """If the parameter is not in the sample set, use the default value.

        :param param_name:                  Parameter name
        :param default_value:               Default value of the parameter
        :param sample_list:                 List of samples

        """

        param_name_list = self.problem_dict.get('names')

        if param_name not in param_name_list:
            return default_value

        else:
            # get index of name in problem dict that corresponds with sample indices
            name_idx = param_name_list.index(param_name)

            return sample_list[name_idx]


class DeltaMomentIndependent:
    """Conduct analysis with the Delta Moment-Independent Measure using Latin Hypercube Sampling."""

    SETTING_OPTIONS = ('Urban', 'Rural')
    EXTENSION_OPTIONS = ('.tif', '.npy', '.csv', '.csv.gz')

    def __init__(self, problem_dict, sample, file_directory, setting, state_name, file_extension, output_file):

        self.problem_dict = problem_dict
        self.sample = sample
        self._file_directory = file_directory
        self._setting = setting
        self._state_name = state_name
        self._file_extension = file_extension
        self.output_file = output_file

    @property
    def file_extension(self):
        """Validate file extension."""

        if self._file_extension in self.EXTENSION_OPTIONS:
            return self._file_extension

        else:
            raise ValueError(f"Provided `file_extension` '{self._file_extension} is not in the acceptable values:  {self.EXTENSION_OPTIONS}")

    @property
    def state_name(self):
        """Target state name."""

        return self._state_name.lower()

    @property
    def setting(self):
        """Either 'Urban' or 'Rural'"""

        if self._setting in self.SETTING_OPTIONS:
            return self._setting

        else:
            raise ValueError(f"Provided `setting` '{self._setting}' is not in the acceptable values:  {self.SETTING_OPTIONS}")

    @property
    def file_directory(self):
        """Validate directory of input files from the model run."""

        if os.path.isdir(self._file_directory):
            return self._file_directory

        else:
            raise NotADirectoryError(f"File directory '{self._file_directory} does not exist.")

    @property
    def file_list(self):
        """Get a list of files to process."""

        files = [os.path.join(self.file_directory, i) for i in os.listdir(self.file_directory) if
                (i.split('_')[0] == self.state_name) and
                (self.file_extension in i) and
                (self.setting in i)]

        return self.validate_list(files)

    @property
    def n_gridcells(self):
        """Get the number of grid cells in an input dataset."""

        # load the first file in the file list
        return self.load_file(self.file_list[0]).shape[0]

    @property
    def data_array(self):
        """Combine output files into a single array.

        :return:                Array; shape = (n_runs, grid cells)

        """

        # shape (n_runs, grid cells)
        arr = np.zeros(shape=(len(self.file_list), self.n_gridcells))

        for index, i in enumerate(self.file_list):
            arr[index, :] = self.load_file(i)

        return arr

    def load_file(self, file):
        """Load data from file."""

        if self.file_extension == '.npy':
            return np.load(file)

        elif self.file_extension == '.csv.gz':
            return pd.read_csv(file, compression='gzip', sep=',')['value'].values

        else:
            raise ValueError(f"Loading '{self.file_extension}' is under development")

    def validate_list(self, in_list):
        """Ensure a list has a length > 0."""

        if len(in_list) > 0:
            return in_list

        else:
            raise ValueError(f"There are no files that match the search criteria of `state_name`='{self.state_name}', `file_extension`='{self.file_extension}', and `setting`='{self.setting}' in the file directory: '{self.file_directory}'")

    def run_analysis(self):
        """Run the sensitivity analysis and write the outputs to a file."""

        pool = Pool(processes=multiprocessing.cpu_count())

        # derive suitability estimates
        results = pool.map(self.delta_gridcell, [i for i in range(self.n_gridcells)])

        # write results to file
        self.write_output(results)

        return results

    def delta_gridcell(self, i):
        """Generate statistics for a gridcell from a n-dim array."""

        out_list = []

        # evaluate
        y = self.data_array[:, i]

        # if all values are the same
        unique_vals = np.unique(y).shape[0]

        if unique_vals > 1:

            # generate the sensitivity indices
            si = delta.analyze(self.problem_dict, self.sample, y, print_to_console=False)

            # write evaluated parameters
            for idx, key in enumerate(self.problem_dict['names']):
                out_list.append(f"{key},{si['delta'][idx]},{si['delta_conf'][idx]},{si['S1'][idx]},{si['S1_conf'][idx]},{i}\n")

        else:

            # write evaluated parameters
            for idx, key in enumerate(self.problem_dict['names']):
                out_list.append(f"{key},{np.nan},{np.nan},{np.nan},{np.nan},{i}\n")

        return out_list

    def write_output(self, result_list):
        """Write output to file."""

        with open(self.output_file, 'w') as out:

            # write header for output file
            out.write('param,delta,delta_conf,S1,S1_conf,gridcell\n')

            for element in result_list:
                for param in element:
                    out.write(param)
