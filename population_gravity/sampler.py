import numpy as np
import pandas as pd


class SampleDataFrame:
    """Construct sampling data frame expected by minimization function.

    :param alpha_range:                 array. evenly spaced alpha values over interval
    :param beta_range:                  array. evenly spaced beta values over interval

    :returns:                           dataframe.  containing alpha and beta parameters and an estimate space holder

    """

    def __init__(self, alpha_range, beta_range):

        self._alpha_range = alpha_range
        self._beta_range = beta_range

    @property
    def beta_length(self):
        """Get the length of the beta range"""
        return self.get_length(self._beta_range)

    @property
    def alpha_length(self):
        """Get the length of the alpha range"""
        return self.get_length(self._alpha_range)

    @staticmethod
    def get_length(arr):
        """Get length of array"""
        return arr.shape[0]

    @property
    def full_length(self):
        """Get the full length of the sample space"""
        return self.alpha_length * self.beta_length

    @property
    def sample_df(self):
        """Build a sample data frame to initialize the values of the minimization"""

        return pd.DataFrame(data={"alpha_param": np.repeat(self._alpha_range, self.beta_length).astype(np.float32),
                                  "beta_param": np.tile(self._beta_range, self.alpha_length).astype(np.float32),
                                  "estimate": np.full(self.full_length, np.nan, dtype=np.float64)})


def equal_interval(alpha_lower, alpha_upper, beta_lower, beta_upper, alpha_intervals, beta_intervals):
    """Evenly spaced numbers over a specified interval.

    :param alpha_lower:                          float. lower alpha boundary
    :param alpha_upper:                          float. upper alpha boundary
    :param beta_lower:                           float. lower beta boundary
    :param beta_upper:                           float. upper beta boundary
    :param alpha_intervals:                      int.   alpha number of intervals
    :param beta_intervals:                       int.   beta number of intervals

    :return:                                     [0] array. evenly spaced alpha values over interval
                                                 [1] array. evenly spaced beta values over interval
                                                 [2] dataframe. holds initial values

    """

    alpha_range = np.linspace(alpha_lower, alpha_upper, alpha_intervals)
    beta_range = np.linspace(beta_lower, beta_upper, beta_intervals)

    return alpha_range, beta_range, SampleDataFrame(alpha_range, beta_range).sample_df
