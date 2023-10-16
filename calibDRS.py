# Define Standard Normal Variate for scatter correction of DRS data

import numpy as np


def snv_drs(drs_input):
    # A standard normal variate (SNV) normalization was used to process spectral data.
    # Using SNV normalization, the mean of each spectrum was set to zero and the
    # standard deviation was set to one.

    # Input: X: (m x n); spectral data where m = number of observations and n = number of features.

    # Output: snv_X: (m x n); snv transformed data.

    snv_output = ((drs_input.T - drs_input.mean(axis=1)) / drs_input.std(axis=1)).T

    return snv_output
