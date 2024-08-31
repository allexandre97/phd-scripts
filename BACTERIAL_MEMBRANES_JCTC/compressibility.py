# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pyedr

import numpy             as np
#import matplotlib.pyplot as plt

from optparse import OptionParser

parser = OptionParser()

parser.add_option('-f', '--edr_file',
                  default = 'POPC/Charmm/EDR/md_0.edr',
                  type    = str,
                  help    = 'input edr file',
                  action  = 'store')

options, arguments = parser.parse_args()

k_B = 1.380649e-23  # Boltzmann constant in J/K

def calculate_expected_value(X):
    """
    Calculate the expected value of an array X using a histogram to estimate the PDF,
    with the number of bins determined by the Freedman-Diaconis rule.
    
    Parameters:
    X (numpy array): Input array.
    
    Returns:
    float: Expected value of X.
    """
    # Calculate the number of bins using the Freedman-Diaconis rule
    iqr = np.subtract(*np.percentile(X, [75, 25]))
    bin_width = 2 * iqr * len(X)**(-1/3)
    num_bins = int(np.ceil((X.max() - X.min()) / bin_width))
    
    # Compute the histogram
    hist, bin_edges = np.histogram(X, bins=num_bins, density=True)
    
    # Calculate the bin centers
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Calculate the expected value
    expected_value = np.sum(hist * bin_centers * np.diff(bin_edges))
    
    return expected_value

def calculate_variance(X):
    """
    Calculate the variance of an array X using a histogram to estimate the PDF,
    with the number of bins determined by the Freedman-Diaconis rule.
    
    Parameters:
    X (numpy array): Input array.
    
    Returns:
    float: Variance of X.
    """
    # Calculate the number of bins using the Freedman-Diaconis rule
    iqr = np.subtract(*np.percentile(X, [75, 25]))
    bin_width = 2 * iqr * len(X)**(-1/3)
    num_bins = int(np.ceil((X.max() - X.min()) / bin_width))
    
    # Compute the histogram
    hist, bin_edges = np.histogram(X, bins=num_bins, density=True)
    
    # Calculate the bin centers
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Calculate the expected value
    expected_value = np.sum(hist * bin_centers * np.diff(bin_edges))
    
    # Calculate the expected value of X^2
    expected_value_X2 = np.sum(hist * bin_centers**2 * np.diff(bin_edges))
    
    # Calculate the variance
    variance = expected_value_X2 - expected_value**2
    
    return variance

if __name__ == "__main__":
    
    edr_dict = pyedr.edr_to_dict(options.edr_file)
    
    #%%
    
    time = edr_dict['Time']
    
    X = edr_dict['Box-X'][time > 400000] * 1e-9 # To metres
    Y = edr_dict['Box-Y'][time > 400000] * 1e-9 # To metres
    T = edr_dict['Temperature'][time > 400000]
    
    A = X*Y
    E_temp   = calculate_expected_value(T)
    E_area   = calculate_expected_value(A)
    Var_area = calculate_variance(A)
    
    K = (k_B * E_temp * E_area)/Var_area # (N * m * T^-1 * T * m^2 ) / (m^4) = N/m 
    
    print(f"Area : {1000*K:.2e} mN · m^-1")