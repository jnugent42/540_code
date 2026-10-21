import uproot
import numpy as np
import pandas as pd
import scipy
import os
import re
import matplotlib.pyplot as plt
import inspect
from concurrent.futures import ThreadPoolExecutor
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from sklearn.neighbors import KernelDensity
from scipy.stats import mode
from iminuit import Minuit
import csv
import pickle
import pylandau
from joblib import Parallel, delayed
from concurrent.futures import ThreadPoolExecutor
import cupy as cp

# Set the working directory
os.chdir("/vols/hyperk/WCTE/WCTE540/nick/data")


def extract_pmt_voltages(file_name, pmt1=True, pmt2=False, pmt3=False, pmt4=False, chunk_size=10000000, num_workers=8):
    """
    Extracts voltage data for specified PMTs from the given ROOT file using chunked reading for optimization.
    
    Args:
        file_name (str): Path to the ROOT file.
        pmt1 (bool): If True, extract data for PMT1.
        pmt2 (bool): If True, extract data for PMT2.
        pmt3 (bool): If True, extract data for PMT3.
        pmt4 (bool): If True, extract data for PMT4.
        chunk_size (int): Number of entries to read at a time (optimize memory use).
        num_workers (int): Number of workers for parallel processing (I/O optimization).
        
    Returns:
        Tuple: Voltage arrays for each selected PMT in millivolts.
               Only the selected PMTs are returned in order.
    """
    voltages = []

    def load_voltages(tree_name):
        all_voltages = []
        for chunk in uproot.iterate([file_name + ":" + tree_name], ["volts"], 
                                    step_size=chunk_size, library="np", num_workers=num_workers):
            volts_gpu = cp.asarray(chunk["volts"])  # Move data to GPU
            all_voltages.append(volts_gpu.astype(cp.float32) * 1000)  # Convert to mV and to float32 on GPU
        return cp.concatenate(all_voltages)

    # Open the ROOT file and check which PMTs to extract
    with uproot.open(file_name) as f:
        with ThreadPoolExecutor(max_workers=num_workers) as executor:
            futures = []
            # Submit tasks for the selected PMTs
            if pmt1 and "PMT1data;1" in f:
                futures.append(executor.submit(load_voltages, "PMT1data;1"))
            if pmt2 and "PMT2data;1" in f:
                futures.append(executor.submit(load_voltages, "PMT2data;1"))
            if pmt3 and "PMT3data;1" in f:
                futures.append(executor.submit(load_voltages, "PMT3data;1"))
            if pmt4 and "PMT4data;1" in f:
                futures.append(executor.submit(load_voltages, "PMT4data;1"))

            # Gather results as they complete
            for future in futures:
                voltages.append(future.result())
    
    return tuple(voltages)

# def extract_pmt_voltages(file_name, pmt1=True, pmt2=False, pmt3=False, pmt4=False, chunk_size=5000000, num_workers=4):
#     """
#     Extracts voltage data for specified PMTs from the given ROOT file using chunked reading for optimization.
    
#     Args:
#         file_name (str): Path to the ROOT file.
#         pmt1 (bool): If True, extract data for PMT1.
#         pmt2 (bool): If True, extract data for PMT2.
#         pmt3 (bool): If True, extract data for PMT3.
#         pmt4 (bool): If True, extract data for PMT4.
#         chunk_size (int): Number of entries to read at a time (optimize memory use).
#         num_workers (int): Number of workers for parallel processing (I/O optimization).
        
#     Returns:
#         Tuple: Voltage arrays for each selected PMT in millivolts.
#                Only the selected PMTs are returned in order.
#     """
#     voltages = []

    # def load_voltages(tree_name):
    #     all_voltages = []
    #     for chunk in uproot.iterate([file_name + ":" + tree_name], ["volts"], 
    #                                 step_size=chunk_size, library="np", num_workers=num_workers):
    #         all_voltages.append(chunk["volts"].astype(np.float32) * 1000)  # Convert to mV and to float32
    #     return np.concatenate(all_voltages)

    # # Open the ROOT file and check which PMTs to extract
    # with uproot.open(file_name) as f:
    #     with ThreadPoolExecutor(max_workers=num_workers) as executor:
    #         futures = []
    #         # Submit tasks for the selected PMTs
    #         if pmt1 and "PMT1data;1" in f:
    #             futures.append(executor.submit(load_voltages, "PMT1data;1"))
    #         if pmt2 and "PMT2data;1" in f:
    #             futures.append(executor.submit(load_voltages, "PMT2data;1"))
    #         if pmt3 and "PMT3data;1" in f:
    #             futures.append(executor.submit(load_voltages, "PMT3data;1"))
    #         if pmt4 and "PMT4data;1" in f:
    #             futures.append(executor.submit(load_voltages, "PMT4data;1"))

    #         # Gather results as they complete
    #         for future in futures:
    #             voltages.append(future.result())
    
    # return tuple(voltages)

def process_window(test_pmt, i, window_size, peak_window_size, noise_start, noise_end):
    window_data = test_pmt[i:i + window_size]
    noise_data = test_pmt[i + noise_start:i + noise_end]
    
    # Use vectorized convolution
    conv_window = np.convolve(window_data, np.ones(peak_window_size) / peak_window_size, mode='valid')
    peak_value = np.max(conv_window)
    
    conv_noise = np.convolve(noise_data, np.ones(peak_window_size) / peak_window_size, mode='valid')
    n, _ = mode(conv_noise)
    
    if np.size(n) != 0:
        return peak_value, n
    else:
        return peak_value, 0

def peak_minus_noise(voltage_data, window_size=2000, peak_window_size=10, noise_start=500, noise_end=2000, n_jobs=-1):
    peaks = []
    noises = []
    
    for test_pmt in voltage_data:
        results = Parallel(n_jobs=n_jobs)(delayed(process_window)(test_pmt, i, window_size, peak_window_size, noise_start, noise_end) for i in range(0, len(test_pmt), window_size))
        
        peak, noise = zip(*results)
        peaks.append(peak)
        noises.append(noise)
    
    return peaks, noises

# def peak_minus_noise(voltage_data, window_size=2000, peak_window_size=10, noise_start=500, noise_end=2000):
#     peaks = []
#     noises = []
#     for test_pmt in voltage_data:
#         peak = []
#         noise = []
#         for i in range(0, len(test_pmt), window_size):
#             window_data = test_pmt[i:i + 200]
#             noise_data = test_pmt[i+noise_start:i+noise_end]
#             peak_value = np.max(np.convolve(window_data, np.ones(peak_window_size)/peak_window_size, mode='valid'))
            
#             mean = np.mean(noise_data)
#             stddev = np.std(noise_data)

#             n ,_ = mode(np.convolve(noise_data, np.ones(peak_window_size)/peak_window_size, mode='valid'))
#             #n = np.median(np.convolve(noise_data, np.ones(peak_window_size)/peak_window_size, mode='valid'))
#             if n!= 0:
#                 peak.append(peak_value)
#             else:
#                 peak.append(1000)
#             noise.append(n)
#         noises.append(noise)
#         peaks.append(peak)
#     return peaks, noises
    


def four_gaussian_with_exp_cons(x, amp0, mean0, stddev0, amp1, mean1, stddev1, amp2, amp3, exp_decay, exp_scale):
    mean2 = mean1 + (mean1 - mean0)
    mean3 = mean2 + (mean2 - mean1)
    stddev2 = np.sqrt(2) * stddev1
    stddev3 = np.sqrt(3) * stddev1
    clipped_x = np.clip(x, None, mean0 + 709 / exp_decay)
    gaussian0 = amp0 * np.exp(-((x - mean0) ** 2) / (2 * stddev0 ** 2))
    gaussian1 = amp1 * np.exp(-((x - mean1) ** 2) / (2 * stddev1 ** 2))
    gaussian2 = amp2 * np.exp(-((x - mean2) ** 2) / (2 * stddev2 ** 2))
    gaussian3 = amp3 * np.exp(-((x - mean3) ** 2) / (2 * stddev3 ** 2))
    exponential = np.where(clipped_x >= mean0, exp_scale * np.exp(-exp_decay * (clipped_x - mean0)), 0)
    exponential = np.where(np.isinf(exponential), 0, exponential)  # Handle overflow by setting to 0
    return gaussian0 + gaussian1 + gaussian2 + gaussian3 + exponential

def nick_multi_gaus(x, amp0, mean0, stddev0, amp1, mean1, stddev1, amp2, a, b, amp3, c, d, exp_decay, exp_scale):
    mean2 = mean1 + ((mean1 - mean0)) * a
    mean3 = mean2 + ((mean2 - mean1)) * c
    stddev2 = np.sqrt(2) * stddev1 * b
    stddev3 = np.sqrt(3) * stddev1 * d
    clipped_x = np.clip(x, None, mean0 + 709 / exp_decay)
    gaussian0 = amp0 * np.exp(-((x - mean0) ** 2) / (2 * stddev0 ** 2))
    gaussian1 = amp1 * np.exp(-((x - mean1) ** 2) / (2 * stddev1 ** 2))
    gaussian2 = amp2 * np.exp(-((x - mean2) ** 2) / (2 * stddev2 ** 2))
    gaussian3 = amp3 * np.exp(-((x - mean3) ** 2) / (2 * stddev3 ** 2))
    exponential = np.where(clipped_x >= mean0, exp_scale * np.exp(-exp_decay * (clipped_x - mean0)), 0)
    exponential = np.where(np.isinf(exponential), 0, exponential)  # Handle overflow by setting to 0
    return gaussian0 + gaussian1 + gaussian2 + gaussian3 + exponential

def four_gaussian_with_exp(x, amp0, mean0, stddev0, amp1, mean1, stddev1, amp2, mean2, stddev2, amp3, mean3, stddev3, exp_decay, exp_scale):
    clipped_x = np.clip(x, None, mean0 + 709 / exp_decay)
    gaussian0 = amp0 * np.exp(-((x - mean0) ** 2) / (2 * stddev0 ** 2))
    gaussian1 = amp1 * np.exp(-((x - mean1) ** 2) / (2 * stddev1 ** 2))
    gaussian2 = amp2 * np.exp(-((x - mean2) ** 2) / (2 * stddev2 ** 2))
    gaussian3 = amp3 * np.exp(-((x - mean3) ** 2) / (2 * stddev3 ** 2))
    exponential = np.where(clipped_x >= mean0, exp_scale * np.exp(-exp_decay * (clipped_x - mean0)), 0)
    exponential = np.where(np.isinf(exponential), 0, exponential)  # Handle overflow by setting to 0
    return gaussian0 + gaussian1 + gaussian2 + gaussian3 + exponential

def multi_gaussian_with_exp_mu3(x, amp0, mean0, stddev0, amp1, mean1, stddev1, amp2, mean2, stddev2, exp_decay, exp_scale):
    clipped_x = np.clip(x, None, mean0 + 709 / exp_decay)
    gaussian0 = amp0 * np.exp(-((x - mean0) ** 2) / (2 * stddev0 ** 2))
    gaussian1 = amp1 * np.exp(-((x - mean1) ** 2) / (2 * stddev1 ** 2))
    gaussian2 = amp2 * np.exp(-((x - mean2) ** 2) / (2 * stddev2 ** 2))
    exponential = np.where(clipped_x >= mean0, exp_scale * np.exp(-exp_decay * (clipped_x - mean0)), 0)
    exponential = np.where(np.isinf(exponential), 0, exponential)  # Handle overflow by setting to 0
    return gaussian0 + gaussian1 + gaussian2 + exponential

def multi_gaussian_with_exp(x, amp0, mean0, stddev0, amp1, mean1, stddev1, amp2, stddev2, exp_decay, exp_scale):
    mean2 = mean1 + (mean1 - mean0)
    clipped_x = np.clip(x, None, mean0 + 709 / exp_decay)
    gaussian0 = amp0 * np.exp(-((x - mean0) ** 2) / (2 * stddev0 ** 2))
    gaussian1 = amp1 * np.exp(-((x - mean1) ** 2) / (2 * stddev1 ** 2))
    gaussian2 = amp2 * np.exp(-((x - mean2) ** 2) / (2 * stddev2 ** 2))
    exponential = np.where(clipped_x >= mean0, exp_scale * np.exp(-exp_decay * (clipped_x - mean0)), 0)
    exponential = np.where(np.isinf(exponential), 0, exponential)  # Handle overflow by setting to 0
    return gaussian0 + gaussian1 + gaussian2 + exponential

def multi_gaussian_with_exp_cons(x, amp0, mean0, stddev0, amp1, mean1, stddev1, amp2, exp_decay, exp_scale):
    mean2 = mean1 + (mean1 - mean0)
    clipped_x = np.clip(x, None, mean0 + 709 / exp_decay)
    stddev2 = np.sqrt(2) * stddev1
    gaussian0 = amp0 * np.exp(-((x - mean0) ** 2) / (2 * stddev0 ** 2))
    gaussian1 = amp1 * np.exp(-((x - mean1) ** 2) / (2 * stddev1 ** 2))
    gaussian2 = amp2 * np.exp(-((x - mean2) ** 2) / (2 * stddev2 ** 2))
    exponential = np.where(clipped_x >= mean0, exp_scale * np.exp(-exp_decay * (clipped_x - mean0)), 0)
    exponential = np.where(np.isinf(exponential), 0, exponential)  # Handle overflow by setting to 0
    return gaussian0 + gaussian1 + gaussian2 + exponential

def multi_gaussian_with_exp_mu3_cons(x, amp0, mean0, stddev0, amp1, mean1, stddev1, amp2, mean2, exp_decay, exp_scale):
    clipped_x = np.clip(x, None, mean0 + 709 / exp_decay)
    stddev2 = np.sqrt(2) * stddev1
    gaussian0 = amp0 * np.exp(-((x - mean0) ** 2) / (2 * stddev0 ** 2))
    gaussian1 = amp1 * np.exp(-((x - mean1) ** 2) / (2 * stddev1 ** 2))
    gaussian2 = amp2 * np.exp(-((x - mean2) ** 2) / (2 * stddev2 ** 2))
    exponential = np.where(clipped_x >= mean0, exp_scale * np.exp(-exp_decay * (clipped_x - mean0)), 0)
    exponential = np.where(np.isinf(exponential), 0, exponential)  # Handle overflow by setting to 0
    return gaussian0 + gaussian1 + gaussian2 + exponential

def multi_gaussian_mu3(x, amp0, mean0, stddev0, amp1, mean1, stddev1, amp2, mean2, stddev2):
    gaussian0 = amp0 * np.exp(-((x - mean0) ** 2) / (2 * stddev0 ** 2))
    gaussian1 = amp1 * np.exp(-((x - mean1) ** 2) / (2 * stddev1 ** 2))
    gaussian2 = amp2 * np.exp(-((x - mean2) ** 2) / (2 * stddev2 ** 2))
    return gaussian0 + gaussian1 + gaussian2

def multi_gaussian(x, amp0, mean0, stddev0, amp1, mean1, stddev1, amp2, stddev2):
    mean2 = mean1 + (mean1 - mean0)
    gaussian0 = amp0 * np.exp(-((x - mean0) ** 2) / (2 * stddev0 ** 2))
    gaussian1 = amp1 * np.exp(-((x - mean1) ** 2) / (2 * stddev1 ** 2))
    gaussian2 = amp2 * np.exp(-((x - mean2) ** 2) / (2 * stddev2 ** 2))
    return gaussian0 + gaussian1 + gaussian2

def multi_gaussian_mu3_cons(x, amp0, mean0, stddev0, amp1, mean1, stddev1, amp2, mean2):
    stddev2 = np.sqrt(2) * stddev1
    gaussian0 = amp0 * np.exp(-((x - mean0) ** 2) / (2 * stddev0 ** 2))
    gaussian1 = amp1 * np.exp(-((x - mean1) ** 2) / (2 * stddev1 ** 2))
    gaussian2 = amp2 * np.exp(-((x - mean2) ** 2) / (2 * stddev2 ** 2))
    return gaussian0 + gaussian1 + gaussian2

def multi_gaussian_cons(x, amp0, mean0, stddev0, amp1, mean1, stddev1, amp2):
    stddev2 = np.sqrt(2) * stddev1
    mean2 = mean1 + (mean1 - mean0)
    gaussian0 = amp0 * np.exp(-((x - mean0) ** 2) / (2 * stddev0 ** 2))
    gaussian1 = amp1 * np.exp(-((x - mean1) ** 2) / (2 * stddev1 ** 2))
    gaussian2 = amp2 * np.exp(-((x - mean2) ** 2) / (2 * stddev2 ** 2))
    return gaussian0 + gaussian1 + gaussian2



#def efficiency_histogram(voltage, initial_guess, cut, pmt_label, fit_function, bins=100, range_=[0, 20], label='', plot=False, midpoint_plot=False, x_label='', bounds=None, method='lm'):
    results = []

    if plot:
        _, ax = plt.subplots(figsize=(8, 6))
        counts, bin_edges, _ = plt.hist(voltage, bins=bins, range=(range_[0], range_[1]), histtype='step', label=label, color='orange')
    if not plot: 
        counts, bin_edges = np.histogram(voltage, bins=bins, range=(range_[0], range_[1]))
    
    bin_midpoints = (bin_edges[:-1] + bin_edges[1:]) / 2
    #Calculate errors

    # Assign errors
    errors = 1 / np.sqrt(np.maximum(counts, 1))


    # kde = KernelDensity(kernel='gaussian', bandwidth=0.3).fit(bin_midpoints.reshape(-1, 1))
    # densities = np.exp(kde.score_samples(bin_midpoints.reshape(-1, 1)))
    # errors = 1 / densities

    cons_bins = np.array(bin_midpoints)[np.array(bin_midpoints) < cut]
    errors = errors[:len(cons_bins)]

    try:
        # Fit the data to the specified model
        popt, _ = curve_fit(fit_function, cons_bins, counts[:len(cons_bins)], p0=initial_guess, bounds=bounds, maxfev=60000, sigma=errors, method=method)

        # Extract the fitted parameters
        cutoff = popt[4] - 1 * popt[5]  # Assuming mean2 and stddev2 are at indices 4 and 5

        # Generate x-values for plotting the fit
        x_fit = np.linspace(min(bin_midpoints), max(bin_midpoints) + 6, 1000)
        y_fit = fit_function(x_fit, *popt)

        # Filter the data: Keep values greater than the cutoff
        kept_data = np.array(voltage)[np.array(voltage) > cutoff]
        efficiency = len(kept_data) / len(voltage) * 100  # Percentage of kept voltage

        if plot:
            # Plot the histogram and the fitted curve
            ax.plot(x_fit, y_fit, label='Fitted Model', color='red')
            # ax.errorbar(bin_midpoints, np.array(counts), yerr=errors, fmt='o', ecolor='r', capthick=2, label='Data')
            ax.axvline(x=cutoff, color='blue', linestyle='--', linewidth=2, label=f'Cutoff at {cutoff:.2f}')
            ax.set_xlabel(x_label, size=15)
            ax.set_yscale('log')
            ax.set_ylabel('Entries', size=15)
            ax.set_ylim(10**2, 10**4)
            if midpoint_plot:
                ax.plot(bin_midpoints, counts, 'o', label='bin midpoints')
                ax.plot(cons_bins, counts[:len(cons_bins)], 'o', label='Midpoints cutoff', color='purple')
            ax.legend()
            ax.set_title(f'PMT {pmt_label}')
        # Store the results
        result = {'PMT': pmt_label, 'Efficiency (%)': efficiency}
        for i, param in enumerate(popt):
            argspec = inspect.getfullargspec(fit_function)
            name = argspec.args[i+1]
            result[name] = param
        results.append(result)

    except RuntimeError as e:
        print(f"Error fitting data for PMT {pmt_label}: {e}")
        cutoff=0
        efficiency=0
        df_results=0
        
    if plot:
        # Adjust layout
        plt.tight_layout()
        plt.show()
        
    # Create a DataFrame for the results
    df_results = pd.DataFrame(results)
    print(df_results.to_markdown())

def landau(x, amp, mean, sigma, c):
    return amp * pylandau.landau(x, mean, sigma) + c

def TTS(voltage_data, window_size=2000, peak_window_size=10):
    means_all = []
    sigmas_all = []
    indices_all = []
    for test_pmt in voltage_data:
        means_1 = []
        sigmas_1 = []
        indices_1 = []
        for i in range(0, len(test_pmt), window_size):
            window_data = test_pmt[i:i + 200]
            peak_value = np.max(np.convolve(window_data, np.ones(peak_window_size)/peak_window_size, mode='valid'))
            index_peak_value = np.argmax(np.convolve(window_data, np.ones(peak_window_size)/peak_window_size, mode='valid'))
            try:
                # Ensure the fitting range is within the bounds of window_data
                fit_range_start = max(0, index_peak_value - 10)
                fit_range_end = min(len(window_data), index_peak_value + 10)
                if fit_range_end - fit_range_start < 3:
                    raise ValueError("Insufficient data points for fitting")
                
                # Ensure the initial guess for peak_value is a single value
                peak_value = peak_value.item() if isinstance(peak_value, np.ndarray) else peak_value
                
                popt, _ = curve_fit(landau, np.arange(fit_range_start, fit_range_end), window_data[fit_range_start:fit_range_end], p0=[peak_value, index_peak_value, 1], maxfev=60000)
                m = popt[1]
                s = popt[2]
                means_1.append(m)
                sigmas_1.append(s)
                indices_1.append(index_peak_value)
            except (RuntimeError, ValueError) as e:
                print(f"Error fitting data at index {i}: {e}")
                means_1.append(np.nan)
                sigmas_1.append(np.nan)
                indices_1.append(np.nan)
        means_all.append(means_1)
        sigmas_all.append(sigmas_1)
        indices_all.append(indices_1)
    return indices_all, means_all, sigmas_all


def efficiency_histogram(data, initial_guess, pmt_label, fit_function, range_=[0, 20], label='', plot=False, bounds=None, step_size=None):

    bin_edges = np.arange(0, range_[1], 0.15625)  
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])  
    bin_contents, _ = np.histogram(data, bins=bin_edges) 
    errors = 1/np.sqrt(np.where(bin_contents == 0, 1, bin_contents)) 

    # Chi2 function for Minuit
    def chi2(*params):
        expected = fit_function(bin_centers, *params)
        return np.sum(((bin_contents - expected) / errors) ** 2) 

    # Initial guesses
    guess = initial_guess

    # Pass the initial guesses into Minuit
    m = Minuit(chi2, *guess)
    m.errordef = Minuit.LEAST_SQUARES  # Chi2 fitting

    # Set parameter limits before fitting
    if bounds:
        m.limits = bounds

    if step_size:
        m.errors = step_size

    m.migrad()  # Perform the fit

    # Extract best-fit parameters
    fit_params = m.values
    cutoff = fit_params[4] - 1 * fit_params[5]
    kept_data = np.array(data)[np.array(data) > cutoff]
    count_rate = len(kept_data) / len(data) * 100 

    # Generate the fitted curve
    x_fit = np.linspace(bin_edges[0], bin_edges[-1], 1000)
    y_fit = fit_function(x_fit, *fit_params)

    # Plot the histogram and fitted function
    if plot:
        plt.figure(figsize=(8, 6))
        plt.plot(x_fit, y_fit, label="Fit", color="red", linewidth=2)
        plt.hist(data, bins=bin_edges, histtype='step', label="Histogram", color="blue")
        plt.axvline(cutoff, color='green', linestyle='--', label='Threshold')
        plt.xlabel("Voltage (mV)", size=20)
        plt.ylabel("Counts", size=20)
        plt.title(label, size=20)
        plt.legend(fontsize=15)
        plt.yscale('log')
        plt.ylim(10**2)
        plt.grid(True)
        plt.text(0.95, 0.95, f'Count Rate: {count_rate:.2f}%', horizontalalignment='right', verticalalignment='top', transform=plt.gca().transAxes, fontsize=15, bbox=dict(facecolor='white', alpha=0.5))  # Add efficiency text to the plot
        plt.show()
        
        result = {'PMT': pmt_label, 'count rate (%)': count_rate}
        for i, param in enumerate(m.values):
            argspec = inspect.getfullargspec(fit_function)
            name = argspec.args[i+1]
            result[name] = param
            result[f'{name}_error'] = m.errors[i]
        
        results = [result]
        df_results = pd.DataFrame(results)

    # Output the results
    print(df_results.to_markdown())
    return(fit_params, cutoff, count_rate)
    


def efficiency_histograms(voltage, initial_guess, cuts, pmt_labels, bins=100, range_=[0, 20], label=''):
    results = []

    fig, axes = plt.subplots(2, 2, figsize=(6, 4))
    axes = axes.flatten()

    for l in range(len(cuts)):
        counts, bin_edges, _ = axes[l].hist(voltage[l], bins=bins, range=(range_[0],range_[1]), histtype='step', label=label, color='orange')
        bin_midpoints = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Initial guesses for the parameters
        initial_guess_f = initial_guess[l]

        cons_bins = np.array(bin_midpoints)[np.array(bin_midpoints) < cuts[l]]

        try:
            # Fit the data to the combined model
            popt, _ = curve_fit(lambda x, amp1, mean1, stddev1, amp2, mean2, stddev2, amp3, stddev3, exp_decay, exp_scale: 
                                multi_gaussian_with_exp(x, amp1, mean1, stddev1, amp2, mean2, stddev2, amp3, stddev3, exp_decay, exp_scale, m=mean1),
                                cons_bins, counts[:len(cons_bins)], p0=initial_guess_f, bounds=(0, np.inf), maxfev=60000)

            # Extract the fitted parameters
            mean2 = popt[4]
            stddev2 = popt[5]
            mean1 = popt[1]

            exp_decay = popt[8]
            exp_scale = popt[9]

            cutoff = mean2 - 1 * stddev2

            # Generate x-values for plotting the fit
            x_fit = np.linspace(min(bin_midpoints), max(bin_midpoints) + 6, 1000)
            y_fit = multi_gaussian_with_exp(x_fit, *popt, m=mean1)

            # Filter the data: Keep values greater than the cutoff
            kept_data = np.array(voltage[l])[np.array(voltage[l]) > cutoff]
            efficiency = len(kept_data) / len(voltage[l]) * 100  # Percentage of kept voltage[l]

            # Plot the histogram and the fitted curve
            axes[l].plot(x_fit, y_fit, label='Fitted Model', color='red')
            axes[l].axvline(x=cutoff, color='blue', linestyle='--', linewidth=2, label=f'Cutoff at {cutoff:.2f}')
            axes[l].set_xlabel('(Peak-noise) mV', size=15)
            axes[l].set_yscale('log')
            axes[l].set_ylabel('Entries', size=15)
            axes[l].set_ylim(10**2, 10**4)
            axes[l].legend()
            axes[l].set_title(f'PMT {pmt_labels[l]}')

            # Store the results
            results.append({
                'PMT': pmt_labels[l],
                'Gaussian 1 Amplitude': popt[0],
                'Gaussian 1 Mean': popt[1],
                'Gaussian 1 Stddev': popt[2],
                'Gaussian 2 Amplitude': popt[3],
                'Gaussian 2 Mean': popt[4],
                'Gaussian 2 Stddev': popt[5],
                'Gaussian 3 Amplitude': popt[6],
                'Gaussian 3 Stddev': popt[7],
                'Exponential Decay Rate': popt[8],
                'Exponential Scale': popt[9],
                'Efficiency (%)': efficiency
            })

        except RuntimeError as e:
            print(f"Error fitting data for PMT {pmt_labels[l]}: {e}")

    # Adjust layout
    plt.tight_layout()
    plt.show()

    # Create a DataFrame for the results
    df_results = pd.DataFrame(results)

    # Display the table
    print(df_results.to_string(index=False))


def calculate_ratios(voltage, window_size=2000, peak_window_size=10, use_rolling_average=True, return_ratios=False, return_indices=False):
    ratios_list = [] if return_ratios else None  # To store the ratios for each dataset
    indices = [] if return_indices else None  # To store the indices where the peak occurs
    average = []  # To store the average or peak values

    for test_pmt in voltage:
        ratios2 = [] if return_ratios else None  # Reset ratios for the current dataset
        test = [] if return_indices else None  # Indices
        test2 = []  # Peak values

        # Loop over the dataset in 2000-sample windows
        for i in range(0, len(test_pmt), window_size):
            # Define the 2000-sample window
            window_data = test_pmt[i:i + 200]

            if use_rolling_average:
                # Calculate rolling average
                rolling_average = np.convolve(window_data, np.ones(peak_window_size) / peak_window_size, mode='valid')
                peak_value = np.max(rolling_average)
                peak_window_start = np.argmax(rolling_average)
            else:
                # Use maximum value as peak
                peak_value = np.max(window_data)
                peak_window_start = np.argmax(window_data)

            test2.append(peak_value)
            if return_indices:
                test.append(peak_window_start)

            # Noise data (fixed section of the window)
            noise_data = test_pmt[i+500:i+2000]
            noise_average ,_ = mode(np.convolve(noise_data, np.ones(peak_window_size)/peak_window_size, mode='valid'))
            noise_peak = np.max(noise_data)

            # Calculate the ratio
            if noise_average != 0:  # Avoid division by zero
                if use_rolling_average:
                    ratio = np.abs((peak_value - noise_average) / noise_average)
                else:
                    ratio = np.abs((peak_value - noise_peak) / noise_peak)
            else:
                ratio = 1000  # Handle cases where noise average is zero

            if return_ratios:
                ratios2.append(ratio)

        average.append(test2)
        if return_indices:
            indices.append(test)
        if return_ratios:
            ratios_list.append(ratios2)  # Store the ratios for the current dataset

    # Conditionally return the arrays based on the flags
    if return_ratios and return_indices:
        return ratios_list, indices, average
    elif return_ratios:
        return ratios_list, average
    elif return_indices:
        return indices, average
    else:
        return average



