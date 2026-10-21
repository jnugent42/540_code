#!/usr/bin/env python
# coding: utf-8

# ## 0. Importing data

# In[3]:


get_ipython().system('pip install cupy-cuda12x')
import cupy as cp
print(cp.cuda.runtime.runtime_version)


# In[ ]:


from lb_functions import *
# Path to the folder where your ROOT files are located
data_folder = "/vols/hyperk/WCTE/WCTE540/nick/data"

# Get a list of all ROOT files in the data folder
root_files = [f for f in os.listdir(data_folder) if f.endswith(".root")]

#Words to extract the specific files needed.
keywords = ['Calibrate', 'LDOn', '2024-09-17']

#Filter out all run 1 files with laser on
run_1 =  sorted([s for s in root_files if all(keyword in s for keyword in keywords)])
print(run_1)

keywords2 = ['Calibrate', 'LDOn', '2024-09-18']

#Filter out all run 1 files with laser on
run_2 =  sorted([s for s in root_files if all(keyword in s for keyword in keywords2)])
print(run_2)


# In[ ]:


voltage_r1 = []
for i in range(len(run_1)):
    pmt1, = extract_pmt_voltages(run_1[i], pmt1=True, chunk_size=2000000, num_workers=4)
    voltage_r1.append(pmt1)


# In[ ]:


voltage_r2 = []
for i in range(len(run_2)):
    pmt1, = extract_pmt_voltages(run_1[i], pmt1=True, chunk_size=2000000, num_workers=4)
    voltage_r2.append(pmt1)


# In[ ]:


kw = ['Ball2', '176.02mA', 'LDON']

test =  sorted([s for s in root_files if all(keyword in s for keyword in kw)])
print(test)

test_ = []
pmt1, = extract_pmt_voltages(test[2], pmt1=True, chunk_size=2000000)
test_.append(pmt1)


# In[ ]:


fitted_params = []
for i in range(0, len(test_[0][0:20000]), 2000):
    means_1 = []
    sigmas_1 = []
    indices_1 = []
    peak_window_size = 10
    window_data = test_[0][i:i + 200]
    conv = np.convolve(window_data, np.ones(peak_window_size)/peak_window_size, mode='valid')
    peak_value = np.max(conv)
    index_peak_value = np.argmax(window_data)
    n ,_ = mode(conv)

    fit_range_start = index_peak_value - 10
    fit_range_end = index_peak_value + 10

    try:
        popt, _ = curve_fit(landau, np.arange(fit_range_start, fit_range_end), window_data[fit_range_start:fit_range_end], p0=[peak_value, index_peak_value, 1, n], maxfev=60000)
        fitted_params.append(popt)
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
        fitted_params.append([np.nan, np.nan, np.nan, np.nan])

plt.plot(test_[0][0:20000])
for j, popt in enumerate(fitted_params):
    if np.isnan(popt[1]):
        continue
    j*=2000
    x = np.arange(j+popt[1] - 10, j+popt[1] + 10, 0.01)
    plt.plot(x, landau(x, popt[0], j+popt[1], popt[2], popt[3]), 'r')

plt.xlabel('Index')
plt.ylabel('Value')
m = 8000
plt.xlim(4100, 4200)
plt.show()


# ## 1. New Technique

# In[ ]:


#Checking waveform
c=0
d=c+50000
peak = np.argmax(voltage_r1[0][c:d])
a = 4
b = c+peak -a
e = b+10


#Finding peaks using new method
x = voltage_r1[0]
peaks, x_ = find_peaks(x, height=0, width=0, distance=1800, prominence=None) #set to 1900
averages = []
for peak in peaks:
    start = peak - 4
    end = peak + 6
    avg = np.mean(voltage_r1[0][start:end])
    averages.append(avg)


# In[ ]:


# # Function to write a list of arrays to a CSV file
# import csv
# output = 'test.csv'
# with open(output,'w') as f:
#     writer = csv.writer(f)
#     writer.writerow(averages)


#Reads function from written csv file, commented out

import csv
from lb_functions import *
os.chdir("/vols/hyperk/WCTE/WCTE540/nick/data")
def read_csv_to_list(filename):
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile)
        return [[float(value) for value in row] for row in reader]

averages = read_csv_to_list('test.csv')[0]


# In[ ]:


#Different initial guesses
cut = 15
amp1 = 4839
mean1 = 1.5
sigma1 = 0.5
amp2 = 931.0
mean2 = 6.5
sigma2 = 1.5
amp3 = 463.0
mean3 = 12.5
sigma3 = 2.5
exp_decay =0.3
exp_scale = 8

#All inital guesses for different functions together
guesses = [[amp1, mean1, sigma1, amp2, mean2, sigma2, amp3, mean3, sigma3, exp_decay, exp_scale],
           [amp1, mean1, sigma1, amp2, mean2, sigma2, amp3, sigma3, exp_decay, exp_scale],
           [amp1, mean1, sigma1, amp2, mean2, sigma2, amp3, exp_decay, exp_scale],
           [amp1, mean1, sigma1, amp2, mean2, sigma2, amp3, mean3, exp_decay, exp_scale],
           [amp1, mean1, sigma1, amp2, mean2, sigma2, amp3, mean3, sigma3],
           [amp1, mean1, sigma1, amp2, mean2, sigma2, amp3, sigma3],
           [amp1, mean1, sigma1, amp2, mean2, sigma2, amp3, mean3],
           [amp1, mean1, sigma1, amp2, mean2, sigma2, amp3]]



functions = [
    multi_gaussian_with_exp_mu3,
    multi_gaussian_with_exp,
    multi_gaussian_with_exp_cons,
    multi_gaussian_with_exp_mu3_cons,
    multi_gaussian_mu3,
    multi_gaussian,
    multi_gaussian_mu3_cons,
    multi_gaussian_cons
]
for i in range(len(functions)):
    efficiency_histogram(averages, guesses[i], cut, pmt_label='B '+functions[i].__name__,
                          fit_function=functions[i], bins=np.arange(0,20,0.15625), range_=[0, 20], label='', plot=True, midpoint_plot=True)


# plt.hist(averages, bins=80, range=(0,25), histtype='step', label='new method')
# #plt.hist(pmn_r1[0], bins=100, range=(0, 30), histtype='step', label='peak-noise')
# plt.xlabel('Average Voltage (mV)', fontsize=20)
# plt.ylabel('Entries', fontsize=20)      
# plt.legend()


# In[ ]:


peaks, data = find_peaks(voltage_r1[0], height=0, width=0, distance=1800, prominence=None)
widths = data['widths']
heights = data['peak_heights']
#peaks, _ = find_peaks(x, height=0, width=0, distance=10, prominence=10)


#Plotting widths and heights, as well as peak positions
# fig, ax1 = plt.subplots(figsize=(12, 6))
# ax1.plot(np.arange(len(x)), x, label ='Laser')
# ax1.plot(peaks, x[peaks], "x", label='peaks') #peaks is the values of the position of peaks and x[peaks] is the height of the peaks

# ax1.vlines(x=peaks, ymin=x[peaks] - _["prominences"],
#            ymax = x[peaks], color = "C1")
# ax1.hlines(y=_["width_heights"], xmin=_["left_ips"],
#            xmax=_["right_ips"], color = "C1")
# #fig.legend(prop={'size': 15})

def gaussian(x, amp, mean, stddev):
    return amp * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))


#plt.xlabel('Index', size=20)

#Fitting guassian to peaks
amps = []
means = []
stds = []
fpvs = []
a1std=[]
for peak in peaks:
    amp = voltage_r1[0][peak]
    mean = peak
    stddev = data["width_heights"][np.where(peaks == peak)[0][0]]/2

    start = max(0, mean - 20)
    end = min(len(voltage_r1[0]), mean + 20)
    try:
        popt, n = curve_fit(gaussian, np.arange(start, end), voltage_r1[0][start:end], p0=[amp, mean, stddev], maxfev=60000)
        amps.append(popt[0])
        means.append(popt[1])
        stds.append(popt[2])
        a1std.append(gaussian(popt[1] + popt[2], *popt))
    except RuntimeError as e:
            print(f"Error fitting peak at index {peak}: {e}")
            fpvs.append(peak)
    #ax1.plot(np.arange(mean-15, mean+15, 0.01), gaussian(np.arange(mean-15, mean+15, 0.01), *popt), label=f'Gaussian fit at {mean}')

#plt.xlim(14100, 14200)
#plt.show()


# In[ ]:


maximum = 35
plt.hist(a1std, bins=np.arange(2,maximum,2*0.15625), range=(0, maximum))
plt.yscale('log')
plt.ylim(10**2)
plt.show()


# In[ ]:


1-0.875


# In[ ]:


fig, ax1 = plt.subplots(figsize=(12,6))
maximum=8
counts_w, bins_w, s_w = ax1.hist(widths, bins=np.arange(1,maximum,4*0.125), range=(0,maximum))
ax1.set_ylabel('Counts', size=20)
ax1.set_xlabel('Width', size=20)
heights = data['peak_heights']


# In[ ]:


np.sort(np.unique(heights))
2.5-0.9375


# In[ ]:


fig, ax1 = plt.subplots(figsize=(12,6))
maxmimum=65
ax1.hist(heights, bins=np.arange(0,maxmimum,1.5625), range=(0,maxmimum))
ax1.set_ylabel('Counts', size=20)
ax1.set_xlabel('Height', size=20)


# In[ ]:


np.diff(np.sort(np.unique(heights * widths_lin)))


# In[ ]:


#Trying to impliment width in different ways
from scipy.optimize import curve_fit

def gaussian(x, amp, mean, stddev):
    return amp * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))

widths_exp = gaussian(widths, (1/12*np.pi*np.var(widths)), np.mean(widths), np.var(widths))+ np.where(widths>np.mean(widths), 0.1*widths,0 )
widths_lin = 0.1*widths+0.1
widths_tv = 1-np.exp(-widths/10)
widths_tanh = np.tanh(widths/10)
widths_sig = 1/(1+np.exp(-widths/10))-0.5

counts_exp, bins_exp = np.histogram(heights * widths_exp, bins=100, range=(0, 100))
counts_lin, bins_lin= np.histogram(heights * widths_lin, bins=100, range=(0, 400))
counts_tv, bins_tv = np.histogram(heights * widths_tv, bins=100, range=(0, 20))
counts_tanh, bins_tanh = np.histogram(heights * widths_tanh, bins=100, range=(0, 20))
counts_sig, bins_sig = np.histogram(heights * widths_sig, bins=100, range=(0, 20))


def gaussian(x, amp, mean, stddev):
    return amp * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))

def three_gaussians(x, amp1, mean1, stddev1, amp2, mean2, stddev2, amp3, stddev3):
    mean3 = mean2 + (mean2-mean1)
    return (gaussian(x, amp1, mean1, stddev1) + 
            gaussian(x, amp2, mean2, stddev2) + 
            gaussian(x, amp3, mean3, stddev3))

def fit_three_gaussians(bin_edges, counts):
    bin_midpoints = (bin_edges[:-1] + bin_edges[1:]) / 2
    initial_guess = [max(counts), 5, 1, 
                     max(counts)/2, 12, 1, 
                     max(counts)/3, 1]
    popt, _ = curve_fit(three_gaussians, bin_midpoints, counts, p0=initial_guess)
    return popt, bin_midpoints

# Fit three gaussians to the bin midpoints of each histogram
# popt_exp, bin_midpoints_exp = fit_three_gaussians(bins_exp, counts_exp)
# popt_lin, bin_midpoints_lin = fit_three_gaussians(bins_lin, counts_lin)
# popt_tv, bin_midpoints_tv = fit_three_gaussians(bins_tv, counts_tv)
# popt_tanh, bin_midpoints_tanh = fit_three_gaussians(bins_tanh, counts_tanh)
# popt_sig, bin_midpoints_sig = fit_three_gaussians(bins_sig, counts_sig)

# # Plot the fitted gaussians
fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.hist(heights * widths_exp, bins=100, range=(0, 30), histtype='step', label='Data')
# ax1.plot(bin_midpoints_exp, three_gaussians(bin_midpoints_exp, *popt_exp), label='Fit')
ax1.set_ylabel('Entries', size=20)
ax1.set_xlabel('Voltage (mV)', size=20)
# ax1.set_yscale('log')
# ax1.set_ylim(10**2)
ax1.legend()
fig.show()

fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.hist(heights * widths_lin, bins=100, range=(0, 30), histtype='step', label='Data')
# ax1.plot(bin_midpoints_lin, three_gaussians(bin_midpoints_lin, *popt_lin), label='Fit')
ax1.set_ylabel('Entries', size=20)
ax1.set_xlabel('Voltage (mV)', size=20)
# ax1.set_yscale('log')
# ax1.set_ylim(8*10**1)
ax1.legend()
fig.show()

fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.hist(heights * widths_tv, bins=100, range=(0, 20), histtype='step', label='Data')
# ax1.plot(bin_midpoints_tv, three_gaussians(bin_midpoints_tv, *popt_tv), label='Fit')
ax1.set_ylabel('Entries', size=20)
ax1.set_xlabel('Voltage (mV)', size=20)
# ax1.set_yscale('log')
# ax1.set_ylim(10**2)
ax1.legend()
fig.show()

fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.hist(heights * widths_tanh, bins=100, range=(0, 20), histtype='step', label='Data')
# ax1.plot(bin_midpoints_tanh, three_gaussians(bin_midpoints_tanh, *popt_tanh), label='Fit')
ax1.set_ylabel('Entries', size=20)
ax1.set_xlabel('Voltage (mV)', size=20)
# ax1.set_yscale('log')
# ax1.set_ylim(10**2)
ax1.legend()
fig.show()

fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.hist(heights * widths_sig, bins=100, range=(0, 6), histtype='step', label='Data')
# ax1.plot(bin_midpoints_sig, three_gaussians(bin_midpoints_sig, *popt_sig), label='Fit')
ax1.set_ylabel('Entries', size=20)
ax1.set_xlabel('Voltage (mV)', size=20)
# ax1.set_yscale('log')
# ax1.set_ylim(10**2)
ax1.legend()
fig.show()


# ## 2. Rolling Average Technique

# ### 2.1 importing data

# In[ ]:


import csv

#Peak minus noise
peaks_r1, noise_r1 = peak_minus_noise(voltage_r1, window_size=2000, peak_window_size=10)
peaks_r2, noise_r2 = peak_minus_noise(voltage_r2, window_size=2000, peak_window_size=10)

# Ensure peaks_r1, noise_r1, peaks_r2, and noise_r2 are arrays
peaks_r1 = np.array(peaks_r1)
noise_r1 = np.array(noise_r1)
peaks_r2 = np.array(peaks_r2)
noise_r2 = np.array(noise_r2)

pmn_r1 = peaks_r1 - noise_r1
pmn_r2 = peaks_r2 - noise_r2

# Function to write a list of arrays to a CSV file
def write_arrays_to_csv(arrays, filename):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for array in arrays:
            writer.writerow(array)

# Write the arrays to CSV files
write_arrays_to_csv(pmn_r1, 'pmn_r1.csv')
write_arrays_to_csv(pmn_r2, 'pmn_r2.csv')
write_arrays_to_csv(peaks_r1, 'peaks_r1.csv')
write_arrays_to_csv(peaks_r2, 'peaks_r2.csv')


# In[ ]:


import csv
from lb_functions import *
os.chdir("/vols/hyperk/WCTE/WCTE540/nick/data")
def read_csv_to_list(filename):
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile)
        return [[float(value) for value in row] for row in reader]

pmn_r1 = read_csv_to_list('pmn_r1.csv')
pmn_r2 = read_csv_to_list('pmn_r2.csv')
peaks_r1 = read_csv_to_list('peaks_r1.csv')
peaks_r2 = read_csv_to_list('peaks_r2.csv')


# ### 2.2 Noise studies

# In[ ]:


def write_arrays_to_csv(arrays, filename):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for array in arrays:
            writer.writerow(array)

# Write the arrays to CSV files
write_arrays_to_csv(noise_r1, 'noises_r1_max.csv')


# In[ ]:


import csv
from lb_functions import *
def read_csv_to_list(filename):
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile)
        return [[float(value) for value in row] for row in reader]

noise_max = read_csv_to_list('noises_r1_max.csv')
noise_average_500_2000 = read_csv_to_list('noises_r1_average_500_2000.csv')
noise_mode = read_csv_to_list('noises_r1_mode.csv')
noise_median = read_csv_to_list('noises_r1_median.csv')


# In[ ]:


d = np.min(np.diff(np.sort(np.unique(np.convolve(voltage_r1[0][500:2000], np.ones(10)/10, mode='valid')))))
u = np.sort(np.unique(np.convolve(voltage_r1[0][500:2000], np.ones(10)/10, mode='valid')))
b=(u[1:]+u[:-1])/2
np.append(b,(u[0]-d))
np.append(b,(u[-1]+d))
plt.hist(np.convolve(voltage_r1[0][500:2000], np.ones(10)/10, mode='valid'), bins=b)
plt.xlabel('Voltage (mV)', size=20)
plt.ylabel('Entries', size=20)
plt.show()


# In[ ]:


np.sort(np.unique(noise_max[0]))


# In[ ]:


np.sort(np.unique(noise_average_500_2000[0]))


# In[ ]:


np.unique(noise_mode[0])


# In[ ]:


np.unique(noise_median[0])


# In[ ]:


plt.figure(figsize=(8, 6))  # Increase the figure size
plt.hist(noise_median[0], bins=np.arange(-0.78125, 0, 0.16525))#, bins=np.arange(-1.09375,0.46875,0.15625), label='Noise')
#plt.yscale('log')
plt.xlabel('Voltage (mV)', size=20)
plt.ylabel('Entries', size=20)

# Calculate mean and standard deviation
mean_noise = np.mean(noise_median[0])
std_noise = np.std(noise_median[0])

# Add legend with mean and standard deviation
plt.legend([f'Noise\n$\mu$: {mean_noise:.2f} mV\n$\sigma$: {std_noise:.2f} mV'])
plt.yscale('log')

plt.show()


# ### 2.3 Fitting

# In[ ]:


pmts = ['B', 'C', 'T', 'V']
for i in range(4):
    plt.figure(figsize=(8, 6))  # Increase the figure size
    plt.hist(pmn_r1[i], bins=np.arange(0, 20, 0.15625), histtype='step', label='Run 1')
    plt.hist(pmn_r2[i], bins=np.arange(0, 20, 0.15625), histtype='step', label='Run 2')
    plt.legend()
    plt.title(f'PMT {pmts[i]} '+ 'peak minus noise', size=20)
    plt.xlabel('Voltage (mV)', size=20)
    plt.ylabel('Entries', size=20)
    plt.yscale('log')
    plt.ylim(10**2)
    plt.show()


# In[ ]:


initial_guess_r1 = [
    [4839, 1.5, 0.5, 931.0, 6.5, 1.5, 463.0, 2.5, 0.3, 8.7],
    [5370.0, 1.5, 0.5, 1193.0, 4.5, 1, 525.0, 2.5, 1, 1],
    [4972.0, 1.5, 0.5, 1135.0, 5.1, 1.5, 545.0, 2.5, 1, 1],
    [5432.0, 1.5, 0.5, 1324.0, 3.7, 1, 500, 3, 1, 2],
]

initial_guess_r2 = [
    [4839, 1.5, 0.5, 944.0, 6.1, 1.5, 463.0, 2.5, 0.3, 8.7],
    [5370.0, 1.5, 0.5, 1193.0, 4.5, 1, 525.0, 2.5, 1, 1],
    [4972.0, 1.5, 0.5, 1212.0, 4.5, 1.5, 545.0, 2.5, 1, 1],
    [5432.0, 1.5, 0.5, 1285.0, 3.5, 1, 500, 3, 1, 2],
]

cut = [15, 12.28, 12, 10]

counts_pmn, bin_edges_pmn= np.histogram(pmn_r1[0], bins=np.arange(0, 20, 0.15625), range=(0,20))
bin_midpoints_pmn = (bin_edges_pmn[:-1] + bin_edges_pmn[1:]) / 2

# sigma1_g_pmn = 0.5
# amp2_g_pmn = np.max(counts_pmn[np.where((bin_midpoints_pmn >= 2.5) & (bin_midpoints_pmn <=7.5))[0]])
# mean2_g_pmn = bin_midpoints_pmn[np.where(counts_pmn == amp2_g_pmn)[0][0]]
# sigma2_g_pmn = 1.5
# amp3_g_pmn = np.max(counts_pmn[np.where(bin_midpoints_pmn >= 8)[0]])
# mean3_g_pmn = bin_midpoints_pmn[np.where(counts_pmn == amp3_g_pmn)[0][0]]
# sigma3_g_pmn = 2.5
# exp_decay_g_pmn = 0.3
# exp_scale_g_pmn = 8.7

# guesses_pmn = [[amp1_g_pmn, mean1_g_pmn, sigma1_g_pmn, amp2_g_pmn, mean2_g_pmn, sigma2_g_pmn, amp3_g_pmn, mean3_g_pmn, sigma3_g_pmn, exp_decay_g_pmn, exp_scale_g_pmn],
#                [amp1_g_pmn, mean1_g_pmn, sigma1_g_pmn, amp2_g_pmn, mean2_g_pmn, sigma2_g_pmn, amp3_g_pmn, sigma3_g_pmn, exp_decay_g_pmn, exp_scale_g_pmn],
#                [amp1_g_pmn, mean1_g_pmn, sigma1_g_pmn, amp2_g_pmn, mean2_g_pmn, sigma2_g_pmn, amp3_g_pmn, exp_decay_g_pmn, exp_scale_g_pmn],
#                [amp1_g_pmn, mean1_g_pmn, sigma1_g_pmn, amp2_g_pmn, mean2_g_pmn, sigma2_g_pmn, amp3_g_pmn, mean3_g_pmn, exp_decay_g_pmn, exp_scale_g_pmn],
#                [amp1_g_pmn, mean1_g_pmn, sigma1_g_pmn, amp2_g_pmn, mean2_g_pmn, sigma2_g_pmn, amp3_g_pmn, mean3_g_pmn, sigma3_g_pmn],
#                [amp1_g_pmn, mean1_g_pmn, sigma1_g_pmn, amp2_g_pmn, mean2_g_pmn, sigma2_g_pmn, amp3_g_pmn, sigma3_g_pmn],
#                [amp1_g_pmn, mean1_g_pmn, sigma1_g_pmn, amp2_g_pmn, mean2_g_pmn, sigma2_g_pmn, amp3_g_pmn, mean3_g_pmn],
#                [amp1_g_pmn, mean1_g_pmn, sigma1_g_pmn, amp2_g_pmn, mean2_g_pmn, sigma2_g_pmn, amp3_g_pmn]]


counts_peaks, bin_edges_peaks = np.histogram(peaks_r1[0], bins=np.arange(0, 22, 0.15625), range=(0,22))
bin_midpoints_peaks = (bin_edges_peaks[:-1] + bin_edges_peaks[1:]) / 2

amp1_g_peaks = np.max(counts_peaks[np.where((bin_midpoints_peaks >= 0) & (bin_midpoints_peaks <=1.5))[0]])
mean1_g_peaks = bin_midpoints_peaks[np.where(counts_peaks == amp1_g_peaks)[0][0]]
sigma1_g_peaks = 0.5
amp2_g_peaks = np.max(counts_peaks[np.where((bin_midpoints_peaks >= 2.5) & (bin_midpoints_peaks <=7.5))[0]])
mean2_g_peaks = bin_midpoints_peaks[np.where(counts_peaks == amp2_g_peaks)[0][0]]
sigma2_g_peaks = 1.5
amp3_g_peaks = np.max(counts_peaks[np.where(bin_midpoints_peaks >= 8)[0]])
mean3_g_peaks = bin_midpoints_peaks[np.where(counts_peaks == amp3_g_peaks)[0][0]]
sigma3_g_peaks = 2.5
exp_decay_g_peaks = 0.3
exp_scale_g_peaks = 8.7

guesses_peaks = [[amp1_g_peaks, mean1_g_peaks, sigma1_g_peaks, amp2_g_peaks, mean2_g_peaks, sigma2_g_peaks, amp3_g_peaks, mean3_g_peaks, sigma3_g_peaks, exp_decay_g_peaks, exp_scale_g_peaks],
                 [amp1_g_peaks, mean1_g_peaks, sigma1_g_peaks, amp2_g_peaks, mean2_g_peaks, sigma2_g_peaks, amp3_g_peaks, sigma3_g_peaks, exp_decay_g_peaks, exp_scale_g_peaks],
                 [amp1_g_peaks, mean1_g_peaks, sigma1_g_peaks, amp2_g_peaks, mean2_g_peaks, sigma2_g_peaks, amp3_g_peaks, exp_decay_g_peaks, exp_scale_g_peaks],
                 [amp1_g_peaks, mean1_g_peaks, sigma1_g_peaks, amp2_g_peaks, mean2_g_peaks, sigma2_g_peaks, amp3_g_peaks, mean3_g_peaks, exp_decay_g_peaks, exp_scale_g_peaks],
                 [amp1_g_peaks, mean1_g_peaks, sigma1_g_peaks, amp2_g_peaks, mean2_g_peaks, sigma2_g_peaks, amp3_g_peaks, mean3_g_peaks, sigma3_g_peaks],
                 [amp1_g_peaks, mean1_g_peaks, sigma1_g_peaks, amp2_g_peaks, mean2_g_peaks, sigma2_g_peaks, amp3_g_peaks, sigma3_g_peaks],
                 [amp1_g_peaks, mean1_g_peaks, sigma1_g_peaks, amp2_g_peaks, mean2_g_peaks, sigma2_g_peaks, amp3_g_peaks, mean3_g_peaks],
                 [amp1_g_peaks, mean1_g_peaks, sigma1_g_peaks, amp2_g_peaks, mean2_g_peaks, sigma2_g_peaks, amp3_g_peaks]]

amp1_l_bound = 0
mean1_l_bound = 0
sigma1_l_bound = 0
amp1_u_bound = 10000
mean1_u_bound = 1.5
sigma1_u_bound = 1

amp2_l_bound = 0
mean2_l_bound = 0
sigma2_l_bound = 0
amp2_u_bound = 5000
mean2_u_bound = 8
sigma2_u_bound = 4

amp3_l_bound = 0
mean3_l_bound = 10
sigma3_l_bound = 0
amp3_u_bound = 1000
mean3_u_bound = 15
sigma3_u_bound = 5

amp4_l_bound = 0
mean4_l_bound = 0
sigma4_l_bound = 0
amp4_u_bound = 500
mean4_u_bound = 20
sigma4_u_bound = 5

exp_decay_l_bound = 0
exp_scale_l_bound = 0
exp_decay_u_bound = 10
exp_scale_u_bound = 10000

amp1_g_pmn = (amp1_l_bound + amp1_u_bound) / 2
mean1_g_pmn = (mean1_l_bound + mean1_u_bound) / 2
sigma1_g_pmn = (sigma1_l_bound + sigma1_u_bound) / 2
amp2_g_pmn = (amp2_l_bound + amp2_u_bound) / 2
mean2_g_pmn = (mean2_l_bound + mean2_u_bound) / 2
sigma2_g_pmn = (sigma2_l_bound + sigma2_u_bound) / 2
amp3_g_pmn = (amp3_l_bound + amp3_u_bound) / 2
mean3_g_pmn = (mean3_l_bound + mean3_u_bound) / 2
sigma3_g_pmn = (sigma3_l_bound + sigma3_u_bound) / 2
exp_decay_g_pmn = (exp_decay_l_bound + exp_decay_u_bound) / 2
exp_scale_g_pmn = (exp_scale_l_bound + exp_scale_u_bound) / 2



guesses_pmn = [[amp1_g_pmn, mean1_g_pmn, sigma1_g_pmn, amp2_g_pmn, mean2_g_pmn, sigma2_g_pmn, amp3_g_pmn, mean3_g_pmn, sigma3_g_pmn, exp_decay_g_pmn, exp_scale_g_pmn],
               [amp1_g_pmn, mean1_g_pmn, sigma1_g_pmn, amp2_g_pmn, mean2_g_pmn, sigma2_g_pmn, amp3_g_pmn, sigma3_g_pmn, exp_decay_g_pmn, exp_scale_g_pmn],
               [amp1_g_pmn, mean1_g_pmn, sigma1_g_pmn, amp2_g_pmn, mean2_g_pmn, sigma2_g_pmn, amp3_g_pmn, exp_decay_g_pmn, exp_scale_g_pmn],
               [amp1_g_pmn, mean1_g_pmn, sigma1_g_pmn, amp2_g_pmn, mean2_g_pmn, sigma2_g_pmn, amp3_g_pmn, mean3_g_pmn, exp_decay_g_pmn, exp_scale_g_pmn],
               [amp1_g_pmn, mean1_g_pmn, sigma1_g_pmn, amp2_g_pmn, mean2_g_pmn, sigma2_g_pmn, amp3_g_pmn, mean3_g_pmn, sigma3_g_pmn],
               [amp1_g_pmn, mean1_g_pmn, sigma1_g_pmn, amp2_g_pmn, mean2_g_pmn, sigma2_g_pmn, amp3_g_pmn, sigma3_g_pmn],
               [amp1_g_pmn, mean1_g_pmn, sigma1_g_pmn, amp2_g_pmn, mean2_g_pmn, sigma2_g_pmn, amp3_g_pmn, mean3_g_pmn],
               [amp1_g_pmn, mean1_g_pmn, sigma1_g_pmn, amp2_g_pmn, mean2_g_pmn, sigma2_g_pmn, amp3_g_pmn]]


bounds_pmn = [
([amp1_l_bound, mean1_l_bound, sigma1_l_bound, amp2_l_bound, mean2_l_bound, sigma2_l_bound, amp3_l_bound, mean3_l_bound, sigma3_l_bound, exp_decay_l_bound, exp_scale_l_bound], 
     [amp1_u_bound, mean1_u_bound, sigma1_u_bound, amp2_u_bound, mean2_u_bound, sigma2_u_bound, amp3_u_bound, mean3_u_bound, sigma3_u_bound, exp_decay_u_bound, exp_scale_u_bound]),
([amp1_l_bound, mean1_l_bound, sigma1_l_bound, amp2_l_bound, mean2_l_bound, sigma2_l_bound, amp3_l_bound, sigma3_l_bound, exp_decay_l_bound, exp_scale_l_bound], 
    [amp1_u_bound, mean1_u_bound, sigma1_u_bound, amp2_u_bound, mean2_u_bound, sigma2_u_bound, amp3_u_bound, sigma3_u_bound, exp_decay_u_bound, exp_scale_u_bound]),
([amp1_l_bound, mean1_l_bound, sigma1_l_bound, amp2_l_bound, mean2_l_bound, sigma2_l_bound, amp3_l_bound, exp_decay_l_bound, exp_scale_l_bound], 
    [amp1_u_bound, mean1_u_bound, sigma1_u_bound, amp2_u_bound, mean2_u_bound, sigma2_u_bound, amp3_u_bound, exp_decay_u_bound, exp_scale_u_bound]),
([amp1_l_bound, mean1_l_bound, sigma1_l_bound, amp2_l_bound, mean2_l_bound, sigma2_l_bound, amp3_l_bound, mean3_l_bound, exp_decay_l_bound, exp_scale_l_bound], 
    [amp1_u_bound, mean1_u_bound, sigma1_u_bound, amp2_u_bound, mean2_u_bound, sigma2_u_bound, amp3_u_bound, mean3_u_bound, exp_decay_u_bound, exp_scale_u_bound]),
([amp1_l_bound, mean1_l_bound, sigma1_l_bound, amp2_l_bound, mean2_l_bound, sigma2_l_bound, amp3_l_bound, mean3_l_bound, sigma3_l_bound], 
    [amp1_u_bound, mean1_u_bound, sigma1_u_bound, amp2_u_bound, mean2_u_bound, sigma2_u_bound, amp3_u_bound, mean3_u_bound, sigma3_u_bound]),
([amp1_l_bound, mean1_l_bound, sigma1_l_bound, amp2_l_bound, mean2_l_bound, sigma2_l_bound, amp3_l_bound, sigma3_l_bound], 
    [amp1_u_bound, mean1_u_bound, sigma1_u_bound, amp2_u_bound, mean2_u_bound, sigma2_u_bound, amp3_u_bound, sigma3_u_bound]),
([amp1_l_bound, mean1_l_bound, sigma1_l_bound, amp2_l_bound, mean2_l_bound, sigma2_l_bound, amp3_l_bound, mean3_l_bound], 
    [amp1_u_bound, mean1_u_bound, sigma1_u_bound, amp2_u_bound, mean2_u_bound, sigma2_u_bound, amp3_u_bound, mean3_u_bound]),
([amp1_l_bound, mean1_l_bound, sigma1_l_bound, amp2_l_bound, mean2_l_bound, sigma2_l_bound, amp3_l_bound], 
    [amp1_u_bound, mean1_u_bound, sigma1_u_bound, amp2_u_bound, mean2_u_bound, sigma2_u_bound, amp3_u_bound])
]

bounds_41 = ([amp1_l_bound, mean1_l_bound, sigma1_l_bound, amp2_l_bound, mean2_l_bound, sigma2_l_bound, amp3_l_bound, amp4_l_bound, exp_decay_l_bound, exp_scale_l_bound], 
             [amp1_u_bound, mean1_u_bound, sigma1_u_bound, amp2_u_bound, mean2_u_bound, sigma2_u_bound, amp3_u_bound, amp4_u_bound, exp_decay_u_bound, exp_scale_u_bound])

bounds_42 = ([amp1_l_bound, mean1_l_bound, sigma1_l_bound, amp2_l_bound, mean2_l_bound, sigma2_l_bound, amp3_l_bound, mean3_l_bound, sigma3_l_bound, amp4_l_bound, mean4_l_bound, sigma4_l_bound, exp_decay_l_bound, exp_scale_l_bound], 
             [amp1_u_bound, mean1_u_bound, sigma1_u_bound, amp2_u_bound, mean2_u_bound, sigma2_u_bound, amp3_u_bound, mean3_u_bound, sigma3_u_bound, amp4_u_bound, mean4_u_bound, sigma4_u_bound, exp_decay_u_bound, exp_scale_u_bound])

functions = [
    multi_gaussian_with_exp_mu3,
    multi_gaussian_with_exp,
    multi_gaussian_with_exp_cons,
    multi_gaussian_with_exp_mu3_cons,
    multi_gaussian_mu3,
    multi_gaussian,
    multi_gaussian_mu3_cons,
    multi_gaussian_cons
]


# In[ ]:


counts_peaks, bin_edges_peaks = np.histogram(pmn_r1[0], bins=np.arange(0, 20, 0.15625), range=(0,20))
bin_midpoints_peaks = (bin_edges_peaks[:-1] + bin_edges_peaks[1:]) / 2
b = np.max(counts_peaks[np.where((bin_midpoints_peaks >= 0) & (bin_midpoints_peaks <=2.5))[0]])
a = bin_midpoints_peaks[np.where(counts_peaks == b)[0][0]]

print(a)
print(b)


# In[ ]:


guesses = [
    [2200, 1.5, 0.5, 931.0, 6.5, 1.5, 463.0, 2.5, 0.3, 150, 15, 2.5, 0.3, 8.7],
    [2000, 1, 0.5, 900, 4.9, 1, 400, 8, 3, 100, 18, 4, 0.3, 100],
    [2000, 1, 0.5, 900, 4.9, 1, 400, 10, 3, 100, 18, 4, 0.4, 100],
    [2000, 1, 0.5, 900, 4.9, 1, 400, 10, 3, 100, 18, 4, 0.4, 100]
]

bounds = [
    [(0, None), (0, 2), (0, 10), (0, None), (5, 7), (0, 10), (0, None), (10, 15), (0, 10), (0, None), (17, 20), (0,10), (0, 1), (0, None) ],
    #[(0, 10000), (0, 1.2), (0, 1), (0, 10000), (5, 7.5), (0, 5), (0, None), (10, 15), (0, 10), (0, None), (15, 20), (0, 10), (0, None), (0, 100)],
    [(0, None), (0, 1.5), (0, 10), (0, None), (2, 6), (0, 10), (0, None), (6, 11), (0, 10), (0, None), (10, 20), (0, 10), (0, None), (0, None)],
    [(0, None), (0, 2), (0, 1), (0, None), (3, 6), (0, 10), (0, None), (6, 11), (0, 10), (0, None), (12.5, 20), (0, 10), (0, None), (0, None)],
    [(0, None), (0, 2), (0, 5), (0, None), (3, 5), (0, 10), (0, None), (6, 8), (0, 10), (0, None), (10, 20), (0, 10), (0, None), (0, None)]
]

labels = ['Peak-noise PMT B', 'Peak-noise PMT C', 'Peak-noise PMT T', 'Peak-noise PMT V']

pmt_labels = ['B', 'C', 'T', 'V']

maximum_pmn = 20

for i in range(4):
    _, cutoff, _ = efficiency_histogram(pmn_r1[i], guesses[i], pmt_label=pmt_labels[i], fit_function=four_gaussian_with_exp, range_=[0, maximum_pmn], label=labels[i], plot=True, bounds=bounds[i])
    print(cutoff)


# In[ ]:


import csv
from lb_functions import *
os.chdir("/vols/hyperk/WCTE/WCTE540/nick/data")
def read_csv_to_list(filename):
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile)
        return [[float(value) for value in row] for row in reader]

pmn_r1 = read_csv_to_list('pmn_r1.csv')
pmn_r2 = read_csv_to_list('pmn_r2.csv')
peaks_r1 = read_csv_to_list('peaks_r1.csv')
peaks_r2 = read_csv_to_list('peaks_r2.csv')


# In[ ]:


guesses = [
    [2200, 1.5, 0.5, 931.0, 6.5, 1.5, 463.0, 1, 1, 150, 1, 1, 0.3, 8.7],
    [2000, 1, 0.5, 900, 4.9, 1.5, 400, 1, 1, 100, 1, 1, 0.3, 100],
    [2000, 1, 0.5, 900, 4.9, 1.5, 400, 1, 1, 100, 1, 1, 0.4, 100],
    [2000, 1, 0.5, 900, 4.9, 1.5, 400, 1, 1, 100, 1, 1, 0.4, 100]
]

bounds = [
    [(0, None), (0, 2), (0, 10), (0, None), (5, 7), (0, 10), (0, None), (0.5, 2), (0.5, 2), (0, None), (0.5, 2), (0.5, 2), (0, 1), (0, None) ],
    [(0, None), (0, 1.5), (0, 10), (0, None), (2, 6), (0, 10), (0, None), (0.5, 2), (0.5, 2), (0, None), (0.5, 2), (0.5, 2), (0, None), (0, None)],
    [(0, None), (0, 2), (0, 10), (0, None), (3, 6), (0, 10), (0, None), (0.5, 2), (0.5, 2), (0, None), (0.5, 2), (0.5, 2), (0, None), (0, None)],
    [(0, None), (0, 2), (0, 10), (0, None), (3, 5), (0, 10), (0, None), (0.5, 2), (0.5, 2), (0, None), (0.5, 2), (0.5, 2), (0, None), (0, None)]
]
step_sizes = [10, 0.1, 0.01, 10, 0.1, 0.01, 10, 0.1, 0.1, 10, 0.1, 0.1, 0.1, 10]

labels = ['PMT B Peak minus Noise', 'PMT C Peak minus Noise', 'PMT T Peak minus Noise', 'PMT V Peak minus Noise']

pmt_labels = ['B', 'C', 'T', 'V']

maximum_pmn = 20

eff_r1 = []
eff_r2 = []
for i in range(4):
    _, cutoff, efficiency_r1 = efficiency_histogram(pmn_r1[i], guesses[i], pmt_label=pmt_labels[i], fit_function=nick_multi_gaus, range_=[0, maximum_pmn], label=labels[i], plot=True, bounds=bounds[i], step_size=step_sizes)
    _, cutoff, efficiency_r2 = efficiency_histogram(pmn_r2[i], guesses[i], pmt_label=pmt_labels[i], fit_function=nick_multi_gaus, range_=[0, maximum_pmn], label=labels[i], plot=True, bounds=bounds[i], step_size=step_sizes)
    eff_r1.append(efficiency_r1)
    eff_r2.append(efficiency_r2)


# In[ ]:


import numpy as np

np.array(eff_r2) / eff_r2[-1]


# In[ ]:


import numpy as np
from iminuit import Minuit
import matplotlib.pyplot as plt

# Define your four_gaussian_with_exp function here

maximum_pmn = 20
maximum_peak = 22
amp4_g_peaks = 150
mean4_g_pmn = 15
sigma4_g_pmn = 2.5

# Assuming pmn_r1[0] is your histogram data
data = pmn_r1[0]
bin_edges = np.arange(0, 20, 0.15625)  # Adjust binning as needed
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])  # Compute bin centers
bin_contents, _ = np.histogram(data, bins=bin_edges)  # Get histogram counts
errors = np.sqrt(bin_contents)  # Poisson uncertainty (sqrt(counts))
errors = np.where(errors == 0, 1, errors)  # Replace zero uncertainties with minimum value
bin_width = bin_edges[1] - bin_edges[0]

# Chi2 function for Minuit
def chi2(amp1, mean1, stddev1, amp2, mean2, stddev2,
         amp3, mean3, stddev3, amp4, mean4, stddev4,
         exp_decay, exp_scale):
    expected = four_gaussian_with_exp(bin_centers, amp1, mean1, stddev1, amp2, mean2, stddev2,
                                      amp3, mean3, stddev3, amp4, mean4, stddev4,
                                      exp_decay, exp_scale)
    return np.sum(((bin_contents - expected) / errors) ** 2)

# Initial guesses
guess_42 = [amp1_g_pmn, mean1_g_pmn, sigma1_g_pmn,
            amp2_g_pmn, mean2_g_pmn, sigma2_g_pmn,
            amp3_g_pmn, mean3_g_pmn, sigma3_g_pmn,
            amp4_g_peaks, mean4_g_pmn, sigma4_g_pmn,
            exp_decay_g_pmn, exp_scale_g_pmn]

# Pass the initial guesses into Minuit
m = Minuit(chi2, *guess_42)
m.errordef = Minuit.LEAST_SQUARES  # Chi2 fitting

# Set parameter limits before fitting
m.limits = [

(0, None), # amp1 >= 0

(0, 2.5), # mean1 within x-range

(0, 10), # stddev1 > 0

(0, None), # amp2 >= 0

(5, 7), # mean2 within x-range

(0, 10), # stddev2 > 0

(0, None), # amp3 >= 0

(10, 15), # mean3 within x-range

(0, 10), # stddev3 > 0

(0, None), # amp4 >= 0

(17, 20), # mean4 within x-range

(0, 10), # stddev4 > 0

(0, None), # exp_decay > 0

(0, None) # exp_scale > 0

]

m.migrad()  # Perform the fit

# Output the results
print("Fit Results:")
print(m.values)  # Fitted parameter values
print("Parameter Errors:")
print(m.errors)  # Parameter uncertainties

# Extract best-fit parameters
fit_params = m.values
cutoff = fit_params['mean2'] - 1 * fit_params['stddev2']
kept_data = np.array(data)[np.array(data) > cutoff]
efficiency = len(kept_data) / len(data) * 100  # Percentage of kept data
print(f"Efficiency: {efficiency:.2f}%")



# Generate the fitted curve
x_fit = np.linspace(bin_edges[0], bin_edges[-1], 1000)
y_fit = four_gaussian_with_exp(x_fit, *fit_params)
y_fit_scaled = y_fit * bin_width

# Plot the histogram and fitted function
plt.figure(figsize=(8, 6))
plt.plot(x_fit, y_fit, label="Fit", color="red", linewidth=2)
plt.hist(data, bins=bin_edges, histtype='step', label="Histogram", color="blue")
plt.axvline(cutoff, color='green', linestyle='--', label='Threshold')
plt.xlabel("Voltage (mV)", size=20)
plt.ylabel("Counts", size=20)
plt.title("PMT B Peak Minus Noise", size=20)
plt.legend(fontsize=15)
plt.yscale('log')
plt.ylim(10**2)
plt.grid(True)

# Add efficiency text to the plot
plt.text(0.95, 0.95, f'Efficiency: {efficiency:.2f}%', horizontalalignment='right', verticalalignment='top', transform=plt.gca().transAxes, fontsize=15, bbox=dict(facecolor='white', alpha=0.5))

plt.show()


# In[ ]:


import numpy as np
from iminuit import Minuit
import matplotlib.pyplot as plt

# Define your four_gaussian_with_exp function here

maximum_pmn = 20
mean4_g_pmn = 15
sigma4_g_pmn = 2.5

# Assuming pmn_r1[0] is your histogram data
data = pmn_r2[1]
bin_edges = np.arange(0, 20, 0.15625)  # Adjust binning as needed
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])  # Compute bin centers
bin_contents, _ = np.histogram(data, bins=bin_edges)  # Get histogram counts
errors = np.sqrt(bin_contents)  # Poisson uncertainty (sqrt(counts))
errors = np.where(errors == 0, 1, errors)  # Replace zero uncertainties with minimum value
bin_width = bin_edges[1] - bin_edges[0]

# Chi2 function for Minuit
def chi2(amp1, mean1, stddev1, amp2, mean2, stddev2,
         amp3, mean3, stddev3, amp4, mean4, stddev4,
         exp_decay, exp_scale):
    expected = four_gaussian_with_exp(bin_centers, amp1, mean1, stddev1, amp2, mean2, stddev2,
                                      amp3, mean3, stddev3, amp4, mean4, stddev4,
                                      exp_decay, exp_scale)
    return np.sum(((bin_contents - expected) / errors) ** 2)

# Initial guesses
guess_42 = [2000, 1, 0.5, #[5370.0, 1.5, 0.5, 1193.0, 4.5, 1, 525.0, 2.5, 1, 1]
            900, 4.9, 1,
            400, 8, 3,
            100, 18, 4,
            0.3, 100]

# Pass the initial guesses into Minuit
m = Minuit(chi2, *guess_42)
m.errordef = Minuit.LEAST_SQUARES  # Chi2 fitting

# Set parameter limits before fitting
m.limits = [
    (0, None),   # amp1 >= 0
    (0, 1.5),     # mean1 within x-range
    (0, 10),     # stddev1 > 0
    (0, None),   # amp2 >= 0
    (2, 6),     # mean2 within x-range
    (0, 10),     # stddev2 > 0
    (0, None),   # amp3 >= 0
    (6, 11),     # mean3 within x-range
    (0, 10),     # stddev3 > 0
    (0, None),   # amp4 >= 0
    (10, 20),     # mean4 within x-range
    (0, 10),     # stddev4 > 0
    (0, None),   # exp_decay > 0
    (0, None)    # exp_scale > 0
]

m.migrad()  # Perform the fit

# Output the results
print("Fit Results:")
print(m.values)  # Fitted parameter values
print("Parameter Errors:")
print(m.errors)  # Parameter uncertainties

# Extract best-fit parameters
fit_params = m.values
cutoff = fit_params['mean2'] - 1 * fit_params['stddev2']
kept_data = np.array(data)[np.array(data) > cutoff]
efficiency = len(kept_data) / len(data) * 100  # Percentage of kept data
print(f"Efficiency: {efficiency:.2f}%")



# Generate the fitted curve
x_fit = np.linspace(bin_edges[0], bin_edges[-1], 1000)
y_fit = four_gaussian_with_exp(x_fit, *fit_params)
y_fit_scaled = y_fit * bin_width

# Plot the histogram and fitted function
plt.figure(figsize=(8, 6))
plt.plot(x_fit, y_fit, label="Fit", color="red", linewidth=2)
plt.hist(data, bins=bin_edges, histtype='step', label="Histogram", color="blue")
plt.axvline(cutoff, color='green', linestyle='--', label='Threshold')
plt.xlabel("Voltage (mV)", size=20)
plt.ylabel("Counts", size=20)
plt.title("PMT C Peak Minus Noise", size=20)
plt.yscale('log')
plt.ylim(10**2)
plt.legend(fontsize=15)
plt.grid(True)

# Add efficiency text to the plot
plt.text(0.95, 0.95, f'Efficiency: {efficiency:.2f}%', horizontalalignment='right', verticalalignment='top', transform=plt.gca().transAxes, fontsize=15, bbox=dict(facecolor='white', alpha=0.5))

plt.show()


# In[ ]:


import numpy as np
from iminuit import Minuit
import matplotlib.pyplot as plt

# Define your four_gaussian_with_exp function here

maximum_pmn = 20
mean4_g_pmn = 15
sigma4_g_pmn = 2.5

# Assuming pmn_r1[0] is your histogram data
data = pmn_r2[2]
bin_edges = np.arange(0, 20, 0.15625)  # Adjust binning as needed
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])  # Compute bin centers
bin_contents, _ = np.histogram(data, bins=bin_edges)  # Get histogram counts
errors = np.sqrt(bin_contents)  # Poisson uncertainty (sqrt(counts))
errors = np.where(errors == 0, 1, errors)  # Replace zero uncertainties with minimum value
bin_width = bin_edges[1] - bin_edges[0]

# Chi2 function for Minuit
def chi2(amp1, mean1, stddev1, amp2, mean2, stddev2,
         amp3, mean3, stddev3, amp4, mean4, stddev4,
         exp_decay, exp_scale):
    expected = four_gaussian_with_exp(bin_centers, amp1, mean1, stddev1, amp2, mean2, stddev2,
                                      amp3, mean3, stddev3, amp4, mean4, stddev4,
                                      exp_decay, exp_scale)
    return np.sum(((bin_contents - expected) / errors) ** 2)

# Initial guesses
guess_42 = [2000, 1, 0.5, #[5370.0, 1.5, 0.5, 1193.0, 4.5, 1, 525.0, 2.5, 1, 1]
            900, 4.9, 1,
            400, 10, 3,
            100, 18, 4,
            0.4, 100]

# Pass the initial guesses into Minuit
m = Minuit(chi2, *guess_42)
m.errordef = Minuit.LEAST_SQUARES  # Chi2 fitting

# Set parameter limits before fitting
m.limits = [
    (0, None),   # amp1 >= 0
    (0, 2),     # mean1 within x-range
    (0, 1),     # stddev1 > 0
    (0, None),   # amp2 >= 0
    (3, 6),     # mean2 within x-range
    (0, 10),     # stddev2 > 0
    (0, None),   # amp3 >= 0
    (6, 11),     # mean3 within x-range
    (0, 10),     # stddev3 > 0
    (0, None),   # amp4 >= 0
    (12.5, 20),     # mean4 within x-range
    (0, 10),     # stddev4 > 0
    (0, None),   # exp_decay > 0
    (0, None)    # exp_scale > 0
]

m.migrad()  # Perform the fit

# Output the results
print("Fit Results:")
print(m.values)  # Fitted parameter values
print("Parameter Errors:")
print(m.errors)  # Parameter uncertainties

# Extract best-fit parameters
fit_params = m.values
cutoff = fit_params['mean2'] - 1 * fit_params['stddev2']
kept_data = np.array(data)[np.array(data) > cutoff]
efficiency = len(kept_data) / len(data) * 100  # Percentage of kept data
print(f"Efficiency: {efficiency:.2f}%")



# Generate the fitted curve
x_fit = np.linspace(bin_edges[0], bin_edges[-1], 1000)
y_fit = four_gaussian_with_exp(x_fit, *fit_params)
y_fit_scaled = y_fit * bin_width

# Plot the histogram and fitted function
plt.figure(figsize=(8, 6))
plt.plot(x_fit, y_fit, label="Fit", color="red", linewidth=2)
plt.hist(data, bins=bin_edges, histtype='step', label="Histogram", color="blue")
plt.axvline(cutoff, color='green', linestyle='--', label='Threshold')
plt.xlabel("Voltage (mV)", size=20)
plt.ylabel("Counts", size=20)
plt.title("PMT T Peak Minus Noise", size=20)
plt.yscale('log')
plt.ylim(10**2)
plt.legend(fontsize=15)
plt.grid(True)

# Add efficiency text to the plot
plt.text(0.95, 0.95, f'Efficiency: {efficiency:.2f}%', horizontalalignment='right', verticalalignment='top', transform=plt.gca().transAxes, fontsize=15, bbox=dict(facecolor='white', alpha=0.5))

plt.show()


# In[ ]:


import numpy as np
from iminuit import Minuit
import matplotlib.pyplot as plt

# Define your four_gaussian_with_exp function here

maximum_pmn = 20
mean4_g_pmn = 15
sigma4_g_pmn = 2.5

# Assuming pmn_r1[0] is your histogram data
data = pmn_r2[3]
bin_edges = np.arange(0, 20, 0.15625)  # Adjust binning as needed
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])  # Compute bin centers
bin_contents, _ = np.histogram(data, bins=bin_edges)  # Get histogram counts
errors = np.sqrt(bin_contents)  # Poisson uncertainty (sqrt(counts))
errors = np.where(errors == 0, 1, errors)  # Replace zero uncertainties with minimum value
bin_width = bin_edges[1] - bin_edges[0]

# Chi2 function for Minuit
def chi2(amp1, mean1, stddev1, amp2, mean2, stddev2,
         amp3, mean3, stddev3, amp4, mean4, stddev4,
         exp_decay, exp_scale):
    expected = four_gaussian_with_exp(bin_centers, amp1, mean1, stddev1, amp2, mean2, stddev2,
                                      amp3, mean3, stddev3, amp4, mean4, stddev4,
                                      exp_decay, exp_scale)
    return np.sum(((bin_contents - expected) / errors) ** 2)

# Initial guesses
guess_42 = [2000, 1, 0.5, #[5370.0, 1.5, 0.5, 1193.0, 4.5, 1, 525.0, 2.5, 1, 1]
            900, 4.9, 1,
            400, 10, 3,
            100, 18, 4,
            0.4, 100]

# Pass the initial guesses into Minuit
m = Minuit(chi2, *guess_42)
m.errordef = Minuit.LEAST_SQUARES  # Chi2 fitting

# Set parameter limits before fitting
m.limits = [
    (0, None),   # amp1 >= 0
    (0, 2),     # mean1 within x-range
    (0, 5),     # stddev1 > 0
    (0, None),   # amp2 >= 0
    (3, 5),     # mean2 within x-range
    (0, 10),     # stddev2 > 0
    (0, None),   # amp3 >= 0
    (6, 8),     # mean3 within x-range
    (0, 10),     # stddev3 > 0
    (0, None),   # amp4 >= 0
    (10, 20),     # mean4 within x-range
    (0, 10),     # stddev4 > 0
    (0, None),   # exp_decay > 0
    (0, None)    # exp_scale > 0
]

m.migrad()  # Perform the fit

# Output the results
print("Fit Results:")
print(m.values)  # Fitted parameter values
print("Parameter Errors:")
print(m.errors)  # Parameter uncertainties

# Extract best-fit parameters
fit_params = m.values
cutoff = fit_params['mean2'] - 1 * fit_params['stddev2']
kept_data = np.array(data)[np.array(data) > cutoff]
efficiency = len(kept_data) / len(data) * 100  # Percentage of kept data
print(f"Efficiency: {efficiency:.2f}%")



# Generate the fitted curve
x_fit = np.linspace(bin_edges[0], bin_edges[-1], 1000)
y_fit = four_gaussian_with_exp(x_fit, *fit_params)
y_fit_scaled = y_fit * bin_width

# Plot the histogram and fitted function
plt.figure(figsize=(8, 6))
plt.plot(x_fit, y_fit, label="Fit", color="red", linewidth=2)
plt.hist(data, bins=bin_edges, histtype='step', label="Histogram", color="blue")
plt.axvline(cutoff, color='green', linestyle='--', label='Threshold')
plt.xlabel("Voltage (mV)", size=20)
plt.ylabel("Counts", size=20)
plt.title("PMT V Peak Minus Noise", size=20)
plt.yscale('log')
plt.ylim(10**2)
plt.legend(fontsize=15)
plt.grid(True)

# Add efficiency text to the plot
plt.text(0.95, 0.95, f'Efficiency: {efficiency:.2f}%', horizontalalignment='right', verticalalignment='top', transform=plt.gca().transAxes, fontsize=15, bbox=dict(facecolor='white', alpha=0.5))

plt.show()


# In[ ]:


peak_ratios_r1, _ = calculate_ratios(voltage_r1, window_size=2000, peak_window_size=10, use_rolling_average=True, return_ratios=True, return_indices=False)
#peak_ratios_r2 = calculate_ratios(voltage_r2, window_size=2000, peak_window_size=10, use_rolling_average=True, return_ratios=False, return_indices=False)


# In[ ]:


#Plotting peak-noise
for i in range(4):
    fig = plt.figure(figsize=(8,6)) 
    ax = plt.gca()
    ax.hist(peak_ratios_r1[i],  bins=np.arange(0, 50, 0.5), range=(0,50), label='Run 1')
    #ax.hist(peak_ratios_r2[i],  bins=np.arange(0, 20, 0.15625), range=(0,20), alpha=0.5, label='Run 2')
    ax.set_xlabel('(Peak-noise)/noise ', size=25)
    ax.set_ylabel('Entries', size=25)
    ax.set_title('Peak (peak- noise)/noise: PMT '+pmts[i], size=25)
    ax.legend(prop={'size': 15})
    ax.set_ylim(5*10**2)
    ax.set_yscale('log')
    plt.show()


# In[ ]:


#Noise code
window_size = 2000
peak_window_size = 10
noises = []
ranges = [[0,100], [500,2000], [1000, 2000], [1500,2000]]

# Loop over the dataset in 2000-sample windows
for test_pmt in voltage_r2:
    fin = []
    for j in range(len(ranges)):

        noise = [] #testing the noise

        # Loop over the dataset in 2000-sample windows
        for i in range(0, len(test_pmt), window_size):
            # Define the 2000-sample window
            window_data = test_pmt[i:i + window_size]

            if len(window_data) < window_size:  # Skip if the window has fewer than 2000 points (end of data)
                continue

            # Noise data (fixed section of the window)
            noise_data = window_data[ranges[j][0]: ranges[j][1]]

            # Calculate the average of the noise
            noise_average = np.average(noise_data)
            noise.append(noise_average)

        fin.append(noise)
    noises.append(fin)


# In[ ]:


#Plotting noise

means = [np.mean(noise) for noise in noises]
rms_values = [np.sqrt(np.mean(np.array(noise)**2)) for noise in noises]

# Plotting noise histograms
fig = plt.figure()   
ax = plt.gca()
ax.set_title('Noise data for PMT V')
labels = ['0-100', '500-2000', '1000-2000', '1500-2000']
ran = (-1, 2)
for i in range(len(noises[0])):
    ax.hist(noises[0][i], bins=50, range=ran, label=labels[i], histtype='step')

# Adding the legend
ax.legend()
ax.set_xlabel('Noise data mV', size=20)
ax.set_ylabel('Entries', size=20)

# Adding the mean and RMS box in the bottom right corner
textstr = '\n'.join((
    r'$\mu=%.2f,\ \mu_{RMS}=%.2f$' % (means[0], rms_values[0]),
    r'$\mu=%.2f,\ \mu_{RMS}=%.2f$' % (means[1], rms_values[1]),
    r'$\mu=%.2f,\ \mu_{RMS}=%.2f$' % (means[2], rms_values[2]),
    r'$\mu=%.2f,\ \mu_{RMS}=%.2f$' % (means[3], rms_values[3])
))

# Text box properties
props = dict(boxstyle='round', facecolor='white', alpha=0.5)

# Place the box in the bottom right corner
ax.text(0.1, 0.5, textstr, transform=ax.transAxes, fontsize=11,
        verticalalignment='bottom', horizontalalignment='left', bbox=props)

# Show the plot
plt.show()


# In[ ]:


#Waveforms where noise = 0

window_size = 2000
peak_window_size = 10

indices_noise_zero = []  # To store the ratios for each dataset

# Loop over the dataset in 2000-sample windows
for test_pmt in voltage_r1:

    ratios2 = []  # Reset ratios for the current dataset

    # Loop over the dataset in 2000-sample windows
    for i in range(0, len(test_pmt), window_size):
        # Define the 2000-sample window
        window_data = test_pmt[i:i + window_size]

        if len(window_data) < window_size:  # Skip if the window has fewer than 2000 points (end of data)
            continue

        # Calculate rolling average
        rolling_average = np.convolve(window_data, np.ones(peak_window_size)/peak_window_size, mode='valid')
        peak_window_start = np.argmax(rolling_average)
        test.append(peak_window_start)

        # Noise data (fixed section of the window)
        noise_data = window_data[0: 100]
        # Calculate the average of the noise
        noise_average = np.average(noise_data)

        # Calculate the ratio
        if noise_average != 0:  # Avoid division by zero
            pass
        else:
            index = i + peak_window_start
        ratios2.append(ratio)

    ratios_list.append(ratios2)  # Store the ratios for the current dataset


# In[ ]:


#Counting how many NAN values are in list for when noise is zero

for i in range(4):
    contains_nan = np.any(np.isnan(ratios_list[i]))

    if contains_nan:
        print("Array contains NaN values.")
        nan_count = np.sum(np.isnan(ratios_list[i]))
        print(f"The array contains {nan_count} NaN values.")
    else:
        print("Array does not contain any NaN values.")


