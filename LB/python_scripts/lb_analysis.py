#!/usr/bin/env python
# coding: utf-8

# In[2]:


from lb_functions import *
# Path to the folder where your ROOT files are located
data_folder = "/vols/hyperk/WCTE/WCTE540/nick/data"

# Get a list of all ROOT files in the data folder
root_files = [f for f in os.listdir(data_folder) if f.endswith(".root")]


# In[3]:


initial_guess_r1 = [
    [4839, 1.5, 0.5, 931.0, 6.5, 1.5, 463.0, 2.5, 0.3, 8.7],
    [5370.0, 1.5, 0.5, 1193.0, 4.5, 1, 525.0, 2.5, 1, 1],
    [4972.0, 1.5, 0.5, 1135.0, 5.1, 1.5, 545.0, 2.5, 1, 1],
    [5432.0, 1.5, 0.5, 1324.0, 3.7, 1, 500, 3, 1, 2],
]

cuts = [15, 12.28, 12, 10]


# In[ ]:


import re
import numpy as np

def analyze_efficiencies(root_files, keywords, plot_title='', pmt_V=False, fit_plots=False, plot_final=False):
    # Filter and sort files based on keywords
    file = sorted([s for s in root_files if all(keyword in s for keyword in keywords)])

    # Compile the regular expression pattern
    pattern_angle = re.compile(r'_([+-]?\d+)deg_')
    pattern_PMT = re.compile(r'_([BCT])V_')

    # Initialize lists for angles and data
    angles = []
    efficiencies = []

    # Extract angles and data in a single loop
    for filename in file:
        match_angle = pattern_angle.search(filename)
        match_PMT = pattern_PMT.search(filename)

        if match_angle and match_PMT:
            angles.append(float(match_angle.group(1)))
            pmt_type = match_PMT.group(1)

            if pmt_type == 'B':
                initial_guess = [4839, 1.5, 0.5, 931.0, 6.5, 1.5, 463.0, 2.5, 0.3, 8.7]  
                cut = 15
                pmt_label = 'B'
                pmt1, = extract_pmt_voltages(filename, pmt1=True, chunk_size=2000000)
                pmn = peak_minus_noise([pmt1], window_size=2000, peak_window_size=10)
                _, cutoff, efficiency = efficiency_histogram(pmn[0], initial_guess=initial_guess, cut=cut, pmt_label=pmt_label, bins=100, range_=[0, 20], label=plot_title, plot=fit_plots)
                efficiencies.append(efficiency)

            elif pmt_type == 'C':
                initial_guess = [5370.0, 1.5, 0.5, 1193.0, 4.5, 1, 525.0, 2.5, 1, 1] 
                cut = 12.28
                pmt_label = 'C'
                pmt1, = extract_pmt_voltages(filename, pmt1=True, chunk_size=2000000)
                pmn = peak_minus_noise([pmt1], window_size=2000, peak_window_size=10)
                _, cutoff, efficiency = efficiency_histogram(pmn[0], initial_guess=initial_guess, cut=cut, pmt_label=pmt_label, bins=100, range_=[0, 20], label=plot_title, plot=fit_plots)
                efficiencies.append(efficiency)

            elif pmt_type == 'T':
                initial_guess = [4972.0, 1.5, 0.5, 1135.0, 5.1, 1.5, 545.0, 2.5, 1, 1] 
                cut = 12
                pmt_label = 'T'
                pmt1, = extract_pmt_voltages(filename, pmt1=True, chunk_size=2000000)
                pmn = peak_minus_noise([pmt1], window_size=2000, peak_window_size=10)
                _, cutoff, efficiency = efficiency_histogram(pmn[0], initial_guess=initial_guess, cut=cut, pmt_label=pmt_label, bins=100, range_=[0, 20], label=plot_title, plot=fit_plots)
                efficiencies.append(efficiency)

    if plot_final:        
        plt.plot(angles, efficiencies, 'o')
        plt.xlabel('Angle ($\circ$)', size=15)
        plt.ylabel('Efficiency', size=15)
        plt.title(plot_title, size=15)
            # if pmt_V:
            #     pmt1, = extract_pmt_voltages(filename, pmt1=False, pmt2=True, chunk_size=2000000)
            #     _, cutoff, efficiency = efficiency_histogram(pmt1, initial_guess, cut, pmt_label, bins=100, range_=[0, 20], label='', plot=True)
            # data_list.append(pmt1)
    return angles, efficiencies


# In[5]:


keywords = ['LDON', '2mmSpacer', 'V']
test =  sorted([s for s in root_files if all(keyword in s for keyword in keywords)])
print(test)


# In[6]:


pmts = ['BV', 'CV', 'TV']
p = ['B', 'C', 'T']
keywords = ['LDON', '3mmSpacer', 'BV']
pmt_dict = {}

for i, pmt in enumerate(pmts):
    keywords = ['LDON', '2mmSpacer', pmt]
    angles, efficiencies = analyze_efficiencies(root_files, keywords, pmt_V=False, fit_plots=False, plot_final=False)
    pmt_dict[pmt] = (angles, efficiencies)
    plt.plot(angles, efficiencies, 'o', label=f'PMT {p[i]}')
    plt.xlabel('Angle (°)', size=15)
    plt.ylabel('Efficiency', size=15)
    plt.title('Efficiency vs Angle for PMT '+p[i], size=15)
    plt.legend()
    plt.show()

