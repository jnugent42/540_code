#!/usr/bin/env python
# coding: utf-8

# In[3]:


import uproot
import pathlib
import numpy as np
import matplotlib.pyplot as plt
import glob
from lb_functions import *


# In[3]:


def get_files(prefix, pmts, phi_deg):
    files = {}
    for l in ['LDON']:
        filename = prefix+f"{phi_deg}deg_{pmts}_{l}_*.root"
        filenames = glob.glob(filename)
        if len(filenames) == 0:
            if phi_deg > 0:
                filename = prefix+f"+{phi_deg}deg_{pmts}_{l}_*.root"
                filenames = glob.glob(filename)
        if len(filenames) == 0:
            raise FileNotFoundError(filename+" was not found.")
        elif len(filenames) > 1:
            raise FileNotFoundError(filename+" matched multiple files.")
        else:
            files[l] = filenames[0]
    return files


# In[4]:


def get_rate(filename, channel, plot=False):
    guesses = [
    [3500, 1.5, 0.5, 931.0, 6.5, 1.5, 463.0, 1, 1, 150, 1, 1, 0.3, 8.7],
    [3500, 1, 0.5, 900, 4.9, 1.5, 400, 1, 1, 100, 1, 1, 0.3, 100],
    [3500, 1, 0.5, 900, 4.9, 1.5, 400, 1, 1, 100, 1, 1, 0.4, 100],
    [3500, 1, 0.5, 900, 4.9, 1.5, 400, 1, 1, 100, 1, 1, 0.4, 100]
    ]

    bounds = [
    [(0, None), (0, 2), (0, 10), (0, None), (5, 7), (0, 10), (0, None), (0.5, 2), (0.5, 2), (0, None), (0.5, 2), (0.5, 2), (0, 1), (0, None) ],
    [(0, None), (0, 1.5), (0, 10), (0, None), (2, 6), (0, 10), (0, None), (0.5, 2), (0.5, 2), (0, None), (0.5, 2), (0.5, 2), (0, None), (0, None)],
    [(0, None), (0, 2), (0, 10), (0, None), (3, 6), (0, 10), (0, None), (0.5, 2), (0.5, 2), (0, None), (0.5, 2), (0.5, 2), (0, None), (0, None)],
    [(0, None), (0, 2), (0, 10), (0, None), (3, 5), (0, 10), (0, None), (0.5, 2), (0.5, 2), (0, None), (0.5, 2), (0.5, 2), (0, None), (0, None)]
    ]
    if channel == 1:
        pattern_PMT = re.compile(r'_([BCT])V_')
        match_PMT = pattern_PMT.search(filename)
        pmt_type = match_PMT.group(1)

        if pmt_type == 'B':
            initial_guess = guesses[0]  
            bounds = bounds[0]
            step_sizes = [10, 0.1, 0.01, 10, 0.1, 0.01, 10, 0.1, 0.1, 10, 0.1, 0.1, 0.1, 10]
            pmt_label = 'B'
            pmt1, = extract_pmt_voltages(filename, pmt1=True, chunk_size=2000000)
            peak, noise = peak_minus_noise([pmt1], window_size=2000, peak_window_size=10, noise_start=500, noise_end=2000)
            pmn = [np.array(peaks) - np.array(noise) for peaks, noise in zip(peak, noise)]
            _, cutoff, efficiency = efficiency_histogram(pmn[0], initial_guess=initial_guess, fit_function=nick_multi_gaus, pmt_label=pmt_label, range_=[0, 20], label='', plot=plot, bounds=bounds, step_size=step_sizes)
            return efficiency

        elif pmt_type == 'C':
            initial_guess = guesses[1]
            bounds = bounds[1]
            step_sizes = [10, 0.1, 0.01, 10, 0.1, 0.01, 10, 0.1, 0.1, 10, 0.1, 0.1, 0.1, 10]
            pmt_label = 'C'
            pmt1, = extract_pmt_voltages(filename, pmt1=True, chunk_size=2000000)
            peak, noise = peak_minus_noise([pmt1], window_size=2000, peak_window_size=10, noise_start=500, noise_end=2000)
            pmn = [np.array(peaks) - np.array(noise) for peaks, noise in zip(peak, noise)]
            _, cutoff, efficiency = efficiency_histogram(pmn[0], initial_guess=initial_guess, fit_function=nick_multi_gaus, pmt_label=pmt_label, range_=[0, 20], label='', plot=plot, bounds=bounds, step_size=step_sizes)
            return efficiency

        elif pmt_type == 'T':
            initial_guess = guesses[2]
            bounds = bounds[2]
            step_sizes = [10, 0.1, 0.01, 10, 0.1, 0.01, 10, 0.1, 0.1, 10, 0.1, 0.1, 0.1, 10]
            pmt_label = 'T'
            pmt1, = extract_pmt_voltages(filename, pmt1=True, chunk_size=2000000)
            peak, noise = peak_minus_noise([pmt1], window_size=2000, peak_window_size=10, noise_start=500, noise_end=2000)
            pmn = [np.array(peaks) - np.array(noise) for peaks, noise in zip(peak, noise)]
            _, cutoff, efficiency = efficiency_histogram(pmn[0], initial_guess=initial_guess, fit_function=nick_multi_gaus, pmt_label=pmt_label, range_=[0, 20], plot=plot, bounds=bounds, step_size=step_sizes)
            return efficiency

    elif channel == 2:
        pmt2, = extract_pmt_voltages(filename, pmt1=False, pmt2=True, chunk_size=2000000)
        pmn = peak_minus_noise([pmt2], window_size=2000, peak_window_size=10, noise_start=500, noise_end=2000)
        initial_guess = guesses[3]
        bounds = bounds[3]
        step_sizes = [10, 0.1, 0.01, 10, 0.1, 0.01, 10, 0.1, 0.1, 10, 0.1, 0.1, 0.1, 10]
        pmt_label = 'V'
        peak, noise = peak_minus_noise([pmt2], window_size=2000, peak_window_size=10, noise_start=500, noise_end=2000)
        pmn = [np.array(peaks) - np.array(noise) for peaks, noise in zip(peak, noise)]
        _, cutoff, efficiency = efficiency_histogram(pmn[0], initial_guess=initial_guess, fit_function=nick_multi_gaus, pmt_label=pmt_label, range_=[0, 20], plot=plot, bounds=bounds, step_size=step_sizes)
        return efficiency


# In[5]:


def get_average_pulse_integral(filename, channel):
    os.chdir("/vols/hyperk/WCTE/WCTE540/nick/Analysis/landau")
    filename ='Landau_Output_'+filename
    file = uproot.open(pathlib.Path(filename))
    sum_charges = file[f"PMT{channel}_Charge"].member("fTsumwx")
    windows = file[f"PMT{channel}_Corrected"].member("fEntries")
    os.chdir("/vols/hyperk/WCTE/WCTE540/nick/data")
    return sum_charges/windows


# In[6]:


eff_corrections = { # From Kieren
    'T': (0.97009099+0.971746971)/2,
    'B': (0.95720313+0.95437616)/2,
    'V': 1.0,
    'C': (1.00553977+0.98350209)/2,
}
distances = { # mm
    'T': 183+15,
    'B': 183+15,
    'V': 185+15,
    'C': 187+15,
}
opening = 9/2 # mm
for p, d in distances.items():
    opening_angle = 2*np.pi*(1-d/np.sqrt(opening*opening+d*d))
    eff_corrections[p] *= opening_angle
v = eff_corrections['V']
for p in eff_corrections.keys():
    eff_corrections[p] /= v
    print(p, eff_corrections[p])


# In[7]:


pmt_channels = {
    'B': 1,
    'C': 1,
    'T': 1,
    'V': 2,
}
colors = {'T': 'purple', 'C': 'orange', 'B': 'green', 'V': 'cyan'}


# # 0 rotations

# In[46]:


title="0 rotations of fine adjust"


# In[ ]:


os.chdir("/vols/hyperk/WCTE/WCTE540/nick/data")

prefix="Ball2_176.02mA_"
rates = {}
phis = np.arange(-180, 135, 45)
for phi in phis:
    rates[phi] = {}
    for pmts in ['BV', 'CV', 'TV']:
        rates[phi][pmts] = {}
        for n, f in get_files(prefix, pmts, phi).items():
            rates[phi][pmts][n] = {}
            for p in pmts:
                rates[phi][pmts][n][p] = get_rate(f, pmt_channels[p], plot=True)
                print(phi, pmts, n, p, rates[phi][pmts][n][p])
            rates[phi][pmts][n]['M'] = get_average_pulse_integral(f, 3)
            print(phi, pmts, n, 'M', rates[phi][pmts][n]['M'])


# In[26]:


# -90 CV LDOFF C is 0.04068, 10x larger than expected, due to mistake of starting before laser fully powered down...
# replace with the LDOFFa value for this run instead, for now. Later, could reprocess only later captures
rates[-90]['CV']['LDOFF']['C'] = rates[-90]['CV']['LDOFFa']['C']


# In[47]:


phis = np.arange(-180, 135, 45)
with open('0rot_data.pkl', 'rb') as f:
    rates = pickle.load(f)


# In[48]:


fig, ax = plt.subplots(figsize=(6, 5))
all_x = []
all_v = []
all_v_dr1 = []
all_v_dr2 = []
correction_gradient = {}
correction_intercept = {}
for pmts in ['TV', 'CV', 'BV']:
    x = [rates[phi][pmts]['LDON']['M'] for phi in phis]
    all_x += x
    all_v += [rates[phi][pmts]['LDON']['V'] for phi in phis]
    #all_v_dr1 += [rates[phi][pmts]['LDOFF']['V'] for phi in phis]
    #all_v_dr2 += [rates[phi][pmts]['LDOFFa']['V'] for phi in phis]
    x = np.array(x)
    y = [rates[phi][pmts]['LDON'][pmts[0]] for phi in phis]
   # y_dr1 = [rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis]
   # y_dr2 = [rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis]
   # y_dr1 = np.array(y_dr1)
    #y_dr2 = np.array(y_dr2)
    #y = y - 0.5*(y_dr1+y_dr2)
    ax.scatter(x, y, label=f"PMT {pmts[0]}", c=colors[pmts[0]])
    a, b = np.polyfit(x, y, 1)
    correction_gradient[pmts[0]] = a
    correction_intercept[pmts[0]] = b
    ax.plot(x, a*x+b, c=colors[pmts[0]])
    print(pmts[0], a, b)
all_x = np.array(all_x)
all_v = np.array(all_v)
all_v_dr1 = np.array(all_v_dr1)
all_v_dr2 = np.array(all_v_dr2)
all_v = all_v #- 0.5*(all_v_dr1 + all_v_dr2)
ax.scatter(all_x, all_v, label=f"PMT V", c=colors['V'])
a, b = np.polyfit(all_x, all_v, 1)
print("V", a, b)
correction_gradient['V'] = a
correction_intercept['V'] = b
ax.plot(all_x, a*all_x+b, c=colors['V'])
ax.legend()
ax.set_xlabel("Average monitor PMT pulse integral")
ax.set_ylabel("Main PMT rate")


# In[49]:


for pmts in ['TV', 'CV', 'BV']:
    v_pmt_rates = np.array([rates[phi][pmts]['LDON']['V'] for phi in phis])
    # v_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF']['V'] for phi in phis])
    # v_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa']['V'] for phi in phis])
    # v_pmt_dr = (v_pmt_dr1+v_pmt_dr2)/2
    o_pmt_rates = np.array([rates[phi][pmts]['LDON'][pmts[0]] for phi in phis])
    # o_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis])
    # o_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis])
    # o_pmt_dr = (o_pmt_dr1+o_pmt_dr2)/2
    m_pmt_q = np.array([rates[phi][pmts]['LDON']['M'] for phi in phis])
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(phis, o_pmt_rates, '.-', c='blue', label=f"{pmts[0]} PMT raw rate")
    # ax.plot(phis, o_pmt_dr1, '.--', c='blue', label=f"{pmts[0]} PMT DR 1")
    # ax.plot(phis, o_pmt_dr1, '.-.', c='blue', label=f"{pmts[0]} PMT DR 2")
    ax.plot(phis, v_pmt_rates, '.-', c='orange', label=f"V PMT raw rate")
    # ax.plot(phis, v_pmt_dr1, '.--', c='orange', label=f"V PMT DR 1")
    # ax.plot(phis, v_pmt_dr2, '.-.', c='orange', label=f"V PMT DR 1")
    #ax.set_xlim(-45, 45)
    ax.set_ylim(50,80)
    ax.set_xlabel("Ball phi [deg]")
    ax.set_ylabel("PMT Count Rate (%)")
    ax.yaxis.label.set_color('blue')
    # ax2 = ax.twinx()
    # ax2.plot(phis, m_pmt_q, '.--', c='black', label=f"Monitor")
    # ax2.set_ylim(1.2, 1.4)
    # ax2.set_ylabel("Monitor average pulse integral")
    ax.legend()
    # ax2.legend(loc="lower right")


# In[50]:


fig, ax = plt.subplots(figsize=(6, 5))
for o in 'TCB':
    pmts = o+'V'
    v_pmt_rates = np.array([rates[phi][pmts]['LDON']['V'] for phi in phis])
    #v_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF']['V'] for phi in phis])
    #v_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa']['V'] for phi in phis])
    #v_pmt_dr = (v_pmt_dr1+v_pmt_dr2)/2
    ax.plot(phis, v_pmt_rates, '.-', label=f"V PMT rate while measuring with {o} PMT")
#ax.set_xlim(-45, 45)
ax.set_ylim(50, 80)
ax.set_xlabel("Ball phi [deg]")
ax.set_ylabel("PMT Count Rate (%)")
ax.legend()


# In[51]:


v_pmt_rates = np.array([rates[phi][pmts]['LDON']['V'] for phi in phis])
# v_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF']['V'] for phi in phis])
# v_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa']['V'] for phi in phis])
# v_pmt_dr = (v_pmt_dr1+v_pmt_dr2)/2
o_pmt_rates = np.array([rates[phi][pmts]['LDON'][pmts[0]] for phi in phis])
# o_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis])
# o_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis])
# o_pmt_dr = (o_pmt_dr1+o_pmt_dr2)/2
m_pmt_q = np.array([rates[phi][pmts]['LDON']['M'] for phi in phis])
v_corrected = v_pmt_rates - correction_gradient['V']*(m_pmt_q - 1.275)
o_corrected = o_pmt_rates - correction_gradient[pmts[0]]*(m_pmt_q - 1.275)


# In[52]:


fig, ax = plt.subplots(1, 3, figsize=(12, 5))
for i, pmts in enumerate(['TV', 'CV', 'BV']):

    ax[i].plot(phis, o_pmt_rates, '.:', c='blue', label=f"{pmts[0]} PMT raw rate")
    ax[i].plot(phis, o_corrected, '.-', c='blue', label=f"{pmts[0]} PMT corrected")
    ax[i].plot(phis, v_pmt_rates, '.:', c='orange', label=f"V PMT raw rate")
    ax[i].plot(phis, v_corrected, '.-', c='orange', label=f"V PMT corrected")
    #ax.set_xlim(-45, 45)
    ax[i].set_ylim(50,80)
    ax[i].set_xlabel("Ball phi [deg]")
    ax[i].set_ylabel("PMT Count Rate (%)")
    ax[i].yaxis.label.set_color('blue')
    #ax2 = ax[i].twinx()
    #ax2.plot(phis, m_pmt_q, '.--', c='black', label=f"Monitor")
    #ax2.set_ylim(1.2, 1.6)
    #ax2.set_ylabel("Monitor average pulse integral")
    ax[i].legend()
    #ax2.legend(loc="lower right")
fig.suptitle(title)
plt.tight_layout()


# In[53]:


fig, ax = plt.subplots(figsize=(6, 5))
all_x = []
all_v = []
all_v_dr1 = []
all_v_dr2 = []
for pmts in ['TV', 'CV', 'BV']:
    x = [rates[phi][pmts]['LDON']['M'] for phi in phis]
    all_x += x
    all_v += [rates[phi][pmts]['LDON']['V'] for phi in phis]
    # all_v_dr1 += [rates[phi][pmts]['LDOFF']['V'] for phi in phis]
    # all_v_dr2 += [rates[phi][pmts]['LDOFFa']['V'] for phi in phis]
    x = np.array(x)
    y = [rates[phi][pmts]['LDON'][pmts[0]] for phi in phis]
    # y_dr1 = [rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis]
    # y_dr2 = [rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis]
    y = y - correction_gradient[pmts[0]]*(x-1.275)
    ax.scatter(x, y, label=f"PMT {pmts[0]}", c=colors[pmts[0]])
    a, b = np.polyfit(x, y, 1)
    ax.plot(x, a*x+b, c=colors[pmts[0]])
    print(pmts[0], a, b)
all_x = np.array(all_x)
all_v = np.array(all_v)
all_v = all_v - correction_gradient['V']*(all_x-1.275)
ax.scatter(all_x, all_v, label=f"PMT V", c=colors['V'])
a, b = np.polyfit(all_x, all_v, 1)
print("V", a, b)
ax.plot(all_x, a*all_x+b, c=colors['V'])
ax.legend()
ax.set_xlabel("Average monitor PMT pulse integral")
ax.set_ylabel("Main PMT rate")


# In[54]:


raw_rates = []
corr_rates = []

fig, ax = plt.subplots(figsize=(6, 5))
for pmts in ['TV', 'CV', 'BV']:
    m = [rates[phi][pmts]['LDON']['M'] for phi in phis]
    v = [rates[phi][pmts]['LDON']['V'] for phi in phis]
    # v_dr1 = [rates[phi][pmts]['LDOFF']['V'] for phi in phis]
    # v_dr2 = [rates[phi][pmts]['LDOFFa']['V'] for phi in phis]
    y = [rates[phi][pmts]['LDON'][pmts[0]] for phi in phis]
    # y_dr1 = [rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis]
    # y_dr2 = [rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis]
    m = np.array(m)
    v = np.array(v)
    # v_dr1 = np.array(v_dr1)
    # v_dr2 = np.array(v_dr2)
    y = np.array(y)
    # y_dr1 = np.array(y_dr1)
    # y_dr2 = np.array(y_dr2)
    y_ = y - correction_gradient[pmts[0]]*(m-1.275)
    v_ = v - correction_gradient['V']*(m-1.275)
    raw_rates.append((np.max(y), np.min(y), np.max(v), np.min(v)))
    corr_rates.append((np.max(y_), np.min(y_), np.max(v_), np.min(v_)))
    print('For PMT ' + pmts + ' the spread in raw count rate is: ' + str(np.max(y)-np.min(y)))
    print('For PMT V of ' + pmts + ' the spread raw in count rate is: ' + str(np.max(v)-np.min(v)))
    print('For PMT ' + pmts + ' the spread in corrected count rate is: ' + str(np.max(y_)-np.min(y_)))
    print('For PMT V of ' + pmts + ' the spread corrected in count rate is: ' + str(np.max(v_)-np.min(v_)))
    ax.plot(phis, y, '.-', label=f"PMT {pmts[0]}", c=colors[pmts[0]])
    ax.plot(phis, v, '.:', label=f"PMT V ({pmts[0]})", c=colors[pmts[0]])
ax.legend(ncol=3)
ax.set_ylim([50,80])
ax.set_xlabel("Ball phi [deg]")
ax.set_ylabel("Corrected PMT rate")
ax.set_title("Corrected by Solid angle & monitor PMT, not relative PMT efficiencies")


# In[55]:


print("The max raw count rate difference is:", np.max(raw_rates)-np.min(raw_rates))
print("The max corrected count rate difference is:", np.max(corr_rates)-np.min(corr_rates))
# print("The average raw count rate difference is:", np.average(np.abs(np.diff(raw_rates))))
# print("The average corrected count rate difference is:", np.average(np.abs(np.diff(corr_rates))))


# In[56]:


rel_rates = []
fig, ax = plt.subplots(figsize=(6, 5))
for pmts in ['TV', 'CV', 'BV']:
    m = [rates[phi][pmts]['LDON']['M'] for phi in phis]
    v = [rates[phi][pmts]['LDON']['V'] for phi in phis]
    # v_dr1 = [rates[phi][pmts]['LDOFF']['V'] for phi in phis]
    # v_dr2 = [rates[phi][pmts]['LDOFFa']['V'] for phi in phis]
    y = [rates[phi][pmts]['LDON'][pmts[0]] for phi in phis]
    # y_dr1 = [rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis]
    # y_dr2 = [rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis]
    m = np.array(m)
    v = np.array(v)
    # v_dr1 = np.array(v_dr1)
    # v_dr2 = np.array(v_dr2)
    y = np.array(y)
    # y_dr1 = np.array(y_dr1)
    # y_dr2 = np.array(y_dr2)
    y = y/eff_corrections[pmts[0]] - correction_gradient[pmts[0]]*(m-1.275)
    v = v  - correction_gradient['V']*(m-1.275)
    ax.plot(phis, (y/v)*100, '.-', label=f"PMT {pmts[0]}", c=colors[pmts[0]])
    rel_rates.append((np.max((y/v)*100), np.min((y/v)*100)))
ax.legend(ncol=3)
ax.set_xlabel("Ball phi [deg]")
ax.set_ylabel("Corrected PMT rate / corrected V PMT rate")
ax.set_title("0 rotations of fine adjust")
ax.set_ylim(80, 100)


# In[57]:


print("The max relative count rate difference is:", np.max(rel_rates)-np.min(rel_rates))
print("The average relative count rate difference is:", 100-np.min(rel_rates))


# # 1 rotation

# In[71]:


title="1 rotation of fine adjust"


# In[78]:


os.chdir("/vols/hyperk/WCTE/WCTE540/nick/data")
prefix="Ball2_1rot_176.02mA_"

rates = {}
phis = np.arange(-180, 180, 90)
for phi in phis:
    rates[phi] = {}
    for pmts in ['BV', 'CV', 'TV']:
        rates[phi][pmts] = {}
        for n, f in get_files(prefix, pmts, phi).items():
            rates[phi][pmts][n] = {}
            for p in pmts:
                rates[phi][pmts][n][p] = get_rate(f, pmt_channels[p], plot=True)
                print(phi, pmts, n, p, rates[phi][pmts][n][p])
            rates[phi][pmts][n]['M'] = get_average_pulse_integral(f, 3)
            print(phi, pmts, n, 'M', rates[phi][pmts][n]['M'])


# In[72]:


with open('1rot_data.pkl', 'rb') as f:
    rates = pickle.load(f)
phis = np.arange(-180, 180, 90)


# In[73]:


fig, ax = plt.subplots(figsize=(6, 5))
all_x = []
all_v = []
all_v_dr1 = []
all_v_dr2 = []
correction_gradient = {}
correction_intercept = {}
for pmts in ['BV', 'CV', 'TV']:
    x = [rates[phi][pmts]['LDON']['M'] for phi in phis]
    all_x += x
    all_v += [rates[phi][pmts]['LDON']['V'] for phi in phis]
    #all_v_dr1 += [rates[phi][pmts]['LDOFF']['V'] for phi in phis]
    #all_v_dr2 += [rates[phi][pmts]['LDOFFa']['V'] for phi in phis]
    x = np.array(x)
    y = [rates[phi][pmts]['LDON'][pmts[0]] for phi in phis]
    #y_dr1 = [rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis]
    #y_dr2 = [rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis]
    #y_dr1 = np.array(y_dr1)
    #y_dr2 = np.array(y_dr2)
    #y = y - 0.5*(y_dr1+y_dr2)
    ax.scatter(x, y, label=f"PMT {pmts[0]}", c=colors[pmts[0]])
    a, b = np.polyfit(x, y, 1)
    correction_gradient[pmts[0]] = a
    correction_intercept[pmts[0]] = b
    ax.plot(x, a*x+b, c=colors[pmts[0]])
    print(pmts[0], a, b)
all_x = np.array(all_x)
all_v = np.array(all_v)
#all_v_dr1 = np.array(all_v_dr1)
#all_v_dr2 = np.array(all_v_dr2)
all_v = all_v #- 0.5*(all_v_dr1 + all_v_dr2)
ax.scatter(all_x, all_v, label=f"PMT V", c=colors['V'])
a, b = np.polyfit(all_x, all_v, 1)
print("V", a, b)
correction_gradient['V'] = a
correction_intercept['V'] = b
ax.plot(all_x, a*all_x+b, c=colors['V'])
ax.legend()
ax.set_xlabel("Average monitor PMT pulse integral")
ax.set_ylabel("PMT Count Rate (%)")


# In[74]:


for pmts in ['BV', 'CV', 'TV']:
  v_pmt_rates = np.array([rates[phi][pmts]['LDON']['V'] for phi in phis])
  #v_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF']['V'] for phi in phis])
  #v_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa']['V'] for phi in phis])
  #v_pmt_dr = (v_pmt_dr1+v_pmt_dr2)/2
  o_pmt_rates = np.array([rates[phi][pmts]['LDON'][pmts[0]] for phi in phis])
  #o_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis])
  #o_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis])
  #o_pmt_dr = (o_pmt_dr1+o_pmt_dr2)/2
  m_pmt_q = np.array([rates[phi][pmts]['LDON']['M'] for phi in phis])
  #v_corrected = v_pmt_rates - v_pmt_dr - correction_gradient['V']*(m_pmt_q - 1.275)
  v_corrected = v_pmt_rates - correction_gradient['V']*(m_pmt_q - 1.275)
  #o_corrected = o_pmt_rates - o_pmt_dr - correction_gradient[pmts[0]]*(m_pmt_q - 1.275)
  o_corrected = o_pmt_rates - correction_gradient[pmts[0]]*(m_pmt_q - 1.275)
  fig, ax = plt.subplots(figsize=(6, 5))
  ax.plot(phis, o_pmt_rates, '.-', c='blue', label=f"{pmts[0]} PMT Efficiency")
  # ax.plot(phis, o_pmt_dr1, '.--', c='blue', label=f"{pmts[0]} PMT DR 1")
  # ax.plot(phis, o_pmt_dr1, '.-.', c='blue', label=f"{pmts[0]} PMT DR 2")
  ax.plot(phis, v_pmt_rates, '.-', c='orange', label=f"V PMT Efficiency")
  #  ax.plot(phis, v_pmt_dr1, '.--', c='orange', label=f"V PMT DR 1")
  #  ax.plot(phis, v_pmt_dr2, '.-.', c='orange', label=f"V PMT DR 1")
  #ax.set_xlim(-45, 45)
  ax.set_ylim(50, 80)
  ax.set_xlabel("Ball phi [deg]")
  ax.set_ylabel("PMT Count Rate (%)")
  ax.yaxis.label.set_color('blue')
  ax2 = ax.twinx()
  ax2.plot(phis, m_pmt_q, '.--', c='black', label=f"Monitor")
  #ax2.set_ylim(1.2, 1.4)
  ax2.set_ylabel("Monitor average pulse integral")
  ax.legend(loc="upper right")
  #ax2.legend(loc="lower right")


# In[75]:


fig, ax = plt.subplots(figsize=(6, 5))
for o in 'BCT':
    pmts = o+'V'
    v_pmt_rates = np.array([rates[phi][pmts]['LDON']['V'] for phi in phis])
    #v_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF']['V'] for phi in phis])
    #v_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa']['V'] for phi in phis])
    #v_pmt_dr = (v_pmt_dr1+v_pmt_dr2)/2
    ax.plot(phis, v_pmt_rates, '.-', label=f"V PMT rate while measuring with {o} PMT")
#ax.set_xlim(-45, 45)
ax.set_ylim(50, 80)
ax.set_xlabel("Ball phi [deg]")
ax.set_ylabel("PMT Count Rate (%)")
ax.legend()


# In[76]:


fig, ax = plt.subplots(1, 3, figsize=(12, 5))
for i, pmts in enumerate(['TV', 'CV', 'BV']):
    v_pmt_rates = np.array([rates[phi][pmts]['LDON']['V'] for phi in phis])
    #v_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF']['V'] for phi in phis])
    #v_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa']['V'] for phi in phis])
    #v_pmt_dr = (v_pmt_dr1+v_pmt_dr2)/2
    o_pmt_rates = np.array([rates[phi][pmts]['LDON'][pmts[0]] for phi in phis])
    #o_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis])
    #o_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis])
    #o_pmt_dr = (o_pmt_dr1+o_pmt_dr2)/2
    m_pmt_q = np.array([rates[phi][pmts]['LDON']['M'] for phi in phis])
    #v_corrected = v_pmt_rates - v_pmt_dr - correction_gradient['V']*(m_pmt_q - 1.275)
    v_corrected = v_pmt_rates  - correction_gradient['V']*(m_pmt_q - 1.275)
    #o_corrected = o_pmt_rates - o_pmt_dr - correction_gradient[pmts[0]]*(m_pmt_q - 1.275)
    o_corrected = o_pmt_rates - correction_gradient[pmts[0]]*(m_pmt_q - 1.275)
    ax[i].plot(phis, o_pmt_rates, '.:', c='blue', label=f"{pmts[0]} PMT raw rate")
    ax[i].plot(phis, o_corrected, '.-', c='blue', label=f"{pmts[0]} PMT corrected")
    ax[i].plot(phis, v_pmt_rates, '.:', c='orange', label=f"V PMT raw rate")
    ax[i].plot(phis, v_corrected, '.-', c='orange', label=f"V PMT corrected")
    #ax.set_xlim(-45, 45)
    ax[i].set_ylim(50,80)
    ax[i].set_xlabel("Ball phi [deg]")
    ax[i].set_ylabel("PMT Count Rate (%)")
    ax[i].yaxis.label.set_color('blue')
    #ax2 = ax[i].twinx()
    #ax2.plot(phis, m_pmt_q, '.--', c='black', label=f"Monitor")
    # ax2.set_ylim(1.2, 1.6)
    #ax2.set_ylabel("Monitor average pulse integral")
    ax[i].legend()
    #ax2.legend(loc="lower right")
fig.suptitle(title)
plt.tight_layout()


# In[77]:


fig, ax = plt.subplots(figsize=(6, 5))
all_x = []
all_v = []
all_v_dr1 = []
all_v_dr2 = []
for pmts in ['BV', 'CV', 'TV']:
    x = [rates[phi][pmts]['LDON']['M'] for phi in phis]
    all_x += x
    all_v += [rates[phi][pmts]['LDON']['V'] for phi in phis]
    #all_v_dr1 += [rates[phi][pmts]['LDOFF']['V'] for phi in phis]
    #all_v_dr2 += [rates[phi][pmts]['LDOFFa']['V'] for phi in phis]
    x = np.array(x)
    y = [rates[phi][pmts]['LDON'][pmts[0]] for phi in phis]
    #y_dr1 = [rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis]
    #y_dr2 = [rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis]
    #y_dr1 = np.array(y_dr1)
    #y_dr2 = np.array(y_dr2)
    #y = y - 0.5*(y_dr1+y_dr2) - correction_gradient[pmts[0]]*(x-1.275)
    y = y - correction_gradient[pmts[0]]*(x-1.275)
    ax.scatter(x, y, label=f"PMT {pmts[0]}", c=colors[pmts[0]])
    a, b = np.polyfit(x, y, 1)
    ax.plot(x, a*x+b, c=colors[pmts[0]])
    print(pmts[0], a, b)
all_x = np.array(all_x)
all_v = np.array(all_v)
#all_v_dr1 = np.array(all_v_dr1)
#all_v_dr2 = np.array(all_v_dr2)
#all_v = all_v - 0.5*(all_v_dr1 + all_v_dr2) - correction_gradient['V']*(all_x-1.275)
all_v = all_v - correction_gradient['V']*(all_x-1.275)
ax.scatter(all_x, all_v, label=f"PMT V", c=colors['V'])
a, b = np.polyfit(all_x, all_v, 1)
print("V", a, b)
ax.plot(all_x, a*all_x+b, c=colors['V'])
ax.legend()
ax.set_xlabel("Average monitor PMT pulse integral")
ax.set_ylabel("PMT Count Rate (%)")


# In[78]:


fig, ax = plt.subplots(figsize=(6, 5))

raw_rates = []
corr_rates = []

for pmts in ['BV', 'CV', 'TV']:
    m = [rates[phi][pmts]['LDON']['M'] for phi in phis]
    v = [rates[phi][pmts]['LDON']['V'] for phi in phis]
    # v_dr1 = [rates[phi][pmts]['LDOFF']['V'] for phi in phis]
    # v_dr2 = [rates[phi][pmts]['LDOFFa']['V'] for phi in phis]
    y = [rates[phi][pmts]['LDON'][pmts[0]] for phi in phis]
    # y_dr1 = [rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis]
    # y_dr2 = [rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis]
    m = np.array(m)
    v = np.array(v)
    #v_dr1 = np.array(v_dr1)
    #v_dr2 = np.array(v_dr2)
    y = np.array(y)
    #y_dr1 = np.array(y_dr1)
   # y_dr2 = np.array(y_dr2)
    # y = y - 0.5*(y_dr1+y_dr2) - correction_gradient[pmts[0]]*(m-1.275)
    y_ = y - correction_gradient[pmts[0]]*(m-1.275)
    # v = v - 0.5*(v_dr1+v_dr2) - correction_gradient['V']*(m-1.275)
    v_ = v - correction_gradient['V']*(m-1.275)
    ax.plot(phis, y_, '.-', label=f"PMT {pmts[0]}", c=colors[pmts[0]])
    ax.plot(phis, v_, '.:', label=f"PMT V ({pmts[0]})", c=colors[pmts[0]])
    raw_rates.append((np.max(y), np.min(y), np.max(v), np.min(v)))
    corr_rates.append((np.max(y_), np.min(y_), np.max(v_), np.min(v_)))
    print('For PMT ' + pmts + ' the spread in raw count rate is: ' + str(np.max(y)-np.min(y)))
    print('For PMT V of ' + pmts + ' the spread raw in count rate is: ' + str(np.max(v)-np.min(v)))
    print('For PMT ' + pmts + ' the spread in corrected count rate is: ' + str(np.max(y_)-np.min(y_)))
    print('For PMT V of ' + pmts + ' the spread corrected in count rate is: ' + str(np.max(v_)-np.min(v_)))


ax.legend(ncol=3)
ax.set_ylim([50,80])
ax.set_xlabel("Ball phi [deg]")
ax.set_ylabel("Corrected PMT count rate")
ax.set_title("Corrected by Solid angle & monitor PMT, not relative PMT efficiencies")


# In[79]:


print("The max raw count rate difference is:", np.max(raw_rates)-np.min(raw_rates))
print("The max corrected count rate difference is:", np.max(corr_rates)-np.min(corr_rates))


# In[80]:


rel_rates = []

fig, ax = plt.subplots(figsize=(6, 5))
for pmts in ['BV', 'CV', 'TV']:
    m = [rates[phi][pmts]['LDON']['M'] for phi in phis]
    v = [rates[phi][pmts]['LDON']['V'] for phi in phis]
    # v_dr1 = [rates[phi][pmts]['LDOFF']['V'] for phi in phis]
    # v_dr2 = [rates[phi][pmts]['LDOFFa']['V'] for phi in phis]
    y = [rates[phi][pmts]['LDON'][pmts[0]] for phi in phis]
    # y_dr1 = [rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis]
    # y_dr2 = [rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis]
    m = np.array(m)
    v = np.array(v)
    # v_dr1 = np.array(v_dr1)
    # v_dr2 = np.array(v_dr2)
    y = np.array(y)
    #y_dr1 = np.array(y_dr1)
    #y_dr2 = np.array(y_dr2)
    #y = y/eff_corrections[pmts[0]] - 0.5*(y_dr1+y_dr2) - correction_gradient[pmts[0]]*(m-1.275)
    y = y/eff_corrections[pmts[0]] - correction_gradient[pmts[0]]*(m-1.275)
    #v = v - 0.5*(v_dr1+v_dr2) - correction_gradient['V']*(m-1.275)
    v = v  - correction_gradient['V']*(m-1.275)
    ax.plot(phis, (y/v)*100, '.-', label=f"PMT {pmts[0]}", c=colors[pmts[0]])
    print('For PMT the ' + pmts + ' the spread raw in count rate is: ' + str(np.max((y/v)*100)-np.min((y/v)*100)))
    rel_rates.append((np.max((y/v)*100), np.min((y/v)*100)))
ax.legend(ncol=3)
ax.set_xlabel("Ball phi [deg]")
ax.set_ylabel("Relative Count Rate")
ax.set_title(title)
ax.set_ylim(80, 110)


# In[81]:


print("The max relative count rate difference is:", np.max(rel_rates)-np.min(rel_rates))
print("The average relative count rate difference is:", 100-np.min(rel_rates))


# # 2 rotation

# In[82]:


title="2 rotations of fine adjust"


# In[14]:


prefix="Ball2_2rot_176.02mA_"
rates = {}
phis = np.array([-180, -135, -90, 0, 90, 180])
for phi in phis:
    rates[phi] = {}
    for pmts in ['BV', 'CV', 'TV']:
        rates[phi][pmts] = {}
        for n, f in get_files(prefix, pmts, phi).items():
            rates[phi][pmts][n] = {}
            for p in pmts:
                rates[phi][pmts][n][p] = get_rate(f, pmt_channels[p], plot=True)
                print(phi, pmts, n, p, rates[phi][pmts][n][p])
            rates[phi][pmts][n]['M'] = get_average_pulse_integral(f, 3)
            print(phi, pmts, n, 'M', rates[phi][pmts][n]['M'])


# In[83]:


phis = np.array([-180, -135, -90, 0, 90, 180])
with open('2rot_data.pkl', 'rb') as f:
    rates = pickle.load(f)


# In[86]:


fig, ax = plt.subplots(figsize=(6, 5))
all_x = []
all_v = []
all_v_dr1 = []
all_v_dr2 = []
correction_gradient = {}
correction_intercept = {}
for pmts in ['BV', 'CV', 'TV']:
    x = [rates[phi][pmts]['LDON']['M'] for phi in phis]
    all_x += x
    all_v += [rates[phi][pmts]['LDON']['V'] for phi in phis]
    # all_v_dr1 += [rates[phi][pmts]['LDOFF']['V'] for phi in phis]
    # all_v_dr2 += [rates[phi][pmts]['LDOFFa']['V'] for phi in phis]
    x = np.array(x)
    y = [rates[phi][pmts]['LDON'][pmts[0]] for phi in phis]
    # y_dr1 = [rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis]
    # y_dr2 = [rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis]
    # y_dr1 = np.array(y_dr1)
    # y_dr2 = np.array(y_dr2)
    # y = y - 0.5*(y_dr1+y_dr2)
    ax.scatter(x, y, label=f"PMT {pmts[0]}", c=colors[pmts[0]])
    a, b = np.polyfit(x, y, 1)
    correction_gradient[pmts[0]] = a
    correction_intercept[pmts[0]] = b
    ax.plot(x, a*x+b, c=colors[pmts[0]])
    print(pmts[0], a, b)
all_x = np.array(all_x)
all_v = np.array(all_v)
# all_v_dr1 = np.array(all_v_dr1)
# all_v_dr2 = np.array(all_v_dr2)
# all_v = all_v - 0.5*(all_v_dr1 + all_v_dr2)
ax.scatter(all_x, all_v, label=f"PMT V", c=colors['V'])
a, b = np.polyfit(all_x, all_v, 1)
print("V", a, b)
correction_gradient['V'] = a
correction_intercept['V'] = b
ax.plot(all_x, a*all_x+b, c=colors['V'])
ax.legend()
ax.set_xlabel("Average monitor PMT pulse integral")
ax.set_ylabel("Main PMT rate")


# In[87]:


for pmts in ['TV', 'CV', 'BV']:
    v_pmt_rates = np.array([rates[phi][pmts]['LDON']['V'] for phi in phis])
    # v_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF']['V'] for phi in phis])
    # v_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa']['V'] for phi in phis])
    # v_pmt_dr = (v_pmt_dr1+v_pmt_dr2)/2
    o_pmt_rates = np.array([rates[phi][pmts]['LDON'][pmts[0]] for phi in phis])
    # o_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis])
    # o_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis])
    # o_pmt_dr = (o_pmt_dr1+o_pmt_dr2)/2
    m_pmt_q = np.array([rates[phi][pmts]['LDON']['M'] for phi in phis])
    v_corrected = v_pmt_rates - correction_gradient['V']*(m_pmt_q - 1.275)
    o_corrected = o_pmt_rates - correction_gradient[pmts[0]]*(m_pmt_q - 1.275)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(phis, o_pmt_rates, '.-', c='blue', label=f"{pmts[0]} PMT raw rate")
    # ax.plot(phis, o_pmt_dr1, '.--', c='blue', label=f"{pmts[0]} PMT DR 1")
    # ax.plot(phis, o_pmt_dr1, '.-.', c='blue', label=f"{pmts[0]} PMT DR 2")
    ax.plot(phis, v_pmt_rates, '.-', c='orange', label=f"V PMT raw rate")
    # ax.plot(phis, v_pmt_dr1, '.--', c='orange', label=f"V PMT DR 1")
    # ax.plot(phis, v_pmt_dr2, '.-.', c='orange', label=f"V PMT DR 1")
    #ax.set_xlim(-45, 45)
    ax.set_ylim(50, 80)
    ax.set_xlabel("Ball phi [deg]")
    ax.set_ylabel("Fractional PMT Count Rate")
    ax.yaxis.label.set_color('blue')
    ax2 = ax.twinx()
    ax2.plot(phis, m_pmt_q, '.--', c='black', label=f"Monitor")
    ax2.set_ylim(1.2, 1.4)
    ax2.set_ylabel("Monitor average pulse integral")
    ax.legend()
    # ax2.legend(loc="lower right")


# In[88]:


fig, ax = plt.subplots(figsize=(6, 5))
for o in 'BCT':
    pmts = o+'V'
    v_pmt_rates = np.array([rates[phi][pmts]['LDON']['V'] for phi in phis])
    #v_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF']['V'] for phi in phis])
    #v_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa']['V'] for phi in phis])
    #v_pmt_dr = (v_pmt_dr1+v_pmt_dr2)/2
    ax.plot(phis, v_pmt_rates,'.-', label=f"V PMT rate while measuring with {o} PMT")
#ax.set_xlim(-45, 45)
ax.set_ylim(50, 80)
ax.set_xlabel("Ball phi [deg]")
ax.set_ylabel("PMT Count Rate (%)")
ax.legend()


# In[89]:


fig, ax = plt.subplots(1, 3, figsize=(12, 5))
for i, pmts in enumerate(['TV', 'CV', 'BV']):
    v_pmt_rates = np.array([rates[phi][pmts]['LDON']['V'] for phi in phis])
    # v_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF']['V'] for phi in phis])
    # v_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa']['V'] for phi in phis])
    # v_pmt_dr = (v_pmt_dr1+v_pmt_dr2)/2
    o_pmt_rates = np.array([rates[phi][pmts]['LDON'][pmts[0]] for phi in phis])
    # o_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis])
    # o_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis])
    # o_pmt_dr = (o_pmt_dr1+o_pmt_dr2)/2
    m_pmt_q = np.array([rates[phi][pmts]['LDON']['M'] for phi in phis])
    v_corrected = v_pmt_rates - correction_gradient['V']*(m_pmt_q - 1.275)
    o_corrected = o_pmt_rates - correction_gradient[pmts[0]]*(m_pmt_q - 1.275)
    ax[i].plot(phis, o_pmt_rates, '.:', c='blue', label=f"{pmts[0]} PMT raw rate")
    ax[i].plot(phis, o_corrected, '.-', c='blue', label=f"{pmts[0]} PMT corrected")
    ax[i].plot(phis, v_pmt_rates, '.:', c='orange', label=f"V PMT raw rate")
    ax[i].plot(phis, v_corrected, '.-', c='orange', label=f"V PMT corrected")
    #ax.set_xlim(-45, 45)
    ax[i].set_ylim(50, 80)
    ax[i].set_xlabel("Ball phi [deg]")
    ax[i].set_ylabel(" PMT Count Rate (%)")
    ax[i].yaxis.label.set_color('blue')
    #ax2 = ax[i].twinx()
    #ax2.plot(phis, m_pmt_q, '.--', c='black', label=f"Monitor")
    #ax2.set_ylim(1.2, 1.6)
    #ax2.set_ylabel("Monitor average pulse integral")
    ax[i].legend()
    #ax2.legend(loc="lower right")
fig.suptitle(title)
plt.tight_layout()


# In[90]:


fig, ax = plt.subplots(figsize=(6, 5))
all_x = []
all_v = []
all_v_dr1 = []
all_v_dr2 = []
for pmts in ['TV', 'CV', 'BV']:
    x = [rates[phi][pmts]['LDON']['M'] for phi in phis]
    all_x += x
    all_v += [rates[phi][pmts]['LDON']['V'] for phi in phis]
    # all_v_dr1 += [rates[phi][pmts]['LDOFF']['V'] for phi in phis]
    # all_v_dr2 += [rates[phi][pmts]['LDOFFa']['V'] for phi in phis]
    x = np.array(x)
    y = [rates[phi][pmts]['LDON'][pmts[0]] for phi in phis]
    # y_dr1 = [rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis]
    # y_dr2 = [rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis]
    # y_dr1 = np.array(y_dr1)
    # y_dr2 = np.array(y_dr2)
    y = y - correction_gradient[pmts[0]]*(x-1.275)
    ax.scatter(x, y, label=f"PMT {pmts[0]}", c=colors[pmts[0]])
    a, b = np.polyfit(x, y, 1)
    ax.plot(x, a*x+b, c=colors[pmts[0]])
    print(pmts[0], a, b)
all_x = np.array(all_x)
all_v = np.array(all_v)
# all_v_dr1 = np.array(all_v_dr1)
# all_v_dr2 = np.array(all_v_dr2)
all_v = all_v - correction_gradient['V']*(all_x-1.275)
ax.scatter(all_x, all_v, label=f"PMT V", c=colors['V'])
a, b = np.polyfit(all_x, all_v, 1)
print("V", a, b)
ax.plot(all_x, a*all_x+b, c=colors['V'])
ax.legend()
ax.set_xlabel("Average monitor PMT pulse integral")
ax.set_ylabel("Main PMT rate")


# In[91]:


fig, ax = plt.subplots(figsize=(6, 5))

raw_rates = []
corr_rates = []


for pmts in ['TV', 'CV', 'BV']:
    m = [rates[phi][pmts]['LDON']['M'] for phi in phis]
    v = [rates[phi][pmts]['LDON']['V'] for phi in phis]
    # v_dr1 = [rates[phi][pmts]['LDOFF']['V'] for phi in phis]
    # v_dr2 = [rates[phi][pmts]['LDOFFa']['V'] for phi in phis]
    y = [rates[phi][pmts]['LDON'][pmts[0]] for phi in phis]
    # y_dr1 = [rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis]
    # y_dr2 = [rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis]
    m = np.array(m)
    v = np.array(v)
    # v_dr1 = np.array(v_dr1)
    # v_dr2 = np.array(v_dr2)
    y = np.array(y)
    # y_dr1 = np.array(y_dr1)
    # y_dr2 = np.array(y_dr2)
    y_ = y - correction_gradient[pmts[0]]*(m-1.275)
    v_ = v - correction_gradient['V']*(m-1.275)
    raw_rates.append((np.max(y), np.min(y), np.max(v), np.min(v)))
    corr_rates.append((np.max(y_), np.min(y_), np.max(v_), np.min(v_)))
    ax.plot(phis, y, '.-', label=f"PMT {pmts[0]}", c=colors[pmts[0]])
    ax.plot(phis, v, '.:', label=f"PMT V ({pmts[0]})", c=colors[pmts[0]])
    print('For PMT ' + pmts + ' the spread in raw count rate is: ' + str(np.max(y)-np.min(y)))
    print('For PMT V of ' + pmts + ' the spread raw in count rate is: ' + str(np.max(v)-np.min(v)))
    print('For PMT ' + pmts + ' the spread in corrected count rate is: ' + str(np.max(y_)-np.min(y_)))
    print('For PMT V of ' + pmts + ' the spread corrected in count rate is: ' + str(np.max(v_)-np.min(v_)))

ax.legend(ncol=3)
ax.set_ylim(50, 80)
ax.set_xlabel("Ball phi [deg]")
ax.set_ylabel("Corrected PMT rate")
ax.set_title("Corrected by Solid angle & monitor PMT, not relative PMT efficiencies")


# In[92]:


print("The max raw count rate difference is:", (np.max(raw_rates)-np.min(raw_rates)))
print("The max corrected count rate difference is:", (np.max(corr_rates)-np.min(corr_rates)))
# print("The average relative count rate difference is:", np.average(np.abs(np.diff(raw_rates))))
# print("The avergae corrected count rate difference is:", np.average(np.abs(np.diff(corr_rates))))


# In[93]:


rel_rates = []

fig, ax = plt.subplots(figsize=(6, 5))
for pmts in ['TV', 'CV', 'BV']:
    m = [rates[phi][pmts]['LDON']['M'] for phi in phis]
    v = [rates[phi][pmts]['LDON']['V'] for phi in phis]
    # v_dr1 = [rates[phi][pmts]['LDOFF']['V'] for phi in phis]
    # v_dr2 = [rates[phi][pmts]['LDOFFa']['V'] for phi in phis]
    y = [rates[phi][pmts]['LDON'][pmts[0]] for phi in phis]
    # y_dr1 = [rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis]
    # y_dr2 = [rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis]
    m = np.array(m)
    v = np.array(v)
    # v_dr1 = np.array(v_dr1)
    # v_dr2 = np.array(v_dr2)
    y = np.array(y)
    # y_dr1 = np.array(y_dr1)
    # y_dr2 = np.array(y_dr2)
    y = y/eff_corrections[pmts[0]] - correction_gradient[pmts[0]]*(m-1.275)
    v = v - correction_gradient['V']*(m-1.275)
    ax.plot(phis, (y/v)*100, '.-', label=f"PMT {pmts[0]}", c=colors[pmts[0]])
    print('For PMT the ' + pmts + ' the spread raw in count rate is: ' + str(np.max((y/v)*100)-np.min((y/v)*100)))
    rel_rates.append((np.max((y/v)*100), np.min((y/v)*100)))
ax.legend(ncol=3)
ax.set_xlabel("Ball phi [deg]")
ax.set_ylabel("Corrected PMT rate / corrected V PMT rate")
ax.set_title(title)
ax.set_ylim(80, 100)


# In[94]:


print("The max relative count rate difference is:", np.max(rel_rates)-np.min(rel_rates))
print("The average relative count rate difference is:", (100-np.min(rel_rates)))


# # 2mm spacer

# In[95]:


title="2mm spacer (approx. 2.5 rotations)"


# In[16]:


os.chdir("/vols/hyperk/WCTE/WCTE540/nick/data")

prefix="Ball2_2mmSpacer_176.02mA_"
rates = {}
phis = np.array([-180, -135, -90, 0, 90, 135, 180])
for phi in phis:
    print(phi)
    rates[phi] = {}
    for pmts in ['BV', 'CV', 'TV']:
        print(pmts)
        rates[phi][pmts] = {}
        for n, f in get_files(prefix, pmts, phi).items():
            rates[phi][pmts][n] = {}
            for p in pmts:
                rates[phi][pmts][n][p] = get_rate(f, pmt_channels[p], plot=True)
                print(phi, pmts, n, p, rates[phi][pmts][n][p])
            rates[phi][pmts][n]['M'] = get_average_pulse_integral(f, 3)
            print(phi, pmts, n, 'M', rates[phi][pmts][n]['M'])


# In[97]:


fig, ax = plt.subplots(figsize=(6, 5))
all_x = []
all_v = []
all_v_dr1 = []
all_v_dr2 = []
correction_gradient = {}
correction_intercept = {}
for pmts in ['BV', 'CV', 'TV']:
    x = [rates[phi][pmts]['LDON']['M'] for phi in phis]
    all_x += x
    all_v += [rates[phi][pmts]['LDON']['V'] for phi in phis]
    # all_v_dr1 += [rates[phi][pmts]['LDOFF']['V'] for phi in phis]
    # all_v_dr2 += [rates[phi][pmts]['LDOFFa']['V'] for phi in phis]
    x = np.array(x)
    y = [rates[phi][pmts]['LDON'][pmts[0]] for phi in phis]
    # y_dr1 = [rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis]
    # y_dr2 = [rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis]
    # y_dr1 = np.array(y_dr1)
    # y_dr2 = np.array(y_dr2)
    # y = y - 0.5*(y_dr1+y_dr2)
    ax.scatter(x, y, label=f"PMT {pmts[0]}", c=colors[pmts[0]])
    a, b = np.polyfit(x, y, 1)
    correction_gradient[pmts[0]] = a
    correction_intercept[pmts[0]] = b
    ax.plot(x, a*x+b, c=colors[pmts[0]])
    print(pmts[0], a, b)
all_x = np.array(all_x)
all_v = np.array(all_v)
# all_v_dr1 = np.array(all_v_dr1)
# all_v_dr2 = np.array(all_v_dr2)
# all_v = all_v - 0.5*(all_v_dr1 + all_v_dr2)
ax.scatter(all_x, all_v, label=f"PMT V", c=colors['V'])
a, b = np.polyfit(all_x, all_v, 1)
print("V", a, b)
correction_gradient['V'] = a
correction_intercept['V'] = b
ax.plot(all_x, a*all_x+b, c=colors['V'])
ax.legend()
ax.set_xlabel("Average monitor PMT pulse integral")
ax.set_ylabel("Main PMT rate")


# In[98]:


phis = np.array([-180, -135, -90, 0, 90, 135, 180])
with open('2mmspacer_data.pkl', 'rb') as f:
    rates = pickle.load(f)


# In[99]:


for pmts in ['TV', 'CV', 'BV']:
    v_pmt_rates = np.array([rates[phi][pmts]['LDON']['V'] for phi in phis])
    # v_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF']['V'] for phi in phis])
    # v_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa']['V'] for phi in phis])
    # v_pmt_dr = (v_pmt_dr1+v_pmt_dr2)/2
    o_pmt_rates = np.array([rates[phi][pmts]['LDON'][pmts[0]] for phi in phis])
    # o_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis])
    # o_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis])
    # o_pmt_dr = (o_pmt_dr1+o_pmt_dr2)/2
    m_pmt_q = np.array([rates[phi][pmts]['LDON']['M'] for phi in phis])
    v_corrected = v_pmt_rates - correction_gradient['V']*(m_pmt_q - 1.275)
    o_corrected = o_pmt_rates - correction_gradient[pmts[0]]*(m_pmt_q - 1.275)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(phis, o_pmt_rates, '.-', c='blue', label=f"{pmts[0]} PMT raw rate")
    #ax.plot(phis, o_pmt_dr1, '.-.', c='blue', label=f"{pmts[0]} PMT DR 2")
    ax.plot(phis, v_pmt_rates, '.-', c='orange', label=f"V PMT raw rate")
    #ax.plot(phis, v_pmt_dr1, '.--', c='orange', label=f"V PMT DR 1")
    #ax.plot(phis, v_pmt_dr2, '.-.', c='orange', label=f"V PMT DR 1")
    #ax.set_xlim(-45, 45)
    ax.set_ylim(50, 80)
    ax.set_xlabel("Ball phi [deg]")
    ax.set_ylabel("PMT Count Rate (%)")
    ax.yaxis.label.set_color('blue')
    # ax2 = ax.twinx()
    # ax2.plot(phis, m_pmt_q, '.--', c='black', label=f"Monitor")
    # ax2.set_ylim(1.2, 1.4)
    # ax2.set_ylabel("Monitor average pulse integral")
    ax.legend()
    # ax2.legend(loc="lower right")


# In[100]:


fig, ax = plt.subplots(figsize=(6, 5))
for o in 'BCT':
    pmts = o+'V'
    v_pmt_rates = np.array([rates[phi][pmts]['LDON']['V'] for phi in phis])
    # v_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF']['V'] for phi in phis])
    # v_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa']['V'] for phi in phis])
    # v_pmt_dr = (v_pmt_dr1+v_pmt_dr2)/2
    ax.plot(phis, v_pmt_rates, '.-', label=f"V PMT rate while measuring with {o} PMT")
#ax.set_xlim(-45, 45)
ax.set_ylim(50, 80)
ax.set_xlabel("Ball phi [deg]")
ax.set_ylabel("PMT Count Rate (%)")
ax.legend()


# In[101]:


fig, ax = plt.subplots(1, 3, figsize=(12, 5))
for i, pmts in enumerate(['TV', 'CV', 'BV']):
    v_pmt_rates = np.array([rates[phi][pmts]['LDON']['V'] for phi in phis])
    # v_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF']['V'] for phi in phis])
    # v_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa']['V'] for phi in phis])
    # v_pmt_dr = (v_pmt_dr1+v_pmt_dr2)/2
    o_pmt_rates = np.array([rates[phi][pmts]['LDON'][pmts[0]] for phi in phis])
    # o_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis])
    # o_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis])
    # o_pmt_dr = (o_pmt_dr1+o_pmt_dr2)/2
    m_pmt_q = np.array([rates[phi][pmts]['LDON']['M'] for phi in phis])
    v_corrected = v_pmt_rates - correction_gradient['V']*(m_pmt_q - 1.275)
    o_corrected = o_pmt_rates - correction_gradient[pmts[0]]*(m_pmt_q - 1.275)
    ax[i].plot(phis, o_pmt_rates, '.:', c='blue', label=f"{pmts[0]} PMT raw rate")
    ax[i].plot(phis, o_corrected, '.-', c='blue', label=f"{pmts[0]} PMT corrected")
    ax[i].plot(phis, v_pmt_rates, '.:', c='orange', label=f"V PMT raw rate")
    ax[i].plot(phis, v_corrected, '.-', c='orange', label=f"V PMT corrected")
    #ax.set_xlim(-45, 45)
    ax[i].set_ylim(50, 80)
    ax[i].set_xlabel("Ball phi [deg]")
    ax[i].set_ylabel("PMT Count Rate (%)")
    ax[i].yaxis.label.set_color('blue')
    ax2 = ax[i].twinx()
    # ax2.plot(phis, m_pmt_q, '.--', c='black', label=f"Monitor")
    # ax2.set_ylim(1.2, 1.6)
    # ax2.set_ylabel("Monitor average pulse integral")
    ax[i].legend()
    # ax2.legend(loc="lower right")
fig.suptitle(title)
plt.tight_layout()


# In[102]:


fig, ax = plt.subplots(figsize=(6, 5))
all_x = []
all_v = []
all_v_dr1 = []
all_v_dr2 = []
for pmts in ['BV', 'CV', 'TV']:
    x = [rates[phi][pmts]['LDON']['M'] for phi in phis]
    all_x += x
    all_v += [rates[phi][pmts]['LDON']['V'] for phi in phis]
    # all_v_dr1 += [rates[phi][pmts]['LDOFF']['V'] for phi in phis]
    # all_v_dr2 += [rates[phi][pmts]['LDOFFa']['V'] for phi in phis]
    x = np.array(x)
    y = [rates[phi][pmts]['LDON'][pmts[0]] for phi in phis]
    # y_dr1 = [rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis]
    # y_dr2 = [rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis]
    # y_dr1 = np.array(y_dr1)
    # y_dr2 = np.array(y_dr2)
    y = y - correction_gradient[pmts[0]]*(x-1.275)
    ax.scatter(x, y, label=f"PMT {pmts[0]}", c=colors[pmts[0]])
    a, b = np.polyfit(x, y, 1)
    ax.plot(x, a*x+b, c=colors[pmts[0]])
    print(pmts[0], a, b)
all_x = np.array(all_x)
all_v = np.array(all_v)
# all_v_dr1 = np.array(all_v_dr1)
# all_v_dr2 = np.array(all_v_dr2)
all_v = all_v  - correction_gradient['V']*(all_x-1.275)
ax.scatter(all_x, all_v, label=f"PMT V", c=colors['V'])
a, b = np.polyfit(all_x, all_v, 1)
print("V", a, b)
ax.plot(all_x, a*all_x+b, c=colors['V'])
ax.legend()
ax.set_xlabel("Average monitor PMT pulse integral")
ax.set_ylabel("Main PMT rate")


# In[103]:


raw_rates = []
corr_rates = []

fig, ax = plt.subplots(figsize=(6, 5))
for pmts in ['TV', 'CV', 'BV']:
    m = [rates[phi][pmts]['LDON']['M'] for phi in phis]
    v = [rates[phi][pmts]['LDON']['V'] for phi in phis]
    # v_dr1 = [rates[phi][pmts]['LDOFF']['V'] for phi in phis]
    # v_dr2 = [rates[phi][pmts]['LDOFFa']['V'] for phi in phis]
    y = [rates[phi][pmts]['LDON'][pmts[0]] for phi in phis]
    # y_dr1 = [rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis]
    # y_dr2 = [rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis]
    m = np.array(m)
    v = np.array(v)
    # v_dr1 = np.array(v_dr1)
    # v_dr2 = np.array(v_dr2)
    y = np.array(y)
    # y_dr1 = np.array(y_dr1)
    # y_dr2 = np.array(y_dr2)
    y_ = y - correction_gradient[pmts[0]]*(m-1.275)
    v_ = v - correction_gradient['V']*(m-1.275)
    raw_rates.append((np.max(y), np.min(y), np.max(v), np.min(v)))
    corr_rates.append((np.max(y_), np.min(y_), np.max(v_), np.min(v_)))
    print('For PMT ' + pmts + ' the spread in raw count rate is: ' + str(np.max(y)-np.min(y)))
    print('For PMT V of ' + pmts + ' the spread raw in count rate is: ' + str(np.max(v)-np.min(v)))
    print('For PMT ' + pmts + ' the spread in corrected count rate is: ' + str(np.max(y_)-np.min(y_)))
    print('For PMT V of ' + pmts + ' the spread corrected in count rate is: ' + str(np.max(v_)-np.min(v_)))
    ax.plot(phis, y, '.-', label=f"PMT {pmts[0]}", c=colors[pmts[0]])
    ax.plot(phis, v, '.:', label=f"PMT V ({pmts[0]})", c=colors[pmts[0]])
ax.legend(ncol=3)
ax.set_ylim([50, 80])
ax.set_xlabel("Ball phi [deg]")
ax.set_ylabel("Corrected PMT rate")
ax.set_title("Corrected by Solid angle & monitor PMT, not relative PMT efficiencies")


# In[104]:


print("The max raw count rate difference is:", (np.max(raw_rates)-np.min(raw_rates)))
print("The max corrected count rate difference is:", (np.max(corr_rates)-np.min(corr_rates)))
# print("The average relative count rate difference is:", np.average(np.abs(np.diff(rel_rates))))
# print("The avergae corrected count rate difference is:", np.average(np.abs(np.diff(corr_rates))))


# In[105]:


rel_rates = []
fig, ax = plt.subplots(figsize=(6, 5))
for pmts in ['BV', 'CV', 'TV']:
    m = [rates[phi][pmts]['LDON']['M'] for phi in phis]
    v = [rates[phi][pmts]['LDON']['V'] for phi in phis]
    # v_dr1 = [rates[phi][pmts]['LDOFF']['V'] for phi in phis]
    # v_dr2 = [rates[phi][pmts]['LDOFFa']['V'] for phi in phis]
    y = [rates[phi][pmts]['LDON'][pmts[0]] for phi in phis]
    # y_dr1 = [rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis]
    # y_dr2 = [rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis]
    m = np.array(m)
    v = np.array(v)
    # v_dr1 = np.array(v_dr1)
    # v_dr2 = np.array(v_dr2)
    y = np.array(y)
    # y_dr1 = np.array(y_dr1)
    # y_dr2 = np.array(y_dr2)
    y = y/eff_corrections[pmts[0]] - correction_gradient[pmts[0]]*(m-1.275)
    v = v  - correction_gradient['V']*(m-1.275)
    print('For PMT the ' + pmts + ' the spread raw in count rate is: ' + str(np.max((y/v)*100)-np.min((y/v)*100)))
    rel_rates.append((np.max((y/v)*100), np.min((y/v)*100)))
    ax.plot(phis, (y/v)*100, '.-', label=f"PMT {pmts[0]}", c=colors[pmts[0]])
ax.legend(ncol=3)
ax.set_xlabel("Ball phi [deg]")
ax.set_ylabel("Corrected PMT rate / corrected V PMT rate")
ax.set_title(title)
ax.set_ylim(80, 100)


# In[106]:


print("The max relative count rate difference is:", np.max(rel_rates)-np.min(rel_rates))
print("The average relative count rate difference is:", 100-np.min(rel_rates))


# # 3mm spacer

# In[107]:


title="3mm Spacer"


# In[28]:


prefix="Ball2_3mmSpacer_176.02mA_"
rates = {}
phis = np.array([-180, -90, 0, 90])
for phi in phis:
    print(phi)
    rates[phi] = {}
    for pmts in ['BV', 'CV', 'TV']:
        print(pmts)
        rates[phi][pmts] = {}
        for n, f in get_files(prefix, pmts, phi).items():
            rates[phi][pmts][n] = {}
            for p in pmts:
                rates[phi][pmts][n][p] = get_rate(f, pmt_channels[p], plot=True)
                print(phi, pmts, n, p, rates[phi][pmts][n][p])
            rates[phi][pmts][n]['M'] = get_average_pulse_integral(f, 3)
            print(phi, pmts, n, 'M', rates[phi][pmts][n]['M'])


# In[109]:


fig, ax = plt.subplots(figsize=(6, 5))
all_x = []
all_v = []
all_v_dr1 = []
all_v_dr2 = []
correction_gradient = {}
correction_intercept = {}
for pmts in ['BV', 'CV', 'TV']:
    x = [rates[phi][pmts]['LDON']['M'] for phi in phis]
    all_x += x
    all_v += [rates[phi][pmts]['LDON']['V'] for phi in phis]
    # all_v_dr1 += [rates[phi][pmts]['LDOFF']['V'] for phi in phis]
    # all_v_dr2 += [rates[phi][pmts]['LDOFFa']['V'] for phi in phis]
    x = np.array(x)
    y = [rates[phi][pmts]['LDON'][pmts[0]] for phi in phis]
    # y_dr1 = [rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis]
    # y_dr2 = [rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis]
    # y_dr1 = np.array(y_dr1)
    # y_dr2 = np.array(y_dr2)
    # y = y - 0.5*(y_dr1+y_dr2)
    ax.scatter(x, y, label=f"PMT {pmts[0]}", c=colors[pmts[0]])
    a, b = np.polyfit(x, y, 1)
    correction_gradient[pmts[0]] = a
    correction_intercept[pmts[0]] = b
    ax.plot(x, a*x+b, c=colors[pmts[0]])
    print(pmts[0], a, b)
all_x = np.array(all_x)
all_v = np.array(all_v)
# all_v_dr1 = np.array(all_v_dr1)
# all_v_dr2 = np.array(all_v_dr2)
# all_v = all_v - 0.5*(all_v_dr1 + all_v_dr2)
ax.scatter(all_x, all_v, label=f"PMT V", c=colors['V'])
a, b = np.polyfit(all_x, all_v, 1)
print("V", a, b)
correction_gradient['V'] = a
correction_intercept['V'] = b
ax.plot(all_x, a*all_x+b, c=colors['V'])
ax.legend()
ax.set_xlabel("Average monitor PMT pulse integral")
ax.set_ylabel("Main PMT rate")


# In[110]:


phis = np.array([-180, -90, 0, 90])
with open('3mmspacer_data.pkl', 'rb') as f:
    rates = pickle.load(f)


# In[111]:


for pmts in ['BV', 'CV', 'TV']:
    v_pmt_rates = np.array([rates[phi][pmts]['LDON']['V'] for phi in phis])
    # v_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF']['V'] for phi in phis])
    # v_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa']['V'] for phi in phis])
    # v_pmt_dr = (v_pmt_dr1+v_pmt_dr2)/2
    o_pmt_rates = np.array([rates[phi][pmts]['LDON'][pmts[0]] for phi in phis])
    # o_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis])
    # o_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis])
    # o_pmt_dr = (o_pmt_dr1+o_pmt_dr2)/2
    m_pmt_q = np.array([rates[phi][pmts]['LDON']['M'] for phi in phis])
    v_corrected = v_pmt_rates - correction_gradient['V']*(m_pmt_q - 1.275)
    o_corrected = o_pmt_rates - correction_gradient[pmts[0]]*(m_pmt_q - 1.275)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(phis, o_pmt_rates, '.-', c='blue', label=f"{pmts[0]} PMT raw rate")
    # ax.plot(phis, o_pmt_dr1, '.--', c='blue', label=f"{pmts[0]} PMT DR 1")
    # ax.plot(phis, o_pmt_dr1, '.-.', c='blue', label=f"{pmts[0]} PMT DR 2")
    ax.plot(phis, v_pmt_rates, '.-', c='orange', label=f"V PMT raw rate")
    # ax.plot(phis, v_pmt_dr1, '.--', c='orange', label=f"V PMT DR 1")
    # ax.plot(phis, v_pmt_dr2, '.-.', c='orange', label=f"V PMT DR 1")
    #ax.set_xlim(-45, 45)
    ax.set_ylim(50,80)
    ax.set_xlabel("Ball phi [deg]")
    ax.set_ylabel("Rate")
    ax.yaxis.label.set_color('blue')
    # ax2 = ax.twinx()
    # ax2.plot(phis, m_pmt_q, '.--', c='black', label=f"Monitor")
    # ax2.set_ylim(1.2, 1.4)
    # ax2.set_ylabel("Monitor average pulse integral")
    ax.legend()
    # ax2.legend(loc="lower right")


# In[112]:


fig, ax = plt.subplots(figsize=(6, 5))
for o in 'BCT':
    pmts = o+'V'
    v_pmt_rates = np.array([rates[phi][pmts]['LDON']['V'] for phi in phis])
    # v_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF']['V'] for phi in phis])
    # v_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa']['V'] for phi in phis])
    # v_pmt_dr = (v_pmt_dr1+v_pmt_dr2)/2
    ax.plot(phis, v_pmt_rates, '.-', label=f"V PMT rate while measuring with {o} PMT")
# ax.set_xlim(-45, 45)
ax.set_ylim(50, 80)
ax.set_xlabel("Ball phi [deg]")
ax.set_ylabel("Rate")
ax.legend()


# In[113]:


fig, ax = plt.subplots(1, 3, figsize=(12, 5))
for i, pmts in enumerate(['TV', 'CV', 'BV']):
    v_pmt_rates = np.array([rates[phi][pmts]['LDON']['V'] for phi in phis])
    # v_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF']['V'] for phi in phis])
    # v_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa']['V'] for phi in phis])
    # v_pmt_dr = (v_pmt_dr1+v_pmt_dr2)/2
    o_pmt_rates = np.array([rates[phi][pmts]['LDON'][pmts[0]] for phi in phis])
    # o_pmt_dr1 = np.array([rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis])
    # o_pmt_dr2 = np.array([rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis])
    # o_pmt_dr = (o_pmt_dr1+o_pmt_dr2)/2
    m_pmt_q = np.array([rates[phi][pmts]['LDON']['M'] for phi in phis])
    v_corrected = v_pmt_rates - correction_gradient['V']*(m_pmt_q - 1.275)
    o_corrected = o_pmt_rates - correction_gradient[pmts[0]]*(m_pmt_q - 1.275)
    ax[i].plot(phis, o_pmt_rates, '.:', c='blue', label=f"{pmts[0]} PMT raw rate")
    ax[i].plot(phis, o_corrected, '.-', c='blue', label=f"{pmts[0]} PMT corrected")
    ax[i].plot(phis, v_pmt_rates, '.:', c='orange', label=f"V PMT raw rate")
    ax[i].plot(phis, v_corrected, '.-', c='orange', label=f"V PMT corrected")
    #ax.set_xlim(-45, 45)
    ax[i].set_ylim(50, 80)
    ax[i].set_xlabel("Ball phi [deg]")
    ax[i].set_ylabel("PMT Count Rate (%)")
    ax[i].yaxis.label.set_color('blue')
    #ax2 = ax[i].twinx()
    #ax2.plot(phis, m_pmt_q, '.--', c='black', label=f"Monitor")
    #ax2.set_ylim(1.2, 1.6)
    #ax2.set_ylabel("Monitor average pulse integral")
    ax[i].legend()
    ##x2.legend(loc="lower right")
fig.suptitle(title)
plt.tight_layout()


# In[114]:


fig, ax = plt.subplots(figsize=(6, 5))
all_x = []
all_v = []
all_v_dr1 = []
all_v_dr2 = []
for pmts in ['BV', 'CV', 'TV']:
    x = [rates[phi][pmts]['LDON']['M'] for phi in phis]
    all_x += x
    all_v += [rates[phi][pmts]['LDON']['V'] for phi in phis]
    # all_v_dr1 += [rates[phi][pmts]['LDOFF']['V'] for phi in phis]
    # all_v_dr2 += [rates[phi][pmts]['LDOFFa']['V'] for phi in phis]
    x = np.array(x)
    y = [rates[phi][pmts]['LDON'][pmts[0]] for phi in phis]
    # y_dr1 = [rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis]
    # y_dr2 = [rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis]
    # y_dr1 = np.array(y_dr1)
    # y_dr2 = np.array(y_dr2)
    y = y - correction_gradient[pmts[0]]*(x-1.275)
    ax.scatter(x, y, label=f"PMT {pmts[0]}", c=colors[pmts[0]])
    a, b = np.polyfit(x, y, 1)
    ax.plot(x, a*x+b, c=colors[pmts[0]])
    print(pmts[0], a, b)
all_x = np.array(all_x)
all_v = np.array(all_v)
# all_v_dr1 = np.array(all_v_dr1)
# all_v_dr2 = np.array(all_v_dr2)
all_v = all_v - correction_gradient['V']*(all_x-1.275)
ax.scatter(all_x, all_v, label=f"PMT V", c=colors['V'])
a, b = np.polyfit(all_x, all_v, 1)
print("V", a, b)
ax.plot(all_x, a*all_x+b, c=colors['V'])
ax.legend()
ax.set_xlabel("Average monitor PMT pulse integral")
ax.set_ylabel("Main PMT rate")


# In[115]:


raw_rates = []
corr_rates = []
fig, ax = plt.subplots(figsize=(6, 5))
for pmts in ['TV', 'CV', 'BV']:
    m = [rates[phi][pmts]['LDON']['M'] for phi in phis]
    v = [rates[phi][pmts]['LDON']['V'] for phi in phis]
    # v_dr1 = [rates[phi][pmts]['LDOFF']['V'] for phi in phis]
    # v_dr2 = [rates[phi][pmts]['LDOFFa']['V'] for phi in phis]
    y = [rates[phi][pmts]['LDON'][pmts[0]] for phi in phis]
    # y_dr1 = [rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis]
    # y_dr2 = [rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis]
    m = np.array(m)
    v = np.array(v)
    # v_dr1 = np.array(v_dr1)
    # v_dr2 = np.array(v_dr2)
    y = np.array(y)
    # y_dr1 = np.array(y_dr1)
    # y_dr2 = np.array(y_dr2)
    y_ = y - correction_gradient[pmts[0]]*(m-1.275)
    v_ = v - correction_gradient['V']*(m-1.275)
    raw_rates.append((np.max(y), np.min(y), np.max(v), np.min(v)))
    corr_rates.append((np.max(y_), np.min(y_), np.max(v_), np.min(v_)))
    print('For PMT ' + pmts + ' the spread in raw count rate is: ' + str(np.max(y)-np.min(y)))
    print('For PMT V of ' + pmts + ' the spread raw in count rate is: ' + str(np.max(v)-np.min(v)))
    print('For PMT ' + pmts + ' the spread in corrected count rate is: ' + str(np.max(y_)-np.min(y_)))
    print('For PMT V of ' + pmts + ' the spread corrected in count rate is: ' + str(np.max(v_)-np.min(v_)))
    ax.plot(phis, y, '.-', label=f"PMT {pmts[0]}", c=colors[pmts[0]])
    ax.plot(phis, v, '.:', label=f"PMT V ({pmts[0]})", c=colors[pmts[0]])
ax.legend(ncol=3)
ax.set_ylim([50,80])
ax.set_xlabel("Ball phi [deg]")
ax.set_ylabel("Corrected PMT rate")
ax.set_title("Corrected by Solid angle & monitor PMT, not relative PMT efficiencies")


# In[116]:


print("The max raw count rate difference is:", np.max(raw_rates)-np.min(raw_rates))
print("The max corrected count rate difference is:", np.max(corr_rates)-np.min(corr_rates))
# print("The average relative count rate difference is:", np.average(np.abs(np.diff(rel_rates))))
# print("The avergae corrected count rate difference is:", np.average(np.abs(np.diff(corr_rates))))


# In[117]:


rel_rates = []
fig, ax = plt.subplots(figsize=(6, 5))
for pmts in ['TV', 'CV', 'BV']:
    m = [rates[phi][pmts]['LDON']['M'] for phi in phis]
    v = [rates[phi][pmts]['LDON']['V'] for phi in phis]
    # v_dr1 = [rates[phi][pmts]['LDOFF']['V'] for phi in phis]
    # v_dr2 = [rates[phi][pmts]['LDOFFa']['V'] for phi in phis]
    y = [rates[phi][pmts]['LDON'][pmts[0]] for phi in phis]
    # y_dr1 = [rates[phi][pmts]['LDOFF'][pmts[0]] for phi in phis]
    # y_dr2 = [rates[phi][pmts]['LDOFFa'][pmts[0]] for phi in phis]
    m = np.array(m)
    v = np.array(v)
    # v_dr1 = np.array(v_dr1)
    # v_dr2 = np.array(v_dr2)
    y = np.array(y)
    # y_dr1 = np.array(y_dr1)
    # y_dr2 = np.array(y_dr2)
    y = y/eff_corrections[pmts[0]] - correction_gradient[pmts[0]]*(m-1.275)
    v = v - correction_gradient['V']*(m-1.275)
    ax.plot(phis, (y/v)*100, '.-', label=f"PMT {pmts[0]}", c=colors[pmts[0]])
    print('For PMT the ' + pmts + ' the spread raw in count rate is: ' + str(np.max((y/v)*100)-np.min((y/v)*100)))
    rel_rates.append((np.max((y/v)*100), np.min((y/v)*100)))
ax.legend(ncol=3)
ax.set_xlabel("Ball phi [deg]")
ax.set_ylabel("Corrected PMT rate / corrected V PMT rate")
ax.set_title(title)
ax.set_ylim(80, 100)


# In[118]:


print("The max relative count rate difference is:", np.max(rel_rates)-np.min(rel_rates))
print("The average relative count rate difference is:", 100-np.min(rel_rates))


# In[119]:


import numpy as np
import matplotlib.pyplot as plt
uncorr_res = [10.42, 11.96, 13.35, 13.84, 13.90]
corr_res = [7.72, 9.70, 10.78, 12.22, 14.43]
rel_res = [5.81, 3.48, 7.80, 6.82, 9.32]
v_res = [9.80, 12.69, 11.62, 12.22, 18.31]
x_labels = ['0mm', '0.8mm', '1.6mm', '2mm', '3mm']
x = np.arange(len(x_labels))

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(x, uncorr_res, 'o-', label='Uncorrected Uniformity Spread')
ax.plot(x, corr_res, 's-', label='Corrected Uniformity Spread')
ax.plot(x, rel_res, 'd-', label='Relative Uniformity Spread')
ax.plot(x, v_res, 'd-', label='Uniformity rel to V spread')

ax.set_xticks(x)
ax.set_xticklabels(x_labels)
ax.set_xlabel('Acrylic Rod Configuration', size=20)
ax.set_ylabel('Uniformity metric', size=20)
ax.legend(prop={'size': 15})
ax.set_title('Resolution Comparison for Different Configurations', size=20)
ax.set_ylim(0,20)
plt.show()


# In[7]:


pmt1 = extract_pmt_voltages('Ball2_2mmSpacer_176.02mA_0deg_BV_LDON_2024-09-26T14:49:18.root', pmt1=True)
peaks, noise = peak_minus_noise(pmt1, window_size=2000)


# In[14]:


pmn = [np.array(peaks) - np.array(noise) for peaks, noise in zip(peaks, noise)]
#plt.hist(peaks, bins=np.arange(0,30,0.15625), histtype='step')
plt.hist(pmn, bins=np.arange(0,30,0.15625), histtype='step')
plt.hist(noise, bins=np.arange(0,30,0.15625))
plt.yscale('log')
plt.ylim(10**2)
plt.show()


# In[1]:


import numpy as np
import matplotlib.pyplot as plt
uncorr_res = [10.42, 11.96, 13.35, 13.84, 13.90]
corr_res = [7.72, 9.70, 10.78, 12.22, 14.43]
rel_res = [5.81, 3.48, 7.80, 6.82, 9.32]
v_res = [9.80, 12.69, 11.62, 12.22, 18.31]
x_labels = ['0mm', '0.8mm', '1.6mm', '2mm', '3mm']
x = np.arange(len(x_labels))

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(x, uncorr_res, 'o-', label='Uncorrected Uniformity Spread')
ax.plot(x, corr_res, 's-', label='Corrected Uniformity Spread')
ax.plot(x, rel_res, 'd-', label='Relative Uniformity Spread')
ax.plot(x, v_res, 'd-', label='Uniformity rel to V spread')

ax.set_xticks(x)
ax.set_xticklabels(x_labels)
ax.set_xlabel('Acrylic Rod Configuration', size=20)
ax.set_ylabel('Uniformity metric', size=20)
ax.legend(prop={'size': 15})
ax.set_title('Resolution Comparison for Different Configurations', size=20)
ax.set_ylim(0,20)
plt.show()

