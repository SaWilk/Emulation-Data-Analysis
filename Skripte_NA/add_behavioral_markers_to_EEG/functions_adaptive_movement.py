import pandas as pd
import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelextrema
import scipy
from scipy import signal

# paths
folder_marker   = Path(r"R:\AG-Beste-Studien\Emulation\05_automagic\Analyzer\export")
path_behav      = Path(r"R:\AG-Beste-Studien\Emulation\04_data\behav")
path_new_marker = Path(r"R:\AG-Beste-Studien\Emulation\05_automagic\Analyzer\new_markers")

# overview of all relevant marker files
"""
all_marker_files = [f for f in os.listdir(folder_marker) if f.endswith(".vmrk") and "no_practice" in f]
marker_files = []
for f in all_marker_files:
    for i in range(7,11):
        if f"Pilot{i}" in f:
            marker_files.append(f)
"""

# overview of all relevant behavioral folders
all_subdirs = [f.name for f in os.scandir(path_behav) if f.is_dir()]
subdirs = all_subdirs

def load_trial(foldername, trialnum):
    # trialnum should start at 1
    all_files = os.listdir(path_behav/foldername)
    filename = [f for f in all_files if f.split("_")[-1].split(".")[0] == str(trialnum)][0]
    return pd.read_csv(path_behav/foldername/filename)
    
def write_markers(marker_label, markers, sbj_ID, task):
    f = open(path_new_marker/f"{sbj_ID}_{task}_adaptive_movement.txt", "w")
    f.writelines("Sampling rate: 256Hz, SamplingInterval: 3.90625ms\n")
    f.writelines("Type, Description, Position, Length, Channel\n")
    for m in markers:
        f.writelines(f"Stimulus, S {marker_label}, {m}, 1, All\n")
    f.close()

def debug_plot_new(behav_trial, EEG_samples_trial_match_for_plot, time_in_samples):
    error_sig   = np.array(behav_trial["error_abs"].tolist())
    pursuit_sig = np.array(behav_trial["purs-y"].tolist())
    error_sig = get_LP_filtered_signal(error_sig)
    pursuit_sig = get_LP_filtered_signal(pursuit_sig)
    plt.figure(figsize=(7,2.5))
    plt.plot(time_in_samples, error_sig)
    plt.plot(time_in_samples, pursuit_sig)
    plt.vlines([v for v in EEG_samples_trial_match_for_plot], min(pursuit_sig), max(pursuit_sig), "red")
    plt.show()

def debug_plot(behavioral_data, local_extrema):
    plt.figure()
    plt.plot(behavioral_data, "black")
    plt.vlines([v for v in local_extrema], min(behavioral_data), max(behavioral_data), "red")
    plt.show()

def load_marker_file(sbj_ID, task):
    all_match_files = [f for f in os.listdir(folder_marker) if sbj_ID in f and f"_{task}_" in f and "removed_practice.Markers" in f]
    fname = all_match_files[0]
    return pd.read_csv(folder_marker/fname, skiprows=2, header=None, delimiter=", ")

def trial_sample_infos_from_marker_df(df):
    df_starts_ends = df[df[1].isin(["S 27", "S 12", "S 16"])]
    start_trials_samp = []
    dur_trials_in_samp = []
    for i in range(len(df_starts_ends)):
        if i == 0:
            continue
        if df_starts_ends[1].tolist()[i] == "S 27":
            start_trials_samp.append(df_starts_ends[2].tolist()[i])
        elif df_starts_ends[1].tolist()[i - 1] == "S 27":
            dur_trials_in_samp.append(df_starts_ends[2].tolist()[i] - df_starts_ends[2].tolist()[i - 1])
    return [start_trials_samp, dur_trials_in_samp]

def get_LP_filtered_signal(task_sig, cutoff_freq=0.05):
    sos = signal.iirfilter(2, cutoff_freq, btype='lowpass', output='sos')
    return scipy.signal.sosfiltfilt(sos, task_sig.tolist())

def get_extrema(behav_trial_data, filter=True, only_max=False):

    # low-pass filter to cushion behavioral irregularities
    if filter:
        filtered_signal = get_LP_filtered_signal(behav_trial_data)
    else:
        filtered_signal = behav_trial_data.tolist()

    # skip the first 10 samples, so the artifactual first extreme is skipped
    maxima = signal.argrelextrema(np.array(filtered_signal[10:]), np.greater)[0].tolist()
    maxima = [m+10 for m in maxima]
    minima = signal.argrelextrema(np.array(filtered_signal[10:]), np.less)[0].tolist()
    minima = [m+10 for m in minima]
    
    max_vals = [filtered_signal[m] for m in maxima]
    min_vals = [filtered_signal[m] for m in minima]
    
    if not only_max:
        # beware: values not sorted!!
        return [maxima+minima, max_vals+min_vals]
    
    return [maxima, max_vals] # in this case they are sorted

def behavioral_events_in_EEG_format(behav_trial, behav_dimension, dur_trial_samp, start_trial_samp, trial_ind, debug=False):
    
    trial_behav_dim = np.array(behav_trial[behav_dimension].tolist())
    
    # derive both indices and values at indices (keep the n largest deviances)
    # explanation of part '"error" in behav_dimension': boolean, whether "only_max" is true (just be for absolute errors)
    samples      = get_extrema(trial_behav_dim, True, "error" in behav_dimension)[0]
    samples_vals = get_extrema(trial_behav_dim, True, "error" in behav_dimension)[1]
    
    # sort by the minimal number (left/right) of ramp-up/ramp-down
    min_ramps = []

    filtered_signal = get_LP_filtered_signal(trial_behav_dim)

    for samp in samples:

        current_samp = samp - 1
        lefts = 0
        left_sign = np.sign(filtered_signal[current_samp + 1] - filtered_signal[current_samp])
        new_left_sign = left_sign
        # go to left as long as same direction
        while current_samp > 0 and (left_sign == new_left_sign) or (new_left_sign == 0):
            lefts += 1
            current_samp -= 1
            new_left_sign = np.sign(filtered_signal[current_samp + 1] - filtered_signal[current_samp])
        
        current_samp = samp + 1
        rights = 0
        right_sign = np.sign(filtered_signal[current_samp - 1] - filtered_signal[current_samp])
        new_right_sign = right_sign
        # go to right as long as same direction
        while current_samp < len(filtered_signal) - 1 and (right_sign == new_right_sign) or (new_right_sign == 0):
            rights += 1
            current_samp += 1
            new_right_sign = np.sign(filtered_signal[current_samp - 1] - filtered_signal[current_samp])

        min_ramps.append(min([lefts, rights]))

    to_sort_df = pd.DataFrame(np.transpose([samples, min_ramps]))
    to_sort_df = to_sort_df[to_sort_df[1] > 7] # ~ 7 -> 100 ms # absolute cutoff criterion
    to_sort_df = to_sort_df.sort_values(1, ascending=False)
    sorted_samples = to_sort_df[0].tolist()

    # find the cutoff for the trial
    samples_perc = []
    samples_abs = []
    if len(sorted_samples) > 0:
        try:
            sorted_samples = sorted_samples[:10]
        except:
            sorted_samples = sorted_samples[:-1]
        for i, lm in enumerate(sorted_samples):
            # or not "error" in behav_dimension
            samples_perc.append(lm/len(trial_behav_dim))
            samples_abs.append(lm)
    
    n_samples_EEG_trial    = dur_trial_samp   #dur_trials_in_samp[trial_ind]
    start_sample_EEG_trial = start_trial_samp #start_trials_samp[trial_ind]
    
    # debug plot (adjust as necessary)

    # currently plots the first trial
    """
    if debug:
        return samples_abs
        if trial_ind == 0:
            debug_plot(trial_behav_dim, samples_abs)
    #"""
    
    return [
        samples_abs, 
        [int(perc*n_samples_EEG_trial) + start_sample_EEG_trial for perc in samples_perc], 
        [(e/len(filtered_signal))*n_samples_EEG_trial + start_sample_EEG_trial for e in range(len(filtered_signal))]
        ]

# integrative function
def behavioral_events_to_EEG(sbj_ID, task, tol, debug=False):
    # writes new marker file

    # 1. load the marker file
    df_marker = load_marker_file(sbj_ID, task)
    
    # 2. derive the start samples and trial durations
    start_trials_samp, dur_trials_in_samp = trial_sample_infos_from_marker_df(df_marker)
            
    # for all trials, read in the behavioral files
    behavioral_folders = [f for f in subdirs if f"{sbj_ID}" in f and f"_{task}" in f][0]
    
    behav_dims = ["error_abs", "purs-y"]
    all_EEG_samples = []#[[] for _ in range(len(behav_dims))]
    
    for trial_ind in range(len(dur_trials_in_samp)):
        behav_trial = load_trial(behavioral_folders, trial_ind + 1)
        dur_trial_samp   = dur_trials_in_samp[trial_ind]
        start_trial_samp = start_trials_samp[trial_ind]
        
        xyz = behavioral_events_in_EEG_format(
                behav_trial,
                behav_dims[0],
                dur_trial_samp,
                start_trial_samp,
                trial_ind,
                False
            )
        EEG_samples_trial_error = xyz[1]
        EEG_samples_trial_error_for_plot = xyz[0]

        xyz = behavioral_events_in_EEG_format(
                behav_trial,
                behav_dims[1],
                dur_trial_samp,
                start_trial_samp,
                trial_ind,
                False
            )
        EEG_samples_trial_behav = xyz[1]
        EEG_samples_trial_behav_for_plot = xyz[0]
        
        time_in_samples = xyz[2]

        # looking for matches with some temporal tolerance
        EEG_samples_trial_match = []
        match_inds = []
        for i, samp in enumerate(EEG_samples_trial_behav):
            for samp_e in EEG_samples_trial_error:
                if samp in range(samp_e - tol, samp_e + tol + 1):
                    EEG_samples_trial_match.append(samp)
                    match_inds.append(i)
                    break
        EEG_samples_trial_match_for_plot = [s for i, s in enumerate(EEG_samples_trial_behav_for_plot) if i in match_inds]

        #EEG_samples_trial_match = [samp for samp in EEG_samples_trial_error if samp in EEG_samples_trial_behav]
        #EEG_samples_trial_match_for_plot = [samp for samp in EEG_samples_trial_error_for_plot if samp in EEG_samples_trial_behav_for_plot]

        # total debug plot
        # if len(EEG_samples_trial_match_for_plot) > 0:
        #     debug_plot_new(behav_trial, EEG_samples_trial_match_for_plot)
        if trial_ind == 0:
            debug_plot_new(behav_trial, EEG_samples_trial_match, time_in_samples)
            #display(xyz[2])

        # append to continuous list of all markers
        all_EEG_samples += EEG_samples_trial_match
        
    # bring into file format readable by BVA
    write_markers(42, all_EEG_samples, sbj_ID, task)

    print("==========================")
    print(len(all_EEG_samples))
    print("==========================")

    """
    for i, behav_dim in enumerate(behav_dims):
        write_markers(behav_dim, 40+i, all_EEG_samples[i], sbj_ID, task)
    """