import pandas as pd
import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelextrema
import scipy
from scipy import signal

# paths
folder_marker = Path(r"R:\AG-Beste-Studien\Tracking Task\pilot_EEG\export")
path_behav = Path(r"R:\AG-Beste-Studien\Tracking Task\pilot_behav\data")
path_new_marker = Path(r"R:\AG-Beste-Studien\Tracking Task\pilot_EEG\new_behavioral_markers")

# overview of all relevant marker files
all_marker_files = [f for f in os.listdir(folder_marker) if f.endswith(".vmrk")]
marker_files = []
for f in all_marker_files:
    for i in range(7,11):
        if f"Pilot{i}" in f:
            marker_files.append(f)

# overview of all relevant behavioral folders
all_subdirs = [f.name for f in os.scandir(path_behav) if f.is_dir()]
subdirs = []
for d in all_subdirs:
    for i in range(7,11):
        if f"Pilot{i}" in d:
            subdirs.append(d)

def load_trial(foldername, trialnum):
    # trialnum should start at 1
    all_files = os.listdir(path_behav/foldername)
    # debug
    """
    filename = ""
    for f in all_files:
        part1 = f.split("_")[-1].split(".")[0]
        part2 = str(trialnum)][0]
        if part1 == part2:
            filename = f
            break
    #"""
    filename = [f for f in all_files if f.split("_")[-1].split(".")[0] == str(trialnum)][0]
    return pd.read_csv(path_behav/foldername/filename)

# write new markers to file
def write_markers_dev(markers, pilot_num, task):
    f = open(path_new_marker/f"Pilot{pilot_num}_{task}_behavioral_markers.txt", "w")
    f.writelines("Sampling rate: 256Hz, SamplingInterval: 3.90625ms\n")
    f.writelines("Type, Description, Position, Length, Channel\n")
    for m in markers:
        f.writelines(f"Stimulus, S 40, {m}, 1, All\n")
    f.close()
    
def write_markers_dirchange(markers, pilot_num, task):
    f = open(path_new_marker/f"Pilot{pilot_num}_{task}_dirchange.txt", "w")
    f.writelines("Sampling rate: 256Hz, SamplingInterval: 3.90625ms\n")
    f.writelines("Type, Description, Position, Length, Channel\n")
    for m in markers:
        f.writelines(f"Stimulus, S 41, {m}, 1, All\n")
    f.close()
    
def write_markers(f_suffix, marker_label, markers, pilot_num, task):
    f = open(path_new_marker/f"Pilot{pilot_num}_{task}_{f_suffix}.txt", "w")
    f.writelines("Sampling rate: 256Hz, SamplingInterval: 3.90625ms\n")
    f.writelines("Type, Description, Position, Length, Channel\n")
    for m in markers:
        f.writelines(f"Stimulus, S {marker_label}, {m}, 1, All\n")
    f.close()

def debug_plot(behavioral_data, local_extrema):
    plt.figure()
    plt.plot(behavioral_data, "black")
    plt.vlines([v for v in local_extrema], min(behavioral_data), max(behavioral_data), "red")
    plt.show()

def load_marker_file(pilot_num, task):
    fname = f"Pilot{pilot_num}_{task}_Topographic Interpolation.Markers"
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

def get_extrema(behav_trial_data, only_max=False):

    # low-pass filter to cushion behavioral irregularities
    # sos = signal.iirfilter(2, 0.05, btype='lowpass', output='sos')
    # filtered_signal = scipy.signal.sosfiltfilt(sos, behav_trial_data.tolist())

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
    samples      = get_extrema(trial_behav_dim, "error" in behav_dimension)[0]
    samples_vals = get_extrema(trial_behav_dim, "error" in behav_dimension)[1]
    
    # find the cutoff for the trial
    dev_samples_vals_sorted = sorted(samples_vals, reverse=True)
    cutoff_dev = dev_samples_vals_sorted[5]
    
    samples_perc = []
    samples_abs = []
    for i, lm in enumerate(samples):
        if samples_vals[i] > cutoff_dev or not "error" in behav_dimension:
            samples_perc.append(lm/len(trial_behav_dim))
            samples_abs.append(lm)
    
    n_samples_EEG_trial    = dur_trial_samp   #dur_trials_in_samp[trial_ind]
    start_sample_EEG_trial = start_trial_samp #start_trials_samp[trial_ind]
    
    # debug plot (adjust as necessary)

    # currently plots the first trial
    if debug:
        if trial_ind == 0:
            debug_plot(trial_behav_dim, samples_abs)
    
    return [int(perc*n_samples_EEG_trial) + start_sample_EEG_trial for perc in samples_perc]
    
    """ # old
    trial_abs_errors = np.array(behav_trial["error_abs"].tolist())
    trial_purs_ys = behav_trial["purs-y"].tolist()

    # derive both indices and values at indices (keep the n largest deviances)
    dev_samples = get_local_min_max(trial_abs_errors)[0]
    dev_samples_vals = get_local_min_max(trial_abs_errors)[1]

    dirchange_samples = new_improved_max_min(trial_purs_ys, False)[0]
    dirchange_samples_vals = new_improved_max_min(trial_purs_ys, False)[1]

    # find the cutoff for the trial
    dev_samples_vals_sorted = sorted(dev_samples_vals, reverse=True)
    cutoff_dev = dev_samples_vals_sorted[5]

    # only include large deviations (as per the previously identified cutoff)
    dev_samples_perc = []
    dev_samples_abs = []
    for i, lm in enumerate(dev_samples):
        if dev_samples_vals[i] > cutoff_dev:
            dev_samples_perc.append(lm/len(trial_abs_errors))
            dev_samples_abs.append(lm)
    dirchange_samples_perc = []
    dirchange_samples_abs = []
    for i, lm in enumerate(dirchange_samples):
        dirchange_samples_perc.append(lm/len(trial_purs_ys))
        dirchange_samples_abs.append(lm)

    #dev_samples_perc = [lm/len(trial_abs_errors) for i, lm in enumerate(dev_samples) if dev_samples_vals[i] > cutoff_dev]

    # adjust these values to EEG sample units
    n_samples_EEG_trial = dur_trials_in_samp[trial_ind]
    start_sample_EEG_trial = start_trials_samp[trial_ind]

    EEG_dev_samples_trial = [int(perc*n_samples_EEG_trial) + start_sample_EEG_trial for perc in dev_samples_perc]
    EEG_dirchange_samples_trial = [int(perc*n_samples_EEG_trial) + start_sample_EEG_trial for perc in dirchange_samples_perc]
    #"""

# integrative function
def behavioral_events_to_EEG(pilot_num, task, debug=False):
    # writes new marker file

    # 1. load the marker file
    df_marker = load_marker_file(pilot_num, task)
    
    # 2. derive the start samples and trial durations
    start_trials_samp, dur_trials_in_samp = trial_sample_infos_from_marker_df(df_marker)
            
    # for all trials, read in the behavioral files
    behavioral_folders = [f for f in subdirs if f"Pilot{pilot_num}" in f and f"-{task}-" in f][0]
    
    behav_dims = ["error_abs", "purs-y"]
    all_EEG_samples = [[] for _ in range(len(behav_dims))]
    
    for trial_ind in range(len(dur_trials_in_samp)):
        behav_trial = load_trial(behavioral_folders, trial_ind + 1)
        dur_trial_samp   = dur_trials_in_samp[trial_ind]
        start_trial_samp = start_trials_samp[trial_ind]
        
        EEG_samples_trial_error = behavioral_events_in_EEG_format(
                behav_trial,
                behav_dims[0],
                dur_trial_samp,
                start_trial_samp,
                trial_ind,
                debug
            )

        EEG_samples_trial_behav = behavioral_events_in_EEG_format(
                behav_trial,
                behav_dims[1],
                dur_trial_samp,
                start_trial_samp,
                trial_ind,
                debug
            )

        EEG_samples_trial_match = [samp for samp in EEG_samples_trial_error if samp in EEG_samples_trial_behav]

        """
        for i, behav_dim in enumerate(behav_dims):
            EEG_samples_trial = behavioral_events_in_EEG_format(
                behav_trial,
                behav_dim,
                dur_trial_samp,
                start_trial_samp,
                trial_ind,
                debug
            )
            all_EEG_samples[i] += EEG_samples_trial
        """
        
    # bring into file format readable by BVA
    for i, behav_dim in enumerate(behav_dims):
        write_markers(behav_dim, 40+i, all_EEG_samples[i], pilot_num, task)