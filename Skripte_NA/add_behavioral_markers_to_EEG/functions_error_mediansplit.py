import pandas as pd
import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelextrema
import scipy
from scipy import signal
from tqdm import tqdm

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
    
def write_markers(marker_label, markers, sbj_ID, task, condition):
    f = open(path_new_marker/f"pp_AM_{sbj_ID}_{task}_{condition}.txt", "w")
    f.writelines("Sampling rate: 256Hz, SamplingInterval: 3.90625ms\n")
    f.writelines("Type, Description, Position, Length, Channel\n")
    for m in markers:
        f.writelines(f"Stimulus, S {marker_label + condition}, {m}, 1, All\n")
    f.close()

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
    return scipy.signal.sosfiltfilt(sos, np.array(task_sig.tolist()))

def get_behavioral_sample_reappear(marker_df, df_trials): #df_trials based on trial_sample_infos_from_marker
    df_reappear = df[df[1].isin(["S 21"])] #df from load_marker_file
    
    for  
    
def mediansplit():

def behavior_to_EEG():

