function EEG = eeg_detrend(EEG)

disp('eeg_detrend(): Detrending EEG data');

EEG.data = reshape(permute(EEG.data, [2 1 3]), [EEG.pnts EEG.nbchan * EEG.trials]);
EEG.data = detrend(EEG.data);
EEG.data = permute(reshape(EEG.data, [EEG.pnts EEG.nbchan EEG.trials]), [2 1 3]);

end