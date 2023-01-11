function [sig_elec, labels] = CBPT_get_sig_elec_from_timewindow(CBPT_output)

% This function extracts all electrodes assigned to significant clusters in
% the significant time window of a CBPT.
% This requires: output of a CBPT, not averaged over time dimension, but
% over frequency dimension; alpha.
% Output is: a list of channels that are 'significant' in the respective
% time window.
% Adriana Böttcher, 2022

% prepare electrode list: 60 zeros, later to be assigned to eleclabels
sig_elec = zeros(60, 1);

% extract mask
mask = squeeze(CBPT_output.mask);

% if electrode has a 1 somewhere, mark as significant
sig_elec = sum(mask, 2) > 0;

% list of electrode labels
labels = {};
for e = 1:length(sig_elec)
    if sig_elec(e) == 1
        labels{end+1} = CBPT_output.label{e};
    end
end
end