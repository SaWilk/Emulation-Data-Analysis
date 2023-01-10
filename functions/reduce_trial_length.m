function shorter_trials = reduce_trial_length(ft_structure, trl_start, trl_end)
% This function reduces the trial length in a ft structure to the wished
% time window, indicated by trl_start and trl_end.
% Adriana Böttcher, 2022

% copy fieldtrip trial structure
shorter_trials = ft_structure;

% find out which time index equals trl_start and trl_end
% start_ind = find(shorter_trials.time{1} == trl_start);
% end_ind = find(shorter_trials.time{1} == trl_end);

% new: determine first value which is bigger than trial start/end
start_ind = find(shorter_trials.time{1} > trl_start, 1, 'first');
end_ind = find(shorter_trials.time{1} > trl_end, 1, 'first');

% loop though trials in data and reduce every trial to trl_start:trl_end
for trl = 1:size(shorter_trials.trial, 2)
    shorter_trials.trial{trl} = shorter_trials.trial{trl}(:, start_ind:end_ind);
    shorter_trials.time{trl} = shorter_trials.time{trl}(:, start_ind:end_ind);
end 
end