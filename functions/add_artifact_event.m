function [new_event] = add_artifact_event(event,trial_event_idx, INVALID_SAMP)
%add_artifact_event Adds the number of events each trial to be considered
%as artifactual
%   event               = event structure from eeglab structures
%   trial_event_size    = indices that give only the events of the current
% trial
%   INVALID_SAMP        = when the invalid samples are taking place
%   new_event           = event structure with additional field

% add field: to be ignored samples
artifact_error_cell = cell(size(trial_event_idx));
artifact_error_log = [];
new_event = event;
% if the event belongs to the first 125 samples of a trial,
% ignore it.
artifact_error_log = ...
    [event(trial_event_idx).("trial_latency")] <= INVALID_SAMP(end);
artifact_error_cell(artifact_error_log) = {'ARTIFACT'};
artifact_error_cell(~artifact_error_log) = {'VALID'};
[new_event(trial_event_idx).("artifact_error")] = ...
    artifact_error_cell{:};

end