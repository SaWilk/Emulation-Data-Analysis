function [event_idx] = get_epoch_center(EEG,ep)
%GET_EPOCH_CENTER Get the event that is central to an EEGLAB epoch
% Author: Saskia Wilken
% Creation Date: 21.09.22
%   [event_idx] = get_epoch_center(EEG, ep)
%   EEG = EEGLAB EEG structure
%   ep = current epoch we are interested in
%   event_idx = index of central event of epoch in event structure
        % get central peak of epoch
        event_idx = EEG.epoch(ep).event([EEG.epoch(ep).eventlatency{:}] == 0);
        % in the rare case that there are multiple events with a latency of
        % 0, determine which one is the peak
        if length(event_idx) > 1
            which_peak_idx = find(strcmp({EEG.event(event_idx).type}, {'S 40'}) );
            which_peak_idx_2 = find(strcmp({EEG.event(event_idx).type}, {'S 50'}) );
            if ~isempty(which_peak_idx)
                event_idx = event_idx(which_peak_idx);
            elseif ~isempty(which_peak_idx_2)
                event_idx = event_idx(which_peak_idx_2);
            end
        end
end