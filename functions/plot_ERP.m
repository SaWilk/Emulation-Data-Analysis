function [plot_handle, min_ind] = plot_ERP(data, times, base_dur, title_str, subtitle_str)
%plot_ERP Plots ERPs from epoched EEGLAB data
%   Detailed explanation goes here
% data = eeg data to be plotted. 
% base_dur = baseline period duration in ms
% times = time vector of the data
% title_str = main title
% subtitle_str = sub title


[~, min_ind] = min(data,[],'all');
[min_ind(1), min_ind(2)] = ind2sub(size(data), min_ind);
[epoch_lims(1), epoch_lims(2)] = bounds(times);
plot_handle = plot(times, data, 'LineWidth', 1.1);
title(title_str);
subtitle(subtitle_str);
hold on
h(1) = xline(0);
h(2) = xline(times(min_ind(2)), 'r');
% h(4:5) = xline([-avg_peak_dist,avg_peak_dist], 'b');
h(3) = xline(epoch_lims(1) + base_dur, 'g');
legend(h, {'peak', 'lowest deflection',  'baseline period'});%, 'average distance of next peak'});
hold off

end