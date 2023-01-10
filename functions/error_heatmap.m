function error_heatmap(track_data, frame_rate)
%ERROR HEATMAP Plots the error per trial segment
%   Takes the track_data structure and frame rate of not upsampled data and
%   returns a plot of the error per
%   trial segment as a heatmap. For now, automatically calculates the error
%   for the first trial segment respectively, from trial start to 2s. 
% Author: Saskia Wilken

clear all_track_error
segment = [1:frame_rate*2]; % how much data to plot
cnt = 1;
for s = 1:length(track_data)
    for task = 1:length(track_data(s).trials)
        for t = 1:length(track_data(s).trials(task).task_a)
            clear tmp
            tmp = track_data(s).trials(task).task_a(t).error;
            all_track_error(:, cnt) = abs(tmp(segment))*1080;
            cnt = cnt + 1;
        end
    end
end

clear freqs
[limits(1), limits(2)] = bounds(all_track_error, 'all');
resolution = 170;
edges = linspace(limits(1),limits(2), resolution);
for tp = 1:size(all_track_error,1)
    freqs(:,tp) = histcounts(all_track_error(tp,:), edges);
end

times = round(linspace(1, 2000, size(all_track_error,1)));
y_limits = find(prctile(edges, 0) <= edges & prctile(edges, 25) > edges);

%% Create Figure

figure()
imagesc(times, edges(y_limits), freqs(y_limits,:), [15, 200])
set(gca, 'YDir','normal','fontsize', 26);
xline(500, 'LineWidth',4, 'LineStyle', '--', 'Color', 'w');
colormap jet
xlabel('time (ms)')
ylabel('error (pixels)')
a = colorbar;
set(get(a,'label'),'string','frequency', 'FontSize', 26)
set(a, 'LineWidth', 2, 'FontSize', 26)
title('Frequency of error magnitude per frame of all trials')


%% Save Figure
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 13 8], 'color','w')
set(gca, 'LineWidth', 2)
box off

    filename = strjoin("Heatmap Error Distribution", "");
print('-djpeg', strjoin([filename, ".jpg"],""), '-r600');
print('-dsvg', strjoin([filename, ".svg"],""), '-vector');

end