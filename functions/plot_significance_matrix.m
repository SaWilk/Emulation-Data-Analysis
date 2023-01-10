function [] = plot_significance_matrix(plot_data, plot_mask, SORTVAR, ...
    eeg_bins_samp, NB_BINS, test_type, CONTRAST, ...
    CLUST_LAB, ALPHA, x_axis_tick_labels, BOOT_REP, BEHAV_ONLY)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if strcmp(SORTVAR, "pursuit_lat")
    ylab = "pursuit latency";
    yunit = "ms";
elseif strcmp(SORTVAR, "epoch_error_pix")
    ylab = "epoch error";
    yunit = "pix";
end

if BEHAV_ONLY
    c_axis_tick_labels = {"-3", "-2.5", "-2", "-1.5", "-1", ...
        "-0.5", "not significant", "0.5", "1", "1.5", "2", "2.5", "3"};
    z_lims = [-3,3];
    step = 0.5;
else
    z_lims = [-2, 2];
    num_ticks = 10;
%      = (z_lims(2)-z_lims(1))/num_ticks;
    c_axis_tick_labels = unique([round(linspace(z_lims(1), 0, num_ticks/2), 2), round(linspace(0, z_lims(2), num_ticks/2), 2) ])
    step = unique(diff(c_axis_tick_labels))
%     {"-0.04", "-0.03", "-0.02", "-0.01" ...
%         "not sign.", "0.01", "0.02", "0.03", "0.04"};
%     
end
% y axis
% NOTE: This doesn't work since the vincentiles were calculated
% per condition but the plot is shown per contrast, thereore
% the vincentile bin limits are in different locations
%             vin_data1 = round(vin_struct.(sortvar{var}).(conds{1}).behav);
%             vin_data2 = round(vin_struct.(sortvar{var}).(conds{2}).behav);
%
%             for bin = 1:NB_BINS
%                 y_axis_tick_labels{bin} = strjoin([num2str(bin), ": ", ...
%                     num2str(vin_data1(bin)), " ", yunit, " ", conds{1}, ", ", ...
%                     num2str(vin_data2(bin)), " ", yunit, " ", conds{2},], "");
%             end % don't need this since I am now using vincentiles that
%             vary between subjects
% plot_data = stats_struct.(sortvar{var}).(CONTRAST).(CHAN).stats.("tvals");
%             plot_data = stats_struct.(sortvar{var}).(CONTRAST).(CLUST_LAB).robust_d;
%             scale_data = stats_struct.(sortvar{var}).(CONTRAST).(strjoin(CLUST, "")).adj_pval;
% tranform to categorical for plotting

% TODO: Loopify
%             scaled_plot_data = nan(size(plot_data));
%             ind = (plot_data <= 0.0001);
%             scaled_plot_data(ind) = 1;
%             ind = (plot_data > 0.0001 & plot_data <= 0.001) ;
%             scaled_plot_data(ind) = 2;
%             ind = (plot_data > 0.001 & plot_data <= 0.01) ;
%             scaled_plot_data(ind) = 3;
%             ind = (plot_data > 0.01 & plot_data <= 0.05) ;
%             scaled_plot_data(ind) = 4;
%             ind = (plot_data > 0.05);
%             scaled_plot_data(ind) = 5;
% maybe don't even need this
plot_data(~plot_mask) = 0;

LINEWIDTH = 1.5;
FONTSIZE = 14;
PLOTSIZE = [0 , 0, 9, 10];

% plot code:
figure()
set(gcf,'color','w')
%             imagesc(1:length(eeg_bins_samp)-1, 1:NB_BINS, scaled_plot_data, [0.999, 5.001])
imagesc(1:length(eeg_bins_samp)-1, 1:NB_BINS, plot_data, z_lims)
C=jet
middle = round(size(C,1)/2);
C(middle-0:middle+1, :) = 0.8;
%             yellowMap = [linspace(0, 1, 256)', linspace(0, 1, 256)', zeros(256, 1)]
%             C = yellowMap(2)
%             colormap(flipud(C))
colormap(C)
%             a = colorbar('Ticks',linspace(1.5, 4.5, 5), ...
%                 'TickLabels',["> 0.0001", "0.0001 - 0.001", "0.001 - 0.01", "0.01 - 0.5", "> 0.5"]);
%             position = get(a,'Position');
%             ylabel(a,'FDR-corrected p-values','Rotation',270, ...
%                 'Position', [position(1)+2.5, position(2)-0.175],'FontSize',10);
a = colorbar
position = get(a,'Position');
ylabel(a,'cohen d','Rotation',270, ...
    'Position', [position(1)+2.5, position(2)+0.3],'FontSize',FONTSIZE);
set(a, 'LineWidth', LINEWIDTH)
a.TickLabels = c_axis_tick_labels;
a.Ticks = z_lims(1):step:z_lims(2);
xline(1.5:length(eeg_bins_samp)+0.5, 'color', 'w', 'LineWidth', LINEWIDTH )
yline(1.5:NB_BINS+0.5, 'color', 'w', 'LineWidth', LINEWIDTH )
if BEHAV_ONLY
    title(strjoin(["Effect sizes of ", test_type, " two-tailed t-test on ", ...
        BOOT_REP, " Bootstrap Repetitions"], ""), 'Interpreter', 'None')
    subtitle(strjoin(["Contrast: ", CONTRAST, ", Behavioral: ", ...
        SORTVAR, " Only Behavioral, Alpha: ", ALPHA], ""), 'Interpreter', 'None')
else
    title(strjoin(["Effect sizes of ", test_type, " two-tailed t-test on ", ...
        BOOT_REP, " Bootstrap Repetitions"], ""), 'Interpreter', 'None')
    subtitle(strjoin(["Contrast: ", CONTRAST, ", Behavioral: ", ...
        SORTVAR, " Cluster: ", CLUST_LAB, " Alpha: ", ALPHA], ""), 'Interpreter', 'None')
end
% axes
if ~BEHAV_ONLY
    xticklabels(x_axis_tick_labels(1:2:end))
    xticks(1:2:length(eeg_bins_samp))
%     xlabel("eeg data bins")
end

%             ylabel(ylab) % we won't specify the individual vincentiles
%             since they are different for each subject
ylabel(ylab)
%             yticklabels(y_axis_tick_labels)

yticklabels(1:2:NB_BINS)
yticks(1:2:NB_BINS)
set(gca,'YDir','normal')

set(gca,'FontSize', FONTSIZE, 'LineWidth', LINEWIDTH)
box off
% copygraphics(gcf, 'BackgroundColor', 'none', 'ContentType', 'vector');
% pause

% save figure
set(gcf,'PaperUnits','centimeters','PaperPosition',PLOTSIZE)
if BEHAV_ONLY
    filename = strjoin(["Hyp Test Matrix of Behavioral Data Only Sorted by ", SORTVAR, ...
        " Contrast ", CONTRAST], "");
else
    filename = strjoin(["Hyp Test Matrix of Cluster ", ...
        CLUST_LAB, "Sorted by ", SORTVAR, ...
        " Contrast ", CONTRAST], "");
end
print('-djpeg', strjoin([filename, ".jpg"],""), '-r600');
print('-dsvg', strjoin([filename, ".svg"],""), '-vector');

end % end of function