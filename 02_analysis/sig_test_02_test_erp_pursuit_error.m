% The neurophysiology of continuous action monitoring
% Saskia Wilken, Adriana Böttcher, Nico Adelhöfer, Markus Raab, Sven
% Hoffmann & Christian Beste

% T-Tests and Plotting of Vincentile Results
% Calculates T-Tests with Bootstrapping, Effect Sizes and plots the results
% in a) singificance matrices for EEG x Behavioral Data bins and b) into
% vincentile per effect size plots for behavioral data only

% Created by: 
% Saskia Wilken, General Psychology: Judgement, Decision Making & Action, 
% University of Hagen
% 02.11.22


%% Empty everything

format compact
format long G
clear
clc
close all


%% Define Paths

% file location
file_path = matlab.desktop.editor.getActiveFilename;
file_path = fileparts(file_path);
addpath(file_path);

%% Load Data

load stats_struct.mat
load vin_struct.mat


%% T-Tests

% set alpha level for t-tests
ALPHA = 0.001;
FDR_RATE = 0.05; % the q value for the FDR function
contrasts = {"all", "constrand", "occvis"};
N = 5000;
df = N-1;

for con = 2:length(contrasts)
    % switch to right eeg conditions
    CONTRAST = contrasts{con};
    conds = {"all", "all"};
    if strcmp(CONTRAST, "occvis")
        conds = {"occluded", "visible"};
    elseif strcmp(CONTRAST, "constrand")
        conds = {"constant", "random1"};
    end
    for ch = conditionxclust{con-1}
        disp("clusters")
        disp([ch, " ch"])

        CLUST = clusters{ch};
        CLUST_LAB = cluster_labels{ch};
        for var = 1:2
            if strcmp(conds{1}, "all")
                error("didn't specify this condition yet.")
            else
                data_1 = [];
                data_2 = [];
                data = {};
                for s = 1:size(ALLEEG, 2)
                    SUB = strjoin(["subject_", s],"");


                    %% Calculate means of each eeg data bin for subject and vincentile

                    for bin = 1:length(eeg_bins_samp)-1
                        data_1(:, bin, s) = mean(vin_struct.(sortvar{var}).( ...
                            conds{1}).(SUB).(CLUST_LAB)(:, eeg_bins_samp(bin):eeg_bins_samp(bin+1)), 2);
                        data_2(:, bin, s) = mean(vin_struct.(sortvar{var}).( ...
                            conds{2}).(SUB).(CLUST_LAB)(:, eeg_bins_samp(bin):eeg_bins_samp(bin+1)), 2);
                    end
                    behav_data_1(:, s) = vin_struct.(sortvar{var}).( ...
                        conds{1}).(SUB).behav;
                    behav_data_2(:, s) = vin_struct.(sortvar{var}).( ...
                        conds{2}).(SUB).behav;
                end
                data{1} = data_1;
                data{2} = data_2;
                behav_data{1} = behav_data_1;
                behav_data{2} = behav_data_2;
            end 


            %% T-tests with NB_BOOT bootstrap samples

            clear stats df pvals surrog h crit_p adj_ci_cvrg adj_p surrog surrog_behav
            [stats, df, pvals, surrog] = statcond(data, 'paired', 'on', 'method', 'bootstrap', ...
                'naccu', BOOT_REP, 'alpha', ALPHA, 'structoutput', 'on');
            [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals, FDR_RATE, 'dep', 'yes');

            stats_struct.(sortvar{var}).(CONTRAST).(CLUST_LAB) = stats;
            stats_struct.(sortvar{var}).(CONTRAST).(CLUST_LAB).("pvals") = pvals;
            stats_struct.(sortvar{var}).(CONTRAST).(CLUST_LAB).("mask_corr") = h;
            stats_struct.(sortvar{var}).(CONTRAST).(CLUST_LAB).("adj_pval") = adj_p;
%             stats_struct.(sortvar{var}).(CONTRAST).behav.("boot_data") = surrog;
            stats_struct.(sortvar{var}).(CONTRAST).(CLUST_LAB).("robust_d") = nan(size(stats.t));
            % t-tests for behavioral data only
            clear stats df pvals  h crit_p adj_ci_cvrg adj_p boot_mean1 boot_mean2 boot_std1 boot_std2
            [stats, df, pvals, surrog_behav] = statcond(behav_data, 'paired', 'on', 'method', 'bootstrap', ...
                'naccu', BOOT_REP, 'alpha', ALPHA, 'structoutput', 'on');
            [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals, FDR_RATE, 'dep', 'yes');
            stats_struct.(sortvar{var}).(CONTRAST).behav = stats;
            stats_struct.(sortvar{var}).(CONTRAST).behav.("pvals") = pvals;
            stats_struct.(sortvar{var}).(CONTRAST).behav.("mask_corr") = h;
            stats_struct.(sortvar{var}).(CONTRAST).behav.("adj_pval") = adj_p;
            stats_struct.(sortvar{var}).(CONTRAST).behav.("boot_data") = surrog_behav;
            stats_struct.(sortvar{var}).(CONTRAST).behav.("robust_d") = nan(size(stats.t));
            
            for vin = 1:size(stats_struct.(sortvar{var}).(CONTRAST).(CLUST_LAB).t, 1)
                for bin = 1:size(stats_struct.(sortvar{var}).(CONTRAST).(CLUST_LAB).t, 2)
                    clear tmp sav boot_mean1 boot_mean2 boot_std1 boot_std2
                    boot_mean1 = bootstrp(BOOT_REP, @mean, data{1}(bin, vin, :));
                    boot_mean2 = bootstrp(BOOT_REP, @mean, data{2}(bin, vin, :));
                    boot_std1 = bootstrp(BOOT_REP, @std, data{1}(bin, vin, :));
                    boot_std2 = bootstrp(BOOT_REP, @std, data{2}(bin, vin, :));
                    stats_struct.(sortvar{var}).(CONTRAST).(CLUST_LAB).("boot_mean1")(bin, vin) = mean(boot_mean1);
                    stats_struct.(sortvar{var}).(CONTRAST).(CLUST_LAB).("boot_mean2")(bin, vin) = mean(boot_mean2);
                    stats_struct.(sortvar{var}).(CONTRAST).(CLUST_LAB).("boot_std1")(bin, vin) = mean(boot_std1);
                    stats_struct.(sortvar{var}).(CONTRAST).(CLUST_LAB).("boot_std2")(bin, vin) = mean(boot_std2);
                    % calculate effect sizes for eeg data tests
                    % pooled variance
                    sav = sqrt((df*mean(boot_std1)^2 + df*mean(boot_std2)^2) / (df*2));
                    % effect size
                    tmp = ((mean(boot_mean1) - mean(boot_mean2))/sav);
                    stats_struct.(sortvar{var}).(CONTRAST).(CLUST_LAB).robust_d(bin, vin) = tmp;
                end
                % calculate effect sizes for behav data tests
                clear tmp boot_mean1 boot_mean2 boot_std1 boot_std2 sav
                boot_mean1 = bootstrp(BOOT_REP,@mean,behav_data{1}(vin, :));
                boot_mean2 = bootstrp(BOOT_REP,@mean,behav_data{2}(vin, :));
                boot_std1 = bootstrp(BOOT_REP,@std,behav_data{1}(vin, :));
                boot_std2 = bootstrp(BOOT_REP,@std,behav_data{2}(vin, :));
                stats_struct.(sortvar{var}).(CONTRAST).behav.("boot_mean1")(vin) = mean(boot_mean1);
                stats_struct.(sortvar{var}).(CONTRAST).behav.("boot_mean2")(vin) = mean(boot_mean2);
                stats_struct.(sortvar{var}).(CONTRAST).behav.("boot_std1")(vin) = mean(boot_std1);
                stats_struct.(sortvar{var}).(CONTRAST).behav.("boot_std2")(vin) = mean(boot_std2);
                % calculate effect size by hand
                % pooled variance
                sav = sqrt((df*mean(boot_std1)^2 + df*mean(boot_std2)^2) / (df*2));
                % effect size
                tmp = ((mean(boot_mean1) - mean(boot_mean2))/sav);
                stats_struct.(sortvar{var}).(CONTRAST).behav.robust_d(vin) = tmp;
            end
        end % sortvar
    end % chans
end % contrasts


%% Plot the results EEG per Behavioral

% image with the RT bins on one axis and the time bin on the other
% Colors are the effect sizes, only significant

% generate axis labels for imagesc
% x axis
for bin = 1:length(eeg_bins_ms)-1
    x_axis_tick_labels{bin} = strjoin([...
        num2str(eeg_bins_ms(bin)), " ms - ", num2str(eeg_bins_ms(bin+1)), " ms"], "");
end

for con = 3
    % switch to right eeg conditions
    CONTRAST = contrasts{con};
    conds = {"all", "all"};
    test_type = "one-sample";
    if strcmp(CONTRAST, "occvis")
        conds = {"occluded", "visible"};
        test_type = "two-sample";
    elseif strcmp(CONTRAST, "constrand")
        conds = {"constant", "random1"};
        test_type = "two-sample";
    end
    for ch = conditionxclust{con-1}
        disp("clusters")
        disp([ch, " ch"])

        CLUST = clusters{ch};
        CLUST_LAB = cluster_labels{ch};

        for var = 1:2
            BEHAV_ONLY = false;

            plot_significance_matrix(stats_struct.(sortvar{var}).(CONTRAST).(CLUST_LAB).robust_d, ...
                stats_struct.(sortvar{var}).(CONTRAST).(CLUST_LAB).mask_corr, ...
                sortvar{var}, ...
                eeg_bins_samp, ...
                NB_BINS, ...
                test_type, ...
                CONTRAST, ...
                CLUST_LAB, ...
                ALPHA, ...
                x_axis_tick_labels, ...
                BOOT_REP, ...
                BEHAV_ONLY)

        end % sortvar
    end % chans
end % contrast


%% Plot Behavioral Data Vincentiles in Combination with Effect Size

FONTSIZE = 26;
LINEWIDTH = 2;

% open text file for saving t-test results
fileID = fopen('t-test_results_behavioral_II.txt','w');

for var = 1:2
    SORTVAR = sortvar{var};
    for con = 2:length(contrasts)
        % switch to right eeg conditions
        CONTRAST = contrasts{con};
        conds = {"all", "all"};
        if strcmp(CONTRAST, "occvis")
            conds = {"occluded", "visible"};
        elseif strcmp(CONTRAST, "constrand")
            conds = {"constant", "random1"};
        end
        for s = 1:size(ALLEEG, 2)
            SUB = strjoin(["subject_", s],"");
            behav_data_1(:, s) = vin_struct.(SORTVAR).( ...
                conds{1}).(SUB).behav;
            behav_data_2(:, s) = vin_struct.(SORTVAR).( ...
                conds{2}).(SUB).behav;
        end % subject
        grand_average_1 = mean(behav_data_1, 2)
        grand_average_2 = mean(behav_data_2, 2)


        %% Plot parameters
        if strcmp(SORTVAR, "pursuit_lat")
            ylab = "pursuit latency";
            yunit = "ms";
            LOWER_EXTRA = 40;
            y_lims = [50, 300];
            y_lims_without_extra = y_lims(1)+100;
        elseif strcmp(SORTVAR, "epoch_error_pix")
            ylab = "epoch error";
            yunit = "pix";
            LOWER_EXTRA = 130;
            y_lims = [-50, 225];
            y_lims_without_extra = y_lims(1)+100;
        end

        UPPER_EXTRA = 0;
        y_ticks = round(linspace(y_lims_without_extra, y_lims(2), 6))

        sign_d = stats_struct.(SORTVAR).(CONTRAST).behav.robust_d;
        tmp = sign_d;
        tmp(~stats_struct.(SORTVAR).(CONTRAST).behav.mask_corr) = nan;
        sign_coord = [(1:NB_BINS)', tmp];
        second_ylim = [0, 10];


        %% Plot code
        figure()
        set(gcf,'color','w')
        % plot data
        plot(1:NB_BINS, grand_average_1, '-o', 'LineWidth', LINEWIDTH*1.5, 'Color', 'k', 'MarkerFaceColor','k')
        hold on
        plot(1:NB_BINS, grand_average_2, '-o', 'LineWidth', LINEWIDTH*1.5, 'Color', 'b', 'MarkerFaceColor','b')
        %         yline(y_ticks(1), 'LineStyle', '--', 'LineWidth', 1)

        yyaxis right
        ba = bar(1:NB_BINS, abs(sign_d), 0.25, ...
            'FaceColor',[0 .5 .5], 'LineStyle', 'None');
        ba.FaceAlpha = 0.5;
        % asterisks
        plot(sign_coord(:,1),abs(sign_coord(:, 2))+0.2, '*k', 'MarkerSize',12, 'LineWidth',LINEWIDTH)
        % axis
        % right y
        ylim(second_ylim);
        if max(sign_d) < 0
            y_labels = string(linspace(0, 2, 3))
            yticks(linspace(0, 2, 3))
            yticklabels(y_labels);
        else
            yticks(linspace(0, 2, 3))
        end
        yl = ylabel("effect size (d)")
        yl.Position = [yl.Position(1), yl.Position(2) - second_ylim(2)*0.33];
        set(gca,'ycolor','k')
        % x axis
        xlabel("vincentiles")
        % left y
        yyaxis left
        ylabel(strjoin([ylab, " (", yunit, ")"], ""))
        ylim(y_lims);
        yticks(y_ticks)
        set(gca,'FontSize', FONTSIZE, 'LineWidth', LINEWIDTH)
        legend({conds{1}, conds{2}}, 'location', 'northwest')
        % title
        title(strjoin(["Vincentiles of ", SORTVAR, " with significant effect sizes"], ""), ...
            'Interpreter', 'None')
        subtitle(strjoin(["Contrast: ", CONTRAST, ", Behavioral: ", ...
            SORTVAR, " Alpha: ", ALPHA], ""), 'Interpreter', 'None')
        box off
        % T-Test Results to file
        fprintf(fileID,'%24s \n', strjoin(["contrast: ", CONTRAST, " sortvar: ", SORTVAR], ""));
        % T-Test Results to file
        fprintf(fileID,'%4s %7s %7s %8s %6s %12s %8s \n','vin', 'mean_1', 'mean_2', 't-value', 'df', 'fdr p-value', 'cohen-d');
        % t-test report
        test_report = [(1:10)', ...
            grand_average_1, ...
            grand_average_2, ...
            stats_struct.(SORTVAR).(CONTRAST).behav.t, ...
            repmat(stats_struct.(SORTVAR).(CONTRAST).behav.df, [10, 1]), ...
            stats_struct.(SORTVAR).(CONTRAST).behav.adj_pval, ...
            stats_struct.(SORTVAR).(CONTRAST).behav.robust_d]';

        fprintf(fileID,'%4.0f %7.2f %7.2f %8.2f %6.2f %8.2f %8.2f \r\n',test_report);

        % save figure
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 8])
            filename = strjoin(["Vincentiles and Effect Size of ", SORTVAR, ...
                " Contrast ", CONTRAST], "");
        print('-djpeg', strjoin([filename, ".jpg"],""), '-r600');
        print('-dsvg', strjoin([filename, ".svg"],""), '-vector');
    end
end % var
fclose(fileID);


%% Plot the Epoch Error on Pursuit Lat Vincentiles

SORTVAR = sortvar{1};
FONTSIZE = 26;
LINEWIDTH = 2;

for con = 2:length(contrasts)
    % switch to right eeg conditions
    CONTRAST = contrasts{con};
    conds = {"all", "all"};
    if strcmp(CONTRAST, "occvis")
        conds = {"occluded", "visible"};
    elseif strcmp(CONTRAST, "constrand")
        conds = {"constant", "random1"};
    end
    for s = 1:size(ALLEEG, 2)
        SUB = strjoin(["subject_", s],"");
        % the untransformed behavioral data by all the epoch that have a
        % valid pursuit latency, now condition 1
        cur_behav1 = vin_struct.(sortvar{1}).(conds{1}).(SUB).unsorted_behav;
        cur_behav2 = vin_struct.(sortvar{2}).(conds{1}).(SUB).unsorted_behav;
        [com_struct.(conds{1}).(SUB).("behav"), ...
            com_struct.(conds{1}).(SUB).("plot_behav"), ...
            com_struct.(conds{1}).(SUB).("index"), ...
            com_struct.(conds{1}).(SUB).("behav_data")] = vincentiles_eeg( ...
            cur_behav1, cur_behav2, NB_BINS); % get pursuit latency sorted by epoch error bins for condition 1
        clear cur_behav1 cur_behav2
        % the untransformed behavioral data by all the epoch that have a
        % valid pursuit latency, now condition two
        cur_behav1 = vin_struct.(sortvar{1}).(conds{2}).(SUB).unsorted_behav;
        cur_behav2 = vin_struct.(sortvar{2}).(conds{2}).(SUB).unsorted_behav;
        [com_struct.(conds{2}).(SUB).("behav"), ...
            com_struct.(conds{2}).(SUB).("plot_behav"), ...
            com_struct.(conds{2}).(SUB).("index"), ...
            com_struct.(conds{2}).(SUB).("behav_data")] = vincentiles_eeg( ...
            cur_behav1, cur_behav2, NB_BINS); % get pursuit latency sorted by epoch error bins
        % the epoch error sorted by pursuit latency vincentile
        behav_data_1(:, s) = com_struct.(conds{1}).(SUB).plot_behav;
        behav_data_2(:, s) = com_struct.(conds{2}).(SUB).plot_behav;
    end % subject
    behav_data{1} = behav_data_1;
    behav_data{2} = behav_data_2;

    % perform statistical test
    [stats, df, pvals, surrog_behav] = statcond(behav_data, 'paired', 'on', 'method', 'bootstrap', ...
        'naccu', BOOT_REP, 'alpha', ALPHA, 'structoutput', 'on');
    [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals, FDR_RATE, 'dep', 'yes');

    for vin = 1:size(behav_data{1}, 1)
        clear tmp boot_mean1 boot_mean2 boot_std1 boot_std2 sav
        boot_mean1 = bootstrp(BOOT_REP,@mean,behav_data{1}(vin, :));
        boot_mean2 = bootstrp(BOOT_REP,@mean,behav_data{2}(vin, :));
        boot_std1 = bootstrp(BOOT_REP,@std,behav_data{1}(vin, :));
        boot_std2 = bootstrp(BOOT_REP,@std,behav_data{2}(vin, :));
        sav = sqrt((df*mean(boot_std1)^2 + df*mean(boot_std2)^2) / (df*2));
        % effect size
        tmp = ((mean(boot_mean1) - mean(boot_mean2))/sav);
        com_struct.(CONTRAST).("cohen_d")(vin) = tmp;
    end
    grand_average_1 = mean(behav_data{1}, 2);
    grand_average_2 = mean(behav_data{2}, 2);

    %% Plot parameters

    ylab = "epoch error";
    yunit = "pix";
    y_lims = [50, 125];
    y_lims_without_extra = y_lims(1)+25;

    UPPER_EXTRA = 0;
    y_ticks = round(linspace(y_lims_without_extra, y_lims(2), 6))

    sign_d = com_struct.(CONTRAST).("cohen_d");
    tmp = sign_d;
    tmp(~h) = nan;
    sign_coord = [(1:NB_BINS)', tmp'];

    second_ylim = [0, 10];


    %% Plot code

    figure()
    set(gcf,'color','w')
    % plot data
    plot(1:NB_BINS, grand_average_1, '-o', 'LineWidth', LINEWIDTH*1.5, 'Color', 'k', 'MarkerFaceColor','k')
    hold on
    plot(1:NB_BINS, grand_average_2, '-o', 'LineWidth', LINEWIDTH*1.5, 'Color', 'b', 'MarkerFaceColor','b')
    %         yline(y_ticks(1), 'LineStyle', '--', 'LineWidth', 1)

    yyaxis right
    ba = bar(1:NB_BINS, abs(sign_d), 0.25, ...
        'FaceColor',[0 .5 .5], 'LineStyle', 'None');
    ba.FaceAlpha = 0.5;
    % asterisks
    plot(sign_coord(:,1),abs(sign_coord(:, 2))+0.2, '*k', 'MarkerSize',14, 'LineWidth',LINEWIDTH)
    % axis
    % right y
    ylim(second_ylim);
    if max(sign_d) < 0
        %             y_labels = string(linspace(0, floor(min(sign_d)), 5));
        y_labels = string(linspace(0, 2, 3))
        %             yticks(linspace(0, abs(floor(min(sign_d))), 5));
        yticks(linspace(0, 2, 3))
        yticklabels(y_labels);
    else
        yticks(linspace(0, 2, 3))
    end
    yl = ylabel("effect size (d)")
    yl.Position = [yl.Position(1), yl.Position(2) - second_ylim(2)*0.33];
    set(gca,'ycolor','k')
    % x axis
    xlabel("pursuit latency (ms)")
    % left y
    yyaxis left
    ylabel(strjoin([ylab, " (", yunit, ")"], ""))
    ylim(y_lims);
    yticks(y_ticks)
    set(gca,'FontSize', FONTSIZE, 'LineWidth', LINEWIDTH)
    legend({conds{1}, conds{2}}, 'location', 'northwest')
    % title
    title(strjoin(["Vincentiles of ", sortvar{2}, " sorted by vincentiles of ", ...
        sortvar{1}, " with significant effect sizes"], ""), ...
        'Interpreter', 'None')
    subtitle(strjoin(["Contrast: ", CONTRAST, ", Behavioral: ", ...
        SORTVAR, " Alpha: ", ALPHA], ""), 'Interpreter', 'None')
    box off

    % save figure
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 8])
    filename = strjoin(["Vincentiles and Effect Size of ", sortvar{2}, ...
        " by ", sortvar{1}, " Contrast ", CONTRAST], "");
    print('-djpeg', strjoin([filename, ".jpg"],""), '-r600');
    print('-dsvg', strjoin([filename, ".svg"],""), '-vector');
end
