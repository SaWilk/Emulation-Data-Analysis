%% Vincentile EEG Data According to Performance Measures of Purusit

% Combines all data in one datasets, vincentiles the eeg data according to
% ehavioral data and plots ERPimages that show how EEG data changes based
% on the different performance measures
% Author: Saskia Wilken
% Creation Date: 18.10.22

% get data paths for parent dirs
filepath_parts = strsplit(file_path, filesep);
parent_dir_2 = filepath_parts(1:end-2);
parent_dir_2 = strjoin(parent_dir_2, filesep);
output_dir = strjoin([parent_dir_2, "Emulation-Data-Output"], filesep);
% output
csd_dir = strjoin([output_dir, "csd_transform"], filesep);


%% Load CDS data

eeglab;
close all

% long epochs
cd(csd_dir);
%list all *.set files in inputpath
file_names = dir('*.set');
%concatenate into one cell array
files2read = {file_names.name};
files2read = files2read(~contains(files2read, 'KMY6K')); % remove subject who
% did not properly track during task B
% load eeg data epoched around peaks
for idx = 1:length(files2read)
    ALLEEG2(idx) = pop_loadset('filename',files2read{idx});
end

ALLEEG = ALLEEG2;
clear ALLEEG2


%% Set Parameters

load clust_chans.mat

contrasts = {"all", "constrand", "occvis"};
conditionxclust = [{[1, 2]}, ...
    {3:length(clust)}];
cluster_labels = {"constrand_neg_3", "constrand_neg_1", "occvis_neg_1", ...
    "occvis_neg_3", "occvis_pos_1", "occvis_pos_2", "occvis_pos_3"};
% NO LONGER USING CHANS
% chans = {"P3", "FC3"}; % Which channel do I want to look at
% the constant - random1 cluster and the occluded - visible cluster
clusters = clust;
NB_BINS = 10; % How many bins to use for vincentiles
EEG_BINS = 10;
% Prepare vincentlie structure
vin_struct = struct("pursuit_lat", [], "epoch_error_pix", []);
sortvar = fieldnames(vin_struct);
vin_struct.time_vec = ALLEEG(1).times;
% Prepare stats structure
stats_struct = struct("pursuit_lat", [], "epoch_error_pix", []);
stats_struct.time_vec = ALLEEG(1).times;
SAMP_DIST = 4;
% I will use the 100 ms to 600 ms time range
EEG_RANGE = [50, 550];
clear tmp tmp2
[~, tmp] = min(abs(vin_struct.time_vec - EEG_RANGE(1)));
[~, tmp2] = min(abs(vin_struct.time_vec - EEG_RANGE(2)));
eeg_bins_samp = round(linspace(tmp, tmp2, EEG_BINS+1));
eeg_bins_ms = round(linspace(EEG_RANGE(1), EEG_RANGE(2), EEG_BINS+1));
BOOT_REP = 5000;
COMPARE = true;


%% Vincentile Data Tracking Data

% An alternative to quantiles are Vincentiles, which are computed by sorting
% the data and splitting them in equi-populated bins (there is the same
% number of observations in each bin). Then the mean is computed for each bin

% vincentiles of single-trial EEG data
% were calculated. Basically, vincentiles are calculated by first sort-
% ing a data vector in increasing order and then separating it into
% the desired number of bins.Subsequently,the data are averaged
% within these bins.

% Step by step:
% Tracking Error and Pursuit Latency were collected from online event triggers in the EEG
% by removing all unimportant event field entries and using EEG.event.[...]

cd(all_subj_sets_dir);
%% Remove all unimportant event field entries

for s = 1:length(ALLEEG)
    ALLEEG(s).data = ALLEEG(s).CSD_data;
    tmpEEG(s) = rmfield(ALLEEG(s), "CSD_data");
end

ALLEEG = tmpEEG;
clear tmpEEG

clear z_ALLEEG TMPEEG z_mean_struct mean_struct CONDEEG CONDEEG2 TMPEEG2
sum(arrayfun(@(s) length(s.epoch), ALLEEG)); % how many epochs we need.
s = 1;
while s <= length(ALLEEG)
    EEG = ALLEEG(s);
    del = 0;
    event_idx = zeros(length(EEG.epoch),1);
    for ep = 1:length(EEG.epoch)
        event_idx(ep) = get_epoch_center(EEG, ep);
        %         if ~(strcmp(EEG.event(ev-del).type, 'S 40') | strcmp(EEG.event(ev-del).type, 'S 50'))
        %             EEG.event(ev-del) = [];
        %             del = del +1;
        %         end
    end
    log_idx = unfind(event_idx, length(EEG.event));
    EEG.event(~log_idx) = [];
    EEG = eeg_checkset(EEG);
    ALLEEG(s) = EEG;
    s = s + 1;
end
% WATCH OUT: FROM NOW ON get_epoch_center won't work anymore since the
% EEG.epoch structure is no longer aligned with the EEG.event structure.


%% Create Condition-Specific ALLEEG Datasets

for s = 1:length(ALLEEG)
    EEG = ALLEEG(s);
    EEG = pop_selectevent( EEG, 'artifact_error',{'VALID'},'deleteevents',...
        'off','deleteepochs','on','invertepochs','off'); % remove the first
    % 500 ms peaks from data because we want to associate eeg with behavioral data
    ALLEEG(s) = EEG;
end

% loop that creates condition-specific datasets with only one event per
% epoch (peak)
for con = 2:length(contrasts) % we are not interested in "all"
    disp("contrast")
    CONTRAST = contrasts{con};
    CONDEEG = ALLEEG;

    switch CONTRAST
        case "occvis"
            % remove occluded VISIBLE
            for s = 1:length(CONDEEG)
                EEG = CONDEEG(s);
                EEG = pop_selectevent( EEG, 'task', {'task_b'},'deleteevents',...
                    'off','deleteepochs','on','invertepochs','off');
                EEG = pop_selectevent( EEG, 'OCCL',{'OFF'},'deleteevents',...
                    'off','deleteepochs','on','invertepochs','off');
                CONDEEG(s) = EEG;
                cond_lab = "visible";
                [CONDEEG(s).event.("subject")] = deal(s);

            end
            CONDEEG2 = ALLEEG;
            % remove visible OCCLUDED
            for s = 1:length(CONDEEG2)
                EEG = CONDEEG2(s);
                EEG = pop_selectevent( EEG, 'OCCL',{'ON'},'deleteevents',...
                    'off','deleteepochs','on','invertepochs','off');
                CONDEEG2(s) = EEG;
                cond_lab2 = "occluded";
                [CONDEEG2(s).event.("subject")] = deal(s);

            end
        case "constrand"
            %     remove random CONSTANT
            CONDEEG = ALLEEG;

            for s = 1:length(CONDEEG)
                EEG = CONDEEG(s);
                EEG = pop_selectevent( EEG, 'task', {'task_a'},'deleteevents',...
                    'off','deleteepochs','on','invertepochs','off');
                EEG = pop_selectevent( EEG, 'TRAJ',{'RANDOM1', 'RANDOM2'},'deleteevents',...
                    'off','deleteepochs','on','invertepochs','off');
                CONDEEG(s) = EEG;
                cond_lab = "constant";
                [CONDEEG(s).event.("subject")] = deal(s);

            end
            CONDEEG2 = ALLEEG;
            % remove constant RANDOM1
            for s = 1:length(CONDEEG2)
                EEG = CONDEEG2(s);
                EEG = pop_selectevent( EEG, 'task', {'task_a'},'deleteevents',...
                    'off','deleteepochs','on','invertepochs','off');
                EEG = pop_selectevent( EEG, 'TRAJ',{'RANDOM2', 'CONST'},'deleteevents',...
                    'off','deleteepochs','on','invertepochs','off');
                CONDEEG2(s) = EEG;
                cond_lab2 = "random1";
                [CONDEEG2(s).event.("subject")] = deal(s);

            end
    end % switch contrast


    %% Put all subjects' data in the same dataset 1

    % initialize dataset
    TMPEEG = CONDEEG(1);
    TMPEEG.data = CONDEEG(1).data;
    % get all data and all event structures
    for s = 2:length(CONDEEG)
        z_tmpEEG = CONDEEG(s);
        z_tmpEEG.data = [];
        z_tmpEEG.data = CONDEEG(s).data;
        % concatenate subject data in one dataset
        TMPEEG.data = cat(3, z_tmpEEG.data, TMPEEG.data);
        % get epoch vectors
        epoch_vec = TMPEEG.event(end).epoch + [z_tmpEEG.event.epoch];
        % get event indices
        event_idx = size(TMPEEG.event,2)+1:size(TMPEEG.event,2)+size(z_tmpEEG.event, 2);
        % add events from previous subject to TMPEEG
        TMPEEG.event(:, event_idx) = z_tmpEEG.event;
        for ev = 1:length(event_idx)
            TMPEEG.event(:, event_idx(ev)).epoch = epoch_vec(ev);
        end
        TMPEEG.epoch(:, size(TMPEEG.epoch,2)+1:...
            size(TMPEEG.epoch,2)+size(z_tmpEEG.epoch, 2)) = z_tmpEEG.epoch;
        epoch_vec = [];
    end
    TMPEEG.trials = length(TMPEEG.epoch);

    clear vec
    vec = [TMPEEG.event.epoch_error]*1080;
    for i = 1:length(vec)
        TMPEEG.event(i).('epoch_error_pix') = vec(i);
    end
    clear CONDEEG


    %% Put all subjects' data in the same dataset 2

    if ~strcmp(CONTRAST, "all")
        TMPEEG2 = CONDEEG2(1);
        TMPEEG2.data = CONDEEG2(1).data;
        % put all epochs of all subjects in one dataset
        % get all data and all event structures
        for s = 2:length(CONDEEG2)
            z_tmpEEG = CONDEEG2(s);
            z_tmpEEG.data = [];
            z_tmpEEG.data = CONDEEG2(s).data;
            % concatenate subject data in one dataset
            TMPEEG2.data = cat(3, z_tmpEEG.data, TMPEEG2.data);
            %      TMPEEG.data(:,:,data_idx) =
            % get epoch vectors
            epoch_vec = TMPEEG2.event(end).epoch + [z_tmpEEG.event.epoch];
            % get event indices
            event_idx = size(TMPEEG2.event,2)+1:size(TMPEEG2.event,2)+size(z_tmpEEG.event, 2);
            % add events from previous subject to TMPEEG
            TMPEEG2.event(:, event_idx) = z_tmpEEG.event;
            for ev = 1:length(event_idx)
                TMPEEG2.event(:, event_idx(ev)).epoch = epoch_vec(ev);
            end
            TMPEEG2.epoch(:, size(TMPEEG2.epoch,2)+1:...
                size(TMPEEG2.epoch,2)+size(z_tmpEEG.epoch, 2)) = z_tmpEEG.epoch;
            epoch_vec = [];
        end
        TMPEEG2.trials = length(TMPEEG2.epoch);
        clear vec
        vec = [TMPEEG2.event.epoch_error]*1080;
        for i = 1:length(vec)
            TMPEEG2.event(i).('epoch_error_pix') = vec(i);
        end
    end

    clear only_data z_tmpEEG vec CONDEEG2


    %% Delete Epochs without valid Pursuit Latency

    TMPEEG_purs = TMPEEG;

    del_idx = [];
    count = 1;
    for ep = 1:length([TMPEEG_purs.event])
        if isempty(TMPEEG_purs.event(ep).pursuit_lat)
            del_idx(count) = TMPEEG_purs.event(ep).epoch;  % get all the epochs without valid pursuit latency entry
            count = count + 1;
        end
    end
    missing_purs_lat = length(del_idx);
    disp(strjoin(["Out of ", size(TMPEEG.data, 3), " epochs, ", ...
        missing_purs_lat, " do not have a valid pursuit latency (", ...
        missing_purs_lat/(size(TMPEEG.data, 3)/100), " %)"], ""));
    % Condition "all" yields:
    % Out of 64737 epochs, 29214 do not have a valid pursuit latency (45.1272 %)
    TMPEEG_purs.data(:, :, del_idx) = []; % note that I copied CSD_data to data
    TMPEEG_purs.event(del_idx) = []; % note that I copied CSD_data to data
    for ev = 1:length([TMPEEG_purs.event.epoch])
        TMPEEG_purs.event(ev).epoch = ev;
    end

    if ~strcmp(CONTRAST, "all")
        TMPEEG2_purs = TMPEEG2;

        del_idx = [];
        count = 1;
        for ep = 1:length([TMPEEG2_purs.event])
            if isempty(TMPEEG2_purs.event(ep).pursuit_lat)
                del_idx(count) = TMPEEG2_purs.event(ep).epoch;  % get all the epochs without valid pursuit latency entry
                count = count + 1;
            end
        end
        TMPEEG2_purs.data(:, :, del_idx) = []; % note that I copied CSD_data to data
        TMPEEG2_purs.event(del_idx) = []; % note that I copied CSD_data to data
        for ev = 1:length([TMPEEG2_purs.event.epoch])
            TMPEEG2_purs.event(ev).epoch = ev;
        end
    end

    % Initially, RTs were sorted and the corresponding
    % trial number was collected, also. Following this, the single-trial
    % data of the ERL channels of interest (i.e.,PO7/PO8) were sorted
    % according to the RTs, using the collected trial numbers. Then, RTs
    % and single-trial EEG data were binned (separately for each chan-
    % nel).% In the present study, 10 bins were calculated. Finally, ERLs
    % were computed for each bin


    %% Determine which variables and conditions and set parameters accordingly

    clear tmp
    [~, tmp] = min(abs(vin_struct.time_vec +200));
    plot_area = tmp:length(vin_struct.time_vec);
    clear tmp

    for ch = conditionxclust{con-1}
        disp("clusters")
        disp([ch, " ch"])

        CLUST = clusters{ch};
        CLUST_LAB = cluster_labels{ch};
        for lab = 1:length(CLUST)
            chan_idx(lab) = find(strcmp({EEG.chanlocs.labels}, CLUST{lab}));
        end

        for var = 1:2
            disp("sortvar")

            if strcmp(sortvar{var}, "pursuit_lat")
                SORTVAR = "pursuit latency";
            elseif strcmp(sortvar{var}, "epoch_error_pix")
                SORTVAR = "epoch error";
            end

            % Get z_lims / is more SET z_lims now
            for c = 1:2
                conds = {"all", "all"};
                if strcmp(CONTRAST, "occvis")
                    conds = {"occluded", "visible"};
                elseif strcmp(CONTRAST, "constrand")
                    conds = {"constant", "random1"};
                end
                EEG_COND = conds{c};
                clear EEG_struct
                if strcmp(SORTVAR, "pursuit latency")
                    EEG_struct = TMPEEG_purs;
                    ylab = "pursuit latency (ms)";
                    if ~strcmp(CONTRAST, "all")
                            if strcmp(EEG_COND, 'visible')
                                EEG_struct = TMPEEG_purs;
                            end
                            if strcmp(EEG_COND, 'occluded')
                                EEG_struct = TMPEEG2_purs;
                            end
                            if strcmp(EEG_COND, 'constant')
                                EEG_struct = TMPEEG_purs;
                            end
                            if strcmp(EEG_COND, 'random1')
                                EEG_struct = TMPEEG2_purs;
                            end
                    end
                elseif strcmp(SORTVAR, "epoch error")
                    EEG_struct = TMPEEG;
                    ylab = "epoch error (pix)";
                    if ~strcmp(CONTRAST, "all")
                            if strcmp(EEG_COND, 'visible')
                                EEG_struct = TMPEEG_purs;
                            end
                            if strcmp(EEG_COND, 'occluded')
                                EEG_struct = TMPEEG2_purs;
                            end
                            if strcmp(EEG_COND, 'constant')
                                EEG_struct = TMPEEG_purs;
                            end
                            if strcmp(EEG_COND, 'random1')
                                EEG_struct = TMPEEG2_purs;
                            end
                    end
                end


                %% Create vincentiles

                clear ind_subjects_behav ind_subjects_vin grand_average_vin grand_average_behav
                for s = 1:max([EEG_struct.event.subject])
                    SUB = strjoin(["subject_", s],"");
                    dat = [];
                    sub_idx = [];
                    sub_idx = find([EEG_struct.event.subject] == s);
                    dat = [EEG_struct.event(sub_idx).(sortvar{var})];
                    vin_struct.(sortvar{var}).(EEG_COND).(SUB).("unsorted_behav") = dat;
                    eegdata = [];
                    % get mean for one cluster
                    eegdata = squeeze(mean(EEG_struct.data(chan_idx,:,sub_idx)));
                    % get vincentiles
                    [vin_struct.(sortvar{var}).(EEG_COND).(SUB).("behav"), ...
                        vin_struct.(sortvar{var}).(EEG_COND).(SUB).(CLUST_LAB), ...
                        vin_struct.(sortvar{var}).(EEG_COND).(SUB).("index"), ...
                        vin_struct.(sortvar{var}).(EEG_COND).(SUB).("behav_data")] = vincentiles_eeg( ...
                        dat, eegdata, NB_BINS);
                    % get data per sunbject for grand average computation
                    ind_subjects_vin(:, :, s) = vin_struct.(sortvar{var}).(EEG_COND).(SUB).(CLUST_LAB);
                    ind_subjects_behav(:, s) = vin_struct.(sortvar{var}).(EEG_COND).(SUB).("behav");
                end
                % get mean of all subjects per vincentile
                grand_average_vin = mean(ind_subjects_vin, 3);
                grand_average_behav = mean(ind_subjects_behav, 2);


                %% Plot Parameters

                % simple switch to make sure every two following conditions get
                % the same z_lims. 
                if COMPARE
                    %                     z_lims(1) = prctile(grand_average_vin, 5, "all");
                    %                     z_lims(2) = prctile(grand_average_vin, 95, "all");
                    z_lims(1) = min(grand_average_vin, [], "all") - (min(grand_average_vin, [], "all")/100)*20;
                    z_lims(2) = max(grand_average_vin, [], "all") + (max(grand_average_vin, [], "all")/100)*20;
                    z_lims(1:2) = [-max(abs(z_lims)), max(abs(z_lims))];
                    COMPARE = false;
                else
                    COMPARE = true;
                end


                %% Plot ERPimages

                figure()
                set(gcf,'color','w')
                imagesc(vin_struct.time_vec(plot_area), 1:NB_BINS, ...
                    grand_average_vin(:, plot_area), z_lims)
                if strcmp(SORTVAR, "pursuit latency")
                    hold on
                    plot(grand_average_behav, 1:NB_BINS, 'LineWidth', 3, 'Color', 'k')
                    hold off
                end
                title(strjoin(["EEG data sorted by ", SORTVAR, ", condition: ", ...
                    EEG_COND], ""), 'Interpreter', 'None')
                subtitle(strjoin(["Cluster ", CLUST_LAB], ""), 'Interpreter', 'None')
                xlabel("time (ms)")
                ylabel(ylab)
                yticks(1:NB_BINS)
                xline(0, 'LineWidth', 3)
                yticklabels({round(grand_average_behav)})
                set(gca,'YDir','normal')
                a = colorbar;
                position = get(a,'Position');
                ylabel(a,'amplitude (μV / mm²)','Rotation',270, ...
                    'Position', [position(1)+2.5, position(2)-0.1],'FontSize',10);
                set(gca,'YDir','normal')
                set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6])
                filename = strjoin(["Imagesc of Clust ", CLUST_LAB, "Sorted by ", SORTVAR, ...
                    " Condition ", EEG_COND, ".jpg"], "");
                print('-djpeg', filename, '-r600');
                filename = strjoin(["Imagesc of Clust ", CLUST_LAB, "Sorted by ", SORTVAR, ...
                    " Condition ", EEG_COND, ".svg"], "");
                print('-dsvg', filename, '-vector');

                % Finally, ERLs
                % were computed for each bin(see Figure 3B) according to the pro-
                % cedure already described (for each ocniditon)
                % where the ERL effects were most prominent,
                % was selected for statistical quantification. More specifically,
                % the time range from 160 to 510 ms was divided into 10 segments
                % (35 ms each), from which the corresponding averages were computed
                % Now I get the data of each trial b< using the vincentile idx for the time
                % bins and I will use all these unaveraged trials to bootstrap each of them
                % a thousand times. From this distribution of simulated sample means, I
                % will derive t-values and correct them for multiple comparisons.

            end % condition
        end % sortvar
    end % channel
end % contrast

save stats_struct.mat stats_struct
save vin_struct.mat vin_struct
