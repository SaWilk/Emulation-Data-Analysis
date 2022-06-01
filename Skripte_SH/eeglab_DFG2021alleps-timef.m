% data        = eeglab2fieldtrip(EEG, ...
%     'raw');
% 
% cfg        = [];
% cfg.metric = 'zvalue';  % use by default zvalue method
% cfg.method = 'summary'; % use by default summary method
% data       = ft_rejectvisual(cfg,data);
% 
% cfg              = [];
% cfg.output       = 'pow';
% cfg.method       = 'mtmconvol';
% cfg.foi          = 3:1:30;
% cfg.channel      = 'FCz';
% cfg.t_ftimwin    = 5./cfg.foi;  % 3 cycles per time window
% %cfg.tapsmofrq  = 0.4 *cfg.foi;
% cfg.taper        = 'hanning';
% cfg.toi          = -1:0.01:3.996;
% cfg.baselinetype = 'relative';
% cfg.baseline=[-1 0];
% TFRhann3 = ft_freqanalysis(cfg, data);
% ft_singleplotTFR(cfg, TFRhann3);
% %ft_multiplotTFR(cfg, TFRhann3);
% 
% 
% 
% [STUDY, ALLEEG] = std_precomp(STUDY, ALLEEG, {'FCz','Cz','Pz','C3','C4','Oz'},'savetrials','on','rmclust',[57 58] ,'interp','on','recompute','on','spec','on','specparams',{'specmode','fft','logtrials','off'},'ersp','on','erspparams',{'cycles',[3 0.5] ,'nfreqs',60,'ntimesout',500,'freqs',[2 30] });
% 


%ersp2read={names.name};

savepath =  '/Users/hoffmann/WORK/EEG-Data/DFG/EEG-Exchange-DFG/TaskA/eeglabdataalleps/';

chan = 'FCz';
for ind = 1:34
    TMPEEG = ALLEEG(ind);
    channum = find(strcmp(chan,{TMPEEG.chanlocs.labels}));
    
    % save data

 

    

    channum = find(strcmp(chan,{TMPEEG.chanlocs.labels}));
    traj = []; 
    for i=1:TMPEEG.trials
        tr = TMPEEG.epoch(i).eventTRAJ{find(cell2mat(TMPEEG.epoch(i).eventlatency)==0)};
        if strmatch(tr,'CONSTANT')
            traj(i,1)=1;
        elseif strmatch(tr,'RANDOM1')
            traj(i,1)=0;

        elseif strmatch(tr,'RANDOM2')
            traj(i,1)=2;
        
        else
            traj(i,1) = NaN;
        end
    end
   
    
    [ERSP.constant{ind}, tmp, tmp, ersptimes, erspfreqs] = newtimef( TMPEEG.data(channum,:, traj == 1), ...
        TMPEEG.pnts,[TMPEEG.xmin*1000  TMPEEG.xmax*1000], TMPEEG.srate,[3   0.5] ,  'elocs', ...
    TMPEEG.chanlocs, 'chaninfo', TMPEEG.chaninfo, 'caption', chan, 'basenorm','on', ...
    'freqs', [2 30],  'nfreqs',60,'plotphase','off', 'plotersp','off','plotitc','off',...
     'padratio', 1,'trialbase','full'); 

    [ERSP.random{ind} tmp, tmp, ersptimes, erspfreqs] = newtimef( TMPEEG.data(channum,:, traj == 0), ...
        TMPEEG.pnts,[TMPEEG.xmin*1000  TMPEEG.xmax*1000], TMPEEG.srate,[3   0.5] ,  'elocs', ...
    TMPEEG.chanlocs, 'chaninfo', TMPEEG.chaninfo, 'caption', chan, 'basenorm','on', ...
    'freqs', [2 30],  'nfreqs',60,'plotphase','off', 'plotersp','off','plotitc','off',...
     'padratio', 1,'trialbase','full'); 

    [ERSP.random2{ind} tmp, tmp, ersptimes, erspfreqs] = newtimef( TMPEEG.data(channum,:, traj == 2), ...
        TMPEEG.pnts,[TMPEEG.xmin*1000  TMPEEG.xmax*1000], TMPEEG.srate,[3   0.5] ,  'elocs', ...
    TMPEEG.chanlocs, 'chaninfo', TMPEEG.chaninfo, 'caption', chan, 'basenorm','on', ...
    'freqs', [2 30],  'nfreqs',60,'plotphase','off', 'plotersp','off','plotitc','off',...
     'padratio', 1,'trialbase','full'); 


    ERSP.subject{ind} = TMPEEG.subject;
    ERSP.TRAJ{ind} = traj;
 end



for id = 1:34
    average_ersp_constant(:,:,id) = ERSP.constant{id};
    average_ersp_random(:,:,id) = ERSP.random{id};
    average_ersp_random2(:,:,id) = ERSP.random{id};
end
 ga_ersps_constant = mean(average_ersp_constant,3);
 ga_ersps_random = mean(average_ersp_random,3);
 ga_ersps_random2 = mean(average_ersp_random2,3);

 figure; 
 sbplot(1,3,1);imagesc(ersptimes,erspfreqs,ga_ersps_random);title([chan '-random']);set(gca,'YDir','normal')
 sbplot(1,3,2);imagesc(ersptimes,erspfreqs,ga_ersps_random2);title([chan '-random2']);set(gca,'YDir','normal')
 sbplot(1,3,3);imagesc(ersptimes,erspfreqs,ga_ersps_constant);title([chan '-constant']);set(gca,'YDir','normal')
 
% 
%  % 
% 
% %  chan = 'Pz';
% % for ind = 1:34
% %     TMPEEG = ALLEEG(ind);
% %     channum = find(strcmp(chan,{TMPEEG.chanlocs.labels}));
% % 
% %     %TMPEEG = iclabel(TMPEEG);
% %     %neurocomps  = find(TMPEEG.etc.ic_classification.ICLabel.classifications(:,1) > 0.5);
% % 
% %     % keep only brain components
% %     %TMPEEG = pop_subcomp( TMPEEG, neurocomps, 0,1);
% %     % save data
% % 
% %     %TMPEEG = pop_saveset(TMPEEG,'filepath',TMPEEG.filepath,'filename',[TMPEEG.setname '-icaclean']);
% % 
% %     % laplacian
% %     %TMPEEG = pop_currentdensity(TMPEEG, 'method','spline');
% % 
% %     %TMPEEG = pop_saveset(TMPEEG,'filepath',TMPEEG.filepath,'filename',[TMPEEG.setname '-CSD']);
% % 
% %     channum = find(strcmp(chan,{TMPEEG.chanlocs.labels}));
% %     traj = []; 
% %     for i=1:TMPEEG.trials
% %         tr = TMPEEG.epoch(i).eventTRAJ{find(cell2mat(TMPEEG.epoch(i).eventlatency)==0)};
% %         if strmatch(tr,'CONSTANT')
% %             traj(i,1)=1;
% %         elseif strmatch(tr,'RANDOM1')
% %             traj(i,1)=0;
% %         else
% %             traj(i,1) = NaN;
% %         end
% %     end
% %    
% %     
% %     [ERSP.constant{ind}, tmp, tmp, ersptimes, erspfreqs] = newtimef( TMPEEG.data(channum,:, traj == 1), ...
% %         TMPEEG.pnts,[TMPEEG.xmin*1000  TMPEEG.xmax*1000], TMPEEG.srate,[3   0.5] ,  'elocs', ...
% %     TMPEEG.chanlocs, 'chaninfo', TMPEEG.chaninfo, 'caption', chan, 'basenorm','on', ...
% %     'freqs', [2 30],  'nfreqs',60,'plotphase','off', 'plotersp','off','plotitc','off',...
% %      'padratio', 1,'trialbase','full'); 
% % 
% %     [ERSP.random{ind} tmp, tmp, ersptimes, erspfreqs] = newtimef( TMPEEG.data(channum,:, traj == 0), ...
% %         TMPEEG.pnts,[TMPEEG.xmin*1000  TMPEEG.xmax*1000], TMPEEG.srate,[3   0.5] ,  'elocs', ...
% %     TMPEEG.chanlocs, 'chaninfo', TMPEEG.chaninfo, 'caption', chan, 'basenorm','on', ...
% %     'freqs', [2 30],  'nfreqs',60,'plotphase','off', 'plotersp','off','plotitc','off',...
% %      'padratio', 1,'trialbase','full'); 
% % 
% %     ERSP.subject{ind} = TMPEEG.subject;
% %     ERSP.TRAJ{ind} = traj;
% % end
% % 
% % 
% % 
% % for id = 1:34
% %     average_ersp_constant(:,:,id) = ERSP.constant{id};
% %     average_ersp_random(:,:,id) = ERSP.random{id};
% % end
% %  ga_ersps_constant = mean(average_ersp_constant,3);
% %  ga_ersps_random = mean(average_ersp_random,3);
% 
%  figure; 
%  sbplot(1,2,1);imagesc(ersptimes,erspfreqs,ga_ersps_random);title([chan '-random']);set(gca,'YDir','normal')
%  sbplot(1,2,2);imagesc(ersptimes,erspfreqs,ga_ersps_constant);title([chan '-constant']);set(gca,'YDir','normal')
% 
% 
%  %Statistics

    stats_ersps = {average_ersp_random,average_ersp_constant};
    [t df pvals surog] = statcond(stats_ersps, 'method',  'bootstrap', 'naccu', 2000);
    [p_fdr, p_masked] = fdr( pvals,0.0001);
    figure; imagesc(ersptimes,erspfreqs,t .* p_masked);title( 't(random-constant)');set(gca,'YDir','normal');

    stats_ersps2 = {average_ersp_random,average_ersp_random2};
    [t2 df2 pvals2 surog2] = statcond(stats_ersps2, 'method',  'bootstrap', 'naccu', 2000);
    [p_fdr2, p_masked2] = fdr( pvals2,0.0001);


    figure; imagesc(ersptimes,erspfreqs,t2 .* p_masked2);title( 't(random-random2)');set(gca,'YDir','normal');
 
%  
