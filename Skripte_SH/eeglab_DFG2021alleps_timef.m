eeglab;close;

savepath =  '/Users/hoffmann/WORK/EEG-Data/DFG/EEG-Exchange-DFG/TaskA/eeglabdataalleps/';
EEG = pop_loadset('filename',{'038GC_runica-eps-icaclean.set','09WM9_runica-eps-icaclean.set','2EFCY_runica-eps-icaclean.set','5J1DZ_runica-eps-icaclean.set','5STJS_runica-eps-icaclean.set','607SE_runica-eps-icaclean.set','69UUS_runica-eps-icaclean.set','6WFER_runica-eps-icaclean.set','8JX88_runica-eps-icaclean.set','916P5_runica-eps-icaclean.set','91L3H_runica-eps-icaclean.set','BIM2G_runica-eps-icaclean.set','CYLBQ_runica-eps-icaclean.set','DW4MP_runica-eps-icaclean.set','HBKTI_runica-eps-icaclean.set','HGMKL_runica-eps-icaclean.set','HVNAK_runica-eps-icaclean.set','JV4UJ_runica-eps-icaclean.set','KL6JX_runica-eps-icaclean.set','KMY6K_runica-eps-icaclean.set','M0P1Q_runica-eps-icaclean.set','NEYON_runica-eps-icaclean.set','OTQYX_runica-eps-icaclean.set','SGUNS_runica-eps-icaclean.set','TSDAX_runica-eps-icaclean.set','UDU06_runica-eps-icaclean.set','V7EJE_runica-eps-icaclean.set','VMMQ7_runica-eps-icaclean.set','WBKH2_runica-eps-icaclean.set','WLOYV_runica-eps-icaclean.set','WM87B_runica-eps-icaclean.set','Y7J68_runica-eps-icaclean.set','YFE0B_runica-eps-icaclean.set','ZV583_runica-eps-icaclean.set'},'filepath','/Users/hoffmann/WORK/EEG-Data/DFG/EEG-Exchange-DFG/TaskA/eeglabdataalleps/');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0); 

chan = 'FCz';
for ind = 1:34
        TMPEEG = ALLEEG(ind);
        TMPEEG = pop_currentdensity(TMPEEG, 'method','spline');
        channum = find(strcmp(chan,{TMPEEG.chanlocs.labels}));
        traj = []; 
        for i = 1:TMPEEG.trials
            tr = TMPEEG.epoch(i).eventTRAJ{find(cell2mat(TMPEEG.epoch(i).eventlatency)==0)};
            if strmatch(tr,'CONSTANT')
                traj(i,1) = 1;
            elseif strmatch(tr,'RANDOM1')
                traj(i,1) = 0;

            elseif strmatch(tr,'RANDOM2')
                traj(i,1) = 2;
        
            else
                traj(i,1) = NaN;
            end
        end
      for channum = 1:60
            
            [erspconstant(:,:,channum), ~, ~, ersptimes, erspfreqs] = newtimef( TMPEEG.data(channum,:, traj == 1), ...
                TMPEEG.pnts,[TMPEEG.xmin * 1000  TMPEEG.xmax*1000], TMPEEG.srate,[3   0.5] ,  'elocs', ...
                TMPEEG.chanlocs, 'chaninfo', TMPEEG.chaninfo, 'caption', chan, 'basenorm','on', ...
                'freqs', [2 30],  'nfreqs',60,'plotphase','off', 'plotersp','off','plotitc','off',...
                'padratio', 1,'trialbase','full'); close;
            
          
            [ersprandom(:,:,channum), ~, ~, ersptimes, erspfreqs] = newtimef( TMPEEG.data(channum,:, traj == 0), ...
                TMPEEG.pnts,[TMPEEG.xmin * 1000  TMPEEG.xmax * 1000], TMPEEG.srate,[3   0.5] ,  'elocs', ...
                TMPEEG.chanlocs, 'chaninfo', TMPEEG.chaninfo, 'caption', chan, 'basenorm','on', ...
                'freqs', [2 30],  'nfreqs',60,'plotphase','off', 'plotersp','off','plotitc','off',...
                 'padratio', 1,'trialbase','full'); close;
            
            figure;
            [ersprandom2(:,:,channum),~, ~, ersptimes, erspfreqs] = newtimef( TMPEEG.data(channum,:, traj == 2), ...
                TMPEEG.pnts,[TMPEEG.xmin * 1000  TMPEEG.xmax *1000], TMPEEG.srate,[3   0.5] ,  'elocs', ...
                TMPEEG.chanlocs, 'chaninfo', TMPEEG.chaninfo, 'caption', chan, 'basenorm','on', ...
                'freqs', [2 30],  'nfreqs',60,'plotphase','off', 'plotersp','off','plotitc','off',...
                'padratio', 1,'trialbase','full'); 
          
      end
  
     % reoder to match tftopo below
    ERSP.constant{ind} = permute(erspconstant,[2,1,3]);
    ERSP.random{ind} = permute(ersprandom,[2,1,3]);
    ERSP.random2{ind} = permute(ersprandom2,[2,1,3]);

    ERSP.subject{ind} = TMPEEG.subject;
    ERSP.chanlocs{ind} = TMPEEG.chanlocs;
    ERSP.times{ind} = ersptimes;
    ERSP.freqs{ind} = erspfreqs;
    ERSP.TRAJ{ind} = traj;
    clear ersprandom2 ersprandom erspconstant;
 end



for id = 1:34
    average_ersp_constant(:,:,:,id) = ERSP.constant{id};
    average_ersp_random(:,:,:,id) = ERSP.random{id};
    average_ersp_random2(:,:,:,id) = ERSP.random2{id};
end
 ga_ersps_constant = mean(average_ersp_constant,4);
 ga_ersps_random = mean(average_ersp_random,4);
 ga_ersps_random2 = mean(average_ersp_random2,4);


save('ERSP-TaskA-CSD.mat','ERSP','average_ersp_constant','average_ersp_random',...
    'average_ersp_random2','ga_ersps_constant','ga_ersps_random', 'ga_ersps_random2')


 chan = 'FCz';
 channr = find(strcmp(chan,{ERSP.chanlocs{1}.labels}));
  figure; 
   sbplot(1,3,1);imagesc(ersptimes,erspfreqs,ga_ersps_random(:,:,channr)');title([chan '-random']);set(gca,'YDir','normal')
   sbplot(1,3,2);imagesc(ersptimes,erspfreqs,ga_ersps_random2(:,:,channr)');title([chan '-random2']);set(gca,'YDir','normal')
   sbplot(1,3,3);imagesc(ersptimes,erspfreqs,ga_ersps_constant(:,:,channr)');title([chan '-constant']);set(gca,'YDir','normal')

%Stats
 stats_ersps = {average_ersp_random,average_ersp_constant};
     [t df pvals surog] = statcond(stats_ersps, 'method',  'bootstrap', 'naccu', 1000);
     [p_fdr, p_masked] = fdr( pvals,0.0001);
  

 stats_ersps3 = {average_ersp_constant,average_ersp_random2};
    [t3 df3 pvals3 surog3] = statcond(stats_ersps3, 'method',  'bootstrap', 'naccu', 1000);
    [p_fdr3, p_masked3] = fdr( pvals3,0.0001);
    figure; imagesc(ersptimes,erspfreqs,t3(:,:,channum) .* p_masked3(:,:,channum));title( 't(random2-constant)');set(gca,'YDir','normal');
%  
 % TFTOPO
   alltimes = zeros([ size(ersptimes) 60]);
    allfreqs = zeros([ size(erspfreqs) 60]);
    
    for ind=1:60
        alltimes (:,:,ind) = ersptimes;
        allfreqs (:,:,ind) = erspfreqs;
    end

    figure;
    %tfrqdata = permute(ga_ersps_random,[2,1,3]);
    tdata = permute(t,[2,1,3]);
    pdata  = permute(p_masked,[2,1,3]) * 1;
    tftopo(tdata ,alltimes(:,:,1),allfreqs(:,:,1),'chanlocs', ERSP.chanlocs{1}, 'timefreqs', ...
    [100 7; 300 2; 500 3], 'limits', [nan nan nan 35 0 8], 'showchan' ,1);

    figure;
    tdata3 = permute(t3,[2,1,3]);
    pdata3  = permute(p_masked3,[2,1,3]) * 1;
    tftopo(tdata3 ,alltimes(:,:,1),allfreqs(:,:,1),'chanlocs', ERSP.chanlocs{1}, 'timefreqs', ...
    [100 7; 300 2; 500 3], 'limits', [nan nan nan 35 0 8], 'showchan' ,1);
  
  
 
% % 
%  %Statistics

%     stats_ersps = {average_ersp_random,average_ersp_constant};
%     [t df pvals surog] = statcond(stats_ersps, 'method',  'bootstrap', 'naccu', 1000);
%     [p_fdr, p_masked] = fdr( pvals,0.0001);
%     figure; imagesc(ersptimes,erspfreqs,t(:,:,2)' .* p_masked(:,:,2)');title( 't(random-constant)');set(gca,'YDir','normal');
% 
%     stats_ersps2 = {average_ersp_random,average_ersp_random2};
%     [t2 df2 pvals2 surog2] = statcond(stats_ersps2, 'method',  'bootstrap', 'naccu', 1000);
%     [p_fdr2, p_masked2] = fdr( pvals2,0.0001);
%     figure; imagesc(ersptimes,erspfreqs,t2 .* p_masked2);title( 't(random-random2)');set(gca,'YDir','normal');
% % 
% % 
    
%  
