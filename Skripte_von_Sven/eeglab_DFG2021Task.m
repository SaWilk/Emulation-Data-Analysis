% %================================================================%
% % eeglab_zuerich2021                                     %
% % Reads in data                        %
% % defines channels, filters, resamples, , cleans data, saves the contineous  data   %       %
% % X=Subject (integer);                       %
% %                                                              %
% % Sven Hoffmann 17.07.2021                                      %
% %================================================================%
clear;
file_path = matlab.desktop.editor.getActiveFilename;
file_dir = fileparts(file_path);
parent_dir = fileparts(file_dir);

eeglab;
close;
%change current dir to data dir
cd(strcat(parent_dir, "\raw_EEG_2021"));

%ldata path for data loading
pathname = strcat(parent_dir, "\raw_EEG_2021");

% data path for data saving
savepath = strcat(parent_dir, "\processed_EEG_2021");  

%list all *.vhdr files in director
filenames = dir('*.vhdr');

%concatenate into one cell array
files2read = {filenames.name};

% define empty BADATA CELL for collecting names of bad data
BADDATA = {};

for ind= 1 : length(filenames)    
    %define condition
           
            
%           %import 
            TMPEEG = pop_loadbv(pathname, files2read{ind});
          
           % add conditional codes for random and constant trajectory
            TMPEEG = addcodingevents(TMPEEG);
            TMPEEG = rmnoneDFG2021(TMPEEG);

            %define setname
            TMPEEG.setname = files2read{ind}(1:end-4);
            
            %store original filename in set
            TMPEEG.comment = files2read{ind};
            TMPEEG.group = 'Control';
            
            
            %define condition
            TMPEEG.condition = 'TaskA';
            
            % define subject
            substring = files2read{ind}(1:end-7);
            TMPEEG.subject = substring;
            
            
            % Edit channel info, re-reference to average reference
            TMPEEG=pop_chanedit(TMPEEG, 'append',60,'changefield',{61,'labels','FPz'});
            TMPEEG=pop_chanedit(TMPEEG, 'setref',{'1:60','FPz'});
            

            % downsampling
            TMPEEG = pop_resample( TMPEEG, 250);
            TMPEEG.setname = files2read{ind}(1:end-7);

            %hp filter
            TMPEEG = pop_eegfiltnew(TMPEEG, 1, []);
     
          
            %remove line noise
           TMPEEG = pop_cleanline(TMPEEG, 'bandwidth',2,'chanlist', 1:TMPEEG.nbchan ,...
               'computepower',1,'linefreqs',[50 100],'normSpectrum',0,'p',0.01,'pad',2,...
               'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,...
               'winsize',4,'winstep',4);
          
            %remove bad channels 
            TMPEEG.oldchanslocs = TMPEEG.chanlocs;
            TMPEEG = clean_rawdata(TMPEEG, 5, -1, 0.85, -1, -1, -1);
            
            
            %lowpass filte5
			TMPEEG = pop_eegfiltnew(TMPEEG, [], 40);

            %interpolate removed channels
            TMPEEG = pop_interp(TMPEEG, TMPEEG.oldchanslocs, 'spherical');
            
            %rereference to average reference
            TMPEEG = pop_reref( TMPEEG, [ ]);
            
            
            % get data for amica ICA
            %if data are epoched:
            %x = double(reshape(TMPEEG.data,TMPEEG.nbchan,TMPEEG.pnts*TMPEEG.trials));
           %if data are contineous: 
            %x = double(TMPEEG.data);
            
            % estimate rank of the data (necessary since data channels have
            % been interpolated, custom function by Sven Hoffmann (part of
            % eeglab) 
            
            %prepare ICA via data subset
            ICA = eeg_regepochs(TMPEEG);

            %detrend eeg data
            ICA = eeg_detrend(ICA);

            %artifact rejection
            ICA = pop_jointprob(ICA,1,1:ICA.nbchan ,5,5,0,1);

            %select only some random trials
            trl = 1:ICA.trials;
            trl = shuffle(trl);
            ICA = pop_select(ICA, 'trial',trl(1:round(length(trl)/2))) ;
            
            %x = double(reshape(TMPEEG.data,TMPEEG.nbchan,TMPEEG.pnts*TMPEEG.trials));
            %if data are contineous: 
            %select random subsample of trials
             x = double(ICA.data);

            %reshape for ICA
             x = reshape(x,size(x,1),size(x,2)*size(x,3));

            %get rank
             rnk = getrank(x);
             
             % now run amica on subset of data 
             %[TMPEEG.icaweights,TMPEEG.icasphere,mods] = runamica15(x,'do_reject',0,'pcakeep',rnk);
             %TMPEEG.icawinv = pinv(TMPEEG.icaweights * TMPEEG.icasphere); 
             %TMPEEG=eeg_checkset(TMPEEG);

             %extende ifomax (more efficeint choice, amica, bowver, yoield
             %more dipolar components)
              ICA = pop_runica(ICA, 'icatype', 'runica', 'extended',1, 'pca',rnk);
              TMPEEG.icaweights = ICA.icaweights;
              TMPEEG.icasphere = ICA.icasphere;
              TMPEEG.icawinv = ICA.icawinv; 
             TMPEEG.icachansind = ICA.icachansind;
             TMPEEG = eeg_checkset(TMPEEG);

             %save contineous ica  data
             TMPEEG = pop_saveset(TMPEEG,'filename',[files2read{ind}(1:end-7) '_runica'], 'filepath',savepath);
            
             %remove ocular artifacts
             %TMPEEG = dfgremICs2021(TMPEEG,false);
            
                      
             % save ocular artifact cleaned data
             %TMPEEG=pop_saveset(TMPEEG,'filename',[files2read{ind}(1:end-7) '_amicaclean'], 'filepath',savepath);
             
             % make epochs for constant and random trajectory
            TMPEEG = pop_epoch( TMPEEG, {  'S 23'  'S 27'  }, [-1  4], 'newname', [TMPEEG.setname '-eps'], 'epochinfo', 'yes');
            TMPEEG = pop_rmbase( TMPEEG, [-1000 0] ,[]);
            TMPEEG = pop_saveset(TMPEEG,'filename',[files2read{ind}(1:end-7) '-eps'], 'filepath',savepath);
           %end 
           pause(0.30);
end % end file loop