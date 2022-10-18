function [vin, vineeg] = vincentiles_eeg(dat,eegdata,nb_bins)

% vincentiles_eeg.m
%
% Computes nb_bins vincentiles for the data vector 'dat' (from condtrials.RT) and sorts and
% averages the eegdata according to the bins. Needs output from
% collectconds_bar13! eegdata shape must (!) be timepoint x trials
% 
% Originally written by:  Trisha van Zandt
% adapted by Sven Hoffmann for EEG-vincentiles 15/07/2013


% sort data and get length and dimensions
[dat, index] = sort(dat);
n = length(dat);
neeg = size(eegdata,1);
dims = size(dat);
dimseeg = size(eegdata);
if dims(2) == 1,
   dat = dat';
end

% sort eegdata
tmpeeg=zeros(dimseeg);
for ind = 1: length(index)
   tmpeeg(:,ind) = eegdata(:,index(ind));
end;


% create temporary arrays with nb_bins replications of each value
copyx = ones(nb_bins,1)*dat;
vinx = reshape(copyx,n,nb_bins);

copyxeeg = repmat(tmpeeg',[1,1,nb_bins]);

vinxeeg = shiftdim(copyxeeg,2);
vinxeeg = reshape(vinxeeg,dimseeg(2),nb_bins,neeg);


% obtain mean values for each (RT EEG)bin
vin = mean(vinx)';
vineeg = squeeze(mean(vinxeeg,1));


