%% Signfiicant test of ERP and Latency resp. Epoch error

% Author: Saskia Wilken
% Creation Date: 18.10.2022

%% Vincentile Data Tracking Data

% An alternative to quantiles are Vincentiles, which are computed by sorting 
% the data and splitting them in equi-populated bins (there is the same 
% number of observations in each bin). Then the mean is computed for each bin
% In the present study, 10 bins were calculated. Finally, ERLs
% were computed for each bin

vincentile_eeg()

%% Vincentile ERP Data in Time Window of Interest (around P2)

% Finally,ERLs
% were computed for each bin(see Figure 3B) according to the pro-
% cedure already described (for each ocniditon)
% wheretheERLeffectsweremostpromi-
% nent, wasselectedforstatisticalquantification.Morespecifically,
% the timerangefrom160to510mswasdividedinto10segments
% (35mseach),fromwhichthecorrespondingaverageswerecom-
% puted
% allaveraged
% segments in the derived timewindow of interest in the RT bins for
% each experimental condition (except the NoGoORI condition),
% were tested against zero via non-parametric bootstrapping (Efron
% and Tibshirani,1993) using1000bootstrapsamplesforeachtest.
% Allderived p-values wereadjustedbymeansoffalsediscoveryrate
% accordingto Benjamini andYekutieli(2001).

bootstrp(nboot,bootfun,d) 

mafdr(PValues)


%% Plot the results

% scimage iwth the RT bins on one axis and the time bin on the other
% Colors are the thresholds of T-values, also binned into significance
% levels



