function EEG = rmnoneDFG2021(EEG)

rmnone = strmatch('none',{EEG.event.TRAJ}');
EEG = pop_editeventvals(EEG,'delete',rmnone);
clear rmnone:;

return