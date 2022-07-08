function EEG = rm_occl_none(EEG)

rmnone = strmatch('none',{EEG.event.TRAJ}');
EEG = pop_editeventvals(EEG,'delete',rmnone);
clear rmnone:;

return