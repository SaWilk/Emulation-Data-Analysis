function [OUTEEG, HRV] = addRpeak2eeglab(INEEG,ecgchannum,srate,thresh)

OUTEEG=INEEG;

x = INEEG.data(ecgchannum,:);

[hrv, R_t, R_amp, R_index, S_t, S_amp]  = rpeakdetect(x',srate,thresh,0);

 
le = length(R_index);
for index = 1 : le
   OUTEEG.event(end+1) = OUTEEG.event(end);
   OUTEEG.event(end).type='R-peak'; % Add event to end of event list
   OUTEEG.event(end).latency =R_index(index); 
end
%OUTEEG = eeg_checkset(OUTEEG, 'eventconsistency'); 
HRV.hrv=hrv;
HRV.R_t= R_t;
HRV.R_amp=R_amp;
HRV.R_index= R_index;
HRV.S_t=S_t;
HRV.S_amp=S_amp;

eeglab redraw;

return
    
    