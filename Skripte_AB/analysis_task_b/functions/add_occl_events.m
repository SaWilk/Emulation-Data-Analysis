function OUTEEG = add_occl_events(INEEG)

TMPEEG=INEEG;

%% all relevant triggers

exp_start       = 11;
fix             = 12;
trial_start_L   = 13;
trial_start_R   = 14;
trial_end       = 15; 
pause_start     = 16;
pause_end       = 17;
exp_end         = 18;

occlusion       = 20;
reappear        = 21; 

start_constant  = 23;
end_constant    = 24; 
start_startvec  = 26;
end_startvec    = 27;

%% look for occlusion/reappear triggers in event types
% then add a variable for occlusion = on/off and delete all other triggers

for ind = 1:length(TMPEEG.event) %loop through all events
    TMPEEG.event(ind).OCCL =   'none'; % initialize occlusion variable as none
    
    if ~isempty(str2double(TMPEEG.event(ind).type(2:end))) %type(2:end) --> number of trigger
        
        if find(occlusion == str2double(TMPEEG.event(ind).type(2:end)))
            TMPEEG.event(ind).OCCL = 'ON';
            
        elseif find(reappear == str2double(TMPEEG.event(ind).type(2:end))) 
            TMPEEG.event(ind).OCCL = 'OFF';
        end
        
    end
end

%% create output 

OUTEEG=TMPEEG;

return