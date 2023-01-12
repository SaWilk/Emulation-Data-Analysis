function OUTEEG = add_occl_events(INEEG)

% Adapted by Saskia Wilken

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

last_off = 1;

for ind = 1:length(TMPEEG.event) %loop through all events
    TMPEEG.event(ind).OCCL =   'none'; % initialize occlusion variable as none
    
    if ~isempty(str2double(TMPEEG.event(ind).type(2:end))) %type(2:end) --> number of trigger

        if find(occlusion == str2double(TMPEEG.event(ind).type(2:end)))
            [TMPEEG.event(last_off:ind-1).OCCL] = deal('OFF'); % add the 
            % label 'OFF' to all the rows in the OCCL field from the last
            % time it was off up to the current index minus one

            last_on = ind; % save the current index as the position in which 
            % occlusion was last switched on

        elseif find(reappear == str2double(TMPEEG.event(ind).type(2:end))) 
            [TMPEEG.event(last_on:ind-1).OCCL] = deal('ON'); % same vice versa

            last_off = ind;              

        end
        
    end
end
% tiny issue: the very last stretch of triggers after the last occlusion is
% not assigned any value. 

%% create output 

OUTEEG=TMPEEG;

return