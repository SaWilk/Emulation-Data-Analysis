function OUTEEG = addcodingevents(INEEG)

TMPEEG=INEEG;

instruction     = 10;
exp_start       = 11;
fix             = 12;
trial_start_L   = 13;
trial_start_R   = 14;
trial_end       = 15;  %superfluous with fix trigger (except for break)
pause_start     = 16;
pause_end       = 17;
exp_end         = 18;
button          = 19;
occlusion       = 20;
reappear        = 21;
fourth_trial    = 22; % TODO what is that???
start_constant  = 23;
end_constant    = 24;
fourth_trial    = 25; % TODO what is that???
start_startvec  = 26;
end_startvec    = 27;
C_too_early     = 28;
C_just_right    = 29;
C_too_late      = 30;


for ind = 1: length(TMPEEG.event)
    TMPEEG.event(ind).TRAJ    =   'none'; 
    if ~isempty(str2double(TMPEEG.event(ind).type(2:end)))
        if find(end_startvec == str2double(TMPEEG.event(ind).type(2:end)))
            TMPEEG.event(ind).TRAJ = 'RANDOM1';
        elseif find(start_constant == str2double(TMPEEG.event(ind).type(2:end))) 
            TMPEEG.event(ind).TRAJ = 'CONSTANT';
        elseif find(end_constant == str2double(TMPEEG.event(ind).type(2:end))) 
             TMPEEG.event(ind).TRAJ = 'RANDOM2';
        end
    end
end
OUTEEG=TMPEEG;
return