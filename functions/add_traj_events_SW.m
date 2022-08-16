function OUTEEG = add_traj_events(INEEG)

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

last_rand2 = 1;


for ind = 1:length(TMPEEG.event) %loop through all events
    TMPEEG.event(ind).TRAJ =   'none'; % initialize occlusion variable as none

    if ~isempty(str2double(TMPEEG.event(ind).type(2:end))) %type(2:end) --> number of trigger

        if find(end_startvec == str2double(TMPEEG.event(ind).type(2:end)))
            %             TMPEEG.event(ind).TRAJ = 'RANDOM1';
            [TMPEEG.event(last_rand2:ind-1).TRAJ] = deal('RANDOM2'); % add the
            % label 'RANDOM1' to all the rows in the OCCL field from the last
            % time it was off up to the current index minus one

            last_rand1 = ind;

        elseif find(start_constant == str2double(TMPEEG.event(ind).type(2:end)))
            %             TMPEEG.event(ind).TRAJ = 'CONST';
            [TMPEEG.event(last_rand1:ind-1).TRAJ] = deal('RANDOM1'); % add the
            % label 'CONST' to all the rows in the OCCL field from the last

            last_const = ind;

        elseif find(end_constant == str2double(TMPEEG.event(ind).type(2:end)))
            [TMPEEG.event(last_const:ind-1).TRAJ] = deal('CONST'); % add the

            last_rand2 = ind;

        end

    end
end

%% create output

OUTEEG=TMPEEG;

return