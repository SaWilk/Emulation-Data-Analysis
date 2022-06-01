function [shift, start_ind, end_ind] = align_const_traj(const_traj, trial_traj)
%November 2021
%Adriana Böttcher
%------
%align constant trajectories of several trials
%compare constant trajectory and trial trajectory at every point
%calculate variance of difference
%should be minimal at constant trajectory

var_diff = [];
for i = 1:(numel(trial_traj) - numel(const_traj))
   trial_part = trial_traj(i:i+numel(const_traj)-1);
   
   diff = trial_part - const_traj;
%    figure;
%    hold on;
%    plot(const_traj, 'b');
%    plot(trial_part, 'r');
%    plot(diff, 'g');
%    
   var_diff(i) = var(diff);
end
shift = find(var_diff == min(var_diff));
start_ind = shift;
end_ind = shift + length(const_traj);

end