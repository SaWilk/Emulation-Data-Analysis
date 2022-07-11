function [shift, start_ind, end_ind] = align_to_const_traj(const_traj, trial_traj)
%08/07/22
%Adriana Böttcher
%------
% This function aligns a trial trajectory to the constant trajectory and
% returns the shift (how many points does the constant trajectory has to be
% shifted until it fits best), the start-index of the constant trajectory
% in the given trial trajectory, and the end-index.
%
% The function shifts the constant trajectory along the trial trajectory
% and calculates the variance of the difference (trial_traj - const_traj) 
% of every point and identifies the shift where the variance of the
% difference is at its minimum.
%
% The function's input:
% const_traj --> the constant trajectory as matrix [x y] read from file
% const_traj_fixed (should be loaded into the workspace before)
% trial_traj --> the matrix of the trial_traj [x y] which the constant traj
% should be aligned to 
%
% important: the constant trajectory and trial trajectory must have the
% same frame rate, if one is interpolated, the other has to be too, since
% the alignment is based on shifting the const traj along the trial traj
% points.
%------

% create empty vector to be filled with the variance of differences at
% every point of the trial traj
var_diff = [];

% loop through the trial traj
for i = 1:(numel(trial_traj) - numel(const_traj))

    %extract trial part as long as the constant traj
    try
        trial_part = trial_traj(i:i+numel(const_traj)-1);
    catch
        disp("trial part could not be extracted");
    end

    %calculate the differences of all points 
    try
        diff = trial_part - const_traj;
    catch
        disp("differences could not be calculated");

    % calculate variance of differences
    try    
        var_diff(i) = var(diff);
    catch
        disp("error calculating variance of differences");
    end
end

% extract shift, start ind and end ind
try
    shift = find(var_diff == min(var_diff));
    start_ind = shift;
    end_ind = shift + length(const_traj);
catch
    disp("error creating output");
end
end