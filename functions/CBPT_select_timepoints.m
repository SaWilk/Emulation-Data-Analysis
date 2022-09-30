function CBPT_reduced = CBPT_select_timepoints(CBPT, timepoints)
% this function reduces the CBPT output (not averaged over time) to
% specific time points of interest.
% time points are specified as indeces in a vector (timepoints)

% copy original output
CBPT_reduced = CBPT;
% reduce to chosen timepoints
CBPT_reduced.time = CBPT_reduced.time(:, timepoints);
if isfield(CBPT, 'prob')
    CBPT_reduced.prob = (CBPT_reduced.prob(:,:,timepoints));
end
if isfield(CBPT, 'posclusterslabelmat')
    CBPT_reduced.posclusterslabelmat = (CBPT_reduced.posclusterslabelmat(:,:,timepoints));
end
if isfield(CBPT, 'negclusterslabelmat')
    CBPT_reduced.negclusterslabelmat = (CBPT_reduced.negclusterslabelmat(:,:,timepoints));
end
if isfield(CBPT, 'cirange')
    CBPT_reduced.cirange = (CBPT_reduced.cirange(:,:,timepoints));
end
if isfield(CBPT, 'mask')
    CBPT_reduced.mask = (CBPT_reduced.mask(:,:,timepoints));
end
if isfield(CBPT, 'stat')
    CBPT_reduced.stat = (CBPT_reduced.stat(:,:,timepoints));
end
if isfield(CBPT, 'ref')
    CBPT_reduced.ref = (CBPT_reduced.ref(:,:,timepoints));
end

end