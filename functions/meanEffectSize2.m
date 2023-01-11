function effect = meanEffectSize2(x, varargin)
%MEANEFFECTSIZE One-sample or two-sample effect size computations
%   EFFECT = MEANEFFECTSIZE(X) computes the mean difference effect size for
%   a single sample X. X is a single or double vector. EFFECT is a table
%   with a row for each effect size computed. It has a column for the value
%   of the effect size, and a column for the confidence intervals for that
%   effect size, if they are computed.
%
%   EFFECT = MEANEFFECTSIZE(X, Y) computes the mean difference effect size
%   for two samples X and Y. X and Y must both be single or double vectors.
%
%   EFFECT = MEANEFFECTSIZE(...,'NAME1',val1,'NAME2',val2,...) specifies one 
%   or more of the following name/value arguments using any of the previous 
%   syntaxes:
%
%       'Alpha'                 - A value between 0 and 1 for a confidence
%                                 level of 100*(1-alpha)%.  The default is 
%                                 0.05 for 95% confidence. 
%
%       'BootstrapOptions'      - Options for computing bootstrap confidence 
%                                 intervals in parallel, specified as a structure.
%                                 Generate the BootstrapOptions structure with 
%                                 statset('bootci'). Computing bootstrap
%                                 confidence intervals in parallel requires 
%                                 Parallel Computing Toolbox.
%                                 MEANEFFECTSIZE uses the following fields:
%                                   'UseParallel'
%                                   'UseSubstreams'
%                                   'Streams'
%                                 For more information on these fields, see
%                                 PARALLELSTATS.
%
%       'ConfidenceIntervalType - A char or string of 'exact', 'bootstrap',
%                                 or 'none', indicating what kind of
%                                 confidence intervals to compute, if any.
%                                 The default is 'exact' when an exact
%                                 formula is known, and 'bootstrap' otherwise.
%
%       'Effect'                - A char, string array, or cell array of 
%                                 effect sizes to compute. The following are 
%                                 supported for single sample input:
%                                   'cohen'
%                                   'meandiff'
%                                   'robustcohen'
%                           
%                                For two sample input, the function
%                                supports the above effect sizes, along
%                                with the following additional effect
%                                sizes:
%                                   'cliff'
%                                   'glass'
%                                   'kstest'
%                                   'mediandiff'
%                                The default is 'meandiff'.
%                                NOTE: 
%                                   When given 'cohen', the function
%                                   computes the unbiased estimate of Cohen's
%                                   d, which is also sometimes called Hedge's
%                                   g.
%                                   A value of 'glass' is not supported for 
%                                   paired data.
%                                 
%       'Mean'                  - A scalar double or single representing a 
%                                 known mean. This is only used for one sample 
%                                 data. The default is 0.
%
%       'NumBootstraps'         - A positive integer indicating the number
%                                 of bootstrap replicas to use when computing 
%                                 the bootstrap confidence intervals. The
%                                 default is 1000.
%
%       'Paired'                - A logical scalar indicating if the two
%                                 samples are paired. The default is false.
%                                 NOTE:
%                                   If 'Paired' is true, then 'VarianceType'
%                                   must be 'equal'.
%                                   Setting 'Effect' to 'glass' is not
%                                   supported for paired data.
%
%       'VarianceType'          - A char or string of 'unequal' or 'equal',
%                                 indicating if the two sample data
%                                 are assumed to have equal or unequal 
%                                 variances. The default is 'equal'.
%
%   Example:  Compare one set of exam grades against a known mean of 70
%   
%   load examgrades
%   effect = meanEffectSize(grades(:,1), 'Mean', 70)
%
%   Example: Compare two sets of stock returns with Cohen's d to see how
%   different they are
%
%   load stockreturns
%   s1 = stocks(:,1);
%   s3 = stocks(:,3);
%   effect = meanEffectSize(s1, s3, 'Effect', 'cohen')
%   % This yields a Cohen's d of .15604. According to [3], we can
%   % conclude there is a small difference between stock 1 and stock 3
%
%   See also GARDNERALTMANPLOT, TTEST, TTEST2

% References:
% [1] Cousineau, Denis & Goulet-Pelletier, Jean-Christophe. (2021). "A study 
%     of confidence intervals for Cohen's d in within-subject designs with 
%     new proposals". The Quantitative Methods for Psychology.
%
% [2] Algina, J., Keselman, H. J., & Penfield, R. D. (2005). "An Alternative 
%     to Cohen's Standardized Mean Difference Effect Size: A Robust Parameter 
%     and Confidence Interval in the Two Independent Groups Case". Psychological 
%     Methods, 10(3), 317–328.
%
% [3] Hess, Melinda & Kromrey, Jeffrey. (2004). "Robust Confidence Intervals 
%     for Effect Sizes: A Comparative Study of Cohen's d and Cliff's Delta 
%     Under Non-normality and Heterogeneous Variances". Paper Presented at
%     the Annual Meeting of the American Educational Research Association. 
%
% [4] Delacre, Marie, et al. (2021). "Why Hedges  G*s Based on the Non-pooled 
%     Standard Deviation Should Be Reported with Welch's T-test."

% Copyright 2021-2022 The MathWorks, Inc.

% Parse data sample args
[varargin{:}] = convertStringsToChars(varargin{:});
if nargin >= 2 && ~ischar(varargin{1})
    % Given multiple inputs, second input is a sample of data
    y = varargin{1};
    varargin = varargin(2:end);
else
    y = [];
end

% Validate inputs
[x, y, nvPairs, effects] = validateAndParseInputs(x, y, varargin);
disp(effects)

% Get effect sizes
precType = internal.stats.dominantType(x,y);
effSizes = zeros(numel(effects), 1 + 2*~strcmpi(nvPairs.CIType, 'none'), 'like', precType);
ciType = nvPairs.CIType;
effectFullName = cell(1, numel(effects));
for i = 1:numel(effects)
    % Get effect size and handle to exact CI function, if it exists
    switch lower(effects{i})
        case 'cohen'
            [effSizes(i,1), df] = cohensD(x, y, nvPairs);
            effFcn = @(samp1, samp2) cohensD(samp1, samp2, nvPairs);
            exactCIFcn = @() cohensDCI(effSizes(i,1), x, y, nvPairs, df);
            defaultCI = 'exact';
            effectFullName{i} = 'CohensD';

        case 'cliff'
            [effSizes(i,1), diffs] = cliffsDelta(x, y, nvPairs);
            effFcn = @(samp1, samp2) cliffsDelta(samp1, samp2, nvPairs);
            exactCIFcn = @() cliffsDeltaCI(effSizes(i,1), diffs, x, y, nvPairs);
            defaultCI = 'exact';
            effectFullName{i} = 'CliffsDelta';

        case 'meandiff'
            [effSizes(i,1), df] = meanDiff(x, y, nvPairs);
            effFcn = @(samp1, samp2) meanDiff(samp1, samp2, nvPairs);
            exactCIFcn = @() meanDiffCI(effSizes(i,1), x, y, nvPairs, df);
            defaultCI = 'exact';
            effectFullName{i} = 'MeanDifference';

        case 'mediandiff'
            effSizes(i,1) = medianDiff(x, y);
            effFcn = @(samp1, samp2) medianDiff(samp1, samp2);
            defaultCI = 'bootstrap';
            effectFullName{i} = 'MedianDifference';

        case 'glass'
            [effSizes(i,1), df] = glasssDelta(x, y);
            effFcn = @(samp1, samp2) glasssDelta(samp1, samp2);
            exactCIFcn = @() glasssDeltaCI(effSizes(i,1), x, y, nvPairs, df);
            defaultCI = 'exact';
            effectFullName{i} = 'GlasssDelta';

        case 'kstest'
            effSizes(i,1) = kstestEff(x, y);
            effFcn = @(samp1, samp2) kstestEff(samp1, samp2);
            defaultCI = 'bootstrap';
            effectFullName{i} = 'KolmogorovSmirnovStatistic';

        case 'robustcohen'
            effSizes(i,1) = robustCohensD(x, y, nvPairs);
            effFcn = @(samp1, samp2) robustCohensD(samp1, samp2, nvPairs);
            defaultCI = 'bootstrap';
            effectFullName{i} = 'RobustCohensD';
    end
    if isempty(nvPairs.CIType)
       ciType = defaultCI;
    end
    
    isEffInfOrNan = (isnan(effSizes(i,1)) || isinf(effSizes(i,1)));
    switch lower(ciType)
        case 'exact'
            % This branch is only ever hit when functions have defined
            % exact CI functions due to how we validate arguments
            if isEffInfOrNan
                ci = [NaN, NaN];
            else
                ci = exactCIFcn();
            end
        case 'bootstrap'
            % General bootstrap procedure is the same for all effect sizes
            if isEffInfOrNan
                ci = [NaN, NaN];
            else
                ci = getBootstrapCI(x, y, nvPairs, effFcn);
            end
        case 'none'
            ci = [];
    end
    if any(isnan(ci))
        generateNaNCIWarning(effects{i}, nvPairs.EqualVar);
    end
    effSizes(i, 2:end) = ci;
end

% Construct data table
if size(effSizes, 2) == 1
    % No CI
    effect = table(effSizes, 'VariableNames', {'Effect'}, 'RowNames', effectFullName);
else
    % CI
    effect = table(effSizes(:,1), effSizes(:,2:end), 'VariableNames', ...
        {'Effect', 'ConfidenceIntervals'}, 'RowNames', effectFullName);
end
end

%%% EFFECT SIZE FUNCTIONS %%%
function [d, df] = cohensD(x, y, nvPairs)
% Compute unbiased Cohen's d, also sometimes called Hedge's g
nx = numel(x);
ny = numel(y);
meanX = mean(x);
if isempty(y)
    % One sample
    otherMean = nvPairs.Mean;
    stddev = std(x);
    df = nx - 1;
else
    % Two sample
    otherMean = mean(y);
    if ~nvPairs.EqualVar || nvPairs.Paired
        varX = var(x);
        varY = var(y);
        stddev = sqrt((varX + varY)/2);
        if nvPairs.Paired
            df = nx-1;
        else
            df = ((nx-1)*(ny-1)*(varX+varY)^2) / ((ny-1)*varX^2 + (nx-1)*varY^2);
        end
    else
        stddev = pooledSTD(x,y);
        df = nx + ny - 2;
    end
end
d = hedgesCorrection(df) * (meanX - otherMean)/stddev;
end

function ci = cohensDCI(d, x, y, nvPairs, df) 
% Compute exact CI for unbiased Cohen's d, also sometimes called Hedge's g
% The exact CIs are constructed using a noncentral t distribution
nx = numel(x);
ny = numel(y);
if isempty(y)
    lambda = sqrt(nx);
else
    if nvPairs.Paired
        % Use MAG method laid out in [1]
        % This method uses the quantiles of the non-central t distribution
        % to form the confidence intervals, unlike the other methods for CI
        % computation for Cohen's d, which find the CI via numerical
        % optimization
        varX = var(x);
        varY = var(y);
        W = sqrt(varX*varY) / ((varX+varY)/2);
        r = sqrt(nx/(2*(1-corr(x,y)*W)));
        lambda = d * hedgesCorrection(nx-1) * r;
        ci = nctinv(nvPairs.CIQuants, nx-1, lambda) / r;
        return
    else
        if nvPairs.EqualVar
            lambda = sqrt(1 / ((1/nx) + (1/ny)));
        else
            varX = var(x);
            varY = var(y);
            lambda = sqrt((nx*ny*(varX + varY)) / (2*(ny*varX + nx*varY)));
        end
    end
end
% When computing the confidence intervals, remove the correction term
% first, then add it back by multiplying it with the bounds
hedgesCorrect = hedgesCorrection(df);
ci = hedgesCorrect * getNCTConfidenceBounds(d*lambda/hedgesCorrect, df, nvPairs.CIQuants) / lambda;
end

function [cliff, diffs] = cliffsDelta(x, y, nvPairs)
% This function computes Cliff's Delta. This effect size is only valid for
% two sample data
diffs = sign(x - y');
nx = numel(x);
ny = numel(y);
if nvPairs.Paired
    % In the paired setup, we can use Cliff's delta to compare with group
    % or between group differences. Here, we opt for between-group differences,
    % since that is more common in literature.
    % In that case, set the difference within groups (on the diagonal) to 0
    % and remove those comparisons from the count of pairs to compare
    pairedInds = 1:nx+1:(nx^2);
    diffs(pairedInds) = 0;
    numPairs = nx * (nx - 1);
else
    numPairs = nx*ny;
end
cliff = sum(diffs, 'all') / numPairs;
end

function ci = cliffsDeltaCI(delta, diffs, x, y, nvPairs)
% This function computes CI for Cliff's Delta.
nx = numel(x);
ny = numel(y);
zcrit = norminv(nvPairs.CIQuants(2));

if nvPairs.Paired
    % Get the marginal terms for the rows and columns, and get the
    % differences relative to the computed delta
    rowMarginals = sum(diffs, 2) / (ny-1);
    colMarginals = sum(diffs)' / (nx-1);
    rowMarginalDiffs = rowMarginals - delta;
    colMarginalDiffs = colMarginals - delta;

    % Get the difference between the computed delta and all elements of the
    % samples, excluding the between-group differences
    allDifs = diffs - delta;
    allDifs(1:nx+1:end) = 0;

    % Construct the confidence intervals
    numer = (nx-1)^2 * (sum(rowMarginalDiffs.^2) + sum(colMarginalDiffs.^2) + ...
        2*rowMarginalDiffs'*colMarginalDiffs) - sum(allDifs.^2, 'all') - sum(allDifs .* allDifs', 'all');
    denom = nx*(nx-1)*(nx-2)*(nx-3);
    varDelta = numer/denom;
    bound = zcrit * sqrt(varDelta);
    ci = [delta - bound, delta + bound];
else
    if all(min(x) > y) || all(min(y) > x)
        % One sample completely dominates the other
        ci = [delta, delta];
    else
        % Get variance estimate for column-wise, row-wise, and element-wise
        % deltas
        colMarginals = sum(diffs)/nx;
        varCols = sum((colMarginals - delta).^2) / (ny-1);
        rowMarginals = sum(diffs, 2)/ny;
        varRows = sum((rowMarginals - delta).^2) / (nx-1);
        varAll = sum((diffs - delta).^2, 'all') / ((nx-1) * (ny-1));

        % Combine above variance estimates to get the estimated variance on
        % delta, use that to construct bounds
        varDelta = ((ny-1)*varRows + (nx-1)*varCols + varAll) / (nx*ny);
        delta2 = delta^2;
        zcrit2 = zcrit^2;
        numer = zcrit * sqrt(varDelta) * sqrt((1-delta2)^2 + (zcrit2 * varDelta));
        denom = 1 - delta2 + (zcrit2 * varDelta);
        deltaDiff = delta - delta^3;
        ci = [deltaDiff - numer, deltaDiff + numer] / denom;
    end
end
% delta is contained within [-1, 1]. Truncate anything outside that range
% to the endpoints
ci(ci < -1) = -1;
ci(ci > 1) = 1;
end

function [md, df] = meanDiff(x, y, nvPairs)
% Computes the mean difference
if isempty(y)
    % One sample
    md = mean(x) - nvPairs.Mean;
    df = numel(x) - 1;
else
    md = mean(x) - mean(y);
    if nvPairs.EqualVar
        df = (numel(x) - 1) + ~nvPairs.Paired*(numel(y) - 1);
    else
        df = welchSatterthwaiteDF(x, y);
    end
end
end

function ci = meanDiffCI(md, x, y, nvPairs, df)
% Computes the exact CI for the mean difference effect size
if isempty(y)
    % One sample
    stddev = (std(x) / sqrt(df+1));
elseif nvPairs.Paired
    stddev = std(x-y) / sqrt(df+1);
else
    nx = numel(x);
    ny = numel(y);
    if ~nvPairs.EqualVar
        stddev = sqrt(var(x)/nx + var(y)/ny);
    else
        stddev = pooledSTD(x, y) * sqrt(1/nx + 1/ny);
    end
end
crit = tinv(nvPairs.CIQuants(2), df) * stddev;
ci = [md - crit, md + crit];
end

function md = medianDiff(x, y)
% Computes the median difference. Only valid for two sample data
md = median(x) - median(y);
end
 
function [gd, df] = glasssDelta(x, y)
% This function computes Glass's Delta. Only valid for unpaired two sample 
% data
df = numel(x) - 1;
gd = hedgesCorrection(df) * (mean(x) - mean(y)) / std(x);
end

function ci = glasssDeltaCI(gd, x, y, nvPairs, df)
% This function computes CI for Glass's Delta. Only valid for unpaired two 
% sample data.
% Equal or unequal variance doesn't impact Glass's delta, since the std of
% only one population is used
nx = numel(x);
ny = numel(y);
lambda = sqrt((1/nx) + var(y)/(ny*var(x)));
hedgesCorrect = hedgesCorrection(df);
ci = hedgesCorrect * getNCTConfidenceBounds(gd/(lambda*hedgesCorrect), df, nvPairs.CIQuants) * lambda;
end

function ks = kstestEff(x, y)
% Computes Kolmogorov-Smirnov test statistic. Only valid for two sample data
[~, ~, ks] = kstest2(x, y);
end

function d = robustCohensD(x, y, nvPairs)
% Computes the robust version of unbiased Cohen's d
% It uses 20% trimmed mean and winsorised variance
pct = 20;
pctl = [pct, 100 - pct];

trimMeanX = trimmean(x, pct);
twentyPctX = prctile(x, pctl);
x(x < twentyPctX(1)) = twentyPctX(1);
x(x > twentyPctX(2)) = twentyPctX(2);
nx = numel(x);
ny = numel(y);
if isempty(y)
    df = nx - 1;
    otherMean = nvPairs.Mean;
    stddev = std(x);
else
    otherMean = trimmean(y, pct);
    twentyPctY = prctile(y, pctl);
    y(y < twentyPctY(1)) = twentyPctY(1);
    y(y > twentyPctY(2)) = twentyPctY(2);

    if ~nvPairs.EqualVar || nvPairs.Paired
        varX = var(x);
        varY = var(y);
        stddev = sqrt((varX + varY)/2);
        if nvPairs.Paired
            df = nx-1;
        else
            df = ((nx-1)*(ny-1)*(varX+varY)^2) / ((ny-1)*varX^2 + (nx-1)*varY^2);
        end
    else
        stddev = pooledSTD(x,y);
        df = nx + ny - 2;
    end
end

% Algina & Keselman [2] suggest multiplying the computed result by .642
% to ensure that if the data are normal and without outliers, then 
% robustCohensD = CohensD
d = .642 * hedgesCorrection(df) * ((trimMeanX - otherMean) / stddev);
end

%%% EFFECT SIZE HELPERS %%%
function J = hedgesCorrection(df)
% Compute the hedge's correction term to make Cohen's d or Glass's delta
% unbiased
if df == 0
    J = NaN;
else
    J = exp(gammaln(df/2) - log(sqrt(df/2)) - gammaln((df-1)/2));
end
end

function ps = pooledSTD(x, y)
% Compute the pooled standard deviation
nx = numel(x);
ny = numel(y);
ps = sqrt(((nx-1)*var(x) + (ny-1)*var(y)) / (nx+ny-2));
end

function df = welchSatterthwaiteDF(x, y)
% Compute the Welch-Satterthwaite estimate of the degrees of freedom
nx = numel(x);
ny = numel(y);
varX = var(x)/nx;
varY = var(y)/ny;
se = sqrt(varX + varY);
df = se^4 / ((varX^2 / (nx-1)) + (varY^2 / (ny-1)));
end

function ci = getBootstrapCI(x, y, nvPairs, effFcn)
% Get bootstrap CI for the given effect
if isempty(y)
    % One sample
    % Signature of all effect functions take in two samples -  add a
    % wrapper around the function to pass in an empty second sample
    oneSampEffFcn = @(x) effFcn(x, []);
    ci = bootci(nvPairs.NumBoot, {oneSampEffFcn, x}, 'Alpha', nvPairs.Alpha, ...
        'Options', nvPairs.BootOpts)';
elseif nvPairs.Paired
    % Paired samples
     ci = bootci(nvPairs.NumBoot, {effFcn, x, y}, 'Alpha', nvPairs.Alpha, ...
        'Options', nvPairs.BootOpts)';
else
    % Two samples, which may be of different lengths
    % bootci only can work on data of the same size, so instead, we stack
    % the x and y together, and give bootstrap a vector of indices to
    % sample from. The function passed to bootbca will translate from the
    % sampled indices to the corresponding values of x and y
    idx = 1:(numel(x)+numel(y));
    bootfun = @(idx) unpairedSamplesBootstrapFcn(idx, x, y, effFcn);
    ci = bootbca(idx, nvPairs.NumBoot, bootfun, nvPairs.Alpha, nvPairs.BootOpts);
end
end

function effect = unpairedSamplesBootstrapFcn(indices, x, y, effectFcn)
% This function is used with bootci for computing confidence intervals when
% the x and y samples are unequal in size
xIdx = (indices <= numel(x));
yIdx = ~xIdx;
currX = x(indices(xIdx));
currY = y(indices(yIdx) - numel(x));
if isempty(currX) || isempty(currY)
    % Sample generated samples that have no data from X or Y
    % We can't reason about what the effect should be, so return NaN
    effect = NaN;
else
    effect = effectFcn(currX, currY);
end
end

function [ci,bstat] = bootbca(dataInds, nboot, bootfun, alpha, bootstrpOptions)
% corrected and accelerated percentile bootstrap CI
effect = bootfun(dataInds);
bstat = bootstrp(nboot, bootfun, dataInds, 'Options', bootstrpOptions);
% Remove any NaN values, which indicate a sample where all x or all y were
% missing
bstat(isnan(bstat)) = [];

% same as bootcper, this is the bias correction
z0 = norminv(mean(bstat < effect,1) + mean(bstat == effect,1)/2);

% apply jackknife
jstat = jackknife(bootfun, dataInds, 'Options', bootstrpOptions);
% Remove NaN for the same reason as above
jstat(isnan(jstat)) = [];
N = size(jstat,1);
weights = repmat(1/N,N,1);

% acceleration finding
mjstat = sum(jstat.*weights,1); % mean along 1st dim.
score = mjstat - jstat; % score function at stat; ignore (N-1) factor because it cancels out in the skew
iszer = all(score==0,1);
skew = sum((score.^3).*weights,1) ./ ...
    (sum((score.^2).*weights,1).^1.5) /sqrt(N); % skewness of the score function
skew(iszer) = 0;
acc = skew/6;  % acceleration

% transform back with bias corrected and acceleration
z_alpha1 = norminv(alpha/2);
z_alpha2 = -z_alpha1;
pct1 = 100*normcdf(z0 +(z0+z_alpha1)./(1-acc.*(z0+z_alpha1)));
pct1(z0==Inf) = 100;
pct1(z0==-Inf) = 0;
pct2 = 100*normcdf(z0 +(z0+z_alpha2)./(1-acc.*(z0+z_alpha2)));
pct2(z0==Inf) = 100;
pct2(z0==-Inf) = 0;

% inverse of ECDF
m = numel(effect);
lower = zeros(1,m);
upper = zeros(1,m);
for i=1:m
    lower(i) = prctile(bstat(:,i),pct2(i),1);
    upper(i) = prctile(bstat(:,i),pct1(i),1);
end
ci = sort([lower;upper],1);
end

function ci = getNCTConfidenceBounds(nctval, df, quants)
% This function finds the confidence intervals for the Cohen's d family of
% effect sizes. Doing do requires numeric optimization to find two 
% non-centrality parameters nct_l and nct_u. With a default Alpha
% (95% confidence intervals), nct_l and nct_u are found such that the input
% nctval has 2.5% and 97.5% cumulative probability, respectively.
ci = [NaN, NaN];

% nctval is theoretically unbounded, but practically we can construct
% bounds where the 0s likely live. Construct a set of bounds that puts the
% nctval in the middle, and sets the bounds to be on a higher order than
% the nctval itself. This should generate a valid bound in most cases
% The bounds are constructed with the larger value as the lower bound due
% to the CDF - the nct param roughly controls where the center of the distribution
% is. bound >> nctval means the center of the distribution is far away from
% nctval, and nctval likely lives in the lower tail. The CDF returns ~0 in 
% this case. bound << nctval means the opposite, so nctval is likely in the 
% upper tail, and the cdf will return ~1
bounds = [nctval + abs(10*nctval), nctval - abs(10*nctval)];
boundOptimFcn = @(nct, q) nctcdf(nctval, df, nct) - q;
if (boundOptimFcn(bounds(1), quants(2)) < 0) && (boundOptimFcn(bounds(2), quants(2)) > 0)
    % The CDF is monotonic, so the one and only zero lives in these bounds
    [ci(1), ~, errL] = fzero(@(nct_l) boundOptimFcn(nct_l, quants(2)), bounds);
else
    % 0 is outside the bounds - use nctval as the start estimate instead of
    % the bounds
    [ci(1), ~, errL] = fzero(@(nct_l) boundOptimFcn(nct_l, quants(2)), nctval);
end

% Do the same as above but for the upper CI bound
if (boundOptimFcn(bounds(1), quants(1)) < 0) && (boundOptimFcn(bounds(2), quants(1)) > 0)
    [ci(2), ~, errU] = fzero(@(nct_u) boundOptimFcn(nct_u, quants(1)), bounds);
else
    [ci(2), ~, errU] = fzero(@(nct_u) boundOptimFcn(nct_u, quants(1)), nctval);
end

if any([errL, errU] < 0)
    error(message('stats:effectsize:NoSolution'));
end
end

function generateNaNCIWarning(effName, eqVar)
% Generate a warning about NaN ci/effect size
switch lower(effName)
    case 'cliff'
        id = 'stats:effectsize:TooFewSamples';
    case 'meandiff'
        % meandiff can generate NaNs if there are too few samples
        % or if the variance is zero, based on if the variances of the
        % samples are equal or not
        if eqVar
            id = 'stats:effectsize:TooFewSamples';
        else
            id = 'stats:effectsize:ZeroVarianceSamples';
        end
    otherwise
        id = 'stats:effectsize:ZeroVarianceSamples';
end
warning(message(id, effName));
end

%%% VALIDATION FUNCTIONS %%%
function [x, y, nvPairs, effects] = validateAndParseInputs(x, y, nvArgs)
% Parses the NV pairs, and validates all given inputs

% Parse NV pairs
names = {'Mean', 'Effect', 'VarianceType', 'ConfidenceIntervalType', 'NumBootstraps', 'BootstrapOptions', 'Paired', 'Alpha'};
vals  = {  0, {'meandiff'}, 'equal',            [],                  1000,       statset('bootci'),   false,    .05};
[knownMean, effects, varType, ciType, numBoot, bootOpts, paired, alpha] = ...
    internal.stats.parseArgs(names, vals, nvArgs{:});

% Validate x
if istall(x)
    error(message('stats:effectsize:TallNotSupported'));
end
validateattributes(x, {'double', 'single'}, {'real', 'vector', 'nonempty', 'nonsparse'}, mfilename, 'x')

% Validate y, if given
if ~isempty(y)
    if istall(y)
        error(message('stats:effectsize:TallNotSupported'));
    end
    validateattributes(y, {'double', 'single'}, {'real', 'vector', 'nonsparse'}, mfilename, 'y')
end

% Validate basic attributes of remaining NV pairs
validateattributes(alpha, {'double', 'single'}, {'real', 'scalar', 'nonempty', 'nonsparse', 'nonnan', '>' 0, '<', 1},...
    mfilename, 'Alpha');
validateattributes(knownMean, {'double', 'single'}, {'real', 'scalar', 'nonempty', 'nonsparse', 'nonnan', 'finite'}, ...
    mfilename, 'Mean');
validateattributes(numBoot, {'numeric'}, {'real', 'scalar', 'nonempty', 'integer', 'nonnan', 'finite', 'positive'}, ...
    mfilename, 'NumBootstraps');
if ~isfloat(numBoot)
    % Given an int type
    numBoot = double(numBoot);
end
validateattributes(bootOpts, {'struct'}, {}, mfilename, 'BootstrapOptions');
paired = internal.stats.parseOnOff(paired, 'Paired');

% Validate varType and ciType, which can only be fixed values
% Convert varType to a logical, indicating equal variance or not
varType = validatestring(varType, {'equal', 'unequal'}, mfilename, 'VarianceType');
varType = strcmpi(varType, 'equal');
if ~isempty(ciType)
    ciType = validatestring(ciType, {'exact', 'bootstrap', 'none'}, mfilename, 'ConfidenceIntervalType');
end

% Validate effects
effectList = {'cohen', 'cliff', 'meandiff', 'mediandiff', ...
    'glass', 'kstest', 'robustcohen'};
validateattributes(effects, {'cell', 'char'}, {'vector'}, mfilename, 'Effect');
if ischar(effects)
    effects = cellstr(effects);
end
for i = 1:numel(effects)
    effects{i} = validatestring(effects{i}, effectList, mfilename, 'Effect');
end
effects = unique(effects, 'stable');

% Check if one sample data is given, but a two sample effect was requested
if isempty(y)
    twoSampEffects = {'cliff', 'mediandiff', 'glass', 'kstest'};
    badEff = find(ismember(effects, twoSampEffects), 1);
    if ~isempty(badEff)
        error(message('stats:effectsize:NotSupportedForOneSample', effects{badEff}));
    end
end

% Not all effect sizes have exact method for CI computation. If 'exact' is 
% requested for an effect that doesn't have an exact method, error
if strcmpi(ciType, 'exact')
    noExactMethod = {'mediandiff', 'kstest', 'robustcohen'};
    noExactSpecified = find(ismember(effects, noExactMethod), 1);
    if ~isempty(noExactSpecified)
        error(message('stats:effectsize:NoExactMethodKnown', effects{noExactSpecified}))
    end
end

% If the samples are paired, we are not able to also support computing effects
% when the variance is unequal
if paired && ~varType
    error(message('stats:effectsize:PairedAndUnequalVarUnsupported'))
end

% Computing Glass's delta with paired data is not supported
if paired && any(strcmpi(effects, 'glass'))
    error(message('stats:effectsize:PairedGlassDeltaNotSupported'))
end

% If samples are paired, then x and y need to be the same length
x = x(:);
y = y(:);
if paired && (numel(x) ~= numel(y))
    error(message('stats:effectsize:PairedDataDiffSizes'))
end

% Remove NaNs. For paired samples, remove the row if either x or y is NaN
yWasPopulated = ~isempty(y);
if paired
    nanVals = (isnan(x) | isnan(y));
    x(nanVals) = [];
    y(nanVals) = [];
else
    x(isnan(x)) = [];
    y(isnan(y)) = [];
end

% After removing all NaN values, if there is no data left, error
if isempty(x) || (yWasPopulated && isempty(y))
    error(message('stats:effectsize:AllNaNData'));
end

quants = [alpha/2, 1 - (alpha/2)];
nvPairs = struct('Alpha', alpha, 'CIQuants', quants, 'Mean', knownMean, 'EqualVar', varType,...
    'CIType', ciType, 'NumBoot', numBoot, 'BootOpts', bootOpts, 'Paired', paired);
end