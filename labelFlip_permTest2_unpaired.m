function [p, obsStat, permStats, info] = labelFlip_permTest2_unpaired(x, y, varargin)
%PERMTEST2_UNPAIRED Two-sample (independent) permutation test between vectors.
%
%   [p, obsStat, permStats, info] = labelFlip_permTest2_unpaired(x, y, 'Name', Value, ...)
%
% OVERVIEW
%   This function implements a two-sample permutation (randomization) test for
%   INDEPENDENT observations. It tests whether two groups x and y differ with
%   respect to a user-chosen statistic (default: mean(x) - mean(y)).
%
%   Null hypothesis (H0):
%       The group labels are exchangeable; equivalently, x and y are drawn from
%       the same distribution (for independent samples).
%
%   Under H0, permuting group labels across all observations is valid.
%
% WHEN THIS IS APPROPRIATE
%   Use this when each data point is an independent experimental unit:
%     - one measurement per animal in each group
%     - independent samples from different subjects
%
% INPUTS
%   x, y : numeric vectors. NaNs are removed.
%
% NAME-VALUE OPTIONS
%   'nPerm'     : number of permutations (default: 10000)
%   'statFun'   : function handle: statFun(x, y) -> scalar
%                default: @(a,b) mean(a) - mean(b)
%   'tail'      : 'two-sided' (default), 'right', or 'left'
%                - two-sided: extreme in magnitude
%                - right    : unusually large (positive direction)
%                - left     : unusually small (negative direction)
%   'seed'      : [] (default) or integer seed for reproducibility
%   'returnNull': true (default) or false; whether to return permStats
%
% OUTPUTS
%   p         : permutation p-value (Monte Carlo estimate with +1 correction)
%   obsStat   : observed statistic using original labels
%   permStats : vector of permuted statistics (empty if returnNull=false)
%   info      : struct with details (nUsed, nPerm, tail, seedUsed, method, etc.)
%
% P-VALUE COMPUTATION (Monte Carlo +1 correction)
%   Let b be the number of permuted statistics at least as extreme as obsStat.
%   Then p = (b + 1) / (nPerm + 1). This avoids returning p = 0 when obsStat is
%   more extreme than all sampled permutations.
%
% EXAMPLES
%   % Default: difference in means (two-sided)
%   [p, obs, null] = labelFlip_permTest2_unpaired(x, y, 'nPerm', 20000);
%
%   % Difference in medians
%   statFun = @(a,b) median(a) - median(b);
%   [p, obs, null] = labelFlip_permTest2_unpaired(x, y, 'statFun', statFun);
%
%   % One-sided (x > y)
%   [p, obs] = labelFlip_permTest2_unpaired(x, y, 'tail', 'right');
%
% SEE ALSO
%   permtest2_paired

% ---------- Parse inputs ----------
validateattributes(x, {'numeric'}, {'vector'}, mfilename, 'x', 1);
validateattributes(y, {'numeric'}, {'vector'}, mfilename, 'y', 2);

x = x(:);
y = y(:);
xx=x(:);
yy=y(:);

% Remove NaNs
x = x(~isnan(x));
y = y(~isnan(y));

nx = numel(x);
ny = numel(y);

if nx < 1 || ny < 1
    error('labelFlip_permTest2_unpaired:EmptyAfterNaNRemoval', ...
        'After removing NaNs, x and y must each contain at least 1 value.');
end

pIn = inputParser;
pIn.FunctionName = mfilename;

addParameter(pIn, 'nPerm', 10000, @(v) validateattributes(v, {'numeric'}, {'scalar','integer','>=',1}));
addParameter(pIn, 'statFun', [], @(f) isempty(f) || isa(f, 'function_handle'));
addParameter(pIn, 'tail', 'two-sided', @(s) ischar(s) || isstring(s));
addParameter(pIn, 'seed', [], @(v) isempty(v) || (isscalar(v) && isnumeric(v) && isfinite(v)));
addParameter(pIn, 'returnNull', true, @(v) islogical(v) && isscalar(v));

parse(pIn, varargin{:});
opt = pIn.Results;

if isempty(opt.statFun)
    opt.statFun = @(a,b) mean(a) - mean(b);
end

tail = lower(string(opt.tail));
if ~any(tail == ["two-sided","right","left"])
    error('labelFlip_permTest2_unpaired:BadTail', 'tail must be ''two-sided'', ''right'', or ''left''.');
end

% Seed handling (reproducibility)
seedUsed = [];
if ~isempty(opt.seed)
    seedUsed = opt.seed;
    rng(opt.seed, 'twister');
end

% ---------- Observed statistic ----------
obsStat = opt.statFun(x, y);
if ~isscalar(obsStat) || ~isfinite(obsStat)
    error('labelFlip_permTest2_unpaired:BadStatFun', ...
        'statFun must return a finite scalar. Got: %s', mat2str(obsStat));
end

% ---------- Permutations (Monte Carlo) ----------
N = nx + ny;
pooled = [x; y];

if opt.returnNull
    permStats = zeros(opt.nPerm, 1);
else
    permStats = [];
end

b = 0; % count of "as or more extreme" permutations

for i = 1:opt.nPerm
    idx = randperm(N);
    xPerm = pooled(idx(1:nx));
    yPerm = pooled(idx(nx+1:end));
    s = opt.statFun(xPerm, yPerm);

    if ~isscalar(s) || ~isfinite(s)
        error('labelFlip_permTest2_unpaired:BadStatDuringPerm', ...
            'statFun returned a non-finite scalar during permutations.');
    end

    if opt.returnNull
        permStats(i) = s;
    end

    switch tail
        case "two-sided"
            if abs(s) >= abs(obsStat)
                b = b + 1;
            end
        case "right"
            if s >= obsStat
                b = b + 1;
            end
        case "left"
            if s <= obsStat
                b = b + 1;
            end
    end
end

% Monte Carlo p-value with +1 correction
p = (b + 1) / (opt.nPerm + 1);

% ---------- Info ----------
info = struct();
info.method      = 'unpaired_label_permutation';
info.nPerm       = opt.nPerm;
info.tail        = char(tail);
info.nX          = nx;
info.nY          = ny;
info.seedUsed    = seedUsed;
info.statFun     = func2str(opt.statFun);
info.nanRemovedX = sum(isnan(xx));
info.nanRemovedY = sum(isnan(yy));
info.countExtreme = b;
info.pValue       = p;

end
