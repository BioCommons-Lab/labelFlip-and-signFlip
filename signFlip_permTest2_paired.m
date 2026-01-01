function [p, obsStat, permStats, info] = signFlip_permTest2_paired(x, y, varargin)
%PERMTEST2_PAIRED Paired permutation test (sign-flip) between matched vectors.
%
%   [p, obsStat, permStats, info] = signFlip_permTest2_paired(x, y, 'Name', Value, ...)
%
% OVERVIEW
%   This function implements a paired permutation test for repeated-measures /
%   matched-pairs data (e.g., before vs after within the same subject).
%
%   It operates on within-pair differences d_i = x_i - y_i and generates the null
%   distribution by randomly flipping the sign of each difference. This is the
%   correct permutation analogue of a paired t-test under minimal assumptions.
%
%   Null hypothesis (H0):
%       Under the null hypothesis, swapping x and y within any pair does not change the joint distribution. 
%       Equivalently, the distribution of within-pair differences is symmetric around zero.
%       (or, more generally, labels within each pair are exchangeable under H0).
%
%   Permutation mechanism:
%       For each pair i, swap labels within that pair (equivalently multiply d_i
%       by +1 or -1). Pairing is preserved; only within-pair assignment changes.
%
% WHEN THIS IS APPROPRIATE
%   Use this when data are paired:
%     - same animal measured in two conditions
%     - same neuron pre vs post manipulation
%     - matched pairs by design
%
% INPUTS
%   x, y : numeric vectors of the same length. NaN pairs are removed (if either
%          member of a pair is NaN, that pair is dropped).
%
% NAME-VALUE OPTIONS
%   'nPerm'     : number of permutations (default: 10000)
%   'statFun'   : function handle: statFun(d) -> scalar, where d = x - y
%                default: @(d) mean(d)
%   'tail'      : 'two-sided' (default), 'right', or 'left'
%                Here obsStat is computed on differences d:
%                - 'right' tests whether x > y on average (obsStat unusually large)
%                - 'left'  tests whether x < y on average
%   'seed'      : [] (default) or integer seed for reproducibility
%   'returnNull': true (default) or false
%
% OUTPUTS
%   p         : permutation p-value (Monte Carlo estimate with +1 correction)
%   obsStat   : observed statistic computed from differences (default mean(d))
%   permStats : vector of permuted statistics (empty if returnNull=false)
%   info      : struct with details (nPairsUsed, nPerm, tail, seedUsed, etc.)
%
% EXAMPLES
%   % Paired before/after
%   [p, obs, null] = signFlip_permTest2_paired(before, after, 'nPerm', 20000);
%
%   % Median of differences (robust)
%   statFun = @(d) median(d);
%   [p, obs] = signFlip_permTest2_paired(before, after, 'statFun', statFun);
%
%   % One-sided: expect increase after treatment (x > y)
%   [p, obs] = signFlip_permTest2_paired(after, before, 'tail', 'right');
%
% SEE ALSO
%   permtest2_unpaired

% ---------- Parse inputs ----------
validateattributes(x, {'numeric'}, {'vector'}, mfilename, 'x', 1);
validateattributes(y, {'numeric'}, {'vector'}, mfilename, 'y', 2);

x = x(:);
y = y(:);

if numel(x) ~= numel(y)
    error('signFlip_permTest2_paired:LengthMismatch', ...
        'For paired tests, x and y must have the same number of elements.');
end

% Remove NaN pairs (drop any pair where either is NaN)
mask = ~(isnan(x) | isnan(y));
nPairsRemoved = sum(~mask);
x = x(mask);
y = y(mask);

nPairs = numel(x);
if nPairs < 1
    error('signFlip_permTest2_paired:EmptyAfterNaNRemoval', ...
        'After removing NaN pairs, at least 1 valid pair is required.');
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
    opt.statFun = @(d) mean(d);
end

tail = lower(string(opt.tail));
if ~any(tail == ["two-sided","right","left"])
    error('signFlip_permTest2_paired:BadTail', 'tail must be ''two-sided'', ''right'', or ''left''.');
end

% Seed handling
seedUsed = [];
if ~isempty(opt.seed)
    seedUsed = opt.seed;
    rng(opt.seed, 'twister');
end

% ---------- Differences and observed statistic ----------
d = x - y;
obsStat = opt.statFun(d);

if ~isscalar(obsStat) || ~isfinite(obsStat)
    error('signFlip_permTest2_paired:BadStatFun', ...
        'statFun must return a finite scalar. Got: %s', mat2str(obsStat));
end

% ---------- Sign-flip permutations (Monte Carlo) ----------
if opt.returnNull
    permStats = zeros(opt.nPerm, 1);
else
    permStats = [];
end

b = 0; % count of "as or more extreme" permutations

for i = 1:opt.nPerm
    % Random sign for each pair: +1 or -1 with equal probability
    signs = (rand(nPairs, 1) > 0.5) * 2 - 1; % gives +1/-1
    dPerm = d .* signs;

    s = opt.statFun(dPerm);

    if ~isscalar(s) || ~isfinite(s)
        error('signFlip_permTest2_paired:BadStatDuringPerm', ...
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
info.method       = 'paired_sign_flip';
info.nPerm        = opt.nPerm;
info.tail         = char(tail);
info.nPairs       = nPairs;
info.seedUsed     = seedUsed;
info.statFun      = func2str(opt.statFun);
info.countExtreme = b;
info.pValue       = p;
info.nanPairsRemoved = nPairsRemoved;

end
