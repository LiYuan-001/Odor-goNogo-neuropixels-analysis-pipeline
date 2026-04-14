function [acg, t] = ACG_2(ts, varargin)
% ACG_2  Memory-efficient autocorrelogram (ACG) computation
%
%   [acg, t] = ACG_2(ts, 'binSize', 0.001, 'duration', 1, 'norm', 'rate')
%
% Inputs:
%   ts        : spike timestamps (sorted, in seconds)
%
% Optional parameters:
%   'binSize' : bin width (s)
%   'duration': half-width of correlogram (s)
%   'norm'    : 'none' | 'rate' | 'prob'
%
% Outputs:
%   acg : autocorrelogram counts (or rate if normalized)
%   t   : time lags (s)
%
% Note:
%   This implementation avoids building an n×n matrix — it scales ~O(n).

% -------------------------------------------------------
p = inputParser;
addParameter(p, 'binSize', 0.001);
addParameter(p, 'duration', 1);
addParameter(p, 'norm', 'none');
parse(p, varargin{:});

binSize = p.Results.binSize;
duration = p.Results.duration;
normType = lower(p.Results.norm);

% -------------------------------------------------------
ts = sort(ts(:));         % ensure column + sorted
nSpikes = numel(ts);
edges = -duration:binSize:duration;
acg = zeros(1, numel(edges)-1);

% -------------------------------------------------------
% Sliding window approach
j1 = 1;
for i = 1:nSpikes
    % move lower index j1 until ts(i) - ts(j1) <= duration
    while j1 < nSpikes && ts(i) - ts(j1) > duration
        j1 = j1 + 1;
    end
    
    % move forward to include all spikes within +duration
    j2 = i + 1;
    while j2 <= nSpikes && ts(j2) - ts(i) <= duration
        dt = ts(j2) - ts(i);
        binIdx = floor((dt + duration) / binSize) + 1;
        if binIdx >= 1 && binIdx <= numel(acg)
            acg(binIdx) = acg(binIdx) + 1;
        end
        j2 = j2 + 1;
    end
    
    % backward direction (negative lags)
    j3 = i - 1;
    while j3 >= 1 && ts(i) - ts(j3) <= duration
        dt = ts(j3) - ts(i);
        binIdx = floor((dt + duration) / binSize) + 1;
        if binIdx >= 1 && binIdx <= numel(acg)
            acg(binIdx) = acg(binIdx) + 1;
        end
        j3 = j3 - 1;
    end
end

t = edges(1:end-1) + binSize/2;

% -------------------------------------------------------
% Normalization
switch normType
    case 'rate'
        acg = acg / (nSpikes * binSize);  % Hz
    case 'prob'
        acg = acg / sum(acg);
end

end