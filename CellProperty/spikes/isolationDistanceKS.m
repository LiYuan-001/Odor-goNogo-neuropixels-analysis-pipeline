function [isoDist, Lratio] = isolationDistanceKS(sc,pc, clusterId, nPCsUse)
% isolationDistanceKS  Isolation Distance (and L-ratio) from Kilosort/Phy files.
%
%   [isoDist, Lratio] = isolationDistanceKS(phyDir, clusterId, 'Name',Value,...)
%
% % Inputs
% sc = readNPY(fullfile(phyDir, 'spike_clusters.npy'));         % [nSpikes, 1]
% pc = readNPY(fullfile(phyDir, 'pc_features.npy'));            % [nSpikes, nPCs,nLocCh]
%   clusterId  : numeric cluster id (matching entries in spike_clusters.npy)
%  nPCsUse: number of PCs per channel to keep

% Outputs
%   isoDist  : Isolation Distance (squared Mahalanobis Nth neighbor)
%   Lratio   : L-ratio (sum of chi2 tail probs / Ncluster), [] if not computed
%
% Notes
%   - Isolation Distance is defined as the Nth smallest squared Mahalanobis
%     distance of NON-cluster spikes to the cluster’s Gaussian, where
%     N = #spikes in the cluster.
%   - If there are < N non-cluster spikes available, isoDist = NaN.

opt.UseZScore = 1;
opt.MaxOther = Inf;
opt.Reg = 1e-6;
% ---- Build 2D feature matrix: [nSpikes x (nLocCh*PCsPerChan)]
X = reshape(pc(:,1:nPCsUse,:), size(pc,1), []);   % flatten local PCs

% (Optional) z-score features to stabilize covariance
if opt.UseZScore
    muX = mean(X,1); sigX = std(X,[],1); sigX(sigX==0)=1;
    X = (X - muX) ./ sigX;
end

inC  = (sc == clusterId);
Xc   = X(inC, :);
Xnot = X(~inC, :);

Nc = size(Xc,1);
D  = size(Xc,2);
isoDist = NaN;  Lratio = NaN;

% sanity checks
if Nc < 10 || Nc <= D
    % too few spikes or too many features → ill-conditioned covariance
    return
end

% Optionally subsample non-cluster spikes for speed
if isfinite(opt.MaxOther) && size(Xnot,1) > opt.MaxOther
    idx = randperm(size(Xnot,1), opt.MaxOther);
    Xnot = Xnot(idx, :);
end

% ---- Fit Gaussian to the cluster
mu = mean(Xc, 1);
C  = cov(Xc);                                % [D x D]
% ridge regularization to avoid singularity
lam = opt.Reg * trace(C) / max(D,1);
C = C + lam * eye(D);

% ---- Mahalanobis distances of NON-cluster spikes to cluster Gaussian
[L, flag] = chol(C, 'lower');
if flag ~= 0
    % covariance still not PD
    return
end
diff = bsxfun(@minus, Xnot, mu);
y = L \ diff';                     % solve L*y = diff'
d2 = sum(y.^2, 1)';                % squared Mahalanobis distances

% Need at least Nc non-cluster spikes
if numel(d2) < Nc
    return
end

% ---- Isolation Distance: Nth smallest D^2
d2s = sort(d2, 'ascend');
isoDist = d2s(Nc);

% ---- L-ratio (optional; useful to return too)
% L-ratio = mean chi-square tail probability of non-cluster spikes
% df = dimensionality D
Lratio = mean(1 - chi2cdf(d2, D));

end
