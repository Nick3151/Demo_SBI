function motif_reliability = compute_motif_reliability(W, X, Xhat, fs, neurons_active_motifs, active_time_motifs, onsets_motifs, k, thresh_onset, tie_tolerance_frames)
% Compute reliability metrics for the k-th motif.
%
% Inputs:
%   W                       N x K x L motif tensor
%   X                       N x T real thresholded dF/F data
%   Xhat                    N x T reconstructed data; kept for call consistency, not used here
%   fs                      sampling rate
%   neurons_active_motifs   cell array; neurons_active_motifs{k} gives active neurons for motif k
%   active_time_motifs      cell array; active_time_motifs{k} gives occurrence times for motif k
%   onsets_motifs           cell array; kept for call consistency, not used here
%   k                       motif index
%   thresh_onset            optional; threshold for determining motif onsets, default 0.05
%   tie_tolerance_frames    optional; timing difference treated as a tie, default 0
%
% Outputs:
%   motif_reliability.S_LOO
%       R x 1 cross-validated whole-motif reproducibility.
%       Computed as cosine similarity between each held-out binary event mask
%       and the leave-one-occurrence-out mean binary event mask.
%
%   motif_reliability.q_active
%       N_active x 1 participation reliability of each motif neuron.
%
%   motif_reliability.Z_active
%       N_active x R logical matrix. Z_active(n,r)=1 means neuron n participated
%       in occurrence r inside its motif support window.
%
%   motif_reliability.K_tau_b
%       R x 1 Kendall tau-b sequence reliability per occurrence.


if nargin < 9 || isempty(thresh_onset)
    thresh_onset = 0.05;
end

if nargin < 10 || isempty(tie_tolerance_frames)
    tie_tolerance_frames = 0;
end

tie_tolerance_frames = max(0, round(tie_tolerance_frames));

[N,K,L] = size(W);

if k < 1 || k > K
    error('k must be between 1 and K.');
end

if size(X,1) ~= N
    error('X must have the same number of rows as W.');
end

Wk = squeeze(W(:,k,:));
neurons_active = neurons_active_motifs{k};
neurons_active = neurons_active(:);

active_times = round(active_time_motifs{k});
active_times = active_times(:);

R = length(active_times);

if R == 0
    error('No motif occurrences found for this motif.');
end

T_total = size(X,2);

if any(active_times < 0) || any(active_times + L > T_total)
    error('Some motif occurrence windows are outside X. Check active_time_motifs{k} and L.');
end

X_tmp = zeros(N, R, L);

for i = 1:R
    X_tmp(:,i,:) = X(:,active_times(i)+(1:L));
end

Wk_active = Wk(neurons_active, :);
X_tmp_active = X_tmp(neurons_active,:,:);
N_active = length(neurons_active);

% Motif-defined onset time for each active neuron.
onsets_motif_neurons = find_onset_neurons(Wk_active, thresh_onset);

% Observed onset of each active neuron in each motif occurrence.
onsets_traces = nan(N_active,R);

for r = 1:R
    onsets_traces(:,r) = find_onset_neurons(squeeze(X_tmp_active(:,r,:)), thresh_onset);
end

Z_active = ~isnan(onsets_traces);

% Same visual sorting idea, but here sorted by motif-defined onset.
sort_key = onsets_motif_neurons;
sort_key(isnan(sort_key)) = inf;
[~,idx_sort] = sort(sort_key);

neurons_active_sorted = neurons_active(idx_sort);
onsets_motif_neurons_sorted = onsets_motif_neurons(idx_sort);
onsets_traces_sorted = onsets_traces(idx_sort,:);

% -------------------------------------------------------------------------
% 1. Cross-validated whole-motif reproducibility
% -------------------------------------------------------------------------
% Uses continuous X_tmp_active, not binarized traces and not weighted by W.
%
% S_LOO(r) compares occurrence r with the leave-one-occurrence-out
% mean real trace profile.
%
% This measures whether real motif occurrences reproduce a common
% continuous spatiotemporal profile.

S_LOO = nan(R,1);

if R >= 2
    for r = 1:R
        ind_train = true(1,R);
        ind_train(r) = false;

        % Held-out real occurrence: N_active x L
        X_r = reshape(X_tmp_active(:,r,:), N_active, L);

        % Leave-one-occurrence-out real template: N_active x L
        X_template = reshape(mean(X_tmp_active(:,ind_train,:), 2), N_active, L);

        % Continuous whole-motif reproducibility
        S_LOO(r) = local_cosine( ...
            double(X_r(:)), ...
            double(X_template(:)));
    end
end

S_LOO_summary = local_summary(S_LOO);

% -------------------------------------------------------------------------
% 2. Reliability of participating neurons
% -------------------------------------------------------------------------
q_active = mean(double(Z_active), 2);
n_occurrences_active = sum(Z_active, 2);

% Occurrence-level completeness is useful to report beside q_active.
C_occurrence = mean(double(Z_active), 1)';

latency_error_frames = onsets_traces - repmat(onsets_motif_neurons, [1,R]);

latency_median_frames = local_nanmedian_rows(latency_error_frames);
latency_jitter_frames = local_nanstd_rows(latency_error_frames);

latency_median_sec = latency_median_frames ./ fs;
latency_jitter_sec = latency_jitter_frames ./ fs;

neuron_reliability_table = table( ...
    neurons_active(:), ...
    onsets_motif_neurons(:), ...
    onsets_motif_neurons(:)./fs, ...
    q_active(:), ...
    n_occurrences_active(:), ...
    latency_median_frames(:), ...
    latency_median_sec(:), ...
    latency_jitter_frames(:), ...
    latency_jitter_sec(:), ...
    'VariableNames', { ...
        'neuron_id', ...
        'motif_onset_frame', ...
        'motif_onset_sec', ...
        'q_active', ...
        'n_occurrences_active', ...
        'latency_median_frames', ...
        'latency_median_sec', ...
        'latency_jitter_frames', ...
        'latency_jitter_sec' ...
    });

% -------------------------------------------------------------------------
% 3. Temporal order reliability
% -------------------------------------------------------------------------
% Uses Kendall's tau-b for all motifs.
%
% For N_active >= 2, K_tau_b(r) is computed for each occurrence r when at
% least two motif neurons participated in that occurrence.
%
% For N_active == 2, K_tau_b(r) reduces to:
%   +1 if the observed order matches the motif-defined order
%   -1 if the observed order is reversed
%   NaN if there is a tie or fewer than two participating neurons
%
% For N_active < 2, all values remain NaN.

K_tau_b = nan(R,1);
n_active_order = zeros(R,1);

for r = 1:R

    ind_use = find( ...
        Z_active(:,r) & ...
        ~isnan(onsets_motif_neurons) & ...
        ~isnan(onsets_traces(:,r)) ...
    );

    n_active_order(r) = length(ind_use);

    if n_active_order(r) >= 2
        K_tau_b(r) = local_kendall_taub( ...
            onsets_motif_neurons(ind_use), ...
            onsets_traces(ind_use,r), ...
            tie_tolerance_frames);
    end
end

K_tau_b_summary = local_summary(K_tau_b);

% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------
motif_reliability = struct;

motif_reliability.k = k;
motif_reliability.fs = fs;
motif_reliability.thresh_onset = thresh_onset;
motif_reliability.tie_tolerance_frames = tie_tolerance_frames;

motif_reliability.N = N;
motif_reliability.K = K;
motif_reliability.L = L;
motif_reliability.R = R;
motif_reliability.N_active = N_active;

motif_reliability.neurons_active = neurons_active;
motif_reliability.neurons_active_sorted = neurons_active_sorted;
motif_reliability.active_times = active_times;

motif_reliability.Wk_active = Wk_active;
motif_reliability.X_tmp_active = X_tmp_active;

motif_reliability.onsets_motif_neurons = onsets_motif_neurons;
motif_reliability.onsets_motif_neurons_sorted = onsets_motif_neurons_sorted;

motif_reliability.onsets_traces = onsets_traces;
motif_reliability.onsets_traces_sorted = onsets_traces_sorted;

motif_reliability.Z_active = Z_active;
motif_reliability.q_active = q_active;
motif_reliability.C_occurrence = C_occurrence;

motif_reliability.S_LOO = S_LOO;
motif_reliability.S_LOO_summary = S_LOO_summary;

motif_reliability.neuron_reliability_table = neuron_reliability_table;

motif_reliability.K_tau_b = K_tau_b;
motif_reliability.K_tau_b_summary = K_tau_b_summary;
motif_reliability.n_active_order = n_active_order;
motif_reliability.K_tau_b_method = ...
    'Kendall tau-b between motif-defined onset order and observed onset order';

end


% =========================================================================
% Local helper functions
% =========================================================================

function c = local_cosine(a, b)

a = a(:);
b = b(:);

den = norm(a) * norm(b);

if den == 0
    c = nan;
else
    c = (a' * b) / den;
end

end


function s = local_summary(x)

x = x(:);
x = x(~isnan(x));

s = struct;
s.n_valid = numel(x);

if isempty(x)
    s.mean = nan;
    s.median = nan;
    s.sd = nan;
    s.iqr = [nan nan];
else
    s.mean = mean(x);
    s.median = median(x);

    if numel(x) >= 2
        s.sd = std(x,0);
    else
        s.sd = nan;
    end

    s.iqr = prctile(x,[25 75]);
end

end


function y = local_nanmedian_rows(A)

n = size(A,1);
y = nan(n,1);

for i = 1:n
    v = A(i,:);
    v = v(~isnan(v));

    if ~isempty(v)
        y(i) = median(v);
    end
end

end


function y = local_nanstd_rows(A)

n = size(A,1);
y = nan(n,1);

for i = 1:n
    v = A(i,:);
    v = v(~isnan(v));

    if numel(v) >= 2
        y(i) = std(v,0);
    end
end

end


function tau_b = local_kendall_taub(x, y, tie_tolerance_frames)

x = x(:);
y = y(:);

valid = ~isnan(x) & ~isnan(y);
x = x(valid);
y = y(valid);

m = length(x);

if m < 2
    tau_b = nan;
    return;
end

C = 0;
D = 0;
n_tie_x = 0;
n_tie_y = 0;

for a = 1:(m-1)
    for b = (a+1):m

        dx = local_sign_tol(x(a) - x(b), tie_tolerance_frames);
        dy = local_sign_tol(y(a) - y(b), tie_tolerance_frames);

        if dx == 0
            n_tie_x = n_tie_x + 1;
        end

        if dy == 0
            n_tie_y = n_tie_y + 1;
        end

        prod_sign = dx * dy;

        if prod_sign > 0
            C = C + 1;
        elseif prod_sign < 0
            D = D + 1;
        end
    end
end

n0 = m * (m - 1) / 2;

den = sqrt((n0 - n_tie_x) * (n0 - n_tie_y));

if den == 0
    tau_b = nan;
else
    tau_b = (C - D) / den;
end

end


function s = local_sign_tol(x, tol)

if abs(x) <= tol
    s = 0;
elseif x > 0
    s = 1;
else
    s = -1;
end

end