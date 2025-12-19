function [mean_nnd, R, p_val] = compute_spatial_clustering(all_coords, idx_encoding, x_range, y_range, n_perm)
% compute_spatial_clustering
% Computes Mean NND and Clark-Evans R index for encoding neurons.
% Performs permutation test for significance of R.
%
% INPUTS:
%   all_coords   - N x 2 matrix of (x, y) neuron coordinates
%   idx_encoding - vector of indices of encoding neurons
%   x_range      - (optional) field of view width (default = 512)
%   y_range      - (optional) field of view height (default = 512)
%   n_perm       - (optional) number of permutations (default = 1000)
%
% OUTPUTS:
%   mean_nnd     - mean nearest neighbor distance among encoding neurons
%   R            - Clark-Evans aggregation index
%   p_val        - p-value from permutation test (two-tailed)

    if nargin < 3, x_range = 512; end
    if nargin < 4, y_range = 512; end
    if nargin < 5, n_perm = 1000; end

    N = size(all_coords, 1);
    num_encoding = numel(idx_encoding);

    % Early exit if too few encoding neurons
    if num_encoding < 5
        warning('Too few encoding neurons (n = %d). Returning NaN.', num_encoding);
        mean_nnd = NaN; R = NaN; p_val = NaN;
        return;
    end

    % --- Observed encoding coordinates
    encoding_coords = all_coords(idx_encoding, :);

    % Pairwise distances and mean NND
    dist_matrix = squareform(pdist(encoding_coords));
    dist_matrix(1:num_encoding+1:end) = Inf;
    nnd = min(dist_matrix, [], 2);
    mean_nnd = mean(nnd);

    % Clark-Evans R index
    area_fov = x_range * y_range;
    density = num_encoding / area_fov;
    expected_nnd = 0.5 / sqrt(density);
    R = mean_nnd / expected_nnd;

    % --- Permutation test
    R_null = zeros(n_perm, 1);
    for i = 1:n_perm
        rand_idx = randperm(N, num_encoding);
        rand_coords = all_coords(rand_idx, :);
        d = squareform(pdist(rand_coords));
        d(1:num_encoding+1:end) = Inf;
        rand_nnd = min(d, [], 2);
        rand_mean_nnd = mean(rand_nnd);
        rand_density = num_encoding / area_fov;
        expected_nnd_perm = 0.5 / sqrt(rand_density);
        R_null(i) = rand_mean_nnd / expected_nnd_perm;
    end

    % Two-tailed p-value
    p_val = mean(abs(R_null - 1) >= abs(R - 1));

    % --- Output
    fprintf('Mean NND: %.4f\n', mean_nnd);
    fprintf('Clark-Evans R: %.4f (p = %.4f)\n', R, p_val);
    if (R < 1) && (p_val<.05)
        disp('→ Encoding neurons are spatially clustered.');
    elseif (R > 1) && (p_val<.05)
        disp('→ Encoding neurons are spatially dispersed.');
    else
        disp('→ Encoding neurons are randomly distributed.');
    end
end