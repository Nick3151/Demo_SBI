function X = build_data_matrix(X_all, neurons, frames, inds)
% BUILD_DATA_MATRIX  Concatenate selected fixed-length trials into one matrix.
%
%   X = build_data_matrix(X_all, neurons, frames, inds)
%
%   Extracts a subset of trials from a session-wide activity matrix and
%   stitches them together, in the order given by INDS, into one continuous
%   matrix.
%
%   INPUTS
%     X_all   - [neurons x (nTrials*frames)] activity matrix for the whole
%               session, with trials stored back-to-back along the columns.
%     neurons - number of rows (neurons) in the output; should match
%               size(X_all,1). Used only to preallocate X (no row selection
%               is performed).
%     frames  - number of frames (columns) per trial.
%     inds    - vector of trial indices (1-based) to extract, in the desired
%               output order.
%
%   OUTPUT
%     X       - [neurons x (length(inds)*frames)] matrix containing the
%               requested trials concatenated along the columns.
X = zeros(neurons, length(inds)*frames);
for i=1:length(inds)
    trial = inds(i);
    X(:,(i-1)*frames+1:i*frames) = X_all(:,(trial-1)*frames+1:trial*frames);
end