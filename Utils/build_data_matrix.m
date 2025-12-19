function X = build_data_matrix(X_all, neurons, frames, inds)
X = zeros(neurons, length(inds)*frames);
for i=1:length(inds)
    trial = inds(i);
    X(:,(i-1)*frames+1:i*frames) = X_all(:,(trial-1)*frames+1:trial*frames);
end