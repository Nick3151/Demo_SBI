function choose_params(input_folder_name)
% We want the temporal loadings H to be sparse, motifs W themselves do
% not need to be saprse, we thus set alpha_W to be 1e-6 and run SBI
% with different lambdas and alpha_Hs.

% choose_params.mat contains the SBI results under each combination of lambda and alpha_H
% Ws: nLambdas x nAlpha_Hs cell, motif profiles W
% Hs: nLambdas x nAlpha_Hs cell, motif loadings H
% recon_errors: nLambdas x nAlpha_Hs array, reconstruction cost
% reg_crosses: nLambdas x nAlpha_Hs array, regularization cost
% sparseness_Hs: nLambdas x nAlpha_Hs cell, proportion of non-zero values in H
% is_significant: nLambdas x nAlpha_Hs cell, if each motif is significant
% is_single: nLambdas x nAlpha_Hs cell, if each motif only contains one active neuron
% pvals: nLambdas x nAlpha_Hs cell, motif p values
% corr_prctiles_train: nLambdas x nAlpha_Hs array, 95th percentile of the correlation structures in the 'motifs'
% corr_prctiles_noise: nLambdas x nAlpha_Hs array, 95th percentile of the correlation structures in the 'noise'

root = pwd;
folder_dir = fullfile(root, 'Data');

addpath(genpath(fullfile(root, 'FlexMF')));
addpath(fullfile(root, 'Utils'))
addpath(fullfile(root, 'TFOCS'))

seed = 1;

file_dFF = 'dFF_preprocessed.mat';
file_beh = 'bh_data_Liang.mat';
file_transients = 'SEUDO_autoclassify.mat';
load(fullfile(folder_dir, input_folder_name, file_dFF));
load(fullfile(folder_dir, input_folder_name, file_beh));
load(fullfile(folder_dir, input_folder_name, file_transients));

%% Preprocess data (Based on SEUDO transients)
[N,T] = size(dFF_norm);
mask = zeros(N,T);
counts = zeros(N,1);
for n=1:N
    transient_times = transientInfo(n).times;
    classification = transientInfo(n).classification;
    counts(n) = sum(classification==1);
    for i=1:size(transient_times, 1)
        mask(n,transient_times(i,1):transient_times(i,2)) = (classification(i)==1);
    end
end
X = dFF_preprocess.*mask;
X(X<0) = 0;

res_thresh = 8;
neurons_good = find(counts>res_thresh & ~[transientInfo.isArtifact]');

X = X(neurons_good, :);
[N,T] = size(X);

trials = length(trials_onset);

%% Split into training and test set, stratified partition
rng(seed) 
cv = cvpartition(groups, "KFold",2);
ind_train = find(training(cv,1));
ind_test = find(test(cv,1));
frames = round(5*fs);
X_train = [];
X_test = [];
for i=1:length(ind_train)
    X_train = [X_train, X(:,trials_onset(ind_train(i))+(0:frames-1))];
end
for i=1:length(ind_test)
    X_test = [X_test, X(:,trials_onset(ind_test(i))+(0:frames-1))];
end
T_train = size(X_train,2);
T_test = size(X_test,2);

%% Run FlexMF under different parameters
parpool(12)
sprintf('Max number of computational thread: %d', maxNumCompThreads)
K = 10;

if isempty(regexp(input_folder_name, '.+ML102623', 'once'))
    L = 100;
else
    L = 50;
end

% Normalize training data
frob_norm = norm(X_train(:));
X_train = X_train/frob_norm*K;

nLambdas = 9;
nAlpha_Hs = 9;
lambdas = sort(logspace(-1,-5,nLambdas), 'ascend'); 
alpha_Hs = sort(logspace(-2,-6,nAlpha_Hs), 'ascend'); 
alpha_W = 1e-6;

recon_errors = zeros(nLambdas, nAlpha_Hs);
reg_crosses = zeros(nLambdas, nAlpha_Hs);
sparseness_Hs = cell(nLambdas, nAlpha_Hs);
is_significant = cell(nLambdas, nAlpha_Hs);
is_single = cell(nLambdas, nAlpha_Hs);
pvals = cell(nLambdas, nAlpha_Hs);
Ws = cell(nLambdas, nAlpha_Hs);
Hs = cell(nLambdas, nAlpha_Hs);
corr_prctiles_train = zeros(nLambdas, nAlpha_Hs);
corr_prctiles_noise = zeros(nLambdas, nAlpha_Hs);

parfor idx=1:nLambdas*nAlpha_Hs
    
    [li,hi] = ind2sub([nLambdas, nAlpha_Hs], idx);
    display(['Testing lambda ' num2str(li) '/' num2str(nLambdas)])
    display(['Testing alpha_H ' num2str(hi) '/' num2str(nAlpha_Hs)])
    lambda = lambdas(li);
    alpha_H = alpha_Hs(hi);

    tic
    [W,H,cost,errors,loadings_train,power_train] = FlexMF(X_train,'K',K, 'L', L, 'maxiter', 50,...
        'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, 'neg_prop', 0, 'showPlot', 0, 'verbal', 0);
    toc
    Ws{idx} = W;
    Hs{idx} = H;

    [recon_errors(idx), reg_crosses(idx), reg_W, reg_H] = helper.get_FlexMF_cost(X_train,W,H);

    sparseness_Hs{idx} = mean(H>1e-3, 2);

    p = .05; % desired p value for factors
    display('Testing significance of factors')
    [pvals{idx},is_significant{idx}, is_single{idx}] = test_significance_new(X_test, W, [], p);

    X_train_motifs = helper.reconstruct(W,H);
    X_train_noise = X_train-X_train_motifs;

    display('Evaluate cross correlations')
    tic
    fs = 30;
    C_train = helper.max_xcorr(X_train, 2*fs);
    C_train_sel = triu(C_train,1);
    corr_prctiles_train(idx) = prctile(C_train_sel(:), 95);
    C_train_noise = helper.max_xcorr(X_train_noise, 2*fs);
    C_train_noise_sel = triu(C_train_noise,1);
    corr_prctiles_noise(idx) = prctile(C_train_noise_sel(:), 95);
    toc
end

output_path = fullfile(folder_dir, input_folder_name);
save(fullfile(output_path, 'choose_params.mat'), "recon_errors", "reg_crosses", "sparseness_Hs", "is_significant", "pvals", "is_single",...
    "Hs", "Ws", "corr_prctiles_train", "corr_prctiles_noise", "lambdas", "alpha_Hs", "X_train", "ind_train")
