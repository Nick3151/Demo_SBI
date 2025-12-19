%% Load data
clear all;
close all

root = pwd;
folder_dir = fullfile(root, 'Data');

addpath(genpath(fullfile(root, 'FlexMF')));
addpath(fullfile(root, 'Utils'))
addpath(fullfile(root, 'TFOCS'))

seed = 1;

input_folder_name = '20240618MH041724';
output_path = fullfile(folder_dir, input_folder_name);

file_dFF = 'dFF_preprocessed.mat';
file_beh = 'bh_data_Liang.mat';
file_transients = 'SEUDO_autoclassify.mat';
load(fullfile(folder_dir, input_folder_name, file_dFF));
load(fullfile(folder_dir, input_folder_name, file_beh));
load(fullfile(folder_dir, input_folder_name, file_transients));

%% Preprocess data (Based on SEUDO transients)
res_thresh = 8;

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

f1 = figure;
SimpleXplot([X_train, X_test], [length(ind_train), length(ind_test)], 'frames', frames); 
set(f1,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

%% Normalize training data
K = 10;
frob_norm = norm(X_train(:));
X_train = X_train/frob_norm*K;

%% Find motifs
load('FlexMF_params_all.mat')
params = parameters{input_folder_name, :};
lambda = params(1);
alpha_W = params(2);
alpha_H = params(3);
L = 100;

tic
[W,H,cost,errors,loadings_train,power_train] = FlexMF(X_train,'K',K, 'L', L, 'maxiter', 50,...
    'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, 'neg_prop', 0, 'showPlot', 1);
toc

p = .05; % desired p value for motifs
display('Testing significance of factors on held-out data')
[pvals,is_significant,is_single] = test_significance_new(X_test, W, [], p);

[max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
indSort = hybrid(:,3);

plotAll = 0;
figure; SimpleWHPlot(W(indSort,:,:), H, 'trials', length(ind_train), 'frames', frames, 'is_significant', is_significant, 'plotAll', plotAll)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot(W(indSort,:,:), H, 'trials', length(ind_train), 'frames', frames, 'Data', X_train(indSort,:), 'is_significant', is_significant, 'plotAll', plotAll)
title('FlexMF factors, with raw data')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

%% refit W to test data
frob_norm = norm(X_test(:));
X_test = X_test/frob_norm*K;

display('Refit W to test data')

[W,H_test,~,~,loadings_test,power_test] = FlexMF(X_test,'K',K, 'L', L, 'maxiter', 1,...
    'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, ...
    'W_fixed', 1, 'W_init', W, 'SortFactors', 0, 'neg_prop', 0, 'showPlot', 1);

figure; SimpleWHPlot(W(indSort,:,:), H_test, 'trials', length(ind_test), 'frames', frames, 'is_significant', is_significant, 'plotAll', plotAll)
title('FlexMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot(W(indSort,:,:), H_test, 'trials', length(ind_test), 'frames', frames, 'Data', X_test(indSort,:), 'is_significant', is_significant, 'plotAll', plotAll)
title('FlexMF factors, with raw data')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

%% Determine motif activation time
[active_time_motifs_all, magnitudes_motifs_all, neurons_active_motifs, onsets_motifs] = analyze_motifs_activation([X_train, X_test], [H,H_test],W, 'thresh_active_freq', .6);
[inds_trials_all, onset_time_trials_all] = match_motif_times_to_trials(active_time_motifs_all, onsets_motifs, 'frames', frames);
motif_freq = cellfun(@(x) length(x), inds_trials_all);
min_occur_times = 4;
num_neurons_active = cellfun(@(x) length(x), neurons_active_motifs);
% Select meaningful motifs (be significant and occurred at least 4 times)
k_sel = find((motif_freq>=min_occur_times) & is_significant' & num_neurons_active>0);

%% Look at motifs on test data
figure; SimpleWHPlot(W(indSort,:,:), H_test, 'trials', length(ind_test), 'frames', frames, 'is_significant', is_significant, 'plotAll', plotAll)
title('FlexMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot(W(indSort,:,:), H_test, 'trials', length(ind_test), 'frames', frames, 'Data', X_test(indSort,:), 'is_significant', is_significant, 'plotAll', plotAll)
title('FlexMF factors, with raw data')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])

%% Compare occurance times against behaviors and stimulus
groups_new = groups([ind_train; ind_test]);
ind_HitLow = find(strcmp(groups_new, 'HitLow'))';
ind_MissLow = find(strcmp(groups_new, 'MissLow'))';
ind_ErrLow = find(strcmp(groups_new, 'ErrLow'))';
ind_HitHigh = find(strcmp(groups_new, 'HitHigh'))';
ind_MissHigh = find(strcmp(groups_new, 'MissHigh'))';
ind_ErrHigh = find(strcmp(groups_new, 'ErrHigh'))';

%% Create table and save
varNames = {'HitLow','MissLow','ErrLow','HitHigh','MissHigh','ErrHigh'};
rowNames = cell(K+1,1);
rowNames{1} = 'Ground Truth';
for k=1:K
    rowNames{k+1} = ['Motif #', num2str(k)];
end

overlap_HitLow = cellfun(@(x) intersect(x,ind_HitLow), inds_trials_all, 'UniformOutput', false);
overlap_HitLow = [{ind_HitLow}; overlap_HitLow];

overlap_MissLow = cellfun(@(x) intersect(x,ind_MissLow), inds_trials_all, 'UniformOutput', false);
overlap_MissLow = [{ind_MissLow}; overlap_MissLow];

overlap_ErrLow = cellfun(@(x) intersect(x,ind_ErrLow), inds_trials_all, 'UniformOutput', false);
overlap_ErrLow = [{ind_ErrLow}; overlap_ErrLow];

overlap_HitHigh = cellfun(@(x) intersect(x,ind_HitHigh), inds_trials_all, 'UniformOutput', false);
overlap_HitHigh = [{ind_HitHigh}; overlap_HitHigh];

overlap_MissHigh = cellfun(@(x) intersect(x,ind_MissHigh), inds_trials_all, 'UniformOutput', false);
overlap_MissHigh = [{ind_MissHigh}; overlap_MissHigh];

overlap_ErrHigh = cellfun(@(x) intersect(x,ind_ErrHigh), inds_trials_all, 'UniformOutput', false);
overlap_ErrHigh = [{ind_ErrHigh}; overlap_ErrHigh];

motifs_trials_match = table(overlap_HitLow, overlap_MissLow, overlap_ErrLow, overlap_HitHigh,...
    overlap_MissHigh, overlap_ErrHigh, 'VariableNames', varNames, 'RowNames', rowNames);

save(fullfile(output_path, 'FlexMF_Motifs.mat'),...
    'W', 'H', 'H_test', 'loadings_train', 'loadings_test', 'power_train', 'power_test', ...
    'is_significant', 'X_train', 'X_test', 'ind_train', 'ind_test', 'indSort', ...
    'motifs_trials_match', 'active_time_motifs_all', 'magnitudes_motifs_all', 'inds_trials_all', ...
    'neurons_active_motifs', 'onset_time_trials_all', 'onsets_motifs', 'idx_neuron', 'neurons_good')