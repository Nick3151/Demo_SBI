% Visualize motif traces
clear all
close all
clc

root = pwd;
folder_dir = fullfile(root, 'Data');
addpath(genpath(fullfile(root, 'CoDybase-MATLAB')))
addpath(genpath(fullfile(root, 'FlexMF')))

input_folder_name = '20240618MH041724';
output_path = fullfile(folder_dir, input_folder_name);
file_dFF = 'dFF_preprocessed.mat';
file_beh = 'bh_data_Liang.mat';
load(fullfile(output_path, 'FlexMF_Motifs.mat'))
load(fullfile(folder_dir, input_folder_name, file_beh));
load(fullfile(folder_dir, input_folder_name, file_dFF));

% trials_onset = trials_onset(1:end-1);
[N,K,L] = size(W);
motif_freq = cellfun(@(x) length(x), inds_trials_all);
num_neurons_active = cellfun(@(x) length(x), neurons_active_motifs);
min_occur_times = 4;
% Select meaningful motifs (be significant and occurred at least 4 times)
k_sel = find((motif_freq>=min_occur_times) & is_significant' & num_neurons_active>0);

%% Motif encoding task variables
[stim_encoding, choice_encoding, joint_encoding, ind_stim, ind_choice, ind_joint] = compute_motifs_encoding_index(motifs_trials_match);
[h_stim, h_choice, h_joint, ind_stim, ind_choice, ind_joint, p_choice, min_p, is_testable] = test_motif_encodings_fisher(motifs_trials_match, k_sel);

%% Visualize motif onsets
figure;
hold on
for i=1:length(k_sel)
    k = k_sel(i);
    Wk = squeeze(W(:,k,:));
    onset_time = onset_time_trials_all{k};
    x = i*ones(1,length(onset_time));
    y{i} = onset_time/fs;
    swarmchart(x,y{i}, 'filled', 'XJitterWidth', .5);
end
rectangle('Position',[.5, 1, length(k_sel), 1], 'FaceColor', [.5, .5, .5, .3], 'EdgeColor', 'none')
yline(0, 'LineWidth',1)
ylim([0,5])
ylabel('onset time(s)')
xlabel('Motif #')
set(gca, 'FontSize', 14, 'XTick', 1:length(k_sel), 'xticklabels', num2str(k_sel))

% [h1,p1,ci,stats] = ttest2(y{1}, y{2}, 'Tail','left');
% [h2,p2,ci,stats] = ttest2(y{1}, y{3}, 'Tail','left');

%% Motif reliability metrics for motif k
X = [X_train, X_test, zeros(N,L)];
Xhat_train = helper.reconstruct(W,H);
Xhat_test = helper.reconstruct(W,H_test);
Xhat = [Xhat_train, Xhat_test, zeros(N,L)];
tie_tol_frames = 3;

k = 1;
motif_reliability = compute_motif_reliability( ...
W, X, Xhat, fs, ...
neurons_active_motifs, ...
active_time_motifs_all, ...
onsets_motifs, ...
k, [], tie_tol_frames);

%% Visualize motifs with time traces
Wk = squeeze(W(:,k,:));
f1 = visualize_motif_traces(Wk, X, Xhat, fs, motif_reliability, true);
