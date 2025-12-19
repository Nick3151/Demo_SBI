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

%% Visualize motifs with time traces
frames = round(5*fs);
X = [X_train, X_test, zeros(N,L)];
Xhat_train = helper.reconstruct(W,H);
Xhat_test = helper.reconstruct(W,H_test);
Xhat = [Xhat_train, Xhat_test, zeros(N,L)];

k = 1;
Wk = squeeze(W(:,k,:));

inds_trials_tmp = inds_trials_all{k};
active_time_motifs_tmp = round(active_time_motifs_all{k});
X_tmp = zeros(N,length(inds_trials_tmp), L);
Xhat_tmp = zeros(N,length(inds_trials_tmp), L);
for i=1:length(inds_trials_tmp)
    X_tmp(:,i,:) = X(:,active_time_motifs_tmp(i)+(1:L));
    Xhat_tmp(:,i,:) = Xhat(:,active_time_motifs_tmp(i)+(1:L));
end
X_tmp_mean = squeeze(mean(X_tmp, 2));
X_tmp_std = squeeze(std(X_tmp, [], 2));

neurons_active = neurons_active_motifs{k};
Wk_active = Wk(neurons_active, :);
Xhat_tmp_active = Xhat_tmp(neurons_active,:,:);
X_mean_active = X_tmp_mean(neurons_active, :);
X_std_active = X_tmp_std(neurons_active, :);
X_tmp_active = X_tmp(neurons_active,:,:);
N_active = length(neurons_active);

figure; 
set(gcf,'position',[200,200,900,800])
dnW = prctile(Wk_active(:),99.9); 
XsEdge = [0, 0, L]/fs; 
YsEdge = [0, dnW*N_active, dnW*N_active];

subplot(141)
for ni =1:N_active
    Xs = [1, 1:L, L]/fs;
    Ys = [dnW*(N_active-ni) dnW*(N_active-ni) + Wk_active(ni, :) dnW*(N_active-ni)];
    patch(Xs,Ys, 'k', 'edgecolor', 'none')
    hold on
end
xline(onsets_motifs{k}/fs, 'r', 'LineWidth',1)

plot(XsEdge,YsEdge, 'LineWidth', 2);
xlim([0 L]/fs);ylim([0 dnW*N_active])
set(gca, 'YTick', [])

subplot(142)
hold on
dnX = prctile(Xhat_tmp_active(:),99.9); 
for ni =1:N_active
    Xs = (1:L)/fs;
    for i=1:length(inds_trials_tmp)
        Ys = dnW*(N_active-ni) + squeeze(Xhat_tmp_active(ni, i, :))/dnX*dnW;
        plot(Xs, Ys, 'Color', [0,0,0,.3], 'LineWidth', 1)
    end

end
plot(XsEdge,YsEdge, 'LineWidth', 2);
xlim([0 L]/fs);ylim([0 dnW*N_active])
set(gca, 'YTick', [])

subplot(143)
hold on
dnX = prctile(X_mean_active(:),99.9); 
for ni =1:N_active
    Xs = (1:L)/fs;
    Ys = dnW*(N_active-ni) + X_mean_active(ni, :)/dnX*dnW;
    Ys_low = Ys - X_std_active(ni, :)/dnX*dnW;
    Ys_high = Ys + X_std_active(ni, :)/dnX*dnW;
    plot(Xs, Ys, 'k', 'LineWidth', 2)
    patch([Xs, flip(Xs)], [Ys_low flip(Ys_high)], 'k','FaceAlpha',.3, 'EdgeColor', 'none')
end
plot(XsEdge,YsEdge, 'LineWidth', 2);
xlim([0 L]/fs);ylim([0 dnW*N_active])
set(gca, 'YTick', [])

subplot(144)
hold on
dnX = prctile(X_tmp_active(:),99.9); 
for ni =1:N_active
    Xs = (1:L)/fs;
    for i=1:length(inds_trials_tmp)
        Ys = dnW*(N_active-ni) + squeeze(X_tmp_active(ni, i, :))/dnX*dnW;
        plot(Xs, Ys, 'Color', [0,0,0,.3], 'LineWidth', 1)
    end

end
plot(XsEdge,YsEdge, 'LineWidth', 2);
xlim([0 L]/fs);ylim([0 dnW*N_active])
set(gca, 'YTick', [])
