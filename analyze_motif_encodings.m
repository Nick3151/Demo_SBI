function analyze_motif_encodings(input_folder_name)
% Analyze motif encodings of stimulus, choice, motion and rewards and save
root = pwd;
folder_dir = fullfile(root, 'Data');

output_path = fullfile(folder_dir, input_folder_name);
addpath(fullfile(root, 'Utils'))
file_dFF = 'dFF_preprocessed.mat';
file_beh = 'bh_data_Liang.mat';
load(fullfile(output_path, 'FlexMF_Motifs.mat'))
load(fullfile(folder_dir, input_folder_name, file_beh));
load(fullfile(folder_dir, input_folder_name, file_dFF));

% trials_onset = trials_onset(1:end-1);
[N,K,L] = size(W);
motif_freq = cellfun(@(x) length(x), inds_trials_all);
min_occur_times = 4;
num_neurons_active = cellfun(@(x) length(x), neurons_active_motifs);
% Select meaningful motifs (be significant and occurred at least 4 times)
k_sel = find((motif_freq>=min_occur_times) & is_significant' & num_neurons_active>1);

%% Motif encoding stimulus and choices
[stim_encoding, choice_encoding, joint_encoding, ~, ~, ~] = compute_motifs_encoding_index(motifs_trials_match);
[h_stim, h_choice, h_joint, ind_stim, ind_choice, ind_joint, p_choice, min_p, is_testable] = test_motif_encodings_fisher(motifs_trials_match, k_sel);

%% load motion, reward signal and licking
T = size(dFF_preprocess,2);

if size(bh_data,2)>length(trials_onset)
    bh_data = bh_data(:,1:length(trials_onset));
end

deg_rotary = deg_rotary*.6;
rot_speed = [0, diff(deg_rotary)];
frames_lick = frames_lick(1<frames_lick & frames_lick<T);
frames_result = frames_result(1<=frames_result & frames_result<=trials_onset(end)+5*fs+1);
frames_result_trials = round(frames_result-trials_onset);
frames_reward_trials = frames_result_trials;
frames_reward_trials(bh_data(3,:)~=1) = -1;
licking = zeros(1, T);
licking(round(frames_lick)) = 1;

%% Plot motion against each motif occurrence
trial_frames = round(5*fs);
ind_all = [ind_train; ind_test];
rot_speed_pad = [zeros(1,L), rot_speed, zeros(1,L)];
rot_degree_pad = [zeros(1,L), deg_rotary, zeros(1,L)];
licking_pad = [zeros(1,L), licking, zeros(1,L)];
h_start = zeros(K,1);
h_stop = zeros(K,1);
h_reward = zeros(K,1);
p_start = ones(K,1);
p_stop = ones(K,1);
p_reward = ones(K,1);
motif_onsets_reward = cell(K,1);
R_index = nan(K,1);
p_cluster = nan(K,1);

for j=1:length(k_sel)
    k = k_sel(j);
    if ~isempty(neurons_active_motifs{k})
        inds_trials_tmp = inds_trials_all{k};
        active_time_motifs_tmp = round(active_time_motifs_all{k});
        onset_time_tmp = onset_time_trials_all{k}';
        frame_reward_tmp = frames_reward_trials(ind_all(inds_trials_tmp))';
        
        rot_speed_motifs = zeros(length(inds_trials_tmp), 2*L);
        rot_speed_trials = zeros(length(inds_trials_tmp), trial_frames);
        rot_degree_motifs = zeros(length(inds_trials_tmp), 2*L);
        rot_degree_trials = zeros(length(inds_trials_tmp), trial_frames);
        licking_motifs = zeros(length(inds_trials_tmp), 2*L);
        licking_trials = zeros(length(inds_trials_tmp), trial_frames);
        
        for i=1:length(inds_trials_tmp)
            rot_speed_motifs(i,:) = rot_speed_pad(trials_onset(ind_all(inds_trials_tmp(i)))+round(onset_time_tmp(i))+(1:2*L));
            rot_speed_trials(i,:) = rot_speed_pad(trials_onset(ind_all(inds_trials_tmp(i)))+L+(1:trial_frames));
            rot_degree_motifs(i,:) = rot_degree_pad(trials_onset(ind_all(inds_trials_tmp(i)))+round(onset_time_tmp(i))+(1:2*L))...
            -rot_degree_pad(trials_onset(ind_all(inds_trials_tmp(i)))+L);
            rot_degree_trials(i,:) = rot_degree_pad(trials_onset(ind_all(inds_trials_tmp(i)))+L+(1:trial_frames))...
            -rot_degree_pad(trials_onset(ind_all(inds_trials_tmp(i)))+L);
            licking_motifs(i,:) = licking_pad(trials_onset(ind_all(inds_trials_tmp(i)))+round(onset_time_tmp(i))+(1:2*L));
            licking_trials(i,:) = licking_pad(trials_onset(ind_all(inds_trials_tmp(i)))+L+(1:trial_frames));
        end
        
        rot_speed_motifs_mean = mean(rot_speed_motifs, 1);
        rot_speed_motifs_std = std(rot_speed_motifs, 1);
        rot_speed_trials_mean = mean(rot_speed_trials, 1);
        rot_speed_trials_std = std(rot_speed_trials, 1);
        
        % Stopping frame of motion in each trial
        [frames_start, frames_stop, ~] = find_start_stop_frame(frames_result_trials(inds_trials_tmp), fs, rot_speed_trials);
        
        % test if motif onsets relative to motion stop/reward is more aligned than
        % motif onsets relative to trial onsets
        [h_start(k),p_start(k)] = vartest2(onset_time_tmp(frames_start>0)-frames_start(frames_start>0),...
            onset_time_tmp(frames_start>0), 'Tail','left', 'Alpha',.05/length(k_sel));
        [h_stop(k),p_stop(k)] = vartest2(onset_time_tmp(frames_stop>0)-frames_stop(frames_stop>0),...
            onset_time_tmp(frames_stop>0), 'Tail','left', 'Alpha',.05/length(k_sel));
        [h_reward(k),p_reward(k)] = vartest2(onset_time_tmp(frame_reward_tmp>0)-frame_reward_tmp(frame_reward_tmp>0),...
            onset_time_tmp(frame_reward_tmp>0), 'Tail','left', 'Alpha',.05/length(k_sel));
        h_start(isnan(h_start)) = 0;
        h_stop(isnan(h_stop)) = 0;
        h_reward(isnan(h_reward)) = 0;
        if h_reward(k)
            motif_onsets_reward{k} = onset_time_tmp(frame_reward_tmp>0)-frame_reward_tmp(frame_reward_tmp>0);
        end
        
        % Distribution of motif neurons, Clark-Evans R index
        neurons_active = neurons_active_motifs{k};
        if length(neurons_active)>=5
            if isempty(regexp(input_folder_name, '.+ML102623', 'once'))
                [mean_nnd, R, p] = compute_spatial_clustering(Locations', neurons_good(neurons_active));
            else
                [mean_nnd, R, p] = compute_spatial_clustering(Locations', neurons_good(neurons_active), 1024, 1024);
            end
            R_index(k) = R;
            p_cluster(k) = p;
        end
    end
end

save(fullfile(output_path, 'Motifs_encodings.mat'), 'k_sel', 'h_stim', 'h_choice', ...
    'h_joint', 'ind_stim', 'ind_choice', 'ind_joint', 'stim_encoding',...
    'choice_encoding', 'joint_encoding', 'h_start', 'h_stop', 'h_reward', ...
    'p_start', 'p_stop', 'p_reward', 'p_choice', 'is_testable', 'motif_onsets_reward', ...
    'R_index', 'p_cluster');