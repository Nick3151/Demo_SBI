function analyze_single_neuron_encoding(input_folder_name, varargin)
p = inputParser;
addOptional(p, 'showplot', 1);
addOptional(p, 'res_thresh', 8);
addOptional(p, 'smooth_win', 5);    % smoothing window length for movement encoding
addOptional(p, 'lags', 15); % maximum time lag for movement encoding
addOptional(p, 'Save', 1);

parse(p, varargin{:});
showplot = p.Results.showplot;
res_thresh = p.Results.res_thresh;
smooth_win = p.Results.smooth_win;
lags = p.Results.lags;
Save = p.Results.Save;

root = pwd;
folder_dir = fullfile(root, 'Data');

addpath(fullfile(root, 'Utils'))

disp(input_folder_name)
output_path = fullfile(folder_dir, input_folder_name);

file_dFF = 'dFF_preprocessed.mat';
file_beh = 'bh_data_Liang.mat';
file_transients = 'SEUDO_autoclassify.mat';
load(fullfile(folder_dir, input_folder_name, file_dFF));
load(fullfile(folder_dir, input_folder_name, file_beh));
load(fullfile(folder_dir, input_folder_name, file_transients));

%% Find good neurons
[N,T] = size(dFF);
counts = zeros(N,1);
for n=1:N
    transient_times = transientInfo(n).times;
    classification = transientInfo(n).classification;
    counts(n) = sum(classification==1);
end

neurons_good = find(counts>res_thresh & ~[transientInfo.isArtifact]');
X = dFF(neurons_good, :);
[N,T] = size(X);

trials = length(trials_onset);
trial_frames = round(5*fs);
% Exclude the trial in the end
if size(bh_data,2)>trials
    bh_data = bh_data(:,1:trials);
end
sound_onsets = bh_data(2,:);
ind_Low = find(bh_data(1,:)==0); % Low freq trials
ind_High = find(bh_data(1,:)==1); % High freq trials

%% load motion, reward signal and licking
deg_rotary = deg_rotary*.6;
rot_speed = [0, diff(deg_rotary)];
frames_lick = frames_lick(1<frames_lick & frames_lick<T);
frames_result = frames_result(1<=frames_result & frames_result<=trials_onset(end)+5*fs+1);
frames_result_trials = round(frames_result-trials_onset);
frames_reward_trials = frames_result_trials;
frames_reward_trials(bh_data(3,:)~=1) = -1;

%% Plot trial motion traces for Low and High frequencies
trial_frames = round(5*fs);
Xs = (1:trial_frames)/fs;

if showplot
    figure;
    hold on
    for i=1:trials
        if strcmp(groups{i}, 'HitLow')
            p_HitLow = plot(Xs, deg_rotary(trials_onset(i)+(1:trial_frames))-deg_rotary(sound_onsets(i)-1), 'b');
        elseif strcmp(groups{i}, 'MissLow')
            p_MissLow = plot(Xs, deg_rotary(trials_onset(i)+(1:trial_frames))-deg_rotary(sound_onsets(i)-1), 'k');
        elseif strcmp(groups{i}, 'ErrLow')
            p_ErrLow = plot(Xs, deg_rotary(trials_onset(i)+(1:trial_frames))-deg_rotary(sound_onsets(i)-1), 'r');
        end
    end
%     legend([p_HitLow, p_MissLow, p_ErrLow], {'HitLow', 'MissLow', 'ErrLow'})
    rectangle('Position',[1, -100, 1, 200], 'FaceColor', [.5, .5, .5, .3], 'EdgeColor', 'none')
    ylim([-100,100])
    
    figure;
    hold on
    for i=1:trials
        if strcmp(groups{i}, 'HitHigh')
            p_HitHigh = plot(Xs, deg_rotary(trials_onset(i)+(1:trial_frames))-deg_rotary(sound_onsets(i)-1), 'b');
        elseif strcmp(groups{i}, 'MissHigh')
            p_MissHigh = plot(Xs, deg_rotary(trials_onset(i)+(1:trial_frames))-deg_rotary(sound_onsets(i)-1), 'k');
        elseif strcmp(groups{i}, 'ErrHigh')
            p_ErrHigh = plot(Xs, deg_rotary(trials_onset(i)+(1:trial_frames))-deg_rotary(sound_onsets(i)-1), 'r');
        end
    end
%     legend([p_HitHigh, p_MissHigh, p_ErrHigh], {'HitHigh', 'MissHigh', 'ErrHigh'})
    rectangle('Position',[1, -100, 1, 200], 'FaceColor', [.5, .5, .5, .3], 'EdgeColor', 'none')
    ylim([-100,100])
end

%% Disrtibution of movement onsets, result time and movement offsets
rot_speed_smooth = filtfilt(ones(1,smooth_win)/smooth_win, 1, rot_speed);
rot_speed_trials = zeros(trials, trial_frames);

for i=1:trials
    rot_speed_trials(i,:) = rot_speed(trials_onset(i)+(1:trial_frames));
end

% Starting and stopping frame of motion in each trial
[frames_start_trials, frames_stop_trials, ~] = find_start_stop_frame(frames_result_trials, fs, rot_speed_trials);

ind_Hit = find(strcmp(groups, 'HitLow') | strcmp(groups, 'HitHigh'));
ind_Err = find(strcmp(groups, 'ErrLow') | strcmp(groups, 'ErrHigh'));

if showplot
    figure
    Xs = [ones(1,length(ind_Hit)) 2*ones(1,length(ind_Hit)) 3*ones(1,length(ind_Hit))];
    swarmchart(Xs, [frames_start_trials(ind_Hit)', frames_result_trials(ind_Hit), frames_stop_trials(ind_Hit)']./fs, 12, 'blue', 'filled')
    hold on
    Xs = [ones(1,length(ind_Err)) 2*ones(1,length(ind_Err)) 3*ones(1,length(ind_Err))];
    swarmchart(Xs, [frames_start_trials(ind_Err)', frames_result_trials(ind_Err), frames_stop_trials(ind_Err)']./fs, 12, 'red', 'filled')
    rectangle('Position',[0, 1, 4, 1], 'FaceColor', [.5, .5, .5, .3], 'EdgeColor', 'none')
    set(gca, 'XTick', [1,2,3], 'XTickLabel', {'Movement Start','Result', 'Movement Stop'})
end

%% Sound responsive neurons
dFF_base_mean = zeros(N,trials); % Baseline dF/F. 1s before sound onsets
dFF_win1_mean = zeros(N,trials); % 0~0.5s after sound onsets
dFF_win2_mean = zeros(N,trials); % 0.5~1s after sound onsets
dFF_win3_mean = zeros(N,trials); % 0~0.5s after sound offsets

for i=1:trials
    dFF_base_mean(:,i) = mean(X(:,sound_onsets(i)+round(-fs:0)),2);
    dFF_win1_mean(:,i) = mean(X(:,sound_onsets(i)+round(0:0.5*fs)),2);
    dFF_win2_mean(:,i) = mean(X(:,sound_onsets(i)+round(0.5*fs:fs)),2);
    dFF_win3_mean(:,i) = mean(X(:,sound_onsets(i)+round(fs:1.5*fs)),2);
end

h_high_win1 = false(N,1);
h_low_win1 = false(N,1);
h_high_win2 = false(N,1);
h_low_win2 = false(N,1);
h_high_win3 = false(N,1);
h_low_win3 = false(N,1);

for n=1:N
    [h_high_win1(n), p_high_win1(n)] = ttest(dFF_base_mean(n,ind_High), dFF_win1_mean(n,ind_High), 'Alpha', 0.05/N, 'Tail','left');
    [h_low_win1(n), p_low_win1(n)] = ttest(dFF_base_mean(n,ind_Low), dFF_win1_mean(n,ind_Low), 'Alpha', 0.05/N, 'Tail','left');
    [h_high_win2(n), p_high_win2(n)] = ttest(dFF_base_mean(n,ind_High), dFF_win2_mean(n,ind_High), 'Alpha', 0.05/N, 'Tail','left');
    [h_low_win2(n), p_low_win2(n)] = ttest(dFF_base_mean(n,ind_Low), dFF_win2_mean(n,ind_Low), 'Alpha', 0.05/N, 'Tail','left');
    [h_high_win3(n), p_high_win3(n)] = ttest(dFF_base_mean(n,ind_High), dFF_win3_mean(n,ind_High), 'Alpha', 0.05/N, 'Tail','left');
    [h_low_win3(n), p_low_win3(n)] = ttest(dFF_base_mean(n,ind_Low), dFF_win3_mean(n,ind_Low), 'Alpha', 0.05/N, 'Tail','left');
end
dFF_high_win1 = mean(dFF_win1_mean(:,ind_High),2)-mean(dFF_base_mean(:,ind_High),2);
dFF_high_win2 = mean(dFF_win2_mean(:,ind_High),2)-mean(dFF_base_mean(:,ind_High),2);
dFF_high_win3 = mean(dFF_win3_mean(:,ind_High),2)-mean(dFF_base_mean(:,ind_High),2);
dFF_low_win1 = mean(dFF_win1_mean(:,ind_Low),2)-mean(dFF_base_mean(:,ind_Low),2);
dFF_low_win2 = mean(dFF_win2_mean(:,ind_Low),2)-mean(dFF_base_mean(:,ind_Low),2);
dFF_low_win3 = mean(dFF_win3_mean(:,ind_Low),2)-mean(dFF_base_mean(:,ind_Low),2);

neurons_response_stim = table(neurons_good, h_low_win1, h_low_win2, h_low_win3, h_high_win1, h_high_win2, h_high_win3);

%% Sound and choice encoding neurons (two-way ANOVA)
dFF_bh_mean = zeros(N,trials);
for i=1:trials
    dFF_bh_mean(:,i) = mean(X(:,sound_onsets(i)+round(0:2*fs)),2);
end

dFF_choice_low = zeros(N,1);
dFF_choice_high = zeros(N,1);
h_choice_low = false(N,1);
h_choice_high = false(N,1);

stimulus = regexp(groups, '(High)|(Low)|', 'match');
stimulus = categorical(cellfun(@(x) char(x), stimulus, 'UniformOutput', false));
choices = regexp(groups, '(Hit)|(Miss)|(Err)', 'match');
choices = categorical(cellfun(@(x) char(x), choices, 'UniformOutput', false));

stimulus_encoding = false(N,1);
choice_encoding = false(N,1);
joint_encoding = false(N,1);

ind_stim = cell(N,1); % Low, High
ind_choice = cell(N,1); % Hit, Miss, Error
ind_joint = cell(N,2); % Choice, Stimulus
alpha = 0.05/N;

is_testable = [sum(strcmp(groups, 'HitLow')), sum(strcmp(groups, 'HitHigh'));
            sum(strcmp(groups, 'MissLow')), sum(strcmp(groups, 'MissHigh'));
            sum(strcmp(groups, 'ErrLow')), sum(strcmp(groups, 'ErrHigh'))]>=5; % At least 5 trials to test choice-encoding

for n=1:N
    % Fit GLME model
    tbl = table(dFF_bh_mean(n,:)'-dFF_base_mean(n,:)', stimulus, choices, 'VariableNames', {'dFF', 'Stimulus', 'Choice'});
%     grp_stats = grpstats(tbl, {'Stimulus', 'Choice'}, 'numel');
    choice_order = {'Hit', 'Miss', 'Err'};
    max_idx_high = 0;
    max_idx_low = 0;

    if sum(is_testable(:,2))>=2
        tbl_high = tbl(tbl.Stimulus=='High',:);   
        % Mean dF/F of different groups
        dFF_group_means_high = zeros(3,1);
        for i=1:3
            dFF_group_means_high(i) = mean(tbl_high{tbl_high.Choice==choice_order{i},'dFF'});
        end
%             dFF_group_means_high(~is_testable(:,1)) = -Inf;
        [dFF_choice_high(n),max_idx_high] = max(dFF_group_means_high);
        
        if sum(is_testable(:,2))==3 % if all choices are testable, use ANOVA to test the effects of choices
            lme_high = fitlme(tbl_high, 'dFF ~ 1+Choice');
            stats_high = anova(lme_high, 'DFMethod', 'satterthwaite');
            if stats_high.pValue(2)<alpha
                [h_choice_high(n), p_choice_high(n)] = ttest(tbl_high{tbl_high.Choice==choice_order{max_idx_high},'dFF'}, 0, 'Alpha', .05, 'Tail','right');
            else
                max_idx_high = 0;
            end
        else    % if only two choices are testable, use two-sample t-test
            idx_choice_tmp = find(is_testable(:,2));
            [h_tmp, p_tmp] = ttest2(tbl_high{tbl_high.Choice==choice_order{idx_choice_tmp(1)},'dFF'},...
                tbl_high{tbl_high.Choice==choice_order{idx_choice_tmp(2)},'dFF'}, 'Alpha', alpha);
            if h_tmp && is_testable(max_idx_high,2)
                [h_choice_high(n), p_choice_high(n)] = ttest(tbl_high{tbl_high.Choice==choice_order{max_idx_high},'dFF'}, 0, 'Alpha', .05, 'Tail','right');
            else
                max_idx_high = 0;
            end
        end
    end
    
    if sum(is_testable(:,1))>=2
        tbl_low = tbl(tbl.Stimulus=='Low',:);   
        % Mean dF/F of different groups
        dFF_group_means_low = zeros(3,1);
        for i=1:3
            dFF_group_means_low(i) = mean(tbl_low{tbl_low.Choice==choice_order{i},'dFF'});
        end
%             dFF_group_means_low(~is_testable(:,1)) = -Inf;
        [dFF_choice_low(n),max_idx_low] = max(dFF_group_means_low);
        
        if sum(is_testable(:,1))==3 % if all choices are testable, use ANOVA to test the effects of choices
            lme_low = fitlme(tbl_low, 'dFF ~ 1+Choice');
            stats_low = anova(lme_low, 'DFMethod', 'satterthwaite');
            if stats_low.pValue(2)<alpha
                [h_choice_low(n), p_choice_low(n)] = ttest(tbl_low{tbl_low.Choice==choice_order{max_idx_low},'dFF'}, 0, 'Alpha', .05, 'Tail','right');
            else
                max_idx_low = 0;
            end
        else    % if only two choices are testable, use two-sample t-test
            idx_choice_tmp = find(is_testable(:,1));
            [h_tmp, p_tmp] = ttest2(tbl_low{tbl_low.Choice==choice_order{idx_choice_tmp(1)},'dFF'},...
                tbl_low{tbl_low.Choice==choice_order{idx_choice_tmp(2)},'dFF'}, 'Alpha', alpha);
            if h_tmp && is_testable(max_idx_low,1)
                [h_choice_low(n), p_choice_low(n)] = ttest(tbl_low{tbl_low.Choice==choice_order{max_idx_low},'dFF'}, 0, 'Alpha', .05, 'Tail','right');
            else
                max_idx_low = 0;
            end
        end
    end
    
    choice_encoding(n) = h_choice_high(n) && h_choice_low(n) && (max_idx_high==max_idx_low);
    joint_encoding(n) = ~choice_encoding(n) && (h_choice_low(n) || h_choice_high(n));
    stimulus_encoding(n) = any(neurons_response_stim{n,2:7}); % && ~joint_encoding(n) && ~choice_encoding(n);
    
    if choice_encoding(n)
        ind_choice{n} = max_idx_high;
    end
    if joint_encoding(n)
        if h_choice_low(n)
            ind_joint{n,1} = max_idx_low;
            ind_joint{n,2} = 1;
        elseif h_choice_high(n)
            ind_joint{n,1} = max_idx_high;
            ind_joint{n,2} = 2;
        end
    end
    if stimulus_encoding(n)
        [~,ind_stim{n}] = max([mean(dFF_bh_mean(n,ind_Low)), mean(dFF_bh_mean(n,ind_High))]);
    end
end

fprintf('Choice Encoding neurons\n')
disp(find(choice_encoding))
fprintf('Joint Choice Encoding neurons\n')
disp(find(joint_encoding))

%% Find reward-encoding neurons
% transient onsets aligned with reward time
alpha = 0.05/N;
ind_HitLow = find(strcmp(groups, 'HitLow'));
ind_HitHigh = find(strcmp(groups, 'HitHigh'));
h_reward_low = false(N,1);
h_reward_high = false(N,1);
for n=1:N
    transient_times = transientInfo(neurons_good(n)).times;
    transient_onsets = transient_times(:,1);
    [ind_onsets_trials, transient_onsets_trials] = match_transient_times_to_trials(transient_onsets, trials_onset);
    ind_reward_trials_low = intersect(ind_onsets_trials, ind_HitLow);  % Overlap between transient onsets trials and low reward trials
    prop_reward_trials_low(n) = length(ind_reward_trials_low)/length(ind_HitLow);  % Proportion of transient onsets trials in low reward trials
    ind_reward_trials_high = intersect(ind_onsets_trials, ind_HitHigh);  % Overlap between transient onsets trials and high reward trials
    prop_reward_trials_high(n) = length(ind_reward_trials_high)/length(ind_HitHigh);  % Proportion of transient onsets trials in high reward trials
    ind_sel_low = ismember(ind_onsets_trials, ind_reward_trials_low);
    ind_sel_high = ismember(ind_onsets_trials, ind_reward_trials_high);
    if sum(ind_sel_low)>=5
        [h_reward_low(n),p_reward_low(n)] = vartest2(transient_onsets_trials(ind_sel_low) - frames_reward_trials(ind_onsets_trials(ind_sel_low)),...
            transient_onsets_trials(ind_sel_low), 'Tail','left', 'Alpha',alpha);
    end
    if sum(ind_sel_high)>=5
        [h_reward_high(n),p_reward_high(n)] = vartest2(transient_onsets_trials(ind_sel_high) - frames_reward_trials(ind_onsets_trials(ind_sel_high)),...
            transient_onsets_trials(ind_sel_high), 'Tail','left', 'Alpha',alpha);
    end
end

fprintf('Low frequency Reward Encoding neurons\n')
disp(find(h_reward_low))
fprintf('High frequency Reward Encoding neurons\n')
disp(find(h_reward_high))

%% Motion encoding neurons
% Smooth dFF
X_smooth = filtfilt(ones(1,smooth_win)/smooth_win, 1, double(X)')';
rot_speed_smooth_ITI = zeros(trials-1, round(3*fs)); % Motion in ITI
rot_speed_smooth_pre_sound = zeros(trials-1, round(fs)); % Motion in 0-1s

for i=2:trials
    rot_speed_smooth_ITI(i-1, :) = rot_speed_smooth(trials_onset(i)+(-round(3*fs):-1));
    rot_speed_smooth_pre_sound(i-1, :) = rot_speed_smooth(trials_onset(i)+(0:round(fs)-1));
end

if showplot
    figure; boxplot([mean(abs(rot_speed_smooth_ITI),2), mean(abs(rot_speed_smooth_pre_sound),2)])
    set(gca, 'XTick', [1,2], 'XTickLabel', {'ITI', 'pre-stimulus'})
    ylabel('Mean Absolute Rotary Speed')
end

% Compare mean rotary speed in ITI and pre-sound period
[h_ITI,p_ITI] = ttest(mean(abs(rot_speed_smooth_ITI),2), mean(abs(rot_speed_smooth_pre_sound),2))

% Cross-correlation between motion and dFF during ITI, for different lags
pvalues_motion = zeros(N,2*lags+1);
xcfs_motion = zeros(N,2*lags+1);
pvalues_motion_abs = zeros(N,2*lags+1);
xcfs_motion_abs = zeros(N,2*lags+1);

load(fullfile(folder_dir, input_folder_name, file_beh), 'fs');
load(fullfile(folder_dir, input_folder_name, file_dFF), 'trials_onset');
parfor n=1:N
    fprintf('n=%d\n',n);
    dFF_smooth_ITI = zeros(trials-1, round(3*fs)); % dFF in ITI
    for i=2:trials
        dFF_smooth_ITI(i-1,:) = X_smooth(n, trials_onset(i)+(-round(3*fs):-1));
    end    
    tic
    [xcfs_motion(n,:),~,pvalues_motion(n,:)] = test_crosscorr_trials(dFF_smooth_ITI, rot_speed_smooth_ITI, 'NumLags', lags);
    [xcfs_motion_abs(n,:),~,pvalues_motion_abs(n,:)] = test_crosscorr_trials(dFF_smooth_ITI, abs(rot_speed_smooth_ITI), 'NumLags', lags);
    toc
end

hs_motion = pvalues_motion<=.05/(N*(2*lags+1));
hs_motion_abs = pvalues_motion_abs<=.05/(N*(2*lags+1));
motion_encoding = any(hs_motion,2);
abs_motion_encoding = any(hs_motion_abs,2);

%% Visualize neurons with max movement coefficients
[max_corrs_motion, neuron_inds_motion] = maxk(max(abs(xcfs_motion),[],2),3);
[max_corrs_motion_abs, neuron_inds_motion_abs] = maxk(max(abs(xcfs_motion_abs),[],2),3);
fprintf('Maximum correlation to movement=%f\n', max_corrs_motion(1));
fprintf('Maximum correlation to absolute movement=%f\n', max_corrs_motion_abs(1));

if showplot
    for i=1:3
        figure;
        plot((-lags:lags)/fs, xcfs_motion(neuron_inds_motion(i),:))
        hold on
        yline(0)
        xline(0)
        ylim([-0.5 0.5])
        xlabel('Time(s)')
    end
    
    for i=1:3
        figure;
        plot((-lags:lags)/fs, xcfs_motion_abs(neuron_inds_motion_abs(i),:))
        hold on
        yline(0)
        xline(0)
        ylim([-0.5 0.5])
        xlabel('Time(s)')
    end
end

%% Distributions of max coefficients to movement
[max_xcfs_motion, time_inds_motion] = max(abs(xcfs_motion),[],2);
[max_xcfs_motion_abs, time_inds_motion_abs] = max(abs(xcfs_motion_abs),[],2);
max_xcf_times_motion = zeros(1,N);
max_xcf_times_abs_motion = zeros(1,N);
times_all = (-lags:lags)/fs;
for n=1:N
    max_xcfs_motion(n) = xcfs_motion(n,time_inds_motion(n));
    max_xcfs_motion_abs(n) = xcfs_motion_abs(n,time_inds_motion_abs(n));
    max_xcf_times_motion(n) = times_all(time_inds_motion(n));
    max_xcf_times_abs_motion(n) = times_all(time_inds_motion_abs(n));
end

% Plot max correlations and their occurring times together
if showplot
    figure;
    ax1 = subplot(211);
    scatter(max_xcf_times_motion(any(hs_motion,2)), max_xcfs_motion(any(hs_motion,2)), 12, "blue", 'filled')
    hold on
    yline(0, 'Color', 'k')
    xlabel('Time of Max Cross Correlation to Motion')
    ylabel('Max Cross Correlation')
    ax2 = subplot(212);
    scatter(max_xcf_times_abs_motion(any(hs_motion_abs,2)), max_xcfs_motion_abs(any(hs_motion_abs,2)), 12, "blue", 'filled')
    hold on
    yline(0, 'Color', 'k')
    xlabel('Time of Max Cross Correlation to Absolute Motion')
    ylabel('Max Cross Correlation')
    linkaxes([ax1, ax2], 'xy')
end

%% Visualize complete single neuron encodings
% Concatenate stimulus encodings, choice encodings, reward encodings and movement encodings
neurons_encoding_all = [h_low_win1, h_low_win2, h_low_win3, h_high_win1, h_high_win2, h_high_win3,...
    h_choice_low, h_choice_high, any(hs_motion,2), any(hs_motion_abs,2)];
neurons_dFF_all = [dFF_low_win1,dFF_low_win2,dFF_low_win3,dFF_high_win1,dFF_high_win2,dFF_high_win3,...
    dFF_choice_low, dFF_choice_high,...
    max(abs(xcfs_motion),[],2), max(abs(xcfs_motion_abs),[],2)];
neurons_significant_dFF_all = abs(neurons_dFF_all).*neurons_encoding_all;

% Sort neurons
% Low, High frequencies
ind_stim_tmp = ind_stim;
ind_stim_tmp(cellfun('isempty', ind_stim_tmp)) = {3};
ind_stim_tmp = cell2mat(ind_stim_tmp);
ind_stim_tmp = 3-ind_stim_tmp;

[~,index_encodings] = sortrows([ind_stim_tmp, neurons_significant_dFF_all], [1,2,5,3,6,4,7,8:11], 'descend');

if showplot
    figure;
    imagesc(neurons_significant_dFF_all(index_encodings,:), [0,1])
    cmap = colormap('winter');
    cmap = flip(cmap);
    cmap(1,:) = [1,1,1];
    colormap(cmap);
    
    hold on
    xline([3.5, 6.5])
    xline(8.5, 'LineWidth',2)
    
    set(gca,'xaxisLocation','top', 'XTick', 1:10, 'XTickLabel', {'Low-win1','Low-win2','Low-win3','High-win1','High-win2','High-win3',...
        'Choice-Low', 'Choice-High', 'Motion', 'Abs Motion'})
    ylabel('Neurons')
    set(gcf, 'Position', [200,200,1200,900])
    colorbar
end

%% Save results
if Save
    save(fullfile(output_path, 'single_neuron_encoding_stats.mat'), 'neurons_response_stim', 'stimulus_encoding', ...
        'choice_encoding', 'joint_encoding', 'motion_encoding', 'abs_motion_encoding', 'h_reward_low', 'h_reward_high',...
        'h_choice_low', 'h_choice_high', 'xcfs_motion', 'xcfs_motion_abs', 'hs_motion', 'hs_motion_abs', 'h_ITI', 'p_ITI',...
        'neurons_significant_dFF_all', 'ind_stim', 'ind_choice', 'ind_joint', 'is_testable',...
        'dFF_high_win1', 'dFF_high_win2', 'dFF_high_win3', 'dFF_low_win1', 'dFF_low_win2', 'dFF_low_win3', 'index_encodings')
end
