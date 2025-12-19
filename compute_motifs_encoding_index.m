function [stim_encoding, choice_encoding, joint_encoding, stim_ind, choice_ind, joint_ind] = compute_motifs_encoding_index(motifs_trials_match)
K = height(motifs_trials_match) - 1;
stim_encoding = zeros(K,1);
choice_encoding = zeros(K,1);
joint_encoding = zeros(K,1);
stim_ind = zeros(K,1);
choice_ind = zeros(K,1);
joint_ind = zeros(K,2);

%% Separate trials based on correct/error
for k=1:K
    trials_nums = zeros(3,2,2);
    trial_nums_raw = cellfun(@(x) length(x), motifs_trials_match.Variables);
    % Hit trials
    trials_nums(1,:,1) = trial_nums_raw(1+k,[1,4]); 
    trials_nums(1,:,2) = (trial_nums_raw(1,[1,4])-trial_nums_raw(1+k,[1,4]));
    % Miss trials
    trials_nums(2,:,1) = trial_nums_raw(1+k,[2,5]); 
    trials_nums(2,:,2) = (trial_nums_raw(1,[2,5])-trial_nums_raw(1+k,[2,5]));
    % Error trials
    trials_nums(3,:,1) = trial_nums_raw(1+k,[3,6]); 
    trials_nums(3,:,2) = (trial_nums_raw(1,[3,6])-trial_nums_raw(1+k,[3,6]));
    pCond_stim = sum(trials_nums(:,:,1),1)./sum(trials_nums,[1,3]);
    [~,stim_ind(k)] = max(pCond_stim);
    pCond_choice = sum(trials_nums(:,:,1),2)./sum(trials_nums,[2,3]);
    pCond_choice(isnan(pCond_choice)) = 0;
    [~,choice_ind(k)] = max(pCond_choice);
    pCond_joint = trials_nums(:,:,1)./sum(trials_nums,3);
    pCond_joint(isnan(pCond_joint)) = 0;
    [~,ind] = max(pCond_joint, [],'all','linear');
    [r,c] = ind2sub([3,2], ind);
    joint_ind(k,:) = [r,c];
    pCond_stim_sorted = sort(pCond_stim, 'descend');
    pCond_choice_sorted = sort(pCond_choice, 'descend');
    pCond_joint_sorted = sort(pCond_joint(:), 'descend');
    stim_encoding(k) = pCond_stim_sorted(1)-pCond_stim_sorted(2);
    choice_encoding(k) = pCond_choice_sorted(1)-pCond_choice_sorted(2);
    joint_encoding(k) = pCond_joint_sorted(1)-pCond_joint_sorted(2);
end