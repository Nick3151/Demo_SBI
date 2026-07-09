function [stim_encoding, choice_encoding, joint_encoding, stim_ind, choice_ind, joint_ind] = compute_motifs_encoding_index(motifs_trials_match)
% COMPUTE_MOTIFS_ENCODING_INDEX  Per-motif selectivity for stimulus/choice.
%
%   [stim_encoding, choice_encoding, joint_encoding, ...
%    stim_ind, choice_ind, joint_ind] = compute_motifs_encoding_index(motifs_trials_match)
%
%   For every motif, compares how often the motif occurs across the six trial
%   types and derives conditional-probability-based encoding scores for the
%   stimulus (Low/High), the animal's choice (Hit/Miss/Error), and their joint
%   combination. Unlike the Fisher-test version, this returns graded selectivity
%   scores rather than significance decisions.
%
%   INPUT
%     motifs_trials_match - table with (K+1) rows and 6 variables ordered
%         {HitLow, MissLow, ErrLow, HitHigh, MissHigh, ErrHigh}. Each cell
%         holds the list of trials of that type. Row 1 is the reference row
%         (all trials); rows 2..K+1 give the trials in which each of the K
%         motifs occurred.
%
%   For each motif a 3x2x2 count array is built as
%     trials_nums(choice, stim, occurrence)
%   where choice = {Hit, Miss, Error}, stim = {Low, High}, and occurrence =
%   {motif present, motif absent}. Conditional probabilities
%   P(motif present | category) are computed by marginalizing over the other
%   factors, and the encoding scores are the gap between the top two.
%
%   OUTPUTS (all length K, one entry per motif)
%     stim_encoding   - selectivity gap (top minus second) of the per-stimulus
%                       conditional probabilities.
%     choice_encoding - selectivity gap of the per-choice conditional probs.
%     joint_encoding  - selectivity gap over the joint choice-by-stimulus cells.
%     stim_ind        - index (1=Low, 2=High) of the preferred stimulus.
%     choice_ind      - index (1=Hit, 2=Miss, 3=Error) of the preferred choice.
%     joint_ind       - [K x 2], the [choice, stimulus] indices of the
%                       preferred joint category.
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