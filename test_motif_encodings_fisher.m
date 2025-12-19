function [h_stim, h_choice, h_joint, ind_stim, ind_choice, ind_joint, p_choice, min_p, is_testable] = test_motif_encodings_fisher(motifs_trials_match, k_sel)
% Fisher's exact test
trial_nums_raw = cellfun(@(x) length(x), motifs_trials_match.Variables);
K = length(trial_nums_raw)-1;

h_stim = zeros(K,1);
h_choice = zeros(K,1);
h_joint = zeros(K,1);

ind_stim = cell(K,1);
ind_choice = cell(K,1);
ind_joint = cell(K,2); % Choice, Stimulus
p_choice = ones(K,2); % p values for choice encoding

n_stim = 2;
n_choice = 3;
min_p = ones(n_choice, n_stim);

for i=1:length(k_sel)
    k = k_sel(i);
    % Hit trials
    trials_nums(1,:,1) = trial_nums_raw(1+k,[1,4]); 
    trials_nums(1,:,2) = (trial_nums_raw(1,[1,4])-trial_nums_raw(1+k,[1,4]));
    % Miss trials
    trials_nums(2,:,1) = trial_nums_raw(1+k,[2,5]); 
    trials_nums(2,:,2) = (trial_nums_raw(1,[2,5])-trial_nums_raw(1+k,[2,5]));
    % Error trials
    trials_nums(3,:,1) = trial_nums_raw(1+k,[3,6]); 
    trials_nums(3,:,2) = (trial_nums_raw(1,[3,6])-trial_nums_raw(1+k,[3,6]));
    
    % test stimulus encoding
    tbl_stim = squeeze(sum(trials_nums,1));
    h_tmp = zeros(1,n_stim);
    for s=1:n_stim
        tbl_stim_test = zeros(2,2);
        tbl_stim_test(1,:) = tbl_stim(s,:);
        tbl_stim_test(2,:) = sum(tbl_stim)-tbl_stim(s,:);
        h_tmp(s)= fishertest(tbl_stim_test, 'Tail','right', 'Alpha',.05/(n_stim*length(k_sel)));
    end
    h_stim(k) = any(h_tmp);
    if h_stim(k)
        ind_stim{k} = find(h_tmp);
    end
    
%     % test choice encoding
%     tbl_choice = squeeze(sum(trials_nums,2));
%     tbl_choice_test = zeros(2,2);
%     tbl_choice_test(1,:) = tbl_choice(ind_choice(k),:);
%     tbl_choice_test(2,:) = sum(tbl_choice,1)-tbl_choice_test(1,:);
%     [h_choice(k),p_choice(k),~] = fishertest(tbl_choice_test, 'Tail','right', 'Alpha',.05/3/length(k_sel));
%     
%     % test joint encoding
%     tbl_joint = zeros(2,2);
%     tbl_joint(1,:) = trials_nums(ind_joint(k,1),ind_joint(k,2),:);
%     tbl_joint(2,:) = squeeze(sum(trials_nums,[1,2]))'-tbl_joint(1,:);
%     [h_joint(k),p_joint(k),~] = fishertest(tbl_joint, 'Tail','right', 'Alpha',.05/3/length(k_sel));

    % choice encoding conditioned on stimulus
    h_cond = zeros(n_choice, n_stim);
    p_cond = ones(n_choice, n_stim);
    for s=1:n_stim
        tbl_cond = squeeze(trials_nums(:,s,:));
        for c=1:n_choice
            tbl_cond_test = zeros(2,2);
            tbl_cond_test(1,:) = tbl_cond(c,:);
            tbl_cond_test(2,:) = sum(tbl_cond)-tbl_cond(c,:);
            [h_cond(c,s), p_cond(c,s),~] = fishertest(tbl_cond_test, 'Tail','right', 'Alpha',.05/(length(k_sel)*n_stim*n_choice));
            min_p(c,s) = 1/nchoosek(sum(tbl_cond_test, 'all'), sum(tbl_cond_test(1,:)));
        end
    end

    h_choice(k) = any(all(h_cond,2));
    if h_choice(k)
        ind_choice{k} = find(all(h_cond,2));
        p_choice(k,:) = p_cond(ind_choice{k},:);
    end
    h_joint(k) = any(any(h_cond,2) & ~all(h_cond,2));
    if h_joint(k)
        ind_joint{k,1} = find(any(h_cond,2) & ~all(h_cond,2));
        ind_joint{k,2} = find(h_cond(ind_joint{k,1},:));
    end

    if h_choice(k) || h_joint(k)
        h_stim(k) = 0;
        ind_stim{k} = 0;
    end
end

if ~isempty(k_sel)
    is_testable = min_p<.05/(length(k_sel)*n_stim*n_choice);
else
    is_testable = zeros(n_choice, n_stim);
end