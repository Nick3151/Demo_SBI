function R = compareMotifsByCalciumThresh(input_folder_name, calcium_threshs, root, min_occur_times)
% Compare motifs from one session across calcium thresholds.
%
% This function uses calcium threshold 2 as the reference.
% Threshold 2 loads FlexMF_Motifs.mat.
% Other thresholds load FlexMF_Motifs_thresh=%1.2f.mat.
%
% Usage:
%   R = compareMotifsByCalciumThresh(input_folder_names{n});
%   R = compareMotifsByCalciumThresh(input_folder_names{n}, [1.5 2 2.5]);

if nargin < 2 || isempty(calcium_threshs)
    calcium_threshs = [1.5 2 2.5];
end
if nargin < 3 || isempty(root)
    root = pwd;
end
if nargin < 4 || isempty(min_occur_times)
    min_occur_times = 4;
end

calcium_threshs = calcium_threshs(:)';
ref_idx = find(abs(calcium_threshs - 2) < 1e-10, 1);
if isempty(ref_idx)
    error('calcium_threshs must include 2, because threshold 2 is the reference.');
end

output_path = fullfile(root, 'Data', input_folder_name);
n_thresh = numel(calcium_threshs);

Ws = cell(1, n_thresh);
Hs = cell(1, n_thresh);
H_tests = cell(1, n_thresh);
neurons_active_ids = cell(1, n_thresh);
motif_ids = cell(1, n_thresh);
motif_freqs = cell(1, n_thresh);
loadings_train = cell(1, n_thresh);
loadings_test = cell(1, n_thresh);
freq_tables = cell(1, n_thresh);
freq_files = cell(1, n_thresh);

for t = 1:n_thresh
    calcium_thresh = calcium_threshs(t);
    S = load(getMotifFile(output_path, calcium_thresh));
    S_beh = load(fullfile(output_path, 'bh_data_Liang.mat'));

    % plot, sorting neurons by latency within each factor
    frames = round(5*S_beh.fs);
    [max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(S.W(:,:,:),1);
    indSort = hybrid(:,3);
    plotAll = 1;
    figure; 
    SimpleWHPlot_patch(S.W(indSort,:,:), S.H, 'trials', length(S.ind_train), 'frames', frames, 'Data', S.X_train(indSort,:), 'is_significant', S.is_significant, 'plotAll', plotAll);
    title('SBI motifs train, with raw data')
    set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
    exportgraphics(gcf, fullfile(output_path, sprintf('%s_motifs_visualization_train_thresh=%1.2f.pdf', input_folder_name, calcium_thresh)), 'ContentType', 'vector')
    figure; 
    SimpleWHPlot_patch(S.W(indSort,:,:), S.H, 'trials', length(S.ind_train), 'frames', frames, 'Data', S.X_train(indSort,:), 'is_significant', S.is_significant, 'plotAll', plotAll);
    title('SBI motifs test, with raw data')
    set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
    exportgraphics(gcf, fullfile(output_path, sprintf('%s_motifs_visualization_test_thresh=%1.2f.pdf', input_folder_name, calcium_thresh)), 'ContentType', 'vector')
    close all
    [Ws{t}, Hs{t}, H_tests{t}, neurons_active_ids{t}, ...
        motif_ids{t}, motif_freqs{t}, freq_tables{t}, loadings_train{t}, loadings_test{t}] = ...
        selectMotifs(S, input_folder_name, calcium_thresh, min_occur_times);

    freq_files{t} = fullfile(output_path, sprintf('%s_motif_freq_thresh_%s.csv', ...
        input_folder_name, threshLabel(calcium_thresh)));
    writetable(freq_tables{t}, freq_files{t});
end

motif_freq_table = vertcat(freq_tables{:});
all_freq_file = fullfile(output_path, sprintf('%s_motif_freq_all_thresh.csv', input_folder_name));
writetable(motif_freq_table, all_freq_file);

[Ws_active_common, neurons_active_common] = alignWsToCommonNeurons(Ws, neurons_active_ids);
similarity_freq_table = makeSimilarityFreqTable(input_folder_name, calcium_threshs, ...
    ref_idx, Ws_active_common, Hs, H_tests, motif_ids, motif_freqs, loadings_train, loadings_test);

similarity_freq_file = fullfile(output_path, ...
    sprintf('%s_similarity_vs_motif_freq_thresh2.csv', input_folder_name));
writetable(similarity_freq_table, similarity_freq_file);

R = struct();
R.input_folder_name = input_folder_name;
R.output_path = output_path;
R.calcium_threshs = calcium_threshs;
R.ref_thresh = 2;
R.min_occur_times = min_occur_times;
R.Ws = Ws;
R.Hs = Hs;
R.H_tests = H_tests;
R.neurons_active_ids = neurons_active_ids;
R.neurons_active_common = neurons_active_common;
R.Ws_active_common = Ws_active_common;
R.motif_ids = motif_ids;
R.motif_freqs = motif_freqs;
R.motif_freq_table = motif_freq_table;
R.similarity_freq_table = similarity_freq_table;
R.freq_files = freq_files;
R.all_freq_file = all_freq_file;
R.similarity_freq_file = similarity_freq_file;
R.plot_files = plotMotifThresholds(R, output_path);
R.save_file = fullfile(output_path, sprintf('%s_motif_compare_by_calcium_thresh.mat', input_folder_name));

save(R.save_file, 'R');
end

function motif_file = getMotifFile(output_path, calcium_thresh)
if abs(calcium_thresh - 2) < 1e-10
    motif_file = fullfile(output_path, 'FlexMF_Motifs.mat');
else
    motif_file = fullfile(output_path, sprintf('FlexMF_Motifs_thresh=%1.2f.mat', calcium_thresh));
end

if ~isfile(motif_file)
    error('Motif file not found: %s', motif_file);
end
end

function [W_active, H_sel, H_test_sel, neurons_active_ids, motif_ids_sel, motif_freq_sel,...
    freq_table, loadings_train_sel, loadings_test_sel] = ...
    selectMotifs(S, input_folder_name, calcium_thresh, min_occur_times)

need_vars = {'W', 'H', 'active_time_motifs_all', 'neurons_active_motifs', ...
    'neurons_good', 'is_significant', 'loadings_train', 'loadings_test'};
for i = 1:numel(need_vars)
    if ~isfield(S, need_vars{i})
        error('Missing variable in motif file: %s', need_vars{i});
    end
end

motif_freq = cellfun(@numel, S.active_time_motifs_all(:));
num_neurons_active = cellfun(@numel, S.neurons_active_motifs(:));
is_sig = logical(S.is_significant(:));
motif_id = (1:numel(motif_freq))';

selected = motif_freq >= min_occur_times & is_sig & num_neurons_active > 1;
k_sel = find(selected);

freq_table = table(repmat({input_folder_name}, numel(motif_id), 1), ...
    repmat(calcium_thresh, numel(motif_id), 1), motif_id, motif_freq, ...
    num_neurons_active, is_sig, selected, ...
    'VariableNames', {'session', 'calcium_thresh', 'motif_id', ...
    'motif_occurrence_freq', 'num_active_neurons', 'is_significant', 'selected'});

W_active = cell(1, numel(k_sel));
neurons_active_ids = cell(1, numel(k_sel));

for i = 1:numel(k_sel)
    k = k_sel(i);
    active_idx = S.neurons_active_motifs{k};
    W_active{i} = reshape(S.W(active_idx, k, :), numel(active_idx), size(S.W, 3));
    neurons_active_ids{i} = S.neurons_good(active_idx);
end

H_sel = S.H(k_sel, :);
if isfield(S, 'H_test')
    H_test_sel = S.H_test(k_sel, :);
else
    H_test_sel = [];
end

motif_ids_sel = k_sel(:);
motif_freq_sel = motif_freq(k_sel);
loadings_train_sel = S.loadings_train(k_sel);
loadings_test_sel = S.loadings_test(k_sel);
end

function [Ws_active_common, neurons_active_common] = alignWsToCommonNeurons(Ws, neurons_active_ids)
n_thresh = numel(Ws);
all_neurons = [];

for t = 1:n_thresh
    for k = 1:numel(neurons_active_ids{t})
        all_neurons = [all_neurons; neurons_active_ids{t}{k}(:)]; %#ok<AGROW>
    end
end

neurons_active_common = unique(all_neurons);
Ws_active_common = cell(1, n_thresh);

for t = 1:n_thresh
    K = numel(Ws{t});
    if K == 0
        Ws_active_common{t} = zeros(numel(neurons_active_common), 0, 0);
        continue
    end

    L = size(Ws{t}{1}, 2);
    W_common = zeros(numel(neurons_active_common), K, L);

    for k = 1:K
        motif_neurons = neurons_active_ids{t}{k}(:);
        [~, loc] = ismember(motif_neurons, neurons_active_common);
        Wk = reshape(Ws{t}{k}, numel(motif_neurons), L);
        W_common(loc, k, :) = reshape(Wk, numel(motif_neurons), 1, L);
    end

    Ws_active_common{t} = W_common;
end
end

function T = makeSimilarityFreqTable(input_folder_name, calcium_threshs, ref_idx, ...
    Ws_active_common, Hs, H_tests, motif_ids, motif_freqs, loadings_train, loadings_test)

parts = {};
W_ref = Ws_active_common{ref_idx};
H_ref = Hs{ref_idx};
H_test_ref = H_tests{ref_idx};

if isempty(W_ref) || size(W_ref, 2) == 0 || isempty(H_ref)
    T = emptySimilarityTable();
    return
end

for t = 1:numel(calcium_threshs)
    if t == ref_idx
        continue
    end

    W_other = Ws_active_common{t};
    H_other = Hs{t};
    H_test_other = H_tests{t};
    if isempty(W_other) || size(W_other, 2) == 0 || isempty(H_other)
        continue
    end

    [coeff_W, coeff_H, ids] = helper.similarity_WH(W_ref, H_ref, W_other, H_other);
    [~, coeff_H_test, ~] = helper.similarity_WH(W_ref, H_test_ref, W_other, H_test_other);
    coeff_W = coeff_W(:);
    coeff_H = coeff_H(:);
    coeff_H_test = coeff_H_test(:);
    ids = ids(:);

    n = numel(coeff_W);
    motif_id = motif_ids{t}(:);
    matched_motif_id_thresh2 = nan(n, 1);
    motif_frequency_thresh2 = nan(n, 1);
    loadings_train_thresh2 = nan(n, 1);
    loadings_test_thresh2 = nan(n, 1);

    valid = ids > 0 & ids <= numel(motif_ids{ref_idx});
    matched_motif_id_thresh2(valid) = motif_ids{ref_idx}(ids(valid));
    motif_frequency_thresh2(valid) = motif_freqs{ref_idx}(ids(valid));
    loadings_train_thresh2(valid) = loadings_train{ref_idx}(ids(valid));
    loadings_test_thresh2(valid) = loadings_test{ref_idx}(ids(valid));

    parts{end+1} = table(repmat({input_folder_name}, n, 1), ... %#ok<AGROW>
        repmat(calcium_threshs(t), n, 1), motif_id, matched_motif_id_thresh2, ...
        motif_frequency_thresh2, coeff_W, coeff_H, coeff_H_test, loadings_train_thresh2, loadings_test_thresh2,...
        'VariableNames', {'session', 'calcium_thresh', 'motif_id', ...
        'matched_motif_id_thresh2', 'motif_frequency_thresh2', 'coeff_W', 'coeff_H', 'coeff_H_test',...
        'loadings_train_thresh2', 'loadings_test_thresh2'});
end

if isempty(parts)
    T = emptySimilarityTable();
else
    T = vertcat(parts{:});
end
end

function T = emptySimilarityTable()
T = table(cell(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    'VariableNames', {'session', 'calcium_thresh', 'motif_id', ...
    'matched_motif_id_thresh2', 'motif_frequency_thresh2', 'coeff_W', 'coeff_H'});
end

function label = threshLabel(calcium_thresh)
label = strrep(sprintf('%1.2f', calcium_thresh), '.', 'p');
end
