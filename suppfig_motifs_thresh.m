%% Supplement figure, impact of calcium threshold on motif detection
clear all
close all
clc
input_folder_names = {'20240618MH041724', '20240602FH040624', ...
    '20240326ML102623', '20250127ML111624', '20250114FL111424'};

nSessions = length(input_folder_names);
tables = cell(nSessions, 1);

for n = 1:nSessions
    input_folder_name = input_folder_names{n};
    disp(['Input Folder Name: ', input_folder_name]);

    R = compareMotifsByCalciumThresh(input_folder_name);
    tables{n} = R.similarity_freq_table;
end

similarity_freq_table = vertcat(tables{:});
plotSimilarityVsMotifFreq(similarity_freq_table)

exportgraphics(gcf, 'Motif_similarity_vs_thresh.pdf', "ContentType","vector")

% Motif similarity for one session
thresholds = unique(similarity_freq_table.calcium_thresh(:))';
for i = 1:numel(thresholds)
    idx = similarity_freq_table.calcium_thresh == thresholds(i) & string(similarity_freq_table.session) == '20240618MH041724';
    coeffs_W = similarity_freq_table.coeff_W(idx);
    coeffs_H = similarity_freq_table.coeff_H(idx);
    coeffs_Htest = similarity_freq_table.coeff_H_test(idx);
    fprintf('Mean Coeffs of W=%0.3f\n', mean(coeffs_W));
    fprintf('Mean Coeffs of H=%0.3f\n', mean(coeffs_H));
    fprintf('Mean Coeffs of Htest=%0.3f\n', mean(coeffs_Htest));
end

% Matching rate
fprintf('Mathicng rate: %0.3f\n', sum(~isnan(similarity_freq_table.matched_motif_id_thresh2))/size(similarity_freq_table,1));

% [S_shape, S_activation, ids, details] = helper.similarity_WH_shape_activation(Ws_active_common{3}, Hs{3}, Ws_active_common{2}, Hs{2});
% [S_shape, S_activation, ids, details] = helper.similarity_WH_shape_activation(Ws_active_common{3}, H_tests{3}, Ws_active_common{2}, H_tests{2});