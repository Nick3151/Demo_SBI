function fig = plotSimilarityVsMotifFreq(T, save_file)
% Plot W and H similarity against motif frequency at calcium threshold 2.
%
% Usage:
%   plotSimilarityVsMotifFreq(similarity_freq_table)
%   plotSimilarityVsMotifFreq(similarity_freq_table, 'similarity_vs_freq.pdf')

fig = figure;
set(fig, 'Units', 'normalized', 'Position', [0.1 0.1 0.9 0.45]);

subplot(1, 3, 1); hold on
plotCoeff(T, 'coeff_W');
xlabel('Motif occurrence frequency at threshold 2');
ylabel('W similarity');
title('W')

subplot(1, 3, 2); hold on
plotCoeff(T, 'coeff_H');
xlabel('Motif occurrence frequency at threshold 2');
ylabel('H similarity');
title('H')

subplot(1, 3, 3); hold on
plotCoeff(T, 'coeff_H_test');
xlabel('Motif occurrence frequency at threshold 2');
ylabel('Htest similarity');
title('Htest')

if nargin > 1 && ~isempty(save_file)
    exportgraphics(fig, save_file, 'ContentType', 'vector');
end
end

function plotCoeff(T, coeff_name)
if isempty(T)
    ylim([0 1]);
    box off
    return
end

thresholds = unique(T.calcium_thresh(:))';
legend_text = cell(1, numel(thresholds));

for i = 1:numel(thresholds)
    idx = T.calcium_thresh == thresholds(i);
    scatter(T.motif_frequency_thresh2(idx), T.(coeff_name)(idx), 35, 'filled');
    legend_text{i} = sprintf('thresh = %g', thresholds(i));
end

ylim([0 1]);
box off
if ~isempty(legend_text)
    legend(legend_text, 'Location', 'best');
end
end
