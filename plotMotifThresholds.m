function fig_files = plotMotifThresholds(R, output_path)
% Plot W, H, and H_test motifs across calcium thresholds for one session.

if nargin < 2 || isempty(output_path)
    output_path = R.output_path;
end

fig_files.W = fullfile(output_path, sprintf('%s_motifs_W_thresh.pdf', R.input_folder_name));
fig_files.H = fullfile(output_path, sprintf('%s_motifs_H_thresh.pdf', R.input_folder_name));
fig_files.H_test = fullfile(output_path, sprintf('%s_motifs_Htest_thresh.pdf', R.input_folder_name));

plotWs(R.Ws_active_common, R.calcium_threshs, fig_files.W);
plotHs(R.Hs, R.calcium_threshs, fig_files.H, 'H');

if any(cellfun(@(x) ~isempty(x), R.H_tests))
    plotHs(R.H_tests, R.calcium_threshs, fig_files.H_test, 'H test');
else
    fig_files.H_test = '';
end
end

function plotWs(Ws_active_common, calcium_threshs, save_file)
colors = motifColors();
n_thresh = numel(Ws_active_common);

figure('Visible', 'off');
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])

for i = 1:n_thresh
    W_active_common = Ws_active_common{i};
    subplot(1, n_thresh, i); hold on
    set(gca, 'ColorOrder', colors);
    title(sprintf('thresh = %g', calcium_threshs(i)))

    if isempty(W_active_common) || size(W_active_common, 2) == 0
        axis off
        continue
    end

    [N, K, L] = size(W_active_common);
    sep = ceil(L * 0.1);
    dnW = prctile(W_active_common(:), 99.9) + eps;
    XsEdge = zeros(3, K);
    YsEdge = [zeros(1, K); dnW * N * ones(2, K)];

    for ki = 1:K
        x0 = (ki - 1) * (L + sep);
        Xs = [x0 + 1, x0 + (1:L), x0 + L];

        for ni = 1:N
            Ys = [dnW * (N - ni), dnW * (N - ni) + squeeze(W_active_common(ni, ki, :))', dnW * (N - ni)];
            patch(Xs, Ys, 'k', 'edgecolor', 'none')
        end

        XsEdge(:, ki) = [x0 + 1; x0 + 1; x0 + L];
    end

    plot(XsEdge, YsEdge, 'LineWidth', 2);
    xlim([1 K * (L + sep)]);
    ylim([0 dnW * N]);
    axis off
end

exportgraphics(gcf, save_file, 'ContentType', 'vector')
close(gcf)
end

function plotHs(Hs, calcium_threshs, save_file, plot_name)
base_colors = motifColors();
n_thresh = numel(Hs);

figure('Visible', 'off');
set(gcf, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])

for i = 1:n_thresh
    H = Hs{i};
    subplot(n_thresh, 1, i); hold on
    title(sprintf('%s, thresh = %g', plot_name, calcium_threshs(i)))

    if isempty(H)
        axis off
        continue
    end

    [K, T] = size(H);
    colors = base_colors(1 + mod(0:K-1, size(base_colors, 1)), :);
    dn = max(H(:)) / 2 + eps;

    for ki = K:-1:1
        k = K - ki + 1;
        Xs = [1, 1:T, T];
        Ys = [dn * ki, dn * ki + H(k, :), dn * ki] - dn / 2;
        patch(Xs, Ys, colors(k, :), 'edgecolor', colors(k, :))
    end

    ylim([0 dn * K + dn * 3 + eps]);
    xlim([0 T + 1]);
    axis off
end

exportgraphics(gcf, save_file, 'ContentType', 'vector')
close(gcf)
end

function colors = motifColors()
colors = [[0 .6 .3]; [.7 0 .7]; [1 .6 0]; [.1 .3 .9]; [1 .1 .1]; ...
    [0 .9 .3]; [.4 .2 .7]; [.7 .2 .1]; [.1 .8 1]; [1 .3 .7]; ...
    [.2 .8 .2]; [.7 .4 1]; [.9 .6 .4]; [0 .6 1]; [1 .1 .3]];
end
