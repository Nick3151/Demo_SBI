function f1 = visualize_motif_traces(Wk, X, Xhat, fs, motif_reliability, sorted)
if nargin<6
    sorted = true;
end

[N,L] = size(Wk);

neurons_active = motif_reliability.neurons_active;
neurons_active_sorted = motif_reliability.neurons_active_sorted;
active_times = motif_reliability.active_times;
onsets_motif_neurons = motif_reliability.onsets_motif_neurons;
onsets_motif_neurons_sorted = motif_reliability.onsets_motif_neurons_sorted;
onsets_traces = motif_reliability.onsets_traces;
onsets_traces_sorted = motif_reliability.onsets_traces_sorted;

R = length(active_times);
X_tmp = zeros(N, R, L);
Xhat_tmp = zeros(N, R, L);
for i=1:R
    X_tmp(:,i,:) = X(:,active_times(i)+(1:L));
    Xhat_tmp(:,i,:) = Xhat(:,active_times(i)+(1:L));
end
X_tmp_mean = squeeze(mean(X_tmp, 2));
X_tmp_std = squeeze(std(X_tmp, [], 2));

Wk_active = Wk(neurons_active, :);
Wk_active_sorted = Wk(neurons_active_sorted, :);
Xhat_tmp_active = Xhat_tmp(neurons_active,:,:);
Xhat_tmp_active_sorted = Xhat_tmp(neurons_active_sorted,:,:);
X_mean_active = X_tmp_mean(neurons_active, :);
X_mean_active_sorted = X_tmp_mean(neurons_active_sorted, :);
X_std_active = X_tmp_std(neurons_active, :);
X_std_active_sorted = X_tmp_std(neurons_active_sorted, :);
X_tmp_active = X_tmp(neurons_active,:,:);
X_tmp_active_sorted = X_tmp(neurons_active_sorted,:,:);
N_active = length(neurons_active);

f1 = figure; 
set(gcf,'position',[200,200,900,800])
dnW = prctile(Wk_active(:),99.9); 
XsEdge = [0, 0, L]/fs; 
YsEdge = [0, dnW*N_active, dnW*N_active];

subplot(141)
for ni =1:N_active
    Xs = [1, 1:L, L]/fs;
    if sorted
        Ys = [dnW*(N_active-ni) dnW*(N_active-ni) + Wk_active_sorted(ni, :) dnW*(N_active-ni)];
    else
        Ys = [dnW*(N_active-ni) dnW*(N_active-ni) + Wk_active(ni, :) dnW*(N_active-ni)];
    end
    patch(Xs,Ys, 'k', 'edgecolor', 'none')
    hold on
end
xline(min(onsets_motif_neurons)/fs, 'r', 'LineWidth',1)

plot(XsEdge,YsEdge, 'LineWidth', 2);
xlim([0 L]/fs);ylim([0 dnW*N_active])
set(gca, 'YTick', [])

subplot(142)
hold on
dnX = prctile(Xhat_tmp_active(:),99.9); 
for ni =1:N_active
    Xs = (1:L)/fs;
    for i=1:R
        if sorted
            Ys = dnW*(N_active-ni) + squeeze(Xhat_tmp_active_sorted(ni, i, :))/dnX*dnW;
        else
            Ys = dnW*(N_active-ni) + squeeze(Xhat_tmp_active(ni, i, :))/dnX*dnW;
        end
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
    if sorted
        Ys = dnW*(N_active-ni) + X_mean_active_sorted(ni, :)/dnX*dnW;
        Ys_low = Ys - X_std_active_sorted(ni, :)/dnX*dnW;
        Ys_high = Ys + X_std_active_sorted(ni, :)/dnX*dnW;
    else
        Ys = dnW*(N_active-ni) + X_mean_active(ni, :)/dnX*dnW;
        Ys_low = Ys - X_std_active(ni, :)/dnX*dnW;
        Ys_high = Ys + X_std_active(ni, :)/dnX*dnW;
    end
    
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
    for i=1:R
        if sorted
            Ys = dnW*(N_active-ni) + squeeze(X_tmp_active_sorted(ni, i, :))/dnX*dnW;
        else
            Ys = dnW*(N_active-ni) + squeeze(X_tmp_active(ni, i, :))/dnX*dnW;
        end
        plot(Xs, Ys, 'Color', [0,0,0,.3], 'LineWidth', 1)
    end

end
plot(XsEdge,YsEdge, 'LineWidth', 2);
xlim([0 L]/fs);ylim([0 dnW*N_active])
set(gca, 'YTick', [])
