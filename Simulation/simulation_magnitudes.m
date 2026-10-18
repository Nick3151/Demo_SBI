% Test the impact of moitf magnitudes on synthetic data
clear all
close all
clc
root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
addpath(genpath(fullfile(root, 'FlexMF')));

%%
Trials = 200;
L = 50; % length of each trial
K = 10;
Nmotifs = 10*ones(K, 1);
Nneurons = 5*ones(K, 1); % the number of neurons in each motif
Dt = 3.*ones(K,1); % gap between each member of the motif

Magnitudes = 1:K; % Motif magnitudes

nSim = 10;
seeds = randperm(1000, nSim);
pvals_SeqNMF = cell(nSim,1);
is_significants_SeqNMF = cell(nSim,1);
loadings_SeqNMF = cell(nSim,1);
W_hats_SeqNMF = cell(nSim,1);
H_hats_SeqNMF = cell(nSim,1);
pvals_SBI = pvals_SeqNMF;
is_significants_SBI = is_significants_SeqNMF;
loadings_SBI = loadings_SeqNMF;
W_hats_SBI = W_hats_SeqNMF;
H_hats_SBI = H_hats_SeqNMF;
Ws = W_hats_SeqNMF;
Hs = H_hats_SeqNMF;
TrainingDatas = W_hats_SeqNMF;

parfor i=1:nSim
    tic
    [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
    'Magnitudes', Magnitudes, 'seed', seeds(i));
    groups = zeros(Trials,1);
    for k=1:K
        groups(motif_ind{k}) = k;
    end
    Ws{i} = W;
    
    % Dimension N*L*Trials
    cv = cvpartition(groups, "KFold",2);
    ind_train = find(training(cv,1));
    X_train = X(:,:,ind_train);
    ind_test = find(test(cv,1));
    X_test = X(:,:,ind_test);
    
    % Dimension N*T
    N = size(W,1);
    TrainingData = zeros(N,cv.TrainSize(1)*L);
    for t=1:cv.TrainSize(1)
        TrainingData(:,(t-1)*L+1:t*L) = squeeze(X_train(:,:,t));
    end
    TestData = zeros(N,cv.TestSize(1)*L);
    for t=1:cv.TestSize(1)
        TestData(:,(t-1)*L+1:t*L) = squeeze(X_test(:,:,t));
    end
    
    % Normalize training data
    frob_norm = norm(TrainingData(:));
    TrainingData = TrainingData/frob_norm*K;

    H_train = H(:,ind_train);
    T = cv.TrainSize(1)*L;
    H_train_full = zeros(K,T);
    H_train_full(:,1:L:T) = H_train;
    Hs{i} = H_train_full;

    lambda = .005;

    % Run SeqNMF with multiplication rule
    [W_hat, H_hat, ~,~,loadings_SeqNMF{i},power]= seqNMF(TrainingData,'K',K,'L',L,...
            'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
    p = .05;
    [pvals_SeqNMF{i},is_significants_SeqNMF{i}] = test_significance_new(TestData, W_hat,0,p);
    W_hats_SeqNMF{i} = W_hat;
    H_hats_SeqNMF{i} = H_hat;
    display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])

    % Run SeqNMF with Bregman Iteration
    alpha_W = 1e-6;
    alpha_H = 1e-2;
    [W_hat, H_hat, ~,~,loadings_SBI{i},power]= FlexMF(TrainingData,'K',K,'L',L,...
            'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
    p = .05;
    [pvals_SBI{i},is_significants_SBI{i}] = test_significance_new(TestData, W_hat,0,p);
    W_hats_SBI{i} = W_hat;
    H_hats_SBI{i} = H_hat;
    TrainingDatas{i} = TrainingData;
    display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
    toc
end

save('simulate_results_magnitudes.mat', "Magnitudes", "H_hats_SBI", 'W_hats_SBI', 'is_significants_SBI', 'pvals_SBI', 'loadings_SBI',...
    'loadings_SeqNMF', 'W_hats_SeqNMF', 'H_hats_SeqNMF', "pvals_SeqNMF", "is_significants_SeqNMF", 'Ws', 'Hs', 'TrainingDatas')

%% Visualization
i = 1;
TrainingData = TrainingDatas{i};
L = 50;
[N,Ttrain] = size(TrainingData);
TrainingTrials = Ttrain/L;
W = Ws{i};
H = Hs{i};

% Split Bregman
is_significant = is_significants_SBI{i};
What = W_hats_SBI{i};
Hhat = H_hats_SBI{i};

plotAll = 1;
figure; SimpleWHPlot_patch(What, Hhat, 'trials', TrainingTrials, 'frames', L, 'is_significant', is_significant, 'plotAll', plotAll); title('SBI reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot_patch(What, Hhat, 'trials', TrainingTrials, 'frames', L, 'is_significant',is_significant, 'Data', TrainingData, 'plotAll', plotAll); title('SBI factors, with raw data')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
exportgraphics(gcf, sprintf('Simulation_SBI_magnitudes.pdf'), "ContentType","vector")

% Compute similarity to ground truth
[coeff_W_SBI, coeff_H_SBI, ids_SBI] = helper.similarity_WH(W, H, What, Hhat);

% SeqNMF
is_significant = is_significants_SeqNMF{i};
What = W_hats_SeqNMF{i};
Hhat = H_hats_SeqNMF{i};

plotAll = 1;
figure; SimpleWHPlot_patch(What, Hhat, 'trials', TrainingTrials, 'frames', L, 'is_significant', is_significant, 'plotAll', plotAll); title('SeqNMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot_patch(What, Hhat, 'trials', TrainingTrials, 'frames', L, 'is_significant', is_significant, 'Data', TrainingData, 'plotAll', plotAll); title('SeqNMF factors, with raw data')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
exportgraphics(gcf, sprintf('Simulation_SeqNMF_magnitudes.pdf'), "ContentType","vector")

% Compute similarity to ground truth
[coeff_W_SeqNMF, coeff_H_SeqNMF, ids_SeqNMF] = helper.similarity_WH(W, H, What, Hhat);

%% Number of detected motifs and correlation to ground truth
K = 10;
nSim = length(Ws);
num_success_SeqNMF = zeros(K,1);
num_success_SBI = zeros(K,1);
coeffs_W_SeqNMF = zeros(nSim, K);
coeffs_W_SBI = zeros(nSim, K);
coeffs_H_SeqNMF = zeros(nSim, K);
coeffs_H_SBI = zeros(nSim, K);

for i=1:nSim
    is_significant_SeqNMF = logical(is_significants_SeqNMF{i});
    is_significant_SBI = logical(is_significants_SBI{i});
    W = Ws{i};
    H = Hs{i};
    W_hat_SeqNMF = W_hats_SeqNMF{i};
    W_hat_SBI = W_hats_SBI{i};
    H_hat_SeqNMF = H_hats_SeqNMF{i};
    H_hat_SBI = H_hats_SBI{i};
    % pair with ground truth factor
    [coeff_W_SeqNMF, coeff_H_SeqNMF, ids_SeqNMF] = helper.similarity_WH(W, H, W_hat_SeqNMF, H_hat_SeqNMF);
    [coeff_W_SBI, coeff_H_SBI, ids_SBI] = helper.similarity_WH(W, H, W_hat_SBI, H_hat_SBI);
    ids_success_SeqNMF = ids_SeqNMF(is_significant_SeqNMF);
    ids_success_SBI = ids_SBI(is_significant_SBI);
    
    ids_success_SeqNMF(~ids_success_SeqNMF) = [];
    ids_success_SBI(~ids_success_SBI) = [];
    num_success_SeqNMF(ids_success_SeqNMF) = num_success_SeqNMF(ids_success_SeqNMF)+1;
    num_success_SBI(ids_success_SBI) = num_success_SBI(ids_success_SBI)+1;

    coeff_W_SeqNMF(~ids_SeqNMF) = [];
    coeff_H_SeqNMF(~ids_SeqNMF) = [];
    ids_SeqNMF(~ids_SeqNMF) = [];
    coeff_W_SBI(~ids_SBI) = [];
    coeff_H_SBI(~ids_SBI) = [];
    ids_SBI(~ids_SBI) = [];

    coeffs_W_SeqNMF(i,ids_SeqNMF) = coeff_W_SeqNMF;
    coeffs_W_SBI(i,ids_SBI) = coeff_W_SBI;
    coeffs_H_SeqNMF(i,ids_SeqNMF) = coeff_H_SeqNMF;
    coeffs_H_SBI(i,ids_SBI) = coeff_H_SBI;
end

%% Plot coefficients
figure('Position', [100, 100, 1500, 800]); 
c1 = [0 0.4549 0.2902];
c2 = [0.9020 0 0.7843];
thresh_x = 3.5;
thresh_y = .8;

ax1 = subplot(211);
hold on
x = repmat(2*Magnitudes, [nSim,1]);
swarmchart(x-.25, coeffs_W_SeqNMF, 12, c1, 'filled', 'XJitterWidth', .2);
swarmchart(x+.25, coeffs_W_SBI, 12, c2, 'filled', 'XJitterWidth', .2);
% title('Correlation Coeff Multiplicative Update Rule')
patch(2*[thresh_x, K+.5, K+.5, thresh_x], [0, 0, 1, 1], [.5 .5 .5], 'EdgeColor', 'none', 'FaceAlpha', .3)
yline(thresh_y, 'LineWidth', 2, 'LineStyle','--');
set(ax1, 'FontSize', 12, 'XTick', 2*Magnitudes, 'XTickLabel', [], 'YTick', [0, .5, thresh_y, 1], 'box','off')

ax2 = subplot(212);
hold on
swarmchart(x-.25, coeffs_H_SeqNMF, 12, c1, 'filled', 'XJitterWidth', .2);
swarmchart(x+.25, coeffs_H_SBI, 12, c2, 'filled', 'XJitterWidth', .2);
% title('Correlation Coeff Split Bregman Iteration')
set(ax2, 'FontSize', 12, 'XTick', 2*(1:K), 'XTickLabel', num2str((1:K)'), 'YTick', [0, .5, thresh_y, 1], 'box','off')
patch(2*[thresh_x, K+.5, K+.5, thresh_x], [0, 0, 1, 1], [.5 .5 .5], 'EdgeColor', 'none', 'FaceAlpha', .3)
yline(thresh_y, 'LineWidth', 2, 'LineStyle','--');
% xtickangle(ax, 90)
xlabel('Magnitude')
set(ax2, 'XLim', 2*[0, K+.5], 'YLim', [0,1])
linkaxes([ax1, ax2], 'xy')
exportgraphics(gcf, sprintf('Motif_coeffs_magnitudes.pdf'), "ContentType","vector")

%% Number significant motifs
figure;
hold on
plot(Magnitudes, num_success_SeqNMF/nSim, 'Color', c1);
plot(Magnitudes, num_success_SBI/nSim, 'Color', c2);
xlabel('Magnitude')
legend('seqNMF', 'SBI', 'Location','best')
title('Proportion detected motif')
exportgraphics(gcf, sprintf('Motif_detected_Magnitude.pdf'), "ContentType","vector")