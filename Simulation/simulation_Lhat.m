% Test the impact of Lhat on synthetic data
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
Nmotifs = 2*(1:K);
Nneurons = 5*ones(K, 1); % the number of neurons in each motif
Dt = 3.*ones(K,1); % gap between each member of the motif

Lhats = 10:10:100; % selection of Lhat

nSim = 10;
seeds = randperm(1000, nSim);
pvals_SeqNMF = cell(nSim,length(Lhats));
is_significants_SeqNMF = cell(nSim,length(Lhats));
loadings_SeqNMF = cell(nSim,length(Lhats));
W_hats_SeqNMF = cell(nSim,length(Lhats));
H_hats_SeqNMF = cell(nSim,length(Lhats));
pvals_SBI = pvals_SeqNMF;
is_significants_SBI = is_significants_SeqNMF;
loadings_SBI = loadings_SeqNMF;
W_hats_SBI = W_hats_SeqNMF;
H_hats_SBI = H_hats_SeqNMF;
Ws = W_hats_SeqNMF;
Hs = H_hats_SeqNMF;
TrainingDatas = W_hats_SeqNMF;

for j=1:length(Lhats)
    Lhat = Lhats(j);
    disp(['Lhat=', num2str(Lhats(j))])
    parfor i=1:nSim
        tic
        [X, W, H, X_hat, motif_ind] = generate_data_trials(Trials, L, Nmotifs, Nneurons, Dt, ...
        'seed', seeds(i));
        groups = zeros(Trials,1);
        for k=1:K
            groups(motif_ind{k}) = k;
        end
        Ws{i,j} = W;
        
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
        Hs{i,j} = H_train_full;
    
        lambda = .01;
    
        % Run SeqNMF with multiplication rule
        [W_hat, H_hat, ~,~,loadings_SeqNMF{i,j},power]= seqNMF(TrainingData,'K',K,'L',Lhat,...
                'lambdaL1W', 0, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
        p = .05;
        [pvals_SeqNMF{i,j},is_significants_SeqNMF{i,j}] = test_significance_new(TestData, W_hat,0,p);
        W_hats_SeqNMF{i,j} = W_hat;
        H_hats_SeqNMF{i,j} = H_hat;
        display(['SeqNMF run ' num2str(i) '/' num2str(nSim)])
    
        % Run SeqNMF with Bregman Iteration
        alpha_W = 1e-6;
        alpha_H = 1e-3;
        [W_hat, H_hat, ~,~,loadings_SBI{i,j},power]= FlexMF(TrainingData,'K',K,'L',Lhat,...
                'lambda', lambda, 'alpha_W', alpha_W, 'alpha_H', alpha_H, 'neg_prop', 0, 'maxiter', 50, 'showPlot', 0, 'verbal', 0); 
        p = .05;
        [pvals_SBI{i,j},is_significants_SBI{i,j}] = test_significance_new(TestData, W_hat,0,p);
        W_hats_SBI{i,j} = W_hat;
        H_hats_SBI{i,j} = H_hat;
        TrainingDatas{i,j} = TrainingData;
        display(['FlexMF run ' num2str(i) '/' num2str(nSim)])
        toc
    end
end

save('simulate_results_Lhat.mat', "Lhats", "H_hats_SBI", 'W_hats_SBI', 'is_significants_SBI', 'pvals_SBI', 'loadings_SBI',...
    'loadings_SeqNMF', 'W_hats_SeqNMF', 'H_hats_SeqNMF', "pvals_SeqNMF", "is_significants_SeqNMF", 'Ws', 'Hs', 'TrainingDatas')

%% Visualization
i = 1;
j = 10; % Selection of Lhat
TrainingData = TrainingDatas{i,j};
L = 50;
Lhat = Lhats(j);
[N,Ttrain] = size(TrainingData);
TrainingTrials = Ttrain/L;
W = Ws{i,j};
H = Hs{i,j};

% Split Bregman
is_significant = is_significants_SBI{i,j};
What = W_hats_SBI{i,j};
Hhat = H_hats_SBI{i,j};

plotAll = 1;
figure; SimpleWHPlot_patch(What, Hhat, 'trials', TrainingTrials, 'frames', L, 'is_significant', is_significant, 'plotAll', plotAll); title('SBI reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot_patch(What, Hhat, 'trials', TrainingTrials, 'frames', L, 'is_significant',is_significant, 'Data', TrainingData, 'plotAll', plotAll); title('SBI factors, with raw data')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
exportgraphics(gcf, sprintf('Simulation_SBI_Lhat=%d.pdf', Lhat), "ContentType","vector")

% Compute similarity to ground truth
[coeff_W_SBI, coeff_H_SBI, ids_SBI] = helper.similarity_WH(W, H, What, Hhat);

% SeqNMF
is_significant = is_significants_SeqNMF{i,j};
What = W_hats_SeqNMF{i,j};
Hhat = H_hats_SeqNMF{i,j};

plotAll = 1;
figure; SimpleWHPlot_patch(What, Hhat, 'trials', TrainingTrials, 'frames', L, 'is_significant', is_significant, 'plotAll', plotAll); title('SeqNMF reconstruction')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
figure; SimpleWHPlot_patch(What, Hhat, 'trials', TrainingTrials, 'frames', L, 'is_significant',is_significant, 'Data', TrainingData, 'plotAll', plotAll); title('SeqNMF factors, with raw data')
set(gcf,'Units','normalized','Position',[0.1 0.1 0.8 0.8])
exportgraphics(gcf, sprintf('Simulation_SeqNMF_Lhat=%d.pdf', Lhat), "ContentType","vector")

% Compute similarity to ground truth
[coeff_W_SeqNMF, coeff_H_SeqNMF, ids_SeqNMF] = helper.similarity_WH(W, H, What, Hhat);

%% Number of detected motifs and correlation to ground truth
K = 10;
[nSim, nCond] = size(Ws);
num_success_SeqNMF = zeros(K,nCond);
num_success_SBI = zeros(K,nCond);
coeffs_W_SeqNMF = zeros(nSim, K, nCond);
coeffs_W_SBI = zeros(nSim, K, nCond);
coeffs_H_SeqNMF = zeros(nSim, K, nCond);
coeffs_H_SBI = zeros(nSim, K, nCond);

for i=1:nSim
    for j=1:nCond
        is_significant_SeqNMF = logical(is_significants_SeqNMF{i,j});
        is_significant_SBI = logical(is_significants_SBI{i,j});
        W = Ws{i,j};
        H = Hs{i,j};
        W_hat_SeqNMF = W_hats_SeqNMF{i,j};
        W_hat_SBI = W_hats_SBI{i,j};
        H_hat_SeqNMF = H_hats_SeqNMF{i,j};
        H_hat_SBI = H_hats_SBI{i,j};
        % pair with ground truth factor
        [coeff_W_SeqNMF, coeff_H_SeqNMF, ids_SeqNMF] = helper.similarity_WH(W, H, W_hat_SeqNMF, H_hat_SeqNMF);
        [coeff_W_SBI, coeff_H_SBI, ids_SBI] = helper.similarity_WH(W, H, W_hat_SBI, H_hat_SBI);
        ids_success_SeqNMF = ids_SeqNMF(is_significant_SeqNMF);
        ids_success_SBI = ids_SBI(is_significant_SBI);
        
        ids_success_SeqNMF(~ids_success_SeqNMF) = [];
        ids_success_SBI(~ids_success_SBI) = [];
        num_success_SeqNMF(ids_success_SeqNMF,j) = num_success_SeqNMF(ids_success_SeqNMF,j)+1;
        num_success_SBI(ids_success_SBI,j) = num_success_SBI(ids_success_SBI,j)+1;

        coeff_W_SeqNMF(~ids_SeqNMF) = [];
        coeff_H_SeqNMF(~ids_SeqNMF) = [];
        ids_SeqNMF(~ids_SeqNMF) = [];
        coeff_W_SBI(~ids_SBI) = [];
        coeff_H_SBI(~ids_SBI) = [];
        ids_SBI(~ids_SBI) = [];
    
        coeffs_W_SeqNMF(i,ids_SeqNMF,j) = coeff_W_SeqNMF;
        coeffs_W_SBI(i,ids_SBI,j) = coeff_W_SBI;
        coeffs_H_SeqNMF(i,ids_SeqNMF,j) = coeff_H_SeqNMF;
        coeffs_H_SBI(i,ids_SBI,j) = coeff_H_SBI;
    end
end

%% Make plots   
f1 = figure;
imagesc([Lhats(1), Lhats(end)], [1,K], num_success_SBI/nSim, [0,1])
xlabel('Lhat')
ylabel('K')
title('SBI %detected motif vs Lhats')
colorbar
exportgraphics(gcf, sprintf('Motif_detected_Lhat_SBI.pdf'), "ContentType","vector")

f2 = figure;
imagesc([Lhats(1), Lhats(end)], [1,K], num_success_SeqNMF/nSim, [0,1])
xlabel('Lhat')
ylabel('K')
title('seqNMF %detected motif vs Lhats')
colorbar
exportgraphics(gcf, sprintf('Motif_detected_Lhat_seqNMF.pdf'), "ContentType","vector")

%% Plot coefficients
j = 1;
f3 = figure('Position', [100, 100, 1500, 800]); 
c1 = [0 0.4549 0.2902];
c2 = [0.9020 0 0.7843];
thresh_x = 3.5;
thresh_y = .8;

ax1 = subplot(211);
hold on
% boxplot(coeffs_W_SeqNMF(:,:,j), 'Colors', c1, 'Positions', 2*(1:K)-.25, 'Jitter', .1, 'Widths', .2);
% boxplot(coeffs_W_FlexMF(:,:,j), 'Colors', c2, 'Positions', 2*(1:K)+.25, 'Jitter', .1, 'Widths', .2);
x = repmat(2*(1:K), [nSim,1]);
swarmchart(x-.25, coeffs_W_SeqNMF(:,:,j), 12, c1, 'filled', 'XJitterWidth', .2);
swarmchart(x+.25, coeffs_W_SBI(:,:,j), 12, c2, 'filled', 'XJitterWidth', .2);
% title('Correlation Coeff Multiplicative Update Rule')
patch(2*[thresh_x, K+.5, K+.5, thresh_x], [0, 0, 1, 1], [.5 .5 .5], 'EdgeColor', 'none', 'FaceAlpha', .3)
yline(thresh_y, 'LineWidth', 2, 'LineStyle','--');
set(ax1, 'FontSize', 12, 'XTick', 2*(1:K), 'XTickLabel', [], 'YTick', [0, .5, thresh_y, 1], 'box','off')

ax2 = subplot(212);
hold on
% boxplot(coeffs_H_SeqNMF(:,:,j), 'Colors', c1, 'Positions', 2*(1:K)-.25, 'Jitter', .1, 'Widths', .2);
% boxplot(coeffs_H_FlexMF(:,:,j), 'Colors', c2, 'Positions', 2*(1:K)+.25, 'Jitter', .1, 'Widths', .2);
swarmchart(x-.25, coeffs_H_SeqNMF(:,:,j), 12, c1, 'filled', 'XJitterWidth', .2);
swarmchart(x+.25, coeffs_H_SBI(:,:,j), 12, c2, 'filled', 'XJitterWidth', .2);
% title('Correlation Coeff Split Bregman Iteration')
set(ax2, 'FontSize', 12, 'XTick', 2*(1:K), 'XTickLabel', num2str((1:K)'), 'YTick', [0, .5, thresh_y, 1], 'box','off')
patch(2*[thresh_x, K+.5, K+.5, thresh_x], [0, 0, 1, 1], [.5 .5 .5], 'EdgeColor', 'none', 'FaceAlpha', .3)
yline(thresh_y, 'LineWidth', 2, 'LineStyle','--');
% xtickangle(ax, 90)
% xlabel('# Motif Occurences')
set(ax2, 'XLim', 2*[0, K+.5], 'YLim', [0,1])
linkaxes([ax1, ax2], 'xy')
exportgraphics(gcf, sprintf('Motif_coeffs_Lhat=%d.pdf', Lhats(j)), "ContentType","vector")