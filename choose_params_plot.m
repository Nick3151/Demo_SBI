%% Demo of choosing hyperparameters for running SBI on real data
% Run choose_params.m first
% Select hyperparameters lambda and alpha_H based on results

clear all;
close all;
root = pwd;
folder_dir = fullfile(root, 'Data');

input_folder_name = '20240618MH041724';
load(fullfile(folder_dir, input_folder_name, 'choose_params.mat'))

% load(fullfile(folder_dir, input_folder_name, 'dFF_preprocessed.mat'));
load(fullfile(folder_dir, input_folder_name, 'bh_data_Liang.mat'));

output_path = fullfile(folder_dir, input_folder_name);

nLambdas = length(lambdas);
nAlpha_Hs = length(alpha_Hs);

%% Plot normalized recon and regularization errors
% The same metrics used for parameter selection in seqNMF
% lambda is usually set between 2*lambda_0 to 5*lambda_0,
% lambda_0 is when the normalized recon cost and reg cost are balanced

f1 = figure;
marginX = 0.1;
marginY = 0.1;

minRs = min(reg_crosses(:)); maxregs= max(reg_crosses(:));
Rs = (reg_crosses-minRs)/(maxregs-minRs);  
minCs =  min(recon_errors(:)); maxCs =  max(recon_errors(:)); 
Cs = (recon_errors-minCs)/(maxCs-minCs); 

for li = 1:nLambdas
    for hi = 1:nAlpha_Hs
        ax = subplot('Position', [marginX+(hi-1)*(1-marginX)/nAlpha_Hs, (nLambdas-li)*(1-marginY)/nLambdas, (1-marginX)/nAlpha_Hs*.8, (1-marginY)/nLambdas*.8]);
        b = bar([Rs(li,hi), Cs(li,hi)],'FaceColor','flat');
        b.CData(2,:) = [1 0 0];
        set(gca,'YLim',[0,1],'XTick',[],'YTick',[], ...
            'Position', [marginX+(hi-1)*(1-marginX)/nAlpha_Hs, (nLambdas-li)*(1-marginY)/nLambdas, (1-marginX)/nAlpha_Hs, (1-marginY)/nLambdas])
    end
end
pos_Lambdas = [zeros(nLambdas,1), (nLambdas-1:-1:0)'*(1-marginY)/nLambdas,... 
    (1-marginX)/nAlpha_Hs*ones(nLambdas,1), (1-marginY)/nLambdas*ones(nLambdas,1)];
pos_Alphas = [marginX+(0:nAlpha_Hs-1)'*(1-marginX)/nAlpha_Hs, (1-marginY)*ones(nAlpha_Hs,1),...
    (1-marginX)/nAlpha_Hs*ones(nAlpha_Hs,1), (1-marginY)/nLambdas*ones(nAlpha_Hs,1)];
for li = 1:nLambdas
    annotation('textbox',pos_Lambdas(li,:),'String',sprintf('\\lambda=%0.3e',lambdas(li)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
for hi = 1:nAlpha_Hs
    annotation('textbox',pos_Alphas(hi,:),'String',sprintf('\\alpha_H=%0.3e',alpha_Hs(hi)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end  
set(f1, 'Position',[100 100 900 900])

%% Plot sparseness and cross correlations between neurons
% Auxiliary metric used in SBI
% Find a point in which both sparseH and normalized corr_prctiles_noise is low

sparseH = cellfun(@mean, sparseness_Hs);
sparseH = sparseH/max(sparseH(:));
f2 = figure;
marginX = 0.1;
marginY = 0.1;
for li = 1:nLambdas
    for hi = 1:nAlpha_Hs
        ax = subplot('Position', [marginX+(hi-1)*(1-marginX)/nAlpha_Hs, (nLambdas-li)*(1-marginY)/nLambdas, (1-marginX)/nAlpha_Hs*.8, (1-marginY)/nLambdas*.8]);
        b = bar([sparseH(li,hi), corr_prctiles_noise(li,hi)/corr_prctiles_train(1,1)],'FaceColor','flat');
        b.CData(2,:) = [1 0 0];
        set(gca,'YLim',[0,1],'XTick',[],'YTick',[], ...
            'Position', [marginX+(hi-1)*(1-marginX)/nAlpha_Hs, (nLambdas-li)*(1-marginY)/nLambdas, (1-marginX)/nAlpha_Hs, (1-marginY)/nLambdas])
    end
end

for li = 1:nLambdas
    annotation('textbox',pos_Lambdas(li,:),'String',sprintf('\\lambda=%0.3e',lambdas(li)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end
for hi = 1:nAlpha_Hs
    annotation('textbox',pos_Alphas(hi,:),'String',sprintf('\\alpha_H=%0.3e',alpha_Hs(hi)),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
end  
set(f2, 'Position',[100 100 900 900])

%% Plot X_hat
% Visualize the motifs and their loadings under your selected parameters

frames = round(5*fs);
li = 6;
hi = 4;
plotAll = 1;
f4 = figure; SimpleWHPlot(Ws{li,hi},Hs{li,hi}, 'trials', length(ind_train), 'frames', frames, 'is_significant', is_significant{li,hi}, 'plotAll', plotAll); 
set(f4,'position',[200,200,1600,900])
gca
title(sprintf('Recon \\lambda=%0.3e, \\alpha_H=%0.3e',lambdas(li),alpha_Hs(hi)))

f5 = figure; SimpleWHPlot(Ws{li,hi},Hs{li,hi}, 'trials', length(ind_train), 'frames', frames, 'is_significant', is_significant{li,hi}, 'Data', X_train, 'plotAll', plotAll); 
set(f5,'position',[200,200,1600,900])
gca
title(sprintf('Raw \\lambda=%0.3e, \\alpha_H=%0.3e',lambdas(li),alpha_Hs(hi)))