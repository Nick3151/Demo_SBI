function SimpleXplot(X, varargin)
% Plot raw data with trials
% trials = [# training trials, # test trials]

p = inputParser;
addOptional(p, 'trials', [])
addOptional(p, 'frames', [])
addOptional(p, 'onsets', [])
addOptional(p, 'neg', 0)
addOptional(p, 'n_plot', [])
addOptional(p, 't_plot', [])

parse(p, varargin{:});
trials = p.Results.trials;
frames = p.Results.frames;
onsets = p.Results.onsets;
neg = p.Results.neg;
n_plot = p.Results.n_plot;
t_plot = p.Results.t_plot;

if isempty(trials)
    training_trials = 0;
    test_trials = 0;
elseif length(trials) == 1
    training_trials = trials;
    test_trials = 0;
else
    training_trials = trials(1);
    test_trials = trials(2);
end

if ~isempty(onsets)
    plotOnsets = 1;
else
    plotOnsets = 0;
end

[N,T] = size(X);
if isempty(n_plot)
    n_plot = 1:N;
end
if isempty(t_plot)
    t_plot = 1:T;
end

X_plot = X(n_plot, t_plot);
maxValue = prctile(abs(X_plot),99,'all')+eps;

if neg
    cmap_red = [ones(128,1),linspace(1,0,128)',linspace(1,0,128)'];
    cmap_blue = [linspace(0,1,128)',linspace(0,1,128)',ones(128,1)];
    cmap = [cmap_blue; cmap_red];
    clims = [-maxValue, maxValue]; 
else
    cmap = flipud(gray); 
    clims = [0, maxValue]; 
end

imagesc([t_plot(1), t_plot(end)], [n_plot(1), n_plot(end)], X_plot, clims);
colormap(cmap)
hold on

if plotOnsets
    for trial=1:training_trials
        xline(onsets(trial), 'LineWidth',2, 'Color', [0.4660, 0.6740, 0.1880])
    end
    for trial=training_trials+1:training_trials+test_trials
        xline(onsets(trial), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560]);
    end
else
    for trial=1:training_trials
        xline(trial*frames, 'LineWidth', 2, 'Color', [0.4660, 0.6740, 0.1880]);
    end
    for trial=training_trials+1:training_trials+test_trials
        xline(trial*frames, 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560]);
    end
end  
ind_ticks = mod(n_plot, 10)==0;
yticks(n_plot(ind_ticks))
yticklabels(cellstr(string(n_plot(ind_ticks))))
xlabel('Frame', 'FontSize', 14)
ylabel('#Neuron', 'FontSize', 14)
set(gca, 'Position', [0.05 0.05 0.9 0.9])

if test_trials
    dim = get(gca, 'Position');
    annotation('textbox',[dim(1), dim(2)+dim(4), dim(3)*training_trials/(training_trials+test_trials) 1-dim(4)-dim(2)],...
        'string', 'Training Trials', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom',...
        'LineStyle', 'none', 'Color', '#77AC30', 'FontWeight', 'bold', 'FontSize', 14)
    annotation('textbox',[dim(1)+dim(3)*training_trials/(training_trials+test_trials), dim(2)+dim(4), dim(3)*test_trials/(training_trials+test_trials) 1-dim(4)-dim(2)],...
        'string', 'Test Trials', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom',...
        'LineStyle', 'none', 'Color', '#7E2F8E', 'FontWeight', 'bold', 'FontSize', 14)
end

end