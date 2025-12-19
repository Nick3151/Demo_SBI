function plot_PSTH(PSTH, varargin)
p = inputParser;
addOptional(p, 'fs', 30)
addOptional(p, 'Onset', [])
addOptional(p, 'Err', 'std')
addOptional(p, 'Color', 'k')

parse(p, varargin{:})
fs = p.Results.fs;
onset = p.Results.Onset;
err = p.Results.Err;
color = p.Results.Color;

[N,T] = size(PSTH);

if N>1
    if strcmp(err, 'std')
        Y1 = mean(PSTH)+std(PSTH);
        Y2 = mean(PSTH)-std(PSTH);
    elseif strcmp(err, 'sem')
        Y1 = mean(PSTH)+std(PSTH)/sqrt(N);
        Y2 = mean(PSTH)-std(PSTH)/sqrt(N);
    elseif strcmp(err, 'ci')
        Y1 = mean(PSTH)+tinv(0.975, N-1)*std(PSTH)/sqrt(N);
        Y2 = mean(PSTH)-tinv(0.975, N-1)*std(PSTH)/sqrt(N);
    else
        error('Err must be sem, std or ci!')
    end
    
    patch([(1:T)/fs, flip((1:T)/fs)], [Y1 flip(Y2)], color, 'FaceAlpha',.3, 'EdgeColor', 'none')
    hold on
    plot((1:T)/fs, mean(PSTH), 'Color', color)
elseif N==1
    plot((1:T)/fs, PSTH, 'Color', color)
end

if ~isempty(onset)
    xline(onset, 'k--')
end