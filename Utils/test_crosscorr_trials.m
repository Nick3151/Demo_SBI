function [xcf_trials_mean, h, p_values] = test_crosscorr_trials(dFF, x, varargin)
% Compute cross correlation between dFF traces and variable x,
% both with dimension Trials * L
% Use permutation test(across trials) to assess the significance of the correlations

p = inputParser;
addOptional(p, 'nNull', 10000);
addOptional(p, 'NumLags', 15);
addOptional(p, 'seed', 0);

parse(p, varargin{:})
nNull = p.Results.nNull;
NumLags = p.Results.NumLags;
seed = p.Results.seed;

rng(seed);
[trials, ~] = size(dFF);
dFF = double(dFF);

xcf_trials = zeros(trials, 2*NumLags+1);
for i=1:trials
    xcf_trials(i,:) = crosscorr(dFF(i,:), x(i,:), 'NumLags', NumLags);
end
xcf_trials(isnan(xcf_trials)) = 0;
xcf_trials_mean = mean(xcf_trials);

p_values = zeros(1,2*NumLags+1);
for n=1:nNull
    xNull = x(randperm(trials),:);
    xcf_trials_Null = zeros(trials, 2*NumLags+1);
    for i=1:trials
        xcf_trials_Null(i,:) = crosscorr(dFF(i,:), xNull(i,:), 'NumLags', NumLags);
    end
    xcf_trials_Null(isnan(xcf_trials_Null)) = 0;
    xcf_trials_Null_mean = mean(xcf_trials_Null);
    p_values = p_values + (xcf_trials_mean.*(xcf_trials_Null_mean-xcf_trials_mean)>0)/nNull;
end

h = p_values<.05/(2*NumLags+1);

