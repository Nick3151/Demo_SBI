function [xcf, h, p_values] = test_crosscorr(dFF, x, varargin)
% Compute cross correlation between dFF trace and variable x,
% Use permutation test to assess the significance of the correlations

p = inputParser;
addOptional(p, 'nNull', 10000);
addOptional(p, 'NumLags', 10);
addOptional(p, 'seed', 0);

parse(p, varargin{:})
nNull = p.Results.nNull;
NumLags = p.Results.NumLags;
seed = p.Results.seed;

rng(seed);
T = length(dFF);
dFF = double(dFF);
xcf = crosscorr(dFF, x, 'NumLags', NumLags);

xcf_Null = zeros(nNull, 2*NumLags+1);
p_values = zeros(1,2*NumLags+1);
for n=1:nNull
    xNull = x(randperm(T));
    xcf_Null(n,:) = crosscorr(dFF, xNull, 'NumLags', NumLags);
    p_values = p_values + (xcf.*(xcf_Null(n,:)-xcf)>0)/nNull;
end

h = p_values<.05/(2*NumLags+1);

