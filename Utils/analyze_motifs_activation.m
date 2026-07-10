function [active_time_motifs, magnitudes_motifs, neurons_active_motifs, onsets_motifs] ...
= analyze_motifs_activation(X, H, W, varargin)
% Compute the active time, magnitudes of motifs, 
% Get the onset time of the first active neuron within motif

p = inputParser;
addOptional(p, 'thresh_motif_magnitude', .05)
addOptional(p, 'thresh_onset', .05)
addOptional(p, 'thresh_active_freq', .7)

parse(p, varargin{:});
thresh_motif_magnitude = p.Results.thresh_motif_magnitude;  % Minimum motif magnitude
thresh_onset = p.Results.thresh_onset;  % activity thresh for motif onset
thresh_active_freq = p.Results.thresh_active_freq;  % Minimum active frequency for single neuron within motif

[K,T] = size(H);
[N,~,L] = size(W);
active_time_motifs = cell(K, 1);
magnitudes_motifs = cell(K, 1);
neurons_active_motifs = cell(K, 1);
onsets_motifs = cell(K,1);

H_filt = filtfilt(ones(1,5)/5,1,double(H)')';
H_pad = zeros(K, T+2);
H_pad(:,2:T+1) = H_filt;
dH = diff(H_pad>0,1,2);
X = [X, zeros(N,L)];

for k=1:K
    times = [];
    magnitudes = [];
    dH_tmp = dH(k,:);
    h_start = find(dH_tmp==1);
    h_end = find(dH_tmp==-1);
    for i=1:length(h_start)
        if h_start(i)==h_end(i)-1
            times = [times, h_start(i)];
            magnitudes = [magnitudes, H(k,h_start(i))];
        else
            times = [times, sum(H(k,h_start(i):h_end(i)-1).*(h_start(i):h_end(i)-1))/sum(H(k,h_start(i):h_end(i)-1))]; % geometric mean
            magnitudes = [magnitudes, sum(H(k,h_start(i):h_end(i)-1))];
        end
    end
    % Drop temporal loadings that are too small
    times(magnitudes<thresh_motif_magnitude) = [];
    magnitudes(magnitudes<thresh_motif_magnitude) = [];
    active_time_motifs{k} = times;
    magnitudes_motifs{k} = magnitudes;

    % Test active neurons
    if isempty(times)
        neurons_active_motifs{k} = [];
        onsets_motifs{k} = NaN;
    else
        X_tmp = zeros(N,length(times), L);
        Wk = squeeze(W(:,k,:));
        for i=1:length(times)
            X_tmp(:,i,:) = X(:,round(times(i))+(1:L));
        end
        [h,p, active_times] = test_significance_neurons(Wk, X_tmp, X);
        neurons_active = find(h & (active_times/length(times)>=thresh_active_freq));
        neurons_active_motifs{k} = neurons_active;
        Wk_active = Wk(neurons_active, :);
    
        if isempty(neurons_active)
            onsets_motifs{k} = NaN;
        else
            onset_time = find_onset_neurons(Wk_active, thresh_onset);
            onsets_motifs{k} = min(onset_time);
        end
    end
end