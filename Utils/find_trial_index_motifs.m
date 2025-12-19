function [offsets_motifs, inds_motifs, magnitudes_motifs] = find_trial_index_motifs(H,frames, thresh)
if isempty(thresh)
    thresh = .05;
end
[K,T] = size(H);
offsets_motifs = cell(K, 1);
inds_motifs = cell(K, 1);
magnitudes_motifs = cell(K, 1);
H_filt = filtfilt(ones(1,5)/5,1,double(H'))';
H_pad = zeros(K, T+2);
H_pad(:,2:T+1) = H_filt;
dH = diff(H_pad>0,1,2);
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
    times(magnitudes<thresh) = [];
    magnitudes(magnitudes<thresh) = [];
    
    inds_motifs{k} = round(times/frames)+1;
    offsets_motifs{k} = times-round(times/frames)*frames;
    magnitudes_motifs{k} = magnitudes;
end
