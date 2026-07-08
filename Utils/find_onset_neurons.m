function onset_time_neurons = find_onset_neurons(Wk, thresh_onset)
% Find onset of motif neurons
[N,L] = size(Wk);
onset_time_neurons = nan(N,1);
for n=1:N
    Wk_tmp = double(Wk(n,:));
%     Wk_tmp_smooth = filtfilt(ones(1,5)/5, 1, Wk_tmp);
    onset_time_tmp = find(Wk_tmp/norm(Wk_tmp)>thresh_onset, 1, 'first');
    if ~isempty(onset_time_tmp)
        onset_time_neurons(n) = onset_time_tmp(1);
    end
end