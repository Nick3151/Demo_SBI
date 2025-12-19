function [inds_trials, onset_time_trials] = match_motif_times_to_trials(active_time_motifs, onsets_motifs, varargin)
p = inputParser;
addOptional(p, 'trials_onset', [])
addOptional(p, 'frames', 150)

parse(p, varargin{:});
trials_onset = p.Results.trials_onset;
frames = p.Results.frames;

K = length(active_time_motifs);
inds_trials = cell(K,1);
onset_time_trials = cell(K,1);
for k=1:K
    time_motif = active_time_motifs{k}+onsets_motifs{k};
    if ~isempty(trials_onset)
        inds = [];
        onsets = [];
        for i = 1:length(time_motif)
            tmp = find(time_motif(i)>=trials_onset);  
            if ~isempty(tmp)
                inds = [inds tmp(end)];
                onsets = [onsets time_motif(i)-trials_onset(tmp(end))];
            end
        end
        inds_trials{k} = inds;
        onset_time_trials{k} = onsets;
    else
        inds_trials{k} = floor(time_motif/frames)+1;
        onset_time_trials{k} = time_motif-floor(time_motif/frames)*frames;
    end
end