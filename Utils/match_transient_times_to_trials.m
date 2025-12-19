function [inds_trials, onset_time_trials] = match_transient_times_to_trials(onset_time, varargin)
% Input: transient onset time within entire session, onset times for each trial
% Output: transient onset time within each trial, index of trials with onsets
p = inputParser;
addOptional(p, 'trials_onset', [])
addOptional(p, 'frames', 150)

parse(p, varargin{:});
trials_onset = p.Results.trials_onset;
frames = p.Results.frames;

n = length(onset_time);
if ~isempty(trials_onset)
    inds_trials = [];
    onset_time_trials = [];
    for i = 1:n
        tmp = find(onset_time(i)>=trials_onset);  
        if ~isempty(tmp)
            inds_trials = [inds_trials tmp(end)];
            onset_time_trials = [onset_time_trials onset_time(i)-trials_onset(tmp(end))];
        end
    end

else
    inds_trials = floor(onset_time/frames)+1;
    onset_time_trials = onset_time-floor(onset_time/frames)*frames;
end
