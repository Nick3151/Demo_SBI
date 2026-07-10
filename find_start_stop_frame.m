function [frames_start, frames_stop, peak_speed] = find_start_stop_frame(frames_result, fs, rot_speed_trials)
% Find the start and stop frame of wheel movement
trials = length(frames_result);
assert(trials==size(rot_speed_trials,1), 'Total trials do not match!')
frames_start = zeros(trials,1);
frames_stop = zeros(trials,1);
peak_speed = zeros(trials,1);

% Smooth data in each trial
smooth_win = 1;
frames = size(rot_speed_trials, 2);
rot_speed_smooth = zeros(trials, frames);
for i=1:trials
    rot_speed_smooth(i,:) = filtfilt(ones(1,smooth_win)/smooth_win, 1, rot_speed_trials(i,:));
end

% Estimate baseline speed std (first second)
rot_speed_baseline = std(abs(rot_speed_trials(:,1:round(fs))), 0, 'all');

for i=1:trials 
    [peak_speed(i),motion_peak_frame] = max(abs(rot_speed_smooth(i,round(fs)+1:min(frames_result(i),frames))));
    motion_peak_frame = motion_peak_frame+round(fs);
    % Find candidate onsets and offsets using smoothed data
    rot_direction_tmp = sign(rot_speed_smooth(i,motion_peak_frame));
    frame_rotary_smooth_tmp = (rot_speed_smooth(i,:).*rot_direction_tmp>0);
    frame_start_candidate = find(diff(frame_rotary_smooth_tmp));
    frame_stop_candidate = find(diff(frame_rotary_smooth_tmp)<0);
%     frame_start_tmp = find(rot_direction_tmp*rot_speed_smooth_sel(i,1:motion_peak_idx_trials(i)-1)<=0 & ...
%         rot_direction_tmp*rot_speed_smooth_sel(i,2:motion_peak_idx_trials(i))>0);
%     frame_stop_tmp = find(rot_direction_tmp*rot_speed_smooth_sel(i,motion_peak_idx_trials(i):size(rot_speed_smooth_sel,2)-1)>0 & ...
%         rot_direction_tmp*rot_speed_smooth_sel(i,motion_peak_idx_trials(i)+1:size(rot_speed_smooth_sel,2))<=0);
%     if ~isempty(frame_start_tmp)
%         frames_start(i) = round(fs) + frame_start_tmp(1);
%         if ~isempty(frame_stop_tmp)
%             frames_stop(i) = round(fs) + motion_peak_idx_trials(i)-1 + frame_stop_tmp(1);
%         end    
%     end
    
    % Make sure movement start before and stop after peak velocity time
    frame_start_candidate(frame_start_candidate>motion_peak_frame) = [];
    frame_stop_candidate(frame_stop_candidate<motion_peak_frame) = [];

    % Determine movement onset and offset based on raw trajectories
    if ~isempty(frame_start_candidate)
        rot_speed_raw_tmp = rot_speed_trials(i,frame_start_candidate(end)+(0:2*smooth_win-1));
        idx_tmp = find(rot_speed_raw_tmp.*rot_direction_tmp>0);
        frames_start(i) = frame_start_candidate(end)+idx_tmp(1)-1;
        if ~isempty(frame_stop_candidate)
            rot_speed_raw_tmp = rot_speed_trials(i,frame_stop_candidate(1)+(0:min(2*smooth_win-1, frames-frame_stop_candidate(1))));
            idx_tmp = find(rot_speed_raw_tmp.*rot_direction_tmp<=0);
            frames_stop(i) = frame_stop_candidate(1)+idx_tmp(1)-1;
        else
            frames_stop(i) = frames_result(i);
        end
    end
end
% frames_start(~frames_stop) = 0;
