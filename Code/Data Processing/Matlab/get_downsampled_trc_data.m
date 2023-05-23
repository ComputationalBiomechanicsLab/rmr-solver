function [markers_downsampled, labels, time_downsampled, downsampled_file_name] = get_downsampled_trc_data(trc_file, path_to_file, downsampling_factor, flag_save_file)
% This function allows to get the downsampled version of a trc file, and
% optionally save a new trc_file with that information.
% INPUTS:
% - trc_file: file containitg the marker trajectories in a .trc format
% - path_to_file: path to the file
% - downsampling_factor: retain one data-point every downsampling_factor
%                        ones (if set to 10, the datapoint will have a
%                        frequency 10 times lower than the original one)
% - flag_save_file: true or false, to also save the downsampled markers.
%
% OUTPUTS:
% - markers_downsampled: marker data after downsampling
% - labels: labels for the markers
% - time_downsampled: time vector corresponding to downsampled trajectories
% - downsampled_file_name: name of the file if saved, 0 otherwise
%
% author: Italo Belli (i.belli@tudelft.nl) 2023

% load the original marker trajectories and data from .trc file
[markers, times, labels, units] = readTRC(fullfile(path_to_file, trc_file));

% check whether the downsampling factor makes sense
num_datapoints = size(times, 1);
if downsampling_factor> num_datapoints
    error("The current value of the downsampling factor is excessive (greater than the datapoints).")
end

% downsample the marker data and the time vector
markers_downsampled = markers(1:downsampling_factor:end, :);
time_downsampled = times(1:downsampling_factor:end, :);

% optionally save a copy of the new marker data
if flag_save_file
    frames = 1:max(size(time_downsampled));
    downsampled_file_name = [extractBefore(trc_file, '.trc'), '_downsampled_', num2str(downsampling_factor), '.trc'];
    writeMarkersToTRC(fullfile(path_to_file, downsampled_file_name), markers_downsampled, labels, 1/(time_downsampled(2)-time_downsampled(1)), frames', time_downsampled, units);
else
    downsampled_file_name = 0;
end
