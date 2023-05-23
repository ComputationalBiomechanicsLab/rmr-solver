function [markers_upsampled, labels, times_upsampled, upsampled_file_name] = get_upsampled_trc_data(trc_file_name, path_trc, upsampling_factor, interpolating_method, flag_save_file, new_file_name)
% This function allows to get the upsampled version of a trc file, and
% optionally save a new trc_file with that information.
% INPUTS:
% - trc_file: file containitg the marker trajectories in a .trc format
% - path_to_file: path to the file
% - upsampling_factor: factor increasing the numerosity of the datapoints
%                      if 10 then 10*number datapoints are returned 
% - interpolating_method: 1 linear interpolation
%                         2 Modified Akima cubic Hermite interpolation
% - flag_save_file: true or false, to also save the upsampled markers.
% - new_file_name: optionally set the name of the file containing the noisy
%                  markers (if 0, then default name is used)
%
% OUTPUTS:
% - markers_downsampled: marker data after downsampling
% - labels: labels for the markers
% - time_downsampled: time vector corresponding to downsampled trajectories
% - downsampled_file_name: name of the file if saved, 0 otherwise
%
% author: Italo Belli (i.belli@tudelft.nl) 2023

% load the marker data from .trc file
[markers, times, labels, units] = readTRC(fullfile(path_trc, trc_file_name));
num_time_points = size(times, 1);

% generate new time vector
start_time = times(1);
end_time = times(end);
times_upsampled = linspace(start_time, end_time, upsampling_factor*num_time_points)';

% generate the new marker data
if interpolating_method == 1 % linear interpolation
    markers_upsampled = interp1(times, markers, times_upsampled);
elseif interpolating_method == 2 % Modified Akima cubic Hermite interpolation
    markers_upsampled = interp1(times, markers, times_upsampled, 'makima');
end

if flag_save_file
    frames = 1:max(size(times_upsampled));
    if ~new_file_name
        upsampled_file_name = [extractBefore(trc_file_name, '.trc'), '_upsampled_', num2str(upsampling_factor) '.trc'];
    else
        upsampled_file_name = new_file_name;
    end
    writeMarkersToTRC(fullfile(path_trc, upsampled_file_name), markers_upsampled, labels, 1/(times_upsampled(2)-times_upsampled(1)), frames', times_upsampled, units);
else
    upsampled_file_name = 0;
end
