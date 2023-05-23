function [markers_noisy, noisy_file_name, max_deviation] = add_equal_white_noise_to_marker(trc_file_name, path_trc, labels_to_modify, mean_wn, sigma_wn, cutoff_freq, flag_save_file, new_file_name)
% This function allows to add a white noise (mean_wn, sigma_wn) to various 
% markers equally, specifying the cutoff frequency of the 2nd order 
% Butterworth filter through which the white noise is filtered before 
% adding it to the markers.
% INPUTS:
% - trc_file_name: name of the file containing marker data in .trc format
% - path_trc: path to the .trc file
% - labels_to_modify: array of strings with the name of the markers to
%                     inject equal white noise onto
% - mean_wn: mean of the white noise [m]
% - sigma_wn: standerd deviation of the white noise [m]
% - cutoff_freq: cutoff frequency of the 2nd order lowpass Butterworth 
%                filter through which the white noise is preprocessed [Hz]
% - flag_save_file: set to true to save the resulting markers in a .trc 
%                   file, false otehrwise
% - new_file_name: optionally set the name of the file containing the noisy
%                  markers (if 0, then default name is used)
%
% OUTPUTS: 
% - markers_noisy: marker data with noise injected on the selected ones
% - noisy_file_name: name of the file in which the noisy markers are saved
%                    (if flag_save_file = true). 0 otherwise
% - max_deviation: maximum deviation/distance in the final position of the
%                  markers that is induced by the noise, with respect to 
%                  the original position
%
% Author: Italo Belli (i.belli@tudelft.nl) 2023

NN = 2;         % filter order

% load the marker data from .trc file
[markers, times, labels, units] = readTRC(fullfile(path_trc, trc_file_name));
num_timeInstants = size(times,1); 

% convert everything to SI units
if strcmp(units, 'mm')
    markers = markers/1000;
    units = 'm';
end

% check whether the required markers are present, and retrieve their index
num_markers_wn = max(size(labels_to_modify));

index_markers_wn = zeros(num_markers_wn, 1);

for index_mrkr=1:num_markers_wn
    index_markers_wn(index_mrkr) = find(labels==labels_to_modify(index_mrkr));
    if index_markers_wn(index_mrkr)==0
        error("Marker not found")
    end
end

% transform the marker indexes in column indexes (each marker trajectory is
% a 3D one, saved in 3 different columns in X,Y,Z order)
index_columns_wn = [];

for index_mrkr=1:num_markers_wn
    index_columns_wn = [index_columns_wn, 3*(index_markers_wn(index_mrkr)-1)+1, 3*(index_markers_wn(index_mrkr)-1)+2, 3*(index_markers_wn(index_mrkr)-1)+3];
end

% generate the white noise signal, and filter it
noise = normrnd(mean_wn, sigma_wn, num_timeInstants, 3);

Fs = 1/(times(2)-times(1));         % Butterworth lowpass filter to filter wn
Wn = cutoff_freq/(Fs/2);            % normalized cut off frequency 
[B,A] = butter(NN,Wn);
noise_filtered = filter(B,A,noise); % noise (after LP filtering)

max_deviation = max(sqrt(sum(noise_filtered.^2, 2)));

% add the same white noise to all of the 3D coordinates of the required
% markers

noise_filtered_matrix = repmat(noise_filtered, 1, num_markers_wn);

markers_noisy = markers;
markers_noisy(:, index_columns_wn) = markers_noisy(:, index_columns_wn) + noise_filtered_matrix;

% optionally save a copy of the new marker data
if flag_save_file
    frames = 1:max(size(times));
    if ~new_file_name
        noisy_file_name = [extractBefore(trc_file_name, '.trc'), '_noisy', '.trc'];
    else
        noisy_file_name = new_file_name;
    end
    writeMarkersToTRC(fullfile(path_trc, noisy_file_name), markers_noisy, labels, 1/(times(2)-times(1)), frames', times, units);
else
    noisy_file_name = 0;
end
