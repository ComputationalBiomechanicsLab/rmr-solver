function [new_trc_file_name, err] = extract_markers_by_name(original_trc_file, original_path, marker_names_string_array, saving_path)
% This function allows to extract from a given trc file (containing
% trajectories for various markers) only the markers which are specified by
% the user. 
% INPUTS:
% - original_trc_file: name of the trc file to be considered;
% - original_path: path to the trc file;
% - marker_names_string_array: an array of strings specifying the markers
%                              to be retained
% - saving_path: the path where to save the reduced new trc file
%
% OUTPUT:
% - new_trc_file: name for the new file generated
% - err: error code when writing the new trc file
%
% author Italo Belli (i.belli@tudelft.nl) 2023

% load the trc file
[markers, times, labels, units] = readTRC(fullfile(original_path, original_trc_file));
rate = 1/(times(2)-times(1));
frames = 1:max(size(times));

% check whether the required markers are present, and retrieve their index
num_markers_to_retain = max(size(marker_names_string_array));

index_markers_to_retain = zeros(num_markers_to_retain, 1);

for index_mrkr=1:num_markers_to_retain
    index_markers_to_retain(index_mrkr) = find(labels==marker_names_string_array(index_mrkr));
    if index_markers_to_retain(index_mrkr)==0
        error("Marker not found")
    end
end

% transform the marker indexes in column indexes (each marker trajectory is
% a 3D one, saved in 3 different columns in X,Y,Z order)
index_columns_to_retain = [];

for index_mrkr=1:num_markers_to_retain
    index_columns_to_retain = [index_columns_to_retain, 3*(index_markers_to_retain(index_mrkr)-1)+1, 3*(index_markers_to_retain(index_mrkr)-1)+2, 3*(index_markers_to_retain(index_mrkr)-1)+3];
end

% retain just the correct columns
markers = markers(:, index_columns_to_retain);
new_labels = labels(index_markers_to_retain);

% generate new file name
new_trc_file_name = extractBefore(original_trc_file, '.trc');
for index_mrkr=1:num_markers_to_retain
    new_trc_file_name = [new_trc_file_name, '_', labels{index_markers_to_retain(index_mrkr)}];
end
new_trc_file_name = [new_trc_file_name, '.trc'];

err = writeMarkersToTRC(fullfile(saving_path, new_trc_file_name), markers, new_labels, rate, frames', times, units);
