function [markers, times, labels, units] = readTRC(trcFile)
% Read marker locations from .trc file
% Optionally provide labels and units for marker locations
% USAGE: [markers, labels, units] = readTRC(trcFile)
% Ajay Seth

markers_struct = importdata(trcFile, '\t', 5);
markers_data = markers_struct.data;
frameNums = markers_data(:,1);
times = markers_data(:,2);
markers = markers_data(:,3:end);

fileProps = textscan(markers_struct.textdata{2}, '%s\t');
fileProps = fileProps{1};
propValues = textscan(markers_struct.textdata{3}, '%s\t');
propValues = propValues{1};

ixnm = strmatch('NumMarkers', fileProps);
ixu = strmatch('Units', fileProps);

numMarkers = str2num(propValues{ixnm});
units = propValues{ixu};

labels = textscan(markers_struct.textdata{4}, '%s', numMarkers + 2);
labels = {labels{1}{3: numMarkers + 2}};

