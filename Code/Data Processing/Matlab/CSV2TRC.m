%% File used to load CSV as table, from which data is extracted and put in .trc format
% Written by Irene Beck (i.l.y.beck@student.tudelft.nl) 25-01-2022.
% NB: HARDCODED
addpath('U:\staff-umbrella\PTbot\Analysis\MATLAB\Code Ajay\trc')
addpath('U:\PTbot\Analysis\MATLAB\Code Ajay\trc')

path = uigetdir();
framerate = 50;

opts = detectImportOptions('MarkerCAL.CSV');
opts.DataLines = 6;
opts.VariableNamesLine = 3;

Markers_CSV = readtable('MarkerCAL.CSV', opts);
markersExp = Markers_CSV.Variables;
labelsExp = Markers_CSV.Properties.VariableNames(3:3:end);

markersExp = markersExp(:,3:end);
timesExp = linspace(0,length(markersExp)/framerate,length(markersExp))';

unitsExp = "mm";
if unitsExp == "mm"
    markersExp = markersExp./1000;
    unitsExp = "m";
end

trcSampledFile = ['Calibration_Irene.trc'];
writeMarkersToTRC_custom(trcSampledFile, markersExp, labelsExp, framerate, [], timesExp, unitsExp, path);   % save a new .trc file