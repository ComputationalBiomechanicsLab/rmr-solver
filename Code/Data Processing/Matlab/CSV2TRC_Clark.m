%% File used to load CSV as table, from which data is extracted and put in .trc format
% Adapted to work on Clark's database
% Written by Irene Beck (i.l.y.beck@student.tudelft.nl) 25-01-2022.
% Modified: Italo Belli (i.belli@tudelft.nl) 26-01-2022
% Modified: Irene Beck 27-01-2022
% NB: HARDCODED

addpath('U:\staff-umbrella\PTbot\Analysis\MATLAB\Code Ajay\trc')
addpath('U:\PTbot\Analysis\MATLAB\Code Ajay\trc')

path = uigetdir();
framerate = 50;

opts = detectImportOptions('MarkerCALi.CSV');
opts.DataLines = 6;
opts.VariableNamesLine = 3;

Markers_CSV = readtable('MarkerCALi.CSV', opts);
markersExp = Markers_CSV.Variables;
labelsExp = Markers_CSV.Properties.VariableNames(3:3:end);

for i=1:length(labelsExp)
    labelsExp{i} = labelsExp{i}(10:end);
end

markersExp = markersExp(:,3:end);
timesExp = linspace(0,length(markersExp)/framerate,length(markersExp))';

unitsExp = "mm";
if unitsExp == "mm"
    markersExp = markersExp./1000;
    unitsExp = "m";
end

trcSampledFile = ['Calibration.trc'];
writeMarkersToTRC_custom(trcSampledFile, markersExp, labelsExp, framerate, [], timesExp, unitsExp, path);   % save a new .trc file