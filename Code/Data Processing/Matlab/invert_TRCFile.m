function invert_TRCFile()
% function for inverting data in a .trc file
% The user can choose which file to invert, and the resulting inverted file
% is generated and saved.
%
% If the initial .trc represent a motion from A to B, the output will be a
% motion from B to A.

[trcFile, trcPath] = uigetfile('*.trc', 'Select TRC file to be inverted');

trcInvertedFileName = [trcPath, trcFile(1:end-4), '_inverted.trc'];

[markersExp, timesExp, labelsExp, unitsExp] = readTRC([trcPath, trcFile]);

markersExp_modified = flip(markersExp, 1);

Rate = 1/(timesExp(2)-timesExp(1));

writeMarkersToTRC(trcInvertedFileName, markersExp_modified, labelsExp, Rate, [], timesExp, unitsExp);

disp(['File saved in ', trcInvertedFileName]);
