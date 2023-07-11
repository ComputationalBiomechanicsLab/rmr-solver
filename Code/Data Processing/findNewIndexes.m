% This script is used to get the indexes that are necessary to run the main
% PlotMuscleResults.m
clear; clc

% set the path current folder to be the one where this script is contained
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

% Check whether the required timing matrix has already been computed
if isfile('..\..\..\..\Results\CMC analysis\timings_matrix_for_plotting_CMC_results.mat')
     load('..\..\..\..\Results\CMC analysis\timings_matrix_for_plotting_CMC_results.mat')
else
     error("File corresponding to the timings of the experiments has not been found. Please regenerate it by running convertIndexes2Time.m")
end


fileName='ThoracoscapularShoulderCMC_states.sto'; % specifiy the file name of the states file 
taskNames={'abd01';'abd02';'abd03';'abd21';'abd22';'abd23';'flx01';'flx02';'flx03';'flx21';'flx22';'flx23';...
    'shrug01';'shrug02';'shrug03';'shrug21';'shrug22';'shrug23'};

for i=1:18
    Dir_file(i)=strcat('../../../../Results/CMC analysis/',taskNames(i),'/',fileName);
end 

file=char(Dir_file);

% define the format of the data in the .sto file resulting from CMC analysis 
format_CMC='%f';
for i=1:86
    format_CMC=[format_CMC ' ' '%f'];
end

% allocate the index_matrix in advance, to be filled with the relevant
% time instants (3 for every task, in the order time_start_movement,
% time_end_movement, time_max_movement
index_matrix = zeros(18, 3);

% Let's now load each file at a time and retrieve its meaningful indexes
for i=1:18
    fileID=fopen(file(i,:));
    CMC_data=textscan(fileID,format_CMC,'headerlines',7);
    CMC_data=cell2mat(CMC_data);
    tol = 3 * (CMC_data(2,1)-CMC_data(1,1));
    for instant=[1,2,3]
        index_selected = find(abs(CMC_data(:, 1)-timings_matrix(i, instant))<tol);
        index_matrix(i, instant) = index_selected(1);
    end
    fclose(fileID);
end


save("index_matrix_CMC_newStudy", 'index_matrix');
fprintf("New indexes successfully saved in %s \n", pwd)