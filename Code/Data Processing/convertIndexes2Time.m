% This script is used to decode the CMC indexes which were used in the
% original code from Seth et al. (2019) and convert those into the
% corresponding timings. Those timings can be used to determine
% automatically the instants of interest in analyzing the solution of the
% CMC when run again on the refined IK results (and to cope more
% effectively with the changes in OpenSim versions used).
%
% Original indexes where defined as follows (in terms on start, end, max of
% each movement):
% 
% indsabd01=3043; indendabd01=10875;indmaxabd01=7792;
% indsabd02=1175; indendabd02=7996;indmaxabd02=5129;
% indsabd03=76; indendabd03=7990;indmaxabd03=4471;
% 
% indsabd21=639; indendabd21=8536;indmaxabd21=4439;
% indsabd22=706; indendabd22=8942;indmaxabd22=4646;
% indsabd23=618; indendabd23=8768;indmaxabd23=4986;
% 
% indsflx01=637; indendflx01=7187;indmaxflx01=3669;
% indsflx02=1001; indendflx02=8037;indmaxflx02=4332;
% indsflx03=1103; indendflx03=8177;indmaxflx03=4170;
% 
% indsflx21=32; indendflx21=6972;indmaxflx21=3671;
% indsflx22=967; indendflx22=9357;indmaxflx22=5187;
% indsflx23=1774; indendflx23=9201;indmaxflx23=4584;
% 
% indsshrug01=836; indendshrug01=2930;indmaxshrug01=1536;
% indsshrug02=634; indendshrug02=2716;indmaxshrug02=1558;
% indsshrug03=662; indendshrug03=2724;indmaxshrug03=1651;
% 
% indsshrug21=544; indendshrug21=2303;indmaxshrug21=1480;
% indsshrug22=477; indendshrug22=2924;indmaxshrug22=1608;
% indsshrug23=739; indendshrug23=2895;indmaxshrug23=1739;
%
% This original values are collected into an indexes_matrix to ease the
% processing and readability of the code
%
% Author: Italo Belli (i.belli@tudelft.nl), 2022


clc
clear all

fileName='ThoracoscapularShoulderCMC_states.sto'; %%%%%%%%%%%specifiy the file name of the states file%%%%%%%%%%%%

% Specify the experiments ("tasks") for which the CMC has been run
taskNames={'abd01';'abd02';'abd03';'abd21';'abd22';'abd23';'flx01';'flx02';'flx03';'flx21';'flx22';'flx23';...
    'shrug01';'shrug02';'shrug03';'shrug21';'shrug22';'shrug23'};

% set the path current folder to be the one where this script is contained
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

for i=1:18
    Dir_file(i)=strcat('../../../../Results/original CMC analysis/',taskNames(i),'/',fileName);
end

file=char(Dir_file);

% define the format of the data in the .sto file resulting from CMC analysis 
format_CMC='%f';
for i=1:86
    format_CMC=[format_CMC ' ' '%f'];
end

% allocate the timings_vector in advance, to be filled with the relevant
% time instants (3 for every task, in the order time_start_movement,
% time_end_movement, time_max_movement
timings_matrix = zeros(18, 3);

%% Transform the original indexes in time instants using the original data
indexes_matrix = [3043, 10875, 7792; ... %abd01
                  1175, 7996,  5129; ... %adb02
                  76,   7990,  4471; ... %abd03
                  639,  8536,  4439; ... %adb21
                  706,  8942,  4646; ... %abd22
                  618,  8768,  4986; ... %abd23
                  637,  7187,  3669; ... %flx01
                  1001, 8037,  4332; ... %flx02
                  1103, 8177,  4170; ... %flx03
                  32,   6972,  3671; ... %flx21
                  967,  9357,  5187; ... %flx22
                  1774, 9201,  4584; ... %flx23
                  836,  2930,  1536; ... %shrug01
                  634,  2716,  1558; ... %shrug02
                  662,  2724,  1651; ... %shrug03
                  544,  2303,  1480; ... %shrug21
                  477,  2924,  1608; ... %shrug22
                  739,  2895,  1739];    %shrug23

% Let's now load each file at a time and convert its meaningful indexes
% into time instants
for i=1:18
    fileID=fopen(file(i,:));
    CMC_data=textscan(fileID,format_CMC,'headerlines',7);
    CMC_data=cell2mat(CMC_data);
    for instant=[1,2,3]
        timings_matrix(i, instant) = CMC_data(indexes_matrix(i,instant), 1);
    end
    fclose(fileID);
end

save("timings_matrix_for_plotting_CMC_results", 'timings_matrix');
fprintf("Timings successfully saved in %s", pwd)
