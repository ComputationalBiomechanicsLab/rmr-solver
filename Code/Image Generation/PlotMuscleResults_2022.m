% Main script to analyze the results of the CMC algorithm and RMR solver, and compare them together. 
% EMG-based activations are plot alongside with CMC results and RMR results
% for the tasks selected by the user through a prompt-window.

clc; clear

% set the path current folder to be the one where this script is contained
mfile_name          = mfilename('fullpath');
[pathstr,~,~]  = fileparts(mfile_name);
cd(pathstr);

addpath(pathstr)
cd '..\..\..'
path_to_repo = pwd;
addpath(path_to_repo)

%% SETTINGS for the study
% specifiy the file name of the states file, output of CMC
fileName='ThoracoscapularShoulderCMC_states.sto';

% Specify the experiments ("tasks") for which the CMC has been run
taskNames={'abd01';'abd02';'abd03';'abd21';'abd22';'abd23';'flx01';'flx02';'flx03';'flx21';'flx22';'flx23';...
           'shrug01';'shrug02';'shrug03';'shrug21';'shrug22';'shrug23'};

% Specify the indexes to be used to compare the CMC results with the EMG
% data
muscle_CMC=[36, 37, 38, 49, 51, 53, 61, 89, 67, 57, 64];  % based on the ordering they have in the columns of ThoracoscapularShoulderCMC_states.sto
muscle_EMG=[6, 5, 7, 2, 4, 3, 12, 9, 8, 10, 13];          % the order is the same as in Seth et al. 2019 (in file: Frame, AD, MD, PD, TS, TM, TI, Infra, SA, LatD, PM_Sternal, PM_Clav, TMajor, Bic, Tric)
muscle_RMR = [1, 2, 3, 12, 13, 14, 18, 51, 23, 16, 21];   % index 51 corresponds to the colums in which the average of serratus anterior activation has been saved during preprocessing

% Maximum Voluntary Contraction values (coming from code in Seth et al, 2019)
MVC=[4.175 0.547 3.202 2.225 1.378 0.706 0.564 1.862 0.173 0.448 1.155];

muscle_Name={'Trapezius, scapula middle';
    'Trapezius, scapula superior';
    'Trapezius, scapula inferior';
    'Deltoid anterior';
    'Deltoid posterior';
    'Deltoid middle';
    'Pectoralis major clavical';
    'Serratus anterior';
    'Infraspinatus (Superior)';
    'Latissimus Dorsi (Medialis)';
    'Teres major'};

% load the idexes corresponding to the different phases of each motion, in the CMC results 
load(fullfile(path_to_repo, 'Results/new CMC analysis/index_matrix_CMC_newStudy.mat'));

%% User selection phase.
% The user selects:
% 1. the task to consider (ABD_0kg, ABD_2kg, ...)
% 2. whether to consider the RMR results with GH constraint
% 3. whether to consider the RMR results without GH constraint

% 1. prompt user with the selction of which task to consider
task_index = listdlg('PromptString',{'Task to consider'}, ...
    'SelectionMode','single','ListString', ...
    {'ABD free', 'ABD 2 kg','FLX free','FLX 2 kg', 'SHRUG free', 'SHRUG 2 kg'});
if task_index==1
    disp("Task considered: ABD_0kg")
    task = 'ABD_0kg';
    
elseif task_index==2
    disp("Task considered: ABD_2kg")
    task = 'ABD_2kg';

elseif task_index==3
    disp("Task considered: FLX_0kg")
    task = 'FLX_0kg';

elseif task_index==4
    disp("Task considered: FLX_2kg")
    task = 'FLX_2kg';

elseif task_index==5
    disp("Task considered: SHRUG_0kg")
    task = 'SHRUG_0kg';

elseif task_index==6
    disp("Task considered: SHRUG_2kg")
    task = 'SHRUG_2kg';
end

% 2. Select which RMR solver results to consider (with GH constraint)
answer_rmr_selection = listdlg('PromptString',{'Plot RMR analysis with GH constraint'}, ...
    'SelectionMode','single','ListString', ...
    {'no', '100 Hz', '10 Hz'});

if answer_rmr_selection==1
    disp("Disregarding RMR results (with GH constraint)")
    path_RMR_GH = 0;
    
elseif answer_rmr_selection==2
    disp("Considering RMR results at 100 Hz (with GH constraint)")
    path_RMR_GH = fullfile(path_to_repo, '\Results\Rapid MRS\100 Hz\with GH');

elseif answer_rmr_selection==3
    disp("Considering RMR results at 10 Hz (with GH constraint)")
    path_RMR_GH = fullfile(path_to_repo, '\Results\Rapid MRS\10 Hz\with GH');
end

% 3. Select which RMR solver results to consider (without GH constraint)
answer_rmr_selection = listdlg('PromptString',{'Plot RMR analysis without GH constraint'}, ...
    'SelectionMode','single','ListString', ...
    {'no', '100 Hz', '10 Hz'});

if answer_rmr_selection==1
    disp("Disregarding RMR results (without GH constraint)")
    path_RMR_noGH = 0;
    
elseif answer_rmr_selection==2
    disp("Considering RMR results at 100 Hz (without GH constraint)")
    path_RMR_noGH = fullfile(path_to_repo, '\Results\Rapid MRS\100 Hz\without GH');

elseif answer_rmr_selection==3
    disp("Considering RMR results at 10 Hz (without GH constraint)")
    path_RMR_noGH = fullfile(path_to_repo, '\Results\Rapid MRS\10 Hz\without GH');
end

disp("Considering CMC results")
% Select the directory containing CMC results
for i=1:18
    Dir_file_CMC(i)=fullfile(path_to_repo, 'Results/new CMC analysis/',taskNames(i),'/',fileName);
end
file=char(Dir_file_CMC);

format_CMC='%f';
for i=1:86
    format_CMC=[format_CMC ' ' '%f'];
end

format_EMG='%f';
for i=1:14
    format_EMG=[format_EMG ' ' '%f'];
end

%%%%%%%%%%%%%%%%%%%%PLOT DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analyze each of the 11 muscles individually
for i=1:11
    if strcmp('ABD_0kg',task)
        fileID=fopen(file(1,:));
        CMCabd01=textscan(fileID,format_CMC,'headerlines',7);
        CMCabd01=cell2mat(CMCabd01);
        fclose(fileID);
        fileID=fopen(fullfile(path_to_repo, 'ExperimentalData/EMG/abd01.exp'));
        EMGabd01=textscan(fileID,format_EMG,'headerlines',9);
        EMGabd01=cell2mat(EMGabd01);
        fclose(fileID);
        
        fileID=fopen(file(2,:));
        CMCabd02=textscan(fileID,format_CMC,'headerlines',7);
        CMCabd02=cell2mat(CMCabd02);
        fclose(fileID);
        fileID=fopen(fullfile(path_to_repo, 'ExperimentalData/EMG/abd02.exp'));
        EMGabd02=textscan(fileID,format_EMG,'headerlines',9);
        EMGabd02=cell2mat(EMGabd02);
        fclose(fileID);
        
        fileID=fopen(file(3,:));
        CMCabd03=textscan(fileID,format_CMC,'headerlines',7);
        CMCabd03=cell2mat(CMCabd03);
        fclose(fileID);
        fileID=fopen(fullfile(path_to_repo, 'ExperimentalData/EMG/abd03.exp'));
        EMGabd03=textscan(fileID,format_EMG,'headerlines',9);
        EMGabd03=cell2mat(EMGabd03);
        fclose(fileID);
        
        % load indexes of relevant phases of the motion from file
        indsabd01=index_matrix(1,1); indendabd01=index_matrix(1,2);indmaxabd01=index_matrix(1,3);
        indsabd02=index_matrix(2,1); indendabd02=index_matrix(2,2);indmaxabd02=index_matrix(2,3);
        indsabd03=index_matrix(3,1); indendabd03=index_matrix(3,2);indmaxabd03=index_matrix(3,3);

        % average the estimates for the serratus muscle across bundles
        % (only for CMC, RMR results have been averaged and saved during
        % pre-processing)
        CMCabd01(:,89)=(CMCabd01(:,40)+CMCabd01(:,41)+CMCabd01(:,42))/3;
        CMCabd02(:,89)=(CMCabd02(:,40)+CMCabd02(:,41)+CMCabd02(:,42))/3;
        CMCabd03(:,89)=(CMCabd03(:,40)+CMCabd03(:,41)+CMCabd03(:,42))/3;

        % plot and perform calculations
        [MAE_abd0kg_CMC, MAE_abd0kg_RMR_GH, MAE_abd0kg_RMR_noGH, corr_abd0kg_CMC, corr_abd0kg_RMR_GH, corr_abd0kg_RMR_noGH]=EMG_CMC_RMR_plot(CMCabd01,CMCabd02,CMCabd03,EMGabd01,EMGabd02,EMGabd03,indsabd01,indendabd01,indmaxabd01,indsabd02,indendabd02,indmaxabd02,indsabd03,indendabd03,indmaxabd03,muscle_CMC(i),muscle_EMG(i),MVC(i),muscle_Name(i), path_RMR_GH, path_RMR_noGH, muscle_RMR(i), task);
        disp('MAE_abd0kg')
        disp(muscle_Name(i))
        disp(strcat('CMC:', {'          '}, string(round(MAE_abd0kg_CMC,2)), ' [MAE] ,    ', string(round(corr_abd0kg_CMC,2)), ' [xcorr]'))
        if path_RMR_GH
            disp(strcat('RMR (GH):', {'     '},  string(round(MAE_abd0kg_RMR_GH,2)), ' [MAE] ,    ', string(round(corr_abd0kg_RMR_GH,2)), ' [xcorr]'))
        end
        if path_RMR_noGH
            disp(strcat('RMR (no GH):', {'  '}, string(round(MAE_abd0kg_RMR_noGH,2)), ' [MAE] ,    ', string(round(corr_abd0kg_RMR_noGH,2)), ' [xcorr]'))
        end
        fprintf('\n \n')
        
    elseif strcmp('ABD_2kg',task)
        fileID=fopen(file(4,:));
        CMCabd21=textscan(fileID,format_CMC,'headerlines',7);
        CMCabd21=cell2mat(CMCabd21);
        fclose(fileID);
        fileID=fopen(fullfile(path_to_repo, 'ExperimentalData/EMG/abd21.exp'));
        EMGabd21=textscan(fileID,format_EMG,'headerlines',9);
        EMGabd21=cell2mat(EMGabd21);
        fclose(fileID);
        
        fileID=fopen(file(5,:));
        CMCabd22=textscan(fileID,format_CMC,'headerlines',7);
        CMCabd22=cell2mat(CMCabd22);
        fclose(fileID);
        fileID=fopen(fullfile(path_to_repo, 'ExperimentalData/EMG/abd22.exp'));
        EMGabd22=textscan(fileID,format_EMG,'headerlines',9);
        EMGabd22=cell2mat(EMGabd22);
        fclose(fileID);
        
        fileID=fopen(file(6,:));
        CMCabd23=textscan(fileID,format_CMC,'headerlines',7);
        CMCabd23=cell2mat(CMCabd23);
        fclose(fileID);
        fileID=fopen(fullfile(path_to_repo, 'ExperimentalData/EMG/abd23.exp'));
        EMGabd23=textscan(fileID,format_EMG,'headerlines',9);
        EMGabd23=cell2mat(EMGabd23);
        fclose(fileID);
        
        % load indexes of relevant phases of the motion from file
        indsabd21=index_matrix(4,1); indendabd21=index_matrix(4,2); indmaxabd21=index_matrix(4,3);
        indsabd22=index_matrix(5,1); indendabd22=index_matrix(5,2); indmaxabd22=index_matrix(5,3);
        indsabd23=index_matrix(6,1); indendabd23=index_matrix(6,2); indmaxabd23=index_matrix(6,3);

        % average the estimates for the serratus muscle across bundles
        % (only for CMC, RMR results have been averaged and saved during
        % pre-processing)
        CMCabd21(:,89)=(CMCabd21(:,40)+CMCabd21(:,41)+CMCabd21(:,42))/3;
        CMCabd22(:,89)=(CMCabd22(:,40)+CMCabd22(:,41)+CMCabd22(:,42))/3;
        CMCabd23(:,89)=(CMCabd23(:,40)+CMCabd23(:,41)+CMCabd23(:,42))/3;

        % plot and perform calculations
        [MAE_abd2kg_CMC, MAE_abd2kg_RMR_GH,MAE_abd2kg_RMR_noGH, corr_abd2kg_CMC, corr_abd2kg_RMR_GH,  corr_abd2kg_RMR_noGH]=EMG_CMC_RMR_plot(CMCabd21,CMCabd22,CMCabd23,EMGabd21,EMGabd22,EMGabd23,indsabd21,indendabd21,indmaxabd21,indsabd22,indendabd22,indmaxabd22,indsabd23,indendabd23,indmaxabd23,muscle_CMC(i),muscle_EMG(i),MVC(i),muscle_Name(i), path_RMR_GH, path_RMR_noGH, muscle_RMR(i), task);
        disp('MAE_abd2kg')
        disp(muscle_Name(i))
        disp(strcat('CMC:', {'          '}, string(round(MAE_abd2kg_CMC,2)), ' [MAE] ,    ', string(round(corr_abd2kg_CMC,2)), ' [xcorr]'))
        if path_RMR_GH
            disp(strcat('RMR (GH):', {'     '},  string(round(MAE_abd2kg_RMR_GH,2)), ' [MAE] ,    ', string(round(corr_abd2kg_RMR_GH,2)), ' [xcorr]'))
        end
        if path_RMR_noGH
            disp(strcat('RMR (no GH):', {'  '}, string(round(MAE_abd2kg_RMR_noGH,2)), ' [MAE] ,    ', string(round(corr_abd2kg_RMR_noGH,2)), ' [xcorr]'))
        end
        fprintf('\n \n')
        
    elseif strcmp('FLX_0kg',task)
        fileID=fopen(file(7,:));
        CMCflx01=textscan(fileID,format_CMC,'headerlines',7);
        CMCflx01=cell2mat(CMCflx01);
        fileID=fopen(fullfile(path_to_repo, 'ExperimentalData/EMG/flx01.exp'));
        EMGflx01=textscan(fileID,format_EMG,'headerlines',9);
        EMGflx01=cell2mat(EMGflx01);
        fclose(fileID);
        
        fileID=fopen(file(8,:));
        CMCflx02=textscan(fileID,format_CMC,'headerlines',7);
        CMCflx02=cell2mat(CMCflx02);
        fclose(fileID);
        fileID=fopen(fullfile(path_to_repo, 'ExperimentalData/EMG/flx02.exp'));
        EMGflx02=textscan(fileID,format_EMG,'headerlines',9);
        EMGflx02=cell2mat(EMGflx02);
        fclose(fileID);
        
        fileID=fopen(file(9,:));
        CMCflx03=textscan(fileID,format_CMC,'headerlines',7);
        CMCflx03=cell2mat(CMCflx03);
        fclose(fileID);
        fileID=fopen(fullfile(path_to_repo, 'ExperimentalData/EMG/flx03.exp'));
        EMGflx03=textscan(fileID,format_EMG,'headerlines',9);
        EMGflx03=cell2mat(EMGflx03);
        fclose(fileID);
        
        % load indexes of relevant phases of the motion from file
        indsflx01=index_matrix(7,1); indendflx01=index_matrix(7,2); indmaxflx01=index_matrix(7,3);
        indsflx02=index_matrix(8,1); indendflx02=index_matrix(8,2); indmaxflx02=index_matrix(8,3);
        indsflx03=index_matrix(9,1); indendflx03=index_matrix(9,2); indmaxflx03=index_matrix(9,3);

        % average the estimates for the serratus muscle across bundles
        % (only for CMC, RMR results have been averaged and saved during
        % pre-processing)
        CMCflx01(:,89)=(CMCflx01(:,40)+CMCflx01(:,41)+CMCflx01(:,42))/3;
        CMCflx02(:,89)=(CMCflx02(:,40)+CMCflx02(:,41)+CMCflx02(:,42))/3;
        CMCflx03(:,89)=(CMCflx03(:,40)+CMCflx03(:,41)+CMCflx03(:,42))/3;

        % plot and perform calculations
        [MAE_flx0kg_CMC, MAE_flx0kg_RMR_GH, MAE_flx0kg_RMR_noGH, corr_flx0kg_CMC, corr_flx0kg_RMR_GH, corr_flx0kg_RMR_noGH]=EMG_CMC_RMR_plot(CMCflx01,CMCflx02,CMCflx03,EMGflx01,EMGflx02,EMGflx03,indsflx01,indendflx01,indmaxflx01,indsflx02,indendflx02,indmaxflx02,indsflx03,indendflx03,indmaxflx03,muscle_CMC(i),muscle_EMG(i),MVC(i),muscle_Name(i), path_RMR_GH, path_RMR_noGH, muscle_RMR(i), task);
        disp('MAE_flx0kg')
        disp(muscle_Name(i))
        disp(strcat('CMC:', {'          '}, string(round(MAE_flx0kg_CMC,2)), ' [MAE] ,    ', string(round(corr_flx0kg_CMC,2)), ' [xcorr]'))
        if path_RMR_GH
            disp(strcat('RMR (GH):', {'     '}, string(round(MAE_flx0kg_RMR_GH,2)),  ' [MAE] ,    ', string(round(corr_flx0kg_RMR_GH,2)), ' [xcorr]'))
        end
        if path_RMR_noGH
            disp(strcat('RMR (no GH):', {'  '}, string(round(MAE_flx0kg_RMR_noGH,2)),  ' [MAE] ,    ', string(round(corr_flx0kg_RMR_noGH,2)), ' [xcorr]'))
        end
        fprintf('\n \n')
        
    elseif strcmp('FLX_2kg',task)
        fileID=fopen(file(10,:));
        CMCflx21=textscan(fileID,format_CMC,'headerlines',7);
        CMCflx21=cell2mat(CMCflx21);
        fclose(fileID);
        fileID=fopen(fullfile(path_to_repo, 'ExperimentalData/EMG/flx21.exp'));
        EMGflx21=textscan(fileID,format_EMG,'headerlines',9);
        EMGflx21=cell2mat(EMGflx21);
        fclose(fileID);
        
        fileID=fopen(file(11,:));
        CMCflx22=textscan(fileID,format_CMC,'headerlines',7);
        CMCflx22=cell2mat(CMCflx22);
        fclose(fileID);
        fileID=fopen(fullfile(path_to_repo, 'ExperimentalData/EMG/flx22.exp'));
        EMGflx22=textscan(fileID,format_EMG,'headerlines',9);
        EMGflx22=cell2mat(EMGflx22);
        fclose(fileID);
        
        fileID=fopen(file(12,:));
        CMCflx23=textscan(fileID,format_CMC,'headerlines',7);
        CMCflx23=cell2mat(CMCflx23);
        fclose(fileID);
        fileID=fopen(fullfile(path_to_repo, 'ExperimentalData/EMG/flx23.exp'));
        EMGflx23=textscan(fileID,format_EMG,'headerlines',9);
        EMGflx23=cell2mat(EMGflx23);
        fclose(fileID);
        
        % load indexes of relevant phases of the motion from file
        indsflx21=index_matrix(10,1); indendflx21=index_matrix(10,2); indmaxflx21=index_matrix(10,3);
        indsflx22=index_matrix(11,1); indendflx22=index_matrix(11,2); indmaxflx22=index_matrix(11,3);
        indsflx23=index_matrix(12,1); indendflx23=index_matrix(12,2); indmaxflx23=index_matrix(12,3);

        % average the estimates for the serratus muscle across bundles
        % (only for CMC, RMR results have been averaged and saved during
        % pre-processing)
        CMCflx21(:,89)=(CMCflx21(:,40)+CMCflx21(:,41)+CMCflx21(:,42))/3;
        CMCflx22(:,89)=(CMCflx22(:,40)+CMCflx22(:,41)+CMCflx22(:,42))/3;
        CMCflx23(:,89)=(CMCflx23(:,40)+CMCflx23(:,41)+CMCflx23(:,42))/3;

        % plot and perform calculations
        [MAE_flx2kg_CMC, MAE_flx2kg_RMR_GH, MAE_flx2kg_RMR_noGH, corr_flx2kg_CMC, corr_flx2kg_RMR_GH, corr_flx2kg_RMR_noGH]=EMG_CMC_RMR_plot(CMCflx21,CMCflx22,CMCflx23,EMGflx21,EMGflx22,EMGflx23,indsflx21,indendflx21,indmaxflx21,indsflx22,indendflx22,indmaxflx22,indsflx23,indendflx23,indmaxflx23,muscle_CMC(i),muscle_EMG(i),MVC(i),muscle_Name(i), path_RMR_GH, path_RMR_noGH, muscle_RMR(i), task);
        disp('MAE_flx2kg')
        disp(muscle_Name(i))
        disp(strcat('CMC:', {'          '}, string(round(MAE_flx2kg_CMC,2)), ' [MAE] ,    ', string(round(corr_flx2kg_CMC,2)), ' [xcorr]'))
        if path_RMR_GH
            disp(strcat('RMR (GH):', {'     '}, string(round(MAE_flx2kg_RMR_GH,2)), ' [MAE] ,    ', string(round(corr_flx2kg_RMR_GH,2)), ' [xcorr]'))
        end
        if path_RMR_noGH
            disp(strcat('RMR (no GH):', {'  '}, string(round(MAE_flx2kg_RMR_noGH,2)), ' [MAE] ,    ', string(round(corr_flx2kg_RMR_noGH,2)), ' [xcorr]'))
        end
        fprintf('\n \n')
        
    elseif strcmp('SHRUG_0kg',task)
        fileID=fopen(file(13,:));
        CMCshrug01=textscan(fileID,format_CMC,'headerlines',7);
        CMCshrug01=cell2mat(CMCshrug01);
        fclose(fileID);
        fileID=fopen(fullfile(path_to_repo, 'ExperimentalData/EMG/shrug01.exp'));
        EMGshrug01=textscan(fileID,format_EMG,'headerlines',9);
        EMGshrug01=cell2mat(EMGshrug01);
        fclose(fileID);
        
        fileID=fopen(file(14,:));
        CMCshrug02=textscan(fileID,format_CMC,'headerlines',7);
        CMCshrug02=cell2mat(CMCshrug02);
        fclose(fileID);
        fileID=fopen(fullfile(path_to_repo, 'ExperimentalData/EMG/shrug02.exp'));
        EMGshrug02=textscan(fileID,format_EMG,'headerlines',9);
        EMGshrug02=cell2mat(EMGshrug02);
        fclose(fileID);
        
        fileID=fopen(file(15,:));
        CMCshrug03=textscan(fileID,format_CMC,'headerlines',7);
        CMCshrug03=cell2mat(CMCshrug03);
        fclose(fileID);
        fileID=fopen(fullfile(path_to_repo, 'ExperimentalData/EMG/shrug03.exp'));
        EMGshrug03=textscan(fileID,format_EMG,'headerlines',9);
        EMGshrug03=cell2mat(EMGshrug03);
        fclose(fileID);
        
        % load indexes of relevant phases of the motion from file
        indsshrug01=index_matrix(13,1); indendshrug01=index_matrix(13,2); indmaxshrug01=index_matrix(13,3);
        indsshrug02=index_matrix(14,1); indendshrug02=index_matrix(14,2); indmaxshrug02=index_matrix(14,3);
        indsshrug03=index_matrix(15,1); indendshrug03=index_matrix(15,2); indmaxshrug03=index_matrix(15,3);

        % average the estimates for the serratus muscle across bundles
        % (only for CMC, RMR results have been averaged and saved during
        % pre-processing)
        CMCshrug01(:,89)=(CMCshrug01(:,40)+CMCshrug01(:,41)+CMCshrug01(:,42))/3;
        CMCshrug02(:,89)=(CMCshrug02(:,40)+CMCshrug02(:,41)+CMCshrug02(:,42))/3;
        CMCshrug03(:,89)=(CMCshrug03(:,40)+CMCshrug03(:,41)+CMCshrug03(:,42))/3;

        % plot and perform calculations
        [MAE_shrug0kg_CMC, MAE_shrug0kg_RMR_GH, MAE_shrug0kg_RMR_noGH, corr_shrug0kg_CMC, corr_shrug0kg_RMR_GH, corr_shrug0kg_RMR_noGH]=EMG_CMC_RMR_plot(CMCshrug01,CMCshrug02,CMCshrug03,EMGshrug01,EMGshrug02,EMGshrug03,indsshrug01,indendshrug01,indmaxshrug01,indsshrug02,indendshrug02,indmaxshrug02,indsshrug03,indendshrug03,indmaxshrug03,muscle_CMC(i),muscle_EMG(i),MVC(i),muscle_Name(i), path_RMR_GH, path_RMR_noGH, muscle_RMR(i), task);
        disp('MAE_shrug0kg')
        disp(muscle_Name(i))
        disp(strcat('CMC:', {'          '}, string(round(MAE_shrug0kg_CMC,2)), ' [MAE] ,    ', string(round(corr_shrug0kg_CMC,2)), ' [xcorr]'))
        if path_RMR_GH
            disp(strcat('RMR (GH):', {'     '}, string(round(MAE_shrug0kg_RMR_GH,2)), ' [MAE] ,    ', string(round(corr_shrug0kg_RMR_GH,2)), ' [xcorr]'))
        end
        if path_RMR_noGH
            disp(strcat('RMR (no GH):', {'  '}, string(round(MAE_shrug0kg_RMR_noGH,2)), ' [MAE] ,    ', string(round(corr_shrug0kg_RMR_noGH,2)), ' [xcorr]'))
        end
        fprintf('\n \n')
        
    elseif strcmp('SHRUG_2kg',task)
        fileID=fopen(file(16,:));
        CMCshrug21=textscan(fileID,format_CMC,'headerlines',7);
        CMCshrug21=cell2mat(CMCshrug21);
        fclose(fileID);
        fileID=fopen(fullfile(path_to_repo, 'ExperimentalData/EMG/shrug21.exp'));
        EMGshrug21=textscan(fileID,format_EMG,'headerlines',9);
        EMGshrug21=cell2mat(EMGshrug21);
        fclose(fileID);
        
        fileID=fopen(file(17,:));
        CMCshrug22=textscan(fileID,format_CMC,'headerlines',7);
        CMCshrug22=cell2mat(CMCshrug22);
        fclose(fileID);
        fileID=fopen(fullfile(path_to_repo, 'ExperimentalData/EMG/shrug22.exp'));
        EMGshrug22=textscan(fileID,format_EMG,'headerlines',9);
        EMGshrug22=cell2mat(EMGshrug22);
        fclose(fileID);
        
        fileID=fopen(file(18,:));
        CMCshrug23=textscan(fileID,format_CMC,'headerlines',7);
        CMCshrug23=cell2mat(CMCshrug23);
        fclose(fileID);
        fileID=fopen(fullfile(path_to_repo, 'ExperimentalData/EMG/shrug23.exp'));
        EMGshrug23=textscan(fileID,format_EMG,'headerlines',9);
        EMGshrug23=cell2mat(EMGshrug23);
        fclose(fileID);
        
        % load indexes of relevant phases of the motion from file
        indsshrug21=index_matrix(16,1); indendshrug21=index_matrix(16,2); indmaxshrug21=index_matrix(16,3);
        indsshrug22=index_matrix(17,1); indendshrug22=index_matrix(17,2); indmaxshrug22=index_matrix(17,3);
        indsshrug23=index_matrix(18,1); indendshrug23=index_matrix(18,2); indmaxshrug23=index_matrix(18,3);

        % average the estimates for the serratus muscle across bundles
        % (only for CMC, RMR results have been averaged and saved during
        % pre-processing)
        CMCshrug21(:,89)=(CMCshrug21(:,40)+CMCshrug21(:,41)+CMCshrug21(:,42))/3;
        CMCshrug22(:,89)=(CMCshrug22(:,40)+CMCshrug22(:,41)+CMCshrug22(:,42))/3;
        CMCshrug23(:,89)=(CMCshrug23(:,40)+CMCshrug23(:,41)+CMCshrug23(:,42))/3;

        % plot and perform calculations
        [MAE_shrug2kg_CMC, MAE_shrug2kg_RMR_GH, MAE_shrug2kg_RMR_noGH, corr_shrug2kg_CMC, corr_shrug2kg_RMR_GH, corr_shrug2kg_RMR_noGH]=EMG_CMC_RMR_plot(CMCshrug21,CMCshrug22,CMCshrug23,EMGshrug21,EMGshrug22,EMGshrug23,indsshrug21,indendshrug21,indmaxshrug21,indsshrug22,indendshrug22,indmaxshrug22,indsshrug23,indendshrug23,indmaxshrug23,muscle_CMC(i),muscle_EMG(i),MVC(i),muscle_Name(i), path_RMR_GH, path_RMR_noGH, muscle_RMR(i), task);
        disp('MAE_shrug2kg')
        disp(muscle_Name(i))
        disp(strcat('CMC:', {'          '}, string(round(MAE_shrug2kg_CMC,2)), ' [MAE] ,    ', string(round(corr_shrug2kg_CMC,2)), ' [xcorr]'))
        if path_RMR_GH
            disp(strcat('RMR (GH):', {'     '}, string(round(MAE_shrug2kg_RMR_GH,2)), ' [MAE] ,    ', string(round(corr_shrug2kg_RMR_GH,2)), ' [xcorr]'))
        end
        if path_RMR_noGH
            disp(strcat('RMR (no GH):', {'  '}, string(round(MAE_shrug2kg_RMR_noGH,2)), ' [MAE] ,    ', string(round(corr_shrug2kg_RMR_noGH,2)), ' [xcorr]'))
        end
        fprintf('\n \n')
    end

end