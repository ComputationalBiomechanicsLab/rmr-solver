% Script to run the Rapid Muscle Redundancy (RMR) solver on user-selected experiments.
% The user is prompted with the selection of the tasks to analyze.
% Within the script, it is possible to adjust the downsampling to be
% applied, and whether the analysis should include the glenohumeral
% constraint or not.
%
% Author: Italo Belli (i.belli@tudelft.nl) 2023

close all; clear; clc; beep off;

% Import the OpenSim libraries.
import org.opensim.modeling.*;

% set the path current folder to be the one where this script is contained
mfile_name          = mfilename('fullpath');
[pathstr,~,~]  = fileparts(mfile_name);
cd(pathstr);

% getting path to other folders in this repo
addpath(pathstr)
cd ..\..\..\
path_to_repo = pwd;
addpath(path_to_repo)
addpath(fullfile(path_to_repo, 'Code\Data Processing\'))

% where you have the experimental files (.trc)
trc_path = fullfile(path_to_repo, 'ExperimentalData\Markers');

% where to save the results
saving_path = fullfile(path_to_repo, 'Personal_Results');


% Select model
modelFile_0kg = append(path_to_repo, '\OpenSim Models\for RMR solver\TSM_subject_noWeight.osim');
model_0kg = Model(modelFile_0kg);

modelFile_2kg = append(path_to_repo, '\OpenSim Models\for RMR solver\TSM_subject_2kgWeight.osim');
model_2kg = Model(modelFile_2kg);

% Select the experimental data to be considered
dataset_considered = 'Seth2019';

[files,path] = uigetfile('*.trc', 'Select the .trc files to analyse', trc_path, 'MultiSelect','on');

if iscell(files)
    num_files = size(files, 2);
else
    num_files = 1;
end

% Set the weight for the various scapula coordinates in IK
% This is to achieve a good agreement between scapula upward rotation and
% shoulder elevation (as reported in the paper)
weight_abd = 0.0001;
weight_elev = 0.0001;
weight_up_rot = 0.0002;
weigth_wing = 0.0001;
weight_coord = [weight_abd, weight_elev, weight_up_rot, weigth_wing];

% Downsampling
time_interval = 1;

% Flags (Select whether to enforce constraints)
dynamic_bounds = true;              % enforcing continuity of the activations from one timestep to the next, to respect first-order dynamics
enforce_GH_constraint = true;       % enforcing directional constraint on the glenohumeral joint force

%% Run Rapid Muscle Redundancy (RMR) solver
% preallocating arrays to hold information about the solutions
optimizationStatus = [];
unfeasibility_flag = [];
tOptim = zeros(num_files,1);
result_file_RMR = {};

for trc_file_index=1:num_files
    fprintf('Running RMR on experiment %i \n', trc_file_index)
    if num_files>1
        experiment = append(path, files(trc_file_index));
        experiment = experiment{1};
        has_2kg_weight = str2num(experiment(end-5));      % based on file name
    else
        experiment = append(path,files);
        has_2kg_weight = str2num(experiment(end-5));      % based on file name
    end
    
    % consider the correct model in the analysis, based on the .trc files
    if has_2kg_weight
        [aux_optimization_status, aux_unfeasibility_flags, tOptim(trc_file_index), aux_result_file] = RMR_analysis(dataset_considered, model_2kg, experiment, 0, weight_coord, time_interval, dynamic_bounds, enforce_GH_constraint, saving_path);
    else
        [aux_optimization_status, aux_unfeasibility_flags, tOptim(trc_file_index), aux_result_file] = RMR_analysis(dataset_considered, model_0kg, experiment, 0, weight_coord, time_interval, dynamic_bounds, enforce_GH_constraint, saving_path);
    end
    optimizationStatus(trc_file_index).experiment = aux_optimization_status;
    result_file_RMR{trc_file_index} = aux_result_file;
    unfeasibility_flag(trc_file_index).experiment = aux_unfeasibility_flags;
end
