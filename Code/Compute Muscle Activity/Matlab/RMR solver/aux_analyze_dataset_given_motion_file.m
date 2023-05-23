
close all; clear; clc; beep off;

% Import the OpenSim libraries.
import org.opensim.modeling.*;

% set the path current folder to be the one where this script is contained
mfile_name          = mfilename('fullpath');
[pathstr,~,~]  = fileparts(mfile_name);
cd(pathstr);

% getting path to other folders in this repo
addpath(pathstr)
cd ..\..\..\..\
path_to_repo = pwd;
addpath(path_to_repo)
addpath(fullfile(path_to_repo, 'Code\Data Processing\Matlab\'))

% where to save the results
saving_path = fullfile(path_to_repo, 'Personal_Results');

% Select model
modelFile_0kg = append(path_to_repo, '\Personal_Results\TSM_Ajay2019_noWeight.osim');
model_0kg = Model(modelFile_0kg);

modelFile_2kg = append(path_to_repo, '\Personal_Results\TSM_Ajay2019_2kgWeight.osim');
model_2kg = Model(modelFile_2kg);

% Select the experimental data to be considered
dataset_considered = 'Seth2019';

[files,path] = uigetfile('*.mot', 'Select the .mot files to analyse', path_to_repo, 'MultiSelect','on');

if iscell(files)
    num_files = size(files, 2);
else
    num_files = 1;
end

% Downsampling
time_interval = 1;

% Flags (Select whether to enforce constraint 3 and 4 from the formulation
% reported in the paper)
dynamic_bounds = true;
enforce_GH_constraint = false;

%% Run Muscle Redundancy Solver (MRS)
% preallocating arrays to hold information about the solutions
optimizationStatus = [];
unfeasibility_flag = [];
tOptim = zeros(num_files,1);
result_file_MRS = {};

for motion_file_index=1:num_files
    fprintf('Running MRS on experiment %i \n', motion_file_index)
    if num_files>1
        experiment = files(motion_file_index);
        experiment = experiment{1};
        has_2kg_weight = str2num(experiment(end-5));      % based on file name
    else
        experiment = files;
        has_2kg_weight = str2num(experiment(end-5));      % based on file name
    end
    
    % TODO: find a better way for this!
    has_2kg_weight = 1;

    if has_2kg_weight
        [aux_optimization_status, aux_unfeasibility_flags, tOptim(motion_file_index), aux_result_file] = MRS_analysis(dataset_considered, model_2kg, [], experiment, [], time_interval, dynamic_bounds, enforce_GH_constraint, saving_path);
    else
        [aux_optimization_status, aux_unfeasibility_flags, tOptim(motion_file_index), aux_result_file] = MRS_analysis(dataset_considered, model_0kg, [], experiment, [], time_interval, dynamic_bounds, enforce_GH_constraint, saving_path);
    end
    optimizationStatus(motion_file_index).experiment = aux_optimization_status;
    result_file_MRS{motion_file_index} = aux_result_file;
    unfeasibility_flag(motion_file_index).experiment = aux_unfeasibility_flags;
    fprintf('\n Solved with %i unfeasible solutions \n \n \n', sum(aux_unfeasibility_flags));
end
