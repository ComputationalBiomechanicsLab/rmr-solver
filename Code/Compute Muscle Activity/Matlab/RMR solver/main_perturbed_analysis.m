% Script to run the Rapid Muscle Redundancy (RMR) solver on user-selected
% experiments, considering multiple models.
% The user is prompted with the selection of the tasks to analyze.
% Within the script, it is possible to adjust the downsampling to be
% applied, and whether the analysis should include the glenohumeral
% constraint or not.
% The same analysis is performed multiple times, if various models are
% selected.
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
cd ..\..\..\..\
path_to_repo = pwd;
addpath(path_to_repo)
addpath(fullfile(path_to_repo, 'Code\Data Processing\Matlab\'))

%% Parameters
% where you have the experimental files (.trc)
trc_path = fullfile(path_to_repo, 'ExperimentalData\Markers');

% where you have the models (.osim)
model_path = fullfile(path_to_repo, 'OpenSim Models\for RMR solver\perturbed_100\');

% where to save the results
saving_path = fullfile(path_to_repo, 'Personal_Results');

% use parallel computing (requires Parallel Computing Toolbox)
parallel_comp = true;

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

% Flags (Select whether to enforce constraint 3 and 4 from the formulation
% reported in the paper)
dynamic_bounds = true;              % enforcing continuity of the activations from one timestep to the next, to respect first-order dynamics
enforce_GH_constraint = true;       % enforcing directional constraint on the glenohumeral joint force

%% Script
% Select models
[model_files, model_path] = uigetfile('*.osim', 'Select the models to consider', model_path, 'MultiSelect','on');

if iscell(model_files)
    num_models = size(model_files, 2);
else
    num_models = 1;
    model_files = {model_files};
end

% Select the experimental data to be considered
dataset_general_name = 'Seth2019';

[files,path] = uigetfile('*.trc', 'Select the .trc files to analyse', trc_path, 'MultiSelect','on');

if iscell(files)
    num_files = size(files, 2);
else
    num_files = 1;
end

%% Run Rapid Muscle Redundancy (RMR) solver
% preallocating arrays to hold information about the solutions
optimizationStatus = [];
unfeasibility_flag = [];
tOptim = zeros(num_files*num_models,1);
result_file_RMR = {};

% adjust name of the results according to choice of GH constraint
if enforce_GH_constraint
    choice_GH = 'withGH';
else
    choice_GH = 'noGH';
end

tic

fprintf('Starting RMR analysis\n\n')

if parallel_comp
    parfor model_index=1:num_models   % for parallel computing (requires multiple cores on your machine)
        fprintf("-model %i\n", model_index)
        model = Model(fullfile(model_path, model_files{model_index}));
        info_numModel_weight = extractBefore(extractAfter(model_files{model_index}, 'subject_'), '.osim');      % based on file name
        dataset_considered = append(dataset_general_name, '_', info_numModel_weight, '_', choice_GH);
    
        for trc_file_index=1:num_files
            fprintf('    experiment %i', trc_file_index)
            if num_files>1
                experiment = append(path, files(trc_file_index));
                experiment = experiment{1};
            else
                experiment = append(path,files);
            end
        
            [aux_optimization_status, aux_unfeasibility_flags, aux_tOptim, aux_result_file] = RMR_analysis(dataset_considered, model, experiment, 0, weight_coord, time_interval, dynamic_bounds, enforce_GH_constraint, saving_path);
            
            fprintf(' ... solved with %i unfeasible solutions \n', sum(aux_unfeasibility_flags));
        end
    end
else
    for model_index=1:num_models   % normal sequential analysis
        fprintf("-model %i\n", model_index)
        model = Model(fullfile(model_path, model_files{model_index}));
        info_numModel_weight = extractBefore(extractAfter(model_files{model_index}, 'subject_'), '.osim');      % based on file name
        dataset_considered = append(dataset_general_name, '_', info_numModel_weight, '_', choice_GH);
    
        for trc_file_index=1:num_files
            fprintf('    experiment %i', trc_file_index)
            if num_files>1
                experiment = append(path, files(trc_file_index));
                experiment = experiment{1};
            else
                experiment = append(path,files);
            end
        
            [aux_optimization_status, aux_unfeasibility_flags, aux_tOptim, aux_result_file] = RMR_analysis(dataset_considered, model, experiment, 0, weight_coord, time_interval, dynamic_bounds, enforce_GH_constraint, saving_path);
            
            tOptim(num_files*(model_index-1)+trc_file_index) = aux_tOptim;
            optimizationStatus(num_files*(model_index-1)+trc_file_index).experiment = aux_optimization_status;
            unfeasibility_flag(num_files*(model_index-1)+trc_file_index).experiment = aux_unfeasibility_flags;
            result_file_RMR{num_files*(model_index-1)+trc_file_index} = aux_result_file;
    
            fprintf(' ... solved with %i unfeasible solutions \n', sum(aux_unfeasibility_flags));
        end
    end
end

t_required = toc;