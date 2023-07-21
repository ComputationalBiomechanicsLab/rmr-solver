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

% select the folder of the Leg6Dof9Musc model (shipped with OpenSim)
% Installation of OpenSim is in the C drive normally
model_path = uigetdir('C:\', 'Select the folder containing the Gait2354 example');

% get the motion file from the Stance experiment
motion_file = fullfile(model_path, 'OutputReference', 'subject01_walk1_ik.mot');

% where to save the results
saving_path = fullfile(path_to_repo, 'Personal_Results');


% Select model
modelFile = fullfile(path_to_repo, 'OpenSim Models', 'Gait2354', 'subject01_simbody_adjusted_forRMR.osim');
model = Model(modelFile);

% Select the experimental data to be considered
dataset_considered = 'GaitExample';

% Downsampling
time_interval = 1;

% Flags (Select whether to enforce constraints)
dynamic_bounds = true;               % enforcing continuity of the activations from one timestep to the next, to respect first-order dynamics
enforce_GH_constraint = false;       % enforcing directional constraint on the glenohumeral joint force
apply_external_force = 0;

%% Generate the external force and add it to the model
force_params =[];
force_params.apply_external_force = apply_external_force;
force_1 = [];
if apply_external_force
    external_force_filename = 'subject01_walk1_grf.xml';             % name of the filename in which the force is going to be stored
    
    data_storage = Storage(fullfile(model_path, 'subject01_walk1_grf.mot'));
    external_force = ExternalForce(data_storage, "ground_force_v", "ground_force_p", "ground_torque_", "calcn_r", "ground", "ground");
    
    external_force.print(external_force_filename)

    % Save external force parameters in structure
    force_1.ef_filename = external_force_filename;
    force_1.ef_storage = data_storage;
    force_1.ef = external_force;
end

force_2 = [];
if apply_external_force
    external_force_filename = 'subject01_walk1_grf.xml';             % name of the filename in which the force is going to be stored
    
    data_storage = Storage(fullfile(model_path, 'subject01_walk1_grf.mot'));
    external_force = ExternalForce(data_storage, "l_ground_force_v", "l_ground_force_p", "l_ground_torque_", "calcn_l", "ground", "ground");
    
    external_force.print(external_force_filename)

    % Save external force parameters in structure
    force_2.ef_filename = external_force_filename;
    force_2.ef_storage = data_storage;
    force_2.ef = external_force;
end

force_params.num_forces = 2;
force_params.force{1} = force_1;
force_params.force{2} = force_2;

%% Run Rapid Muscle Redundancy (RMR) solver
disp('Running RMR')

[optimization_status, unfeasibility_flags, tOptim, result_file] = RMR_analysis(dataset_considered, model, 0, motion_file, [], time_interval, dynamic_bounds, enforce_GH_constraint, force_1, saving_path);

fprintf('\n Solved with %i unfeasible solutions \n \n \n', sum(unfeasibility_flags));
