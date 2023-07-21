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
model_path = uigetdir('C:\', 'Select the folder containing the Leg6Dof9Musc model');

% get the motion file from the Stance experiment
motion_file = fullfile(model_path, 'Stance', 'leg69_IK_stance.mot');

% where to save the results
saving_path = fullfile(path_to_repo, 'Personal_Results');


% Select model
modelFile = fullfile(path_to_repo, 'OpenSim Models', 'Leg6Dof9Musc', 'leg6dof9musc_rmr_Millard2012EqMusc.osim');
model = Model(modelFile);

% Select the experimental data to be considered
dataset_considered = 'Stance';

% Downsampling
time_interval = 1;

% Flags (Select whether to enforce constraints)
dynamic_bounds = true;               % enforcing continuity of the activations from one timestep to the next, to respect first-order dynamics
enforce_GH_constraint = false;       % enforcing directional constraint on the glenohumeral joint force
apply_external_force = 1;

%% Generate the external force and add it to the model
force_params =[];
force_params.apply_external_force = apply_external_force;
force_1 = [];
if apply_external_force
    external_force_filename = 'leg69_right_grf.xml';             % name of the filename in which the force is going to be stored
    
    data_storage = Storage(fullfile(model_path, 'Stance', 'leg69_stance_grf.mot'));
    external_force = ExternalForce(data_storage, "ground_force_v", "ground_force_p", "ground_torque_", "calcn_r", "ground", "ground");
    
    external_force.print(external_force_filename)

    % Save external force parameters in structure
    force_1.ef_filename = external_force_filename;
    force_1.ef_storage = data_storage;
    force_1.ef = external_force;
end

force_params.num_forces = 1;
force_params.forces{1} = force_1;

%% Run Rapid Muscle Redundancy (RMR) solver
disp('Running RMR')

[optimization_status, unfeasibility_flags, tOptim, result_file] = RMR_analysis_extForces(dataset_considered, model, 0, motion_file, [], time_interval, dynamic_bounds, enforce_GH_constraint, force_params, saving_path);

