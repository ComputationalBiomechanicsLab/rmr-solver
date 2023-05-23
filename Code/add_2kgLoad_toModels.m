% This code is used to add a 2-kg load in hand to selected thoracoscapular
% shoulder models. The original model is modified adding the required mass
% on the hand body, and modifying its inertia accordingly.

%% Import the OpenSim libraries.
import org.opensim.modeling.*;

% set the path current folder to be the one where this script is contained
mfile_name          = mfilename('fullpath');
[pathstr,~,~]  = fileparts(mfile_name);
cd(pathstr);

% getting path to other folders in this repo
addpath(pathstr)
cd ..\
path_to_repo = pwd;
addpath(path_to_repo)
addpath(fullfile(path_to_repo, 'Code\Data Processing\Matlab\'))

%% Parameters
% where to save the results
saving_path = fullfile(path_to_repo, 'Personal_Results');

% mass to add to the hand [kg]
mass = 2; 

% resulting inertia of the hand [kg*m^2]
inertia = Vec6(0.003, 0.00088545, 0.003, 0, 0, 0);   % the hand has been approximated as a parallelepiped

% body to be modified
body = 'hand';

% common name of the new models
common_name = 'TSM_subject_2kgWeight';

%% Script
% select file name of the original model
[unloadedFileNames, path_unloaded_model] = uigetfile('*.osim', 'Select the OpenSim models to add the load on', path_to_repo, 'MultiSelect','on');

num_models  = size(unloadedFileNames,2);
if num_models==1
    unloadedFileNames = {unloadedFileNames};
end

tic
% loop through all the models and modify the required body properties
for index_model = 1:num_models
    model = Model(fullfile(path_unloaded_model, unloadedFileNames{index_model}));
    body_in_model = model.getBodySet().get(body);
    original_mass = body_in_model.getMass();

    % set new mass
    body_in_model.setMass(original_mass+mass);

    % set new inertia
    body_in_model.set_inertia(inertia);
    
    % save the new model
    loaded_model_name = append(common_name, string(index_model), '.osim');
    model.print(fullfile(saving_path, loaded_model_name));
end
toc