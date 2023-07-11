% This script generates a set of models. Each is identical to the model
% specified by the user, while the markers have been randomly perturbed from
% their original locations, up to a specified maximum distance.
% The user can modify the "Parameters" section to adjust the scripts to
% their needs.
%
% This code has been developed based on the freely available codes at 
% https://simtk.org/projects/quant_uncertain
%
% In our pipeline, we first ran this script with the load-free model as an
% input, and then used the "add_2kgLoad_toModels.m" to programmatically add
% the load to all of the models.

clear; clc;

% Import the OpenSim libraries.
import org.opensim.modeling.*;

% set the path current folder to be the one where this script is contained
mfile_name          = mfilename('fullpath');
[pathstr,~,~]  = fileparts(mfile_name);
cd(pathstr);

% getting path to other folders in this repo
addpath(pathstr)
cd ..\..\

path_to_repo = pwd;
addpath(path_to_repo)
addpath(fullfile(path_to_repo, 'Code\Data Processing\'))

%% Parameters
% select file name of the original model
[inputFileName, path_model] = uigetfile('*.osim', 'Select the OpenSim model to perturb', path_to_repo, 'MultiSelect','off');

% Number of models to generate (set to 100 in the paper's pipeline)
numModels = 100; % 100

% Maximum distance by which each marker will be perturbed, in meters
maxRadius = 0.01; % 0.01

% Specify markers that should not be perturbed
unperturbed_markers = ["Glenoid_Center", "HumHead_Center", "Glenoid_Edge"];  % we retain these markers as they are virtual markers defining the glenohumeral constraint

% location and name of the IK setup file used to assess validity of the 
% perturbed models
path_IK_setup = fullfile(path_to_repo, 'ExperimentalData/IK setup files/');
IK_setup_file = 'IKSetup_2019.xml';

% maximum value of RMSE allowable for a perturbed model to be considered valid
% The rmse_on_data is calculated as an average of the RMSEs resulting from
% running IK on each of the individual experimental recordings
max_rmse_on_data = 0.01;

% dataset of .trc files used to validate the perturbed models with IK.
% The perturbed model will be accepted only if the average RMSE resulting
% from IK is lower or equal to max_rmse_on_data
trc_dataset_path = fullfile(path_to_repo, 'ExperimentalData/Markers');

% set a very high duration of the experiments, so we don't need to load all
% the files to get the real value out of them (IK would work up until there
% is data)
end_time = 100;  % in seconds
start_time = 0;  % time when the IK is started, also in seconds

% Flag to print each model number during the loop (one line printed per model)
printModelInfo = true;

% Flag to report detailed information about the perturbation applied to each
% marker (for debugging purposes). WARNING: clutters the console if many markers
% or models are being processed
printMarkerDetail = false;

% where to save the results
saving_path = fullfile(path_to_repo, 'Personal_Results');

%% Script ---------------------------------------------------------------
% load model and connect all components
model = Model(fullfile(path_model, inputFileName));
model.finalizeConnections();

% PRECOMPUTING USEFUL VARIABLES FOR INVERSE KINEMATICS
% load the names of the trc files composing the dataset, so that we can
% then run IK to check whether the perturbed models we obtain are
% performing well enough
trc_dataset = dir(fullfile(trc_dataset_path, '*.trc'));
trc_file_names = {trc_dataset.name}';
dim_dataset = size(trc_dataset, 1);

% Set the weight for the various scapula coordinates in IK
% This is to achieve a good agreement between scapula upward rotation and
% shoulder elevation (as reported in the paper)
weight_abd = 0.0001;
weight_elev = 0.0001;
weight_up_rot = 0.0002;
weigth_wing = 0.0001;
weight_coord = [weight_abd, weight_elev, weight_up_rot, weigth_wing];

% getting the values of default scapula coordinate 
% we get the values of the coordinates describing the scapula position from 
% the general model in default pose
scapula_abd = model.getJointSet().get(2).get_coordinates(0);
scapula_ele = model.getJointSet().get(2).get_coordinates(1);
scapula_urt = model.getJointSet().get(2).get_coordinates(2);
scapula_wng = model.getJointSet().get(2).get_coordinates(3);

default_sa = scapula_abd.get_default_value();
default_se = scapula_ele.get_default_value();
default_su = scapula_urt.get_default_value();
default_sw = scapula_wng.get_default_value();

% GENERATE PERTURBED MODELS
for iter=1:numModels
    attempt = 1;
    aux_iter = iter;
    while aux_iter==iter
        if printModelInfo
            fprintf("-model %i (attempt %i)\n", iter, attempt)
        end
    
        modelPerturbed = Model(fullfile(path_model, inputFileName));
        markerset = modelPerturbed.updMarkerSet();
        numMarkers = markerset.getSize();
    
        for mrk_id=0:numMarkers-1
            % load the marker considered at this iteration
            marker=markerset.get(mrk_id);
    
            % move just the right markers
            if ~ismember(string(marker.getName()), unperturbed_markers)
    
                % get original location
                oriloc = marker.get_location().getAsMat;
                if printMarkerDetail
                    fprintf('  -> marker %s \n', string(marker.getName()))
                    fprintf('     maximum radius: %4f [m] \n', maxRadius)
                    fprintf('     original location [%4f, %4f, %4f] \n', oriloc(1), oriloc(2), oriloc(3))
                end
                
                % Perturb the marker location within a sphere of radius
                % maxRadius. The location is adjusted until is satisfying the
                % maximum allowed perturbation bounds
                dloc = maxRadius * ones(3,1);
                while norm(dloc)>maxRadius
                    dloc = [-maxRadius + 2*maxRadius*rand(1); ...
                            -maxRadius + 2*maxRadius*rand(1); ...
                            -maxRadius + 2*maxRadius*rand(1)];
                end
    
                % update the marker location to be the new one
                newloc = [oriloc(1) + dloc(1), oriloc(2) + dloc(2), oriloc(3) + dloc(3)];
                marker.set_location(Vec3(newloc(1), newloc(2), newloc(3)));
    
                delta_loc = norm(dloc);
                if printMarkerDetail
                    fprintf('     new location: [%4f, %4f, %4f] \n', newloc(1), newloc(2), newloc(3))
                    fprintf('     distance moved: %4f \n', delta_loc)
                end
            end
        end
    
        % connect all the components in the perturbed model
        modelPerturbed.finalizeConnections();
        
        % RUN IK ON THE WHOLE DATASET
        % setting up a common IK tool
        ikTool = InverseKinematicsTool(fullfile(path_IK_setup, IK_setup_file));
        ikTool.setModel(modelPerturbed);
        ikTool.setStartTime(start_time);
        ikTool.setEndTime(end_time);
    
        % save the marker errors corresponding to the IK solutions
        ikTool.set_report_errors(1)
    
        % set the results directory (for the _ik_marker_errors.sto)
        ikTool.set_results_directory(saving_path)
    
        % set the reference values for the scapula coordinates (last 4 tasks)
        num_IK_tasks = ikTool.getIKTaskSet.getSize();
        
        % set the weight of each coordinate in the tracking tasks
        ikTool.getIKTaskSet.get(num_IK_tasks-4).setWeight(weight_coord(1));
        ikTool.getIKTaskSet.get(num_IK_tasks-3).setWeight(weight_coord(2));
        ikTool.getIKTaskSet.get(num_IK_tasks-2).setWeight(weight_coord(3));
        ikTool.getIKTaskSet.get(num_IK_tasks-1).setWeight(weight_coord(4));
        
        % set also the values here
        IKCoordinateTask.safeDownCast(ikTool.getIKTaskSet.get(num_IK_tasks-4)).setValue(default_sa);
        IKCoordinateTask.safeDownCast(ikTool.getIKTaskSet.get(num_IK_tasks-3)).setValue(default_se);
        IKCoordinateTask.safeDownCast(ikTool.getIKTaskSet.get(num_IK_tasks-2)).setValue(default_su);
        IKCoordinateTask.safeDownCast(ikTool.getIKTaskSet.get(num_IK_tasks-1)).setValue(default_sw);
    
        % check if the model produces realistic RMSE when the IK is run
        % on all the dataset that we have
        rmse_on_data = zeros(dim_dataset, 1);
        
        for trc_file_index=1:dim_dataset
            
            if printModelInfo
                % print the experiment considered
                fprintf('     ... running IK and checking error for %s   ', trc_file_names{trc_file_index})
            end
            % name of the motion file that will be produced
            motion_file_name = append(trc_file_names(trc_file_index), '.mot');
    
            % set the .trc file to consider
            ikTool.setMarkerDataFileName(fullfile(trc_dataset_path, trc_file_names{trc_file_index}));
    
            % set the name and location of the file where the resulting
            % joint angles will be saved
            ikTool.setOutputMotionFileName(fullfile(saving_path, motion_file_name));
    
            % print the complete autogenerated IK setup file
            ikTool.print(fullfile(saving_path, 'RMR_autogenerated_IK_setup.xml'));
    
            % run IK
            ikTool.run();
    
            % load the resulting RMSE from output file of IK
            [labels, data] = readStoFile(fullfile(saving_path, '_ik_marker_errors.sto'));
            
            % retain the 'marker_error_RMS', 3rd column of the data
            rmse_on_trial = data(:,3);
            
            % save the RMSE for this particular experimental trial
            rmse_on_data(trc_file_index) = mean(rmse_on_trial);
            
            if printModelInfo
                % print the error value 
                fprintf('         (error is %4f [m]) \n', mean(rmse_on_trial))
            end
        end
    
        % check whether the perturbed model produced satisfactory
        % results in IK (in terms of cumulative RMSE). If not, discard
        % the model and generate a new one.
        if mean(rmse_on_data)<=max_rmse_on_data
            outputFileName = append(inputFileName(1:end-5), string(iter), '.osim');
            modelPerturbed.print(fullfile(saving_path, outputFileName));
            
            % we can proceed to the generation of the next model
            aux_iter = aux_iter +1;
        else
            % repeat the generation of this model, as the IK error was too
            % high
            attempt = attempt+1;
        end
        
        if printModelInfo
            % print info regarding the observed RMSE
            fprintf('     Cumulative RMSE = %f \n\n', mean(rmse_on_data))
        end
    end

    if printModelInfo
        fprintf('\n\n')
    end
end

