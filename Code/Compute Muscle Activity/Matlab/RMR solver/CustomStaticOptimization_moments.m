% CustomStaticOptimization
% ------------------------
%   This script provides a framework for you to build your own custom code to 
%   solve the static optimization problem. Visit the companion Confluence page
%   for more information and suggestions on how to complete this code: 
%   simtk-confluence.stanford.edu/display/OpenSim/Custom+Static+Optimization+in+MATLAB

%-----------------------------------------------------------------------%
% The OpenSim API is a toolkit for musculoskeletal modeling and         %
% simulation. See http://opensim.stanford.edu and the NOTICE file       %
% for more information. OpenSim is developed at Stanford University     %
% and supported by the US National Institutes of Health (U54 GM072970,  %
% R24 HD065690) and by DARPA through the Warrior Web program.           %
%                                                                       %
% Copyright (c) 2020 Stanford University and the Authors                %
% Author(s): Nick Bianco                                                %
%                                                                       %
% Licensed under the Apache License, Version 2.0 (the "License");       %
% you may not use this file except in compliance with the License.      %
% You may obtain a copy of the License at                               %
% http://www.apache.org/licenses/LICENSE-2.0.                           %
%                                                                       %
% Unless required by applicable law or agreed to in writing, software   %
% distributed under the License is distributed on an "AS IS" BASIS,     %
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or       %
% implied. See the License for the specific language governing          %
% permissions and limitations under the License.                        %
%-----------------------------------------------------------------------%

close all; clear all; clc; beep off;

%% Import the OpenSim libraries.
import org.opensim.modeling.*;

%% SETTINGS from the user (check the every time!)

subject_considered = 'ACC';

force_value = 60;
force_direction = '+X';

motion_file_name = 'SO_motion.mot';
externalLoadFile = 'SO_external_load';
genForcesFile = 'SO_generalized_forces.sto';
print_flag = true;
withviz = false;

weight_coord = 8e-3;    % weight used in the tracking of the scapula

%% Set the correct paths
% set the path current folder to be the one where this script is contained
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

root_dir = fullfile('');

% getting path to other folders in this repo
addpath(pathstr)
addpath('..\..\Data Processing\Matlab\')
cd ..\..\..\
path_to_repo = pwd;
addpath(path_to_repo)

% cd to Personal Results to have all the results saved there
cd([path_to_repo, '\Personal_Results']);

%% Load model and get initial state 
% Select the model you want to use, among the ones available in the
% OpenSim Models folder of this repository
[modelFile, modelPath] = uigetfile('*.osim', 'Select model to be processed (should have GLENOID markers already!)');
model = Model([modelPath,modelFile]);
state = model.initSystem();

%% Select the trc file to be considered from Clark's dataset

[trc_file_used, trcPath] = uigetfile('*.trc', 'Select processed TRC file in PTbot\DataShare\dataset_Clark\MoCap\TRC files');

[markersExp, timesExp, labelsExp, unitsExp] = readTRC([trcPath, trc_file_used]);
start_time = timesExp(1);
end_time =  timesExp(end);

frequency_trc_data = 1/(timesExp(2)-timesExp(1));

%% Getting quantities about GelnoHumeral joint
[maxAngle, ~] = get_glenoid_status(model, state);

% get the glenohumeral joint
alljoints = model.getJointSet;
glen = alljoints.get('GlenoHumeral');

%% getting the values of default scapula coordinate 
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

%% Creating external load to apply
% getting markers on the hand body
mcp5_traj = markersExp(:,52:54);
mcp2_traj = markersExp(:,55:57);

% defining the application point as the center of the two markers on the
% hand
application_point = (mcp2_traj + mcp5_traj)/2;
body = 'hand';
generate_external_loads(force_value, force_direction, application_point, body, frequency_trc_data, 'Handle_Load', externalLoadFile);

%% performing IK
% perform IK on the basis of marker data to retrieve the motion file for
% the coordinates of the model
ikSetupFile = [path_to_repo,'' ...
        '\Code\Scale Model\Utilities\setup_files\Clark\IKSetup_Clark.xml'];

ikTool = InverseKinematicsTool(ikSetupFile);
ikTool.setMarkerDataFileName([trcPath, trc_file_used]);
ikTool.setOutputMotionFileName([path_to_repo, '\Personal_Results\', motion_file_name]);
ikTool.set_report_marker_locations(1);
ikTool.setStartTime(start_time);
ikTool.setEndTime(end_time);
ikTool.print('SO_autogenerated_IK_setup.xml');
ikTool.setModel(model);

% set the reference values for the scapula coordinates (last 4 tasks)
num_IK_tasks = ikTool.getIKTaskSet.getSize();

% set the weight of each coordinate in the tracking tasks
ikTool.getIKTaskSet.get(num_IK_tasks-4).setWeight(weight_coord);
ikTool.getIKTaskSet.get(num_IK_tasks-3).setWeight(weight_coord);
ikTool.getIKTaskSet.get(num_IK_tasks-2).setWeight(weight_coord);
ikTool.getIKTaskSet.get(num_IK_tasks-1).setWeight(weight_coord);

% set also the values here
IKCoordinateTask.safeDownCast(ikTool.getIKTaskSet.get(num_IK_tasks-4)).setValue(default_sa);
IKCoordinateTask.safeDownCast(ikTool.getIKTaskSet.get(num_IK_tasks-3)).setValue(default_se);
IKCoordinateTask.safeDownCast(ikTool.getIKTaskSet.get(num_IK_tasks-2)).setValue(default_su);
IKCoordinateTask.safeDownCast(ikTool.getIKTaskSet.get(num_IK_tasks-1)).setValue(default_sw);

ikTool.run();

%% Use OpenSim tools to generate data for optimization.
% Use the AnalyzeTool to compute coordinate speeds.
fprintf('Running kinematics analysis...\n\n');
analyze = AnalyzeTool(model);
analyze.setName('analyze');
analyze.setCoordinatesFileName(motion_file_name);
analyze.loadStatesFromFile(state);
analyze.setStartTime(start_time);
analyze.setFinalTime(end_time);
analysisSet = analyze.getAnalysisSet();
kinematicsAnalysis = Kinematics();
kinematicsAnalysis.setInDegrees(false);
analysisSet.cloneAndAppend(kinematicsAnalysis);
analyze.addAnalysisSetToModel();
analyze.run();
analyze.print('SO_analyzeTool.xml');

%% Run inverse dynamics tool to compute joint moments.

% create an external load with custom function
% TODO

% set up and run inverse dynamics tool to get the generalized
% forces at every joint
fprintf('Running inverse dynamics...\n\n');
idtool = InverseDynamicsTool();        
idtool.setModel(model);
idtool.setStartTime(start_time);
idtool.setEndTime(end_time);
idtool.setCoordinatesFileName(motion_file_name);
idtool.setOutputGenForceFileName(genForcesFile);
idtool.setExternalLoadsFileName([externalLoadFile, '.xml']);
excludedForces = ArrayStr();
excludedForces.append('muscles');
idtool.setExcludedForces(excludedForces);
idtool.run();
idtool.print('SO_id_settings.xml');

%% Load data into MATLAB arrays.
% Use the loadFilterCropArray() function provided by OpenSim Tutorial to load the 
% coordinate kinematic and generalized force data into MATLAB arrays. This 
% function also filters and crops the loaded array based on its two input 
% arguments (more details in loadFilterCropArray.m).
lowpassFreq = 6.0; % Hz
timeRange = [start_time end_time];
[coordinates, coordNames, time] = ...
        loadFilterCropArray('analyze_Kinematics_q.sto', lowpassFreq, timeRange);
[speeds, speedNames, ~] = ...
        loadFilterCropArray('analyze_Kinematics_u.sto', lowpassFreq, timeRange);

[genForces, forceNames, ~] = ...
        loadFilterCropArray(genForcesFile, lowpassFreq, timeRange);

%% Remove unneeded generalized force data columns
% We create an array to remove all generalized forces whose names contain each substring. 
% 'ground' and 'elbow' are removed, since the model does not have enough
% muscle to match the required generalized forces.
forcesToRemove = {'ground', 'elbow'};

% we build an array of indicies to remove all columns at once (see below)
colsToRemove = [];
for i = 1:length(coordNames)
    for j = 1:length(forcesToRemove)
        if contains(coordNames{i}, forcesToRemove{j})
            colsToRemove = [colsToRemove i];
        end
    end
end

% Store an array contain all the original coordinate names of the model
coordNamesAll = coordNames;

% Remove columns and associated coordinate labels.
genForces(:, colsToRemove) = [];  
coordNames(colsToRemove) = [];

%% Store max isometric force values and disable muscle dynamics
muscles = model.getMuscles();

Fmax = zeros(muscles.getSize(), 1);

for i = 1:muscles.getSize()
   % Downcast base muscle to Millard2012EquilibriumMuscle
   muscle = Millard2012EquilibriumMuscle.safeDownCast(muscles.get(i-1));
   Fmax(i) = muscle.get_max_isometric_force();
   muscle.set_ignore_tendon_compliance(true);
   muscle.set_ignore_activation_dynamics(true);
end

if (withviz == true)
    model.setUseVisualizer(true);
end

% Update the system to include any muscle modeling changes
state = model.initSystem();

%% Get coordinate actuators 
% (since ID does not consider the presence of the AC
% constraint, we need this)
allActs = model.getActuators; 
num_acts = getSize(allActs); 
act_Names = cell(num_acts,1); 
acts = cell(num_acts,1);

for i = 0:num_acts-1
    act_Names{i+1} = char(allActs.get(i)) ;
    acts(i+1) = ScalarActuator.safeDownCast(allActs.get(i));
end

%% Perform static optimization.
% We use FMINCON to solve the static optimization problem at selected time points. 
% We set the 'timeInterval' variable to select the time points to be included in the 
% optimization. For example, if set to 10, every 10th time point is selected. A 
% time interval of 1 will select all available time points.
% Clark's marker data is recorded at 50 Hz. 
timeInterval = 10;

% Update data arrays based on the time interval.
N = size(coordinates, 1);
coordinates = coordinates(1:timeInterval:N, :);
speeds = speeds(1:timeInterval:N, :);
genForces = genForces(1:timeInterval:N, :);
numTimePoints = size(coordinates, 1);
% numTimePoints = 50/timeInterval * 2;      % considering just the first 2 seconds of the experiment

% Create the FMINCON options structure.
options = optimoptions('fmincon','Display','notify-detailed', ...
     'TolCon',1e-4,'TolFun',1e-12,'TolX',1e-8,'MaxFunEvals',100000, ...
     'MaxIter',5000,'Algorithm','active-set');
 
% Construct initial guess and bounds arrays
numCoords = length(coordNames);
numMuscles = muscles.getSize();
lb = [zeros(1,numMuscles), -100*ones(1,numCoords)];
ub = [ones(1,numMuscles), 100*ones(1,numCoords)];
x0 = zeros(1,numMuscles+numCoords);

% We define the activation squared cost as a MATLAB anonymous function
w = [ones(1,numMuscles), 100*ones(1,numCoords)];
cost =@(x) sum(w.*(x.^2));

% Pre-allocate arrays to hold muscle moment arms and multipliers for each time 
% point.
fl = zeros(1, numMuscles);
fv = zeros(1, numMuscles);
fp = zeros(1, numMuscles);
cosPenn = zeros(1, numMuscles);
r = zeros(length(coordNames), numMuscles);

% Pre-allocate an array to hold the solution for all time points.
xsol = zeros(numTimePoints, length(x0));

% store also the computed moments
Moments_SO = zeros(size(genForces));

coords = model.getCoordinateSet();
for i = 1:numTimePoints
    fprintf('Optimizing...time step %i/%i \n',  i, numTimePoints);
    
    % at first, let the actuators use their computed force
    for index = 1:num_acts
         acts{index}.overrideActuation(state, false); 
    end
    
    % Loop through model coordinates to set coordinate values and speeds. We set
    % all coordinates to make sure we have the correct kinematic state when 
    % compute muscle multipliers and moment arms.
    for j = 1:length(coordNamesAll)
        coord = coords.get(coordNamesAll{j});
        coord.setValue(state, coordinates(i,j), false);
        coord.setSpeedValue(state, speeds(i,j));
    end
    
    % Populate the muscle multiplier and moment arm arrays. We must realize the
    % system to the Velocity stage. Here, we only need to loop through the
    % coordinates associated with the generalized forces 
    model.realizeVelocity(state);

    % equilibrate the muscles to make them start in the correct state
    model.equilibrateMuscles(state);

    for j = 1:length(coordNames)
        coord = coords.get(coordNames{j});
        for k = 1:numMuscles
            muscle = Millard2012EquilibriumMuscle.safeDownCast(muscles.get(k-1));
            r(j,k) = muscle.computeMomentArm(state, coord);
            if j == 1
                fl(k) = muscle.getActiveForceLengthMultiplier(state);            
                fv(k) = muscle.getForceVelocityMultiplier(state);
                fp(k) = muscle.getPassiveForceMultiplier(state);
                cosPenn(k) = muscle.getCosPennationAngle(state);
            end
        end 
    end
     
    % Create the linear equality constraint arrays, Aeq*x = Beq.
    % 
    %  Aeq: [numCoords x numMuscles+numCoords]
    %    x: [numMuscles+numCoords x 1]
    %  Beq: [numCoords x 1]
   
    % Compute coefficient matrix (Aeq).
    Amusc = bsxfun(@times, r, fl.*fv.*Fmax'.*cosPenn);
    Ares = eye(numCoords);
    Aeq = [Amusc Ares];
    
    % Compute right-hand side (Beq).
    M = genForces(i, :)';
    Mpass = r*(Fmax.*fp'.*cosPenn');
    Beq = M - Mpass;

    % store the values of active and passive maximum force in the current
    % configuration
    AMuscForce = (fl.*fv.*Fmax'.*cosPenn)'; 
    PMuscForce = (Fmax.*fp'.*cosPenn'); 

    % create a struct containing relevant information to be passed to the
    % function evaluating the nonlinear joint reaction constraint
    % It is much faster to pass these as parameters than to evaluate inside
    % the constraint expression every time
    params.nMuscles = numMuscles; 
    params.state = state;
    params.model = model; 
    params.maxAngle = maxAngle;
    params.AMuscForce = AMuscForce;
    params.PMuscForce = PMuscForce;
    params.acts = acts;
    params.num_acts = num_acts;
    params.glen = glen;
    
    % Call FMINCON to solve the problem
    x = fmincon(cost, x0, [], [], Aeq, Beq, lb, ub, @(x)jntrxncon(x,params), options);
%     x = fmincon(cost, x0, [], [], Aeq, Beq, lb, ub, [], options);   
         
    % Store solution and set guess for next time point
    xsol(i, :) = x;
    x0 = x;   

    % calculate and store the generalized forces produced by SO
    Moments_SO(i, :) = Mpass + Aeq*x';

    if (withviz == true)
        model.getVisualizer.show(state);
    end
end

%% Plot results
% According to the value of the 'print_flag', we print the muscle
% activations that are found

if print_flag
    figure;
    title("Muscle Activations")
    muscleNames = ArrayStr();
    muscles.getNames(muscleNames);
    pgc = linspace(0, 100, size(xsol,1));
    for i = 1:numMuscles
       subplot(5,8,i)
       plot(pgc,xsol(:,i),'b-')
       ylim([0 1])
       muscName = muscleNames.get(i-1).toCharArray';
       title(muscName(1:end), 'interpreter', 'none')
    end
    
    % Plot reserve actuator excitations.
    figure;
    title("Reserve actuators")
    side = ceil(sqrt(length(coordNames)));
    for i = 1:length(coordNames)
        subplot(side,side,i)
        plot(xsol(:,numMuscles+i), 'linewidth', 2);
        title(coordNames{i});
    end
end

%% SAVING THE RESULTS TO FILE
name_file = append('muscle_activations_', subject_considered, '_experiment_', trc_file_used(1));
muscle_order = "";
for i = 1:numMuscles
    muscle_order= [muscle_order, string(muscleNames.get(i-1).toCharArray')];
end

for i=1:length(coordNames)
    muscle_order= [muscle_order, string(coordNames{i})];
end

muscle_order = muscle_order(2:end);

frequency_solution = 50/timeInterval;

save(name_file, 'xsol', 'muscle_order', 'frequency_solution');