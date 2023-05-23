% Custom Static Optimization written on the basis of the official OpenSim example.
% Starting from experimental marker data (in .trc format) the optimal
% muscle activations are found that can reproduce the motion, solving:
%
%   min   sum (w_i * a_i^2)+ sum (w_j * c_j^2)
%   a,c    i                  j                 
%
%   s.t.        0<=a_i<=1            for muscle activations
%              -l<=c_j<=l            for coordinateActuators controls (if present)
%        acc_{j,FD} = acc_{j,data}   constraint on accelerations
%
%
% The code is written specifically to consider a thoracoscapular shoulder
% model that has been already scaled to the biometrics of the subject of
% interest, and the data are loaded assuming to be considering a particular
% dataset (Waterloo's one). However, this script can be generalized to
% consider other models and data without changing its main structure. 
%
% Please check the SETTINGS section, where you are asked if you want to
% include an external force (if so, specify the magnitude, the direction
% with respect to the subject, the body on which it is applied).
%
% author: Italo Belli (i.belli@tudelft.nl) May 2022

close all; clear all; clc; beep off;

%% Import the OpenSim libraries.
import org.opensim.modeling.*;

%% SETTINGS from the user (check the every time!)
subject_considered = 'FLX01_';

apply_external_force = false;
force_value = 20;
force_direction = 'push';
body = 'hand';

motion_file_name = 'flx01_1.mot';
print_flag = true;
withviz = false;

weight_coord = 0.0001; 
% weight_coord = 8e-3;    % weight used in the tracking of the scapula

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
% [modelFile, modelPath] = uigetfile('*.osim', 'Select model to be processed (should have GLENOID markers already!)');
% model = Model([modelPath,modelFile]);
modelFile = 'C:\Users\italobelli\Desktop\GitHub\PTbot_officialCBL\Personal_Results\TSM_2019_glen_noLocked.osim';
model = Model(modelFile);

%% Select the trc file to be considered from Clark's dataset
% [trc_file_used, trcPath] = uigetfile('*.trc', 'Select processed TRC file in PTbot\DataShare\dataset_Clark\MoCap\TRC files');
trcPath = 'C:\Users\italobelli\Desktop\comparing TSM models\FrontiersPubMaterials-latest\ThoracoscapularShoulderPaperMaterials\ThoracoscapularShoulderPaperMaterials\ExperimentalData\Markers\';
trc_file_used = 'FLX01.trc';

[markersExp, timesExp, labelsExp, unitsExp] = readTRC([trcPath, trc_file_used]);
start_time = timesExp(1);
end_time =  timesExp(end);

frequency_trc_data = 1/(timesExp(2)-timesExp(1));

%% Getting quantities about GelnoHumeral joint
model_temp = model.clone();    % create a temporary copy of the model, to be used in the tool

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

%% performing IK
% perform IK on the basis of marker data to retrieve the motion file for
% the coordinates of the model

ikSetupFile = [path_to_repo,'' ...
        '\Code\Scale Model\Utilities\setup_files\Standford2019_experiments\IKSetup_2019.xml'];

ikTool = InverseKinematicsTool(ikSetupFile);
ikTool.setMarkerDataFileName([trcPath, trc_file_used]);
ikTool.setOutputMotionFileName([path_to_repo, '\Personal_Results\', motion_file_name]);
ikTool.set_report_marker_locations(1);
ikTool.setStartTime(start_time);
ikTool.setEndTime(end_time);
ikTool.setModel(model_temp);

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
ikTool.print('SO_autogenerated_IK_setup.xml');

ikTool.run();

model_temp.delete();          % delete the temporary copy. It is just useful to avoid initializing/modifying the main model at this stage   

%% getting the kinematic data that we need
% Use the loadFilterCropArray() function provided by OpenSim Tutorial to load the 
% coordinate kinematic and generalized force data into MATLAB arrays. This 
% function also filters and crops the loaded array based on its two input 
% arguments (more details in loadFilterCropArray.m).
lowpassFreq = 3.0; % Hz
timeRange = [start_time end_time];

% get the coordinates from the output of the IK in rad
[coordinates, coordNames, time] = loadFilterCropArray(motion_file_name, lowpassFreq, timeRange);
coordinates = deg2rad(coordinates);
coordNamesAll = coordNames; % TODO: remove this if it works

% get the velocities for each joint in rad/s
time_step = timesExp(2)-timesExp(1);
speeds = zeros(size(coordinates));
for i=1:size(coordNames,1)
    speeds(:,i) = gradient(coordinates(:,i), time_step);
end
speedNames = coordNames;

% get the accelerations for each coordinate in rad/s^2
accelerations = zeros(size(speeds));
for i=1:size(coordNames,1)
    accelerations(:,i) = gradient(speeds(:,i), time_step);
end
accNames = speedNames;

% check the values of joint states, speeds and accelerations
figure
for i=1:16
subplot(4,4,i)
hold on
plot(coordinates(:,i))
plot(speeds(:,i))
plot(accelerations(:,i))
title(coordNames{i});
hold off
grid on
end
legend("coords", "speeds", "accs")

%% Apply the external force to the model (if needed)
if apply_external_force
    % decode the direction of the force
    if strcmp(force_direction, 'right')
        force_direction_OS = '-X';
    elseif strcmp(force_direction, 'left')
        force_direction_OS = '+X';
    elseif strcmp(force_direction, 'up')
        force_direction_OS = '-Y';
    elseif strcmp(force_direction, 'down')
        force_direction_OS = '+Y';
    elseif strcmp(force_direction, 'push')
        force_direction_OS = '+Z';
    elseif strcmp(force_direction, 'pull')
        force_direction_OS = '-Z';
    end
    
    % getting markers on the hand body
    mcp5_traj = markersExp(:,52:54);
    mcp2_traj = markersExp(:,55:57);
    
    % defining the application point as the center of the two markers on the
    % hand
    application_point = (mcp2_traj + mcp5_traj)/2;
 
    % generate the external force file
    external_force = generate_external_force(force_value, force_direction_OS, application_point, body, frequency_trc_data, 'TSM_external_force');
 
    % apply the external force to the model
    model.addForce(external_force);
    
    % reset the data source file for the force applied
    external_force_storage = Storage('TSM_external_force.mot', false);
    n_forces = model.getForceSet().getSize();
    force_in_model = model.getForceSet().get(n_forces-1);
    ExternalForce.safeDownCast(force_in_model).setDataSource(external_force_storage);
    
    model.finalizeConnections();
end

%% Store max isometric force values and disable muscle dynamics
muscles = model.getMuscles();
numMuscles = muscles.getSize();

for i = 1:numMuscles
   % Downcast base muscle to Millard2012EquilibriumMuscle
   muscle = Millard2012EquilibriumMuscle.safeDownCast(muscles.get(i-1));
   muscleNames{i} = char(muscle);
   muscle.set_ignore_tendon_compliance(true);
   muscle.set_ignore_activation_dynamics(true);
end

if (withviz == true)
    model.setUseVisualizer(true);
end

% Update the system to include any muscle modeling changes
state = model.initSystem();

%% Get coordinate actuators 
allActs = model.getActuators; 
num_acts = getSize(allActs); 

% get all actuators and override actuation for the muscles only
for i = 1:num_acts
    act_Names{i} = char(allActs.get(i-1)) ;
    acts(i) = ScalarActuator.safeDownCast(allActs.get(i-1));
    if i<=numMuscles
        acts(i).overrideActuation(state, true);
    end
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
accelerations = accelerations(1:timeInterval:N, :);
numTimePoints = size(coordinates, 1);
% numTimePoints = 50/timeInterval * 2;      % considering just the first 2 seconds of the experiment

% Create the FMINCON options structure.
options = optimoptions('fmincon','Display','none', ...
     'TolCon',1e-3,'TolFun',1e-3,'TolX',1e-6,'MaxFunEvals',100000, ...
     'MaxIter',10000,'Algorithm','sqp', 'DiffMinChange', 1.0e-2);
 
% Construct initial guess and bounds arrays
numCoords = length(coordNames);
numCoordActs = num_acts-numMuscles;

lb = [zeros(1,numMuscles), -600*ones(1,numCoordActs)];
ub = [ones(1,numMuscles), 600*ones(1,numCoordActs)];
x0 = [0.1* ones(1,numMuscles), zeros(1,numCoordActs)];

% We define the activation squared cost as a MATLAB anonymous function
epsilon = 0;
w = [ones(1,numMuscles), epsilon*ones(1,8), 10*ones(1,9)];     % the cost function is written such that it allows the use of coord acts for the underactuated coordinates
cost =@(x) sum(w.*(x.^2));

% Pre-allocate arrays to hold muscle values for each time point.
fl = zeros(1, numMuscles);
fv = zeros(1, numMuscles);
fp = zeros(1, numMuscles);
cosPenn = zeros(1, numMuscles);
Fmax = zeros(1, numMuscles);

% Pre-allocate an array to hold the solution for all time points.
xsol = zeros(numTimePoints, length(x0));
simulatedAccelerations = zeros(numTimePoints, length(coordNames));

coords = model.getCoordinateSet();

tic
for i = 1:numTimePoints
    fprintf('Optimizing...time step %i/%i \n',  i, numTimePoints);
    
    % Loop through model coordinates to set coordinate values and speeds. We set
    % all coordinates to make sure we have the correct kinematic state when 
    % compute muscle multipliers and moment arms.
    for j = 1:length(coordNamesAll)
        coord = coords.get(coordNamesAll{j});
        coord.setValue(state, coordinates(i,j), false);
        coord.setSpeedValue(state, speeds(i,j));
    end

    % equilibrate the muscles to make them start in the correct state
    model.equilibrateMuscles(state);
    
    % Populate the muscle multiplier and moment arm arrays. We must realize the
    % system to the Velocity stage. Here, we only need to loop through the
    % coordinates associated with the generalized forces 
    model.realizeVelocity(state);

    % populate muscle multipliers
    for k = 1:numMuscles
        muscle = Millard2012EquilibriumMuscle.safeDownCast(muscles.get(k-1));
        fl(k) = muscle.getActiveForceLengthMultiplier(state);            
        fv(k) = muscle.getForceVelocityMultiplier(state);
        fp(k) = muscle.getPassiveForceMultiplier(state);
        cosPenn(k) = muscle.getCosPennationAngle(state);
        Fmax(k) = muscle.getMaxIsometricForce();
    end 
    
    % store the values of active and passive maximum force in the current
    % configuration
    AMuscForce = (fl.*fv.*Fmax.*cosPenn)'; 
    PMuscForce = (Fmax.*fp.*cosPenn)'; 
    % I guess these values should be equivalent to
    % muscle.getActiveFiberForce(state)
    % muscle.getPassiveFiberForce(state) or
    % muscle.getPassiveFiberElasticForce(state) - this one matches!
    % but actually that is not the case -> ask Ajay!

    % create a struct containing relevant information to be passed to the
    % function evaluating the nonlinear joint reaction constraint
    % It is much faster to pass these as parameters than to evaluate inside
    % the constraint expression every time

    params.model = model;
    params.acts = acts;
    params.nActs = num_acts-numMuscles;
    params.state = state;
    params.coords = coords;
    params.coordNames = coordNames;
    params.AccelData = accelerations(i, :);
    params.nMuscles = numMuscles;
    params.muscles = muscles;
    params.AMuscForce = AMuscForce;
    params.PMuscForce = PMuscForce;

    %% Get model controls
    model_controls = model.getControls(state);

    params.model_controls = model_controls;
    
    % Call FMINCON to solve the problem
    [x,fval,exitflag,output] = fmincon(cost, x0, [], [], [], [], lb, ub, @(x)accelCons(x, params, 1), options);
    if exitflag ==0
        [x,fval,exitflag,output] = fmincon(cost, x, [], [], [], [], lb, ub, @(x)accelCons(x, params, 1), options);
        exitflag
    end

    % get best feasible point, if different from what returned by fmincon
    if size(output.bestfeasible,1)>0
        x = output.bestfeasible.x;
    end
         
    % Store solution and set guess for next time point
    xsol(i, :) = x;
    x0 = x;
    x0(end-8:end) = 0; % setting to 0 the initial guess for the redundant actuators. This does not really change the activation results

    % update cost expression to include changes in intial guess
    cost =@(x) sum(w.*(x.^2))+(x-x0)*(x-x0)';

    % Retrieve the optimal accelerations
    % override the muscle forces
    MuscleForces = AMuscForce.*xsol(i, 1:numMuscles)' + PMuscForce;
    for k = 1:numMuscles
        muscle = muscles.get(k-1);
        muscle.setOverrideActuation(state, MuscleForces(k));
    end

    % override the coordinate actuator forces
    for k = numMuscles+1:num_acts
        acts(k).setControls(Vector(1, xsol(i, k)), model_controls);
    end

    model.realizeVelocity(state);
    model.setControls(state, model_controls);
        
    model.realizeAcceleration(state);

    if (withviz == true)
        model.getVisualizer.show(state);
    end
    
    for j = 1:length(coordNames)
        coord = coords.get(coordNames{j});
        simulatedAccelerations(i, j) = coord.getAccelerationValue(state);
    end
end
tOptim = toc;

%% Plot results
% According to the value of the 'print_flag', we print the muscle
% activations that are found

if print_flag
    % plot muscle activations
    figure;
    title("Muscle Activations")
    muscleNames = ArrayStr();
    muscles.getNames(muscleNames);
    pgc = linspace(0, 100, size(xsol,1));
    for i = 1:numMuscles
       subplot(5,8,i)
       hold on
       plot(pgc,xsol(:,i),'b-')
       ylim([0 1])
       muscName = muscleNames.get(i-1).toCharArray';
       title(muscName(1:end), 'interpreter', 'none')
       hold off
    end
    legend("muscle activation")
    
    % Plot reserve actuator excitations.
    figure;
    title("Reserve actuators")
    side = ceil(sqrt(numCoordActs));
    for i = 1:numCoordActs
        subplot(side,side,i)
        plot(xsol(:,numMuscles+i), 'linewidth', 2);
        title(char(acts(numMuscles+i)));
    end
    legend("reserve act value")

    % plot accelerations
    figure;
    title("Accelerations")
    side = ceil(sqrt(length(coordNames)));
    for i = 1:length(coordNames)
        subplot(side,side,i)
        hold on
        plot(accelerations(:, i), 'linewidth', 1.5);
        plot(simulatedAccelerations(:, i), 'linewidth', 1.5);
        xlabel("samples")
        ylabel("[]/s^2")
        grid on
        title(coordNames{i});
        hold off
    end
    legend("measured", "simulated")

    % plot the constraint violation on the accelerations per coordinate
    violation = abs(accelerations-simulatedAccelerations);

    figure
    for i = 1:length(coordNames)
        subplot(side,side,i)
        hold on
        plot(violation(:,i), 'linewidth', 1.5);
        xlabel("samples")
        ylabel("[]/s^2")
        grid on
        title(coordNames{i});
        hold off
    end
    legend("acc violation")

    % plot the constraint violation per timestep
    violation_t = sum(violation, 2);
    figure
    hold on
    scatter(1:numTimePoints ,violation_t, 'filled')
    plot(1:numTimePoints, violation_t, 'blue')
    xlabel("samples")
    ylabel("const violation")
    grid on
    title("Cumulative constraint violation per time-step")
    hold off

end

%% inspect the importance of each term in x in the cost function
cost = xsol*w';
cost_elements = w.* xsol;


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

% rescale the frequancy of the solution knowing that Cark's one is 50 Hz
frequency_solution = 50/timeInterval;

save(name_file, 'xsol', 'muscle_order', 'frequency_solution');
