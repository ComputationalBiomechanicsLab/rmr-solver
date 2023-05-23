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
subject_considered = 'Ajay2019';

apply_external_force = false;
force_value = 0;
force_direction = 'up';         % direction of the force that the person will need to generate to compensate the external one
body = 'hand';

motion_file_name = 'abd21.mot';
print_flag = true;
withviz = false;

weight_coord = 0.0001; 

dynamic_activation_bounds = true;

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
modelFile = 'C:\Users\italobelli\Desktop\GitHub\PTbot_officialCBL\Personal_Results\TSM_Ajay2019_2kgWeight.osim';
model = Model(modelFile);

%% Select the trc file to be considered from Clark's dataset
% [trc_file_used, trcPath] = uigetfile('*.trc', 'Select processed TRC file in PTbot\DataShare\dataset_Clark\MoCap\TRC files');
trcPath = 'C:\Users\italobelli\Desktop\Ajay old studies\FrontiersPubMaterials-latest\ThoracoscapularShoulderPaperMaterials\ThoracoscapularShoulderPaperMaterials\ExperimentalData\Markers\';
trc_file_used = 'abd21.trc';

[markersExp, timesExp, labelsExp, unitsExp] = readTRC([trcPath, trc_file_used]);
start_time = timesExp(1);
end_time =  timesExp(end);

if strcmp(unitsExp, 'mm')
    markersExp = markersExp/1000;
    unitsExp = 'm';
end

frequency_trc_data = 1/(timesExp(2)-timesExp(1));

%% Getting quantities about GelnoHumeral joint
model_temp = model.clone();    % create a temporary copy of the model, to be used in the tool
state = model_temp.initSystem();
[maxAngle, ~] = get_glenoid_status(model_temp, state);

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
[coordinates, coordNames, ~] = loadFilterCropArray(motion_file_name, lowpassFreq, timeRange);
coordinates(:, 1:3) = deg2rad(coordinates(:, 1:3));
coordinates(:, 7:end) = deg2rad(coordinates(:, 7:end));

% get the velocities for each joint in rad/s
time_step_trc = timesExp(2)-timesExp(1);
speeds = zeros(size(coordinates));
for i=1:size(coordNames,1)
    speeds(:,i) = gradient(coordinates(:,i), time_step_trc);
end
speedNames = coordNames;

% get the accelerations for each coordinate in rad/s^2
accelerations = zeros(size(speeds));
for i=1:size(coordNames,1)
    accelerations(:,i) = gradient(speeds(:,i), time_step_trc);
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
        force_direction_OS = '-Z';
    elseif strcmp(force_direction, 'left')
        force_direction_OS = '+Z';
    elseif strcmp(force_direction, 'up')
        force_direction_OS = '-Y';
    elseif strcmp(force_direction, 'down')
        force_direction_OS = '+Y';
    elseif strcmp(force_direction, 'push')
        force_direction_OS = '-X';
    elseif strcmp(force_direction, 'pull')
        force_direction_OS = '+X';
    end
    
    % getting markers for application point (change it to the marker you
    % want!!)
    application_point = markersExp(:,4:6);
 
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
muscleNames = cell(numMuscles,1);

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
acts = cell(num_acts,1);

% get all actuators and override actuation for the muscles only
for i = 1:num_acts
    act_Names{i} = char(allActs.get(i-1)) ;
    acts(i) = ScalarActuator.safeDownCast(allActs.get(i-1));
    if i<=numMuscles
        acts{i}.overrideActuation(state, true);
    end
end

%% Perform static optimization.
% We use FMINCON to solve the static optimization problem at selected time points. 
% We set the 'timeInterval' variable to select the time points to be included in the 
% optimization. For example, if set to 10, every 10th time point is selected. A 
% time interval of 1 will select all available time points. 
timeInterval = 10;
time_step_SO = time_step_trc * timeInterval;

% Update data arrays based on the time interval.
N = size(coordinates, 1);
coordinates = coordinates(1:timeInterval:N, :);
speeds = speeds(1:timeInterval:N, :);
accelerations = accelerations(1:timeInterval:N, :);
numTimePoints = size(coordinates, 1);
% numTimePoints = 50/timeInterval * 2;      % considering just the first 2 seconds of the experiment

% Create the FMINCON options structure.
options = optimoptions('fmincon','Display','iter', ...
     'TolCon',1e-2,'TolFun',1e-3,'TolX',1e-3,'MaxFunEvals',100000, ...
     'MaxIter',10000,'Algorithm','sqp', 'StepTolerance', 1e-8); %, 'DiffMinChange', 1.0e-3);
 
% Construct initial guess and bounds arrays
numCoords = length(coordNames);
numCoordActs = num_acts-numMuscles;
k = 6;
t_act = 0.01;           % activation time constant for muscles
t_deact = 0.04;         % deactivation time constant

lb = [zeros(1,numMuscles), -k*ones(1,numCoordActs)];
ub = [ones(1,numMuscles), k*ones(1,numCoordActs)];
x0 = [0.1* ones(1,numMuscles), zeros(1,numCoordActs)];

% We define the activation squared cost as a MATLAB anonymous function
epsilon = 0;
w = [ones(1,numMuscles), epsilon*ones(1,8), 10*ones(1,9)];     % the cost function is written such that it allows the use of coord acts for the underactuated coordinates
cost =@(x) sum(w.*(x.^2));

% Pre-allocate arrays to be filled in the optimization loop
fl = zeros(1, numMuscles);
fv = zeros(1, numMuscles);
fp = zeros(1, numMuscles);
cosPenn = zeros(1, numMuscles);
Fmax = zeros(1, numMuscles);
A_eq_acc = zeros(numCoords,num_acts);
A_eq_force = zeros(3, num_acts);
xsol = zeros(numTimePoints, length(x0));
simulatedAccelerations = zeros(numTimePoints, length(coordNames));
optimizationStatus = cell(numTimePoints,1);
norm_fv_in_ground = zeros(numTimePoints, 3);
norm_fv_rotated = zeros(numTimePoints, 3);
rel_angle = zeros(numTimePoints);

% get model quantities we still need
coords = model.getCoordinateSet();

tic
% enter in the optimization loop
for i = 1:numTimePoints
    fprintf('Optimizing...time step %i/%i \n',  i, numTimePoints);
    
    % Loop through model coordinates to set coordinate values and speeds. We set
    % all coordinates to make sure we have the correct kinematic state when 
    % compute muscle multipliers and moment arms.
    for j = 1:length(coordNames)
        coord = coords.get(coordNames{j});
        coord.setValue(state, coordinates(i,j), false);
        coord.setSpeedValue(state, speeds(i,j));
    end

    % equilibrate the muscles to make them start in the correct state
    model.realizeVelocity(state);
    model.equilibrateMuscles(state);           % assuming fiber velocity is the same as the overall muscle lengthening velocity
   
    modelControls = model.getControls(state);

    % Populate the muscle multiplier arrays. To do this, we must realize the
    % system to the Velocity stage
    for k = 1:numMuscles
        muscle = Millard2012EquilibriumMuscle.safeDownCast(muscles.get(k-1));
        fl(k) = muscle.getActiveForceLengthMultiplier(state);            
        fv(k) = muscle.getForceVelocityMultiplier(state);
        fp(k) = muscle.getPassiveForceMultiplier(state);
        cosPenn(k) = muscle.getCosPennationAngle(state);
        Fmax(k) = muscle.getMaxIsometricForce();
    end 

    % get the vector Vec_H2GC between humeral head and the glenoid center
    % it is expressed in the ground frame
    [~, Vec_H2GC] = get_glenoid_status(model, state);
    
    % store the values of active and passive maximum force in the current
    % configuration
    AMuscForce = (fl.*fv.*Fmax.*cosPenn)'; 
    PMuscForce = (Fmax.*fp.*cosPenn)'; 

    % create a struct containing relevant information to be passed to the
    % function simulating the accelerations induced in the model
    params_indAcc.model = model;
    params_indAcc.state = state;    
    params_indAcc.AMuscForce = AMuscForce;
    params_indAcc.PMuscForce = PMuscForce;
    params_indAcc.coords = coords;
    params_indAcc.coordNames = coordNames;
    params_indAcc.acts = acts;
    params_indAcc.muscles = muscles;
    params_indAcc.useMuscles = 1;
    params_indAcc.useControls = 1;
    params_indAcc.modelControls = modelControls;

    % evaluate the matrices that allow the acceleration constraint to be
    % written as a linear one
    q_ddot_0 = findInducedAccelerations(zeros(1,num_acts),params_indAcc);
    delQ_delX = eye(num_acts);
    
    for k = 1:num_acts
        incrementalForceAccel = findInducedAccelerations(delQ_delX(k,:),params_indAcc);
        kthColumn =  incrementalForceAccel - q_ddot_0;
        A_eq_acc(:,k) = kthColumn;
    end

    Beq = accelerations(i,:)' - q_ddot_0;

    % evaluate the matrices that allow the reaction force to be expressed
    % as linear
    params_force.model = model;
    params_force.state = state;
    params_force.acts = acts;
    params_force.muscles = muscles;
    params_force.AMuscForce = AMuscForce;
    params_force.PMuscForce = PMuscForce;
    params_force.glen = glen;
    params_force.useControls = 1;
    params_force.modelControls = modelControls;

    % COmpute linearized expression for the joint reaction force, and
    % equivalent force matrix to be used in the optimization
    [F_r0, ~] = findReactionForceAndMomentGH(zeros(1,num_acts), params_force);

    for k = 1:num_acts
        [F_rk, ~] = findReactionForceAndMomentGH(delQ_delX(k,:),params_force);
        kthColumn =  F_rk - F_r0;
        A_eq_force(:,k) = kthColumn;
    end

    % Call FMINCON to solve the problem
    [x,fval,exitflag,output] = fmincon(cost, x0, [], [], A_eq_acc, Beq, lb, ub, @(x)jntrxncon_linForce(x, Vec_H2GC, maxAngle, A_eq_force, F_r0), options);
    if exitflag ==0
        [x,fval,exitflag,output] = fmincon(cost, x, [], [], A_eq_acc, Beq, lb, ub, @(x)jntrxncon_linForce(x, Vec_H2GC, maxAngle, A_eq_force, F_r0), options);    % call the solver again, starting from current x, in case the maximum iterations are exceeded
        exitflag
    end

    optimizationStatus{i} = output;

    % get best feasible point, if different from what returned by fmincon
    if size(output.bestfeasible,1)>0
        x = output.bestfeasible.x;
    end
         
    % Store solution and set guess for next time point
    xsol(i, :) = x;
%     x0 = x;
%     x0(end-8:end) = 0; % setting to 0 the initial guess for the redundant actuators. This does not really change the activation results, but looks reasonable to avoid drifting of the solution

    % update cost expression to include changes in initial guess
    % NOTE: this makes sense if the time grid at which SO is performed is
    % fine enough, otherwise it might be not so correct
%     cost =@(x) sum(w.*(x.^2)) +(x(1:numMuscles)-x0(1:numMuscles))*(x(1:numMuscles)-x0(1:numMuscles))';

    % dynamically update the upper and lower bounds for the activations
    if dynamic_activation_bounds
        for k = 1:numMuscles
            lb(k) = max(x(k) - x(k) * (0.5 + 1.5 * x(k)) * time_step_SO /t_deact, 0);
            ub(k) = min (x(k) + (1-x(k)) * time_step_SO / (t_act * (0.5 + 1.5*x(k))), 1);
        end
    end

    % Retrieve the optimal accelerations
    simulatedAccelerations(i,:) = findInducedAccelerations(x,params_indAcc);

    % retrieve the position of the joint reaction force on the approximated
    % glenoid
    % computing the reaction force vector at the given joint
    force_vec = A_eq_force * x' + F_r0;
    
    % evaluate the relative angle between the reaction force and Vec_H2GC
    cosTheta = max(min(dot(Vec_H2GC,force_vec)/(norm(Vec_H2GC)*norm(force_vec)),1),-1);
    rel_angle(i) = real(acosd(cosTheta));

    % evaluate the position on the glenoid where reaction force is exerted
    norm_Vec_H2GC = Vec_H2GC/norm(Vec_H2GC);
    norm_fv_in_ground(i,:) = force_vec/norm(force_vec);

    beta_angle = atan(norm_Vec_H2GC(3)/norm_Vec_H2GC(1));
    alpha_angle = atan(norm_Vec_H2GC(3)/(sin(beta_angle)*norm_Vec_H2GC(2)));

    Ry = [cos(beta_angle) 0 sin(beta_angle); 0 1 0; -sin(beta_angle) 0 cos(beta_angle)];
    Rz = [cos(alpha_angle) -sin(alpha_angle) 0; sin(alpha_angle) cos(alpha_angle) 0; 0 0 1];

    norm_fv_rotated(i,:) = Rz*Ry*norm_fv_in_ground(i,:)';

    if (withviz == true)
        model.getVisualizer.show(state);
    end
end

tOptim = toc;

%% Plot results
% According to the value of the 'print_flag'
if print_flag
    % plot muscle activations
    f1 = figure;
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
    f1.WindowState = 'maximized';
    saveas(f1, 'MuscleActivations.png')
    
    % Plot reserve actuator excitations.
    f2 = figure;
    title("Reserve actuators")
    side = ceil(sqrt(numCoordActs));
    for i = 1:numCoordActs
        subplot(side,side,i)
        hold on
        plot(pgc, xsol(:,numMuscles+i), 'linewidth', 2);
        title(char(acts{numMuscles+i}));
        hold off
    end
    legend("reserve act value")
    f2.WindowState = 'maximized';
    saveas(f2, 'ReserveActuators.png')

    % plot accelerations
    f3 = figure;
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
    f3.WindowState = 'maximized';
    saveas(f3, 'AccelerationsMatching.png')

    % plot the constraint violation on the accelerations per coordinate
    violation = abs(accelerations-simulatedAccelerations);

    f4 = figure;
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
    f4.WindowState = 'maximized';
    saveas(f4, 'AccViolation.png')

    % plot the constraint violation per timestep
    violation_t = sum(violation, 2);
    f5 = figure;
    hold on
    scatter(1:numTimePoints ,violation_t, 'filled')
    plot(1:numTimePoints, violation_t, 'blue')
    xlabel("samples")
    ylabel("const violation")
    grid on
    title("Cumulative constraint violation per time-step")
    hold off
    f5.WindowState = 'maximized';
    saveas(f5, 'CumulativeAccViolation.png')

    % plot the position of the GH force on the glenoid
    radius = sind(maxAngle);
    p=nsidedpoly(1000, 'Center', [0,0], 'Radius', radius);
    c = linspace(0,timesExp(end),length(norm_fv_rotated));
    f6 = figure;
    hold on
    plot(p, 'FaceColor', 'r')
    for i=1:numTimePoints
        scatter(-norm_fv_rotated(i,3), -norm_fv_rotated(i,1), [], c(i), 'filled')
    end
    hcb = colorbar;
    h = gca;
    set(h, "XTickLabel", [])
    set(h, "YTickLabel", [])
    xlabel("back                                                                       front")   % corresponding roughly to OpenSim X axis (horizontal pointing forward)
    ylabel("down                                                                       up")      % corresponding to OpenSim Y axis (vertical pointing upwards)
    colorTitleHandle = get(hcb,'Title');
    titleString = 'time [s]';
    set(colorTitleHandle ,'String',titleString);
    hold off
    saveas(f6, 'CoPGH.png')

end

%% SAVING THE RESULTS TO FILE
name_trc = extractBefore(trc_file_used, '.');
name_file = append('muscle_activations_', subject_considered, '_experiment_', name_trc);
muscle_order = "";
for i = 1:numMuscles
    muscle_order= [muscle_order, string(muscleNames.get(i-1).toCharArray')];
end

for i=1:length(coordNames)
    muscle_order= [muscle_order, string(coordNames{i})];
end

muscle_order = muscle_order(2:end);

% rescale the frequancy of the solution knowing the freq of the data
frequency_solution = frequency_trc_data/timeInterval;

save(name_file, 'xsol', 'muscle_order', 'frequency_solution');
