% This script allows the user to rescale a (ThoracoScapular) OpenSim
% model to match the biological metrics of a subject, starting from marker
% data (in .trc format) provided by the user.
% Depending on the model considered, the subject and the specific
% experiments used, it is necessary to change (or create) the setup files
% for the tools that are used in it (ScaleTool and IKTool).
%
% authors: Italo Belli (i.belli@tudelft.nl) and Irene Beck
% (I.L.Y.Beck@student.tudelft.nl) 2022
%
% written on the basis of code by Joran de Vet, 2020

clc; clear all; close all;

% set the path current folder to be the one where this script is contained
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

root_dir = fullfile('');

% getting path to other folders in this repo
addpath('..\..\Data Processing\Matlab\')
cd ..\..\..\
path_to_repo = pwd;
cd(pathstr);

%% Import OpenSim libraries
import org.opensim.modeling.*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print_flag = true;                               % determine whether to print info about optimization at run-time

% Generic model that needs to be scaled
% Select the model you want to use, among the ones available in the
% OpenSim Models folder of this repository
[modelFile, modelPath] = uigetfile('*.osim', 'Select Generic model that needs to be scaled');

model = Model([modelPath,modelFile]);

% Directory containing experimental data of the subject of interest
% Select: U:\PTbot\Analysis\MATLAB\
dirName = uigetdir(modelPath,'Select base MATLAB folder accessing PTbot\Analysis\MATLAB on shared drive');
PTbot_path = [dirName, '\..\..\'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preliminary Stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% task-name and combined .trc file-name for InverseKinematics Tool
task_name = 'all_trials';
trcSampledFile = [modelFile(1:end-5),'_all_sampled.trc'];

% .trc file sub-directory
trcDir = [dirName,'.\Code Ajay\trc_for_scaling\'];
trcFiles = dir([trcDir '*.trc']);

if size(trcFiles,1)==0
    error_msg = ['Unable to load .trc files from ', trcDir];
    error(error_msg)
end

% Set output directory to 'Personal_Results', to save results in local
current_dir = pwd;
cd([path_to_repo, '\Personal_Results']);
addpath(current_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 1: Combine and down-sample experimental data in trcFiles into single .trc file used in optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize structures
lastTime = 0;
deltaT = 0;
markersSamp = [];
timesSamp = [];

% Down-sample setting
DownSampFactor = 5;                                                                        % downsampling factor used

for I=1:length(trcFiles)
    % Read in the experimental marker data
    expFileName = [trcDir, trcFiles(I).name];
    [markersExp, timesExp, labelsExp, unitsExp] = readTRC(expFileName);
    markersExp = markersExp(:,1:end-6);                                                     % NOTE THAT WE EXCLUDE THE LAST 2 MARKERS
    labelsExp = labelsExp(1:end-2);
    framesSamp = [1:DownSampFactor:length(timesExp)]';                                      % select the frames to keep after downsampling
    if unitsExp == "mm"                                                                     % Convert experimental marker coordinates to meters if in mm
        markersExp = markersExp./1000;
        unitsExp = "m";
    end
    markersSamp = [markersSamp; markersExp(framesSamp,:)];                                  % Down-sample experimental data to speedup fitting 
    timesSamp = [timesSamp; timesExp(framesSamp,:)+lastTime+deltaT];                        % create a single time vector
    deltaT = timesSamp(2)-timesSamp(1);
    lastTime = timesSamp(end);
end

Rate = 1/(deltaT);                                                                          % get the rate of the new data
 
writeMarkersToTRC(trcSampledFile, markersSamp, labelsExp, Rate, [], timesSamp, unitsExp);   % save a new .trc file (the one that will be used in the IK)

%% Pre-scale the model
% Classic OpenSim scaling that can be done in GUI
% This is not an optimization technique, only average marker distances in a
% static pose are used to provide a better initial guess for the later
% stages of the scaling

% This part has to be checked carefully when a new model is selected
disp('')
disp('-------------------------------------------')
disp('Have you checked the setup files carefully?')
disp('-------------------------------------------')
disp('')

scaleSetupFile= [path_to_repo,'' ...
        '\Code\Scale Model\Utilities\setup_files\Clark\PreScale_setup_Clark.xml']; 

scaleTool= ScaleTool(scaleSetupFile);

model_name = model.getName().toCharArray()';
scaleTool.getGenericModelMaker().setModelFileName([modelPath, model_name '.osim']);

scaleTool.getGenericModelMaker().setMarkerSetFileName([path_to_repo, '' ...
    '\Code\Scale Model\Utilities\marker_set\TSM_marker_set_Clark.xml']);

scaleTool.getModelScaler().setMarkerFileName([PTbot_path, '' ...
    '\DataShare\dataset_Clark\MoCap\TRC files\ACC\5.trc']);                                 % remember to check in the .xml the time range considered (default: 5-6.5) 

scaleTool.getMarkerPlacer().setMarkerFileName([PTbot_path, '' ...
    '\DataShare\dataset_Clark\MoCap\TRC files\ACC\5.trc'])

scaleTool.getModelScaler().setOutputModelFileName([path_to_repo, '\Personal_Results\', model_name '_preScaled.osim']);
scaleTool.getMarkerPlacer().setOutputModelFileName([path_to_repo, '\Personal_Results\', model_name '_preScaled.osim']);
scaleTool.setPathToSubject(root_dir);
scaleTool.print('autogenerated_scale_setup.xml');
scaleTool.run();

% Get the model name and initialize it
model_preScaled = Model(char(scaleTool.getMarkerPlacer().getOutputModelFileName));
model_preScaled.initSystem();

%% Set up InverseKinematics- and Scale Tools for optimization
% These might be generic-model-dependent (because of marker differences)
ikSetupFile = [path_to_repo,'' ...
        '\Code\Scale Model\Utilities\setup_files\Clark\IKSetup_Clark.xml'];

ikTool = InverseKinematicsTool(ikSetupFile);
ikTool.setMarkerDataFileName([path_to_repo, '\Personal_Results\', trcSampledFile]);
ikTool.setOutputMotionFileName([path_to_repo, '\Personal_Results\', task_name '_IK_fitted.mot']);
ikTool.set_report_marker_locations(1);
ikTool.setEndTime(lastTime+deltaT);                                                         % deltaT: necessary as sampled .trc file starts from deltaT and not from 0
ikTool.print('autogenerated_IK_setup.xml');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 2: Optimize the Scapulothoracic Joint to best fit Marker Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We optimize the following thoracoscapular joint parameters: 
% position of the center in parent frame 
% radius_x, radius_y, radius_z(=radius_x)
% rotation_z
%
% Eventually, we are left with 6 optimization variables to find.

% retrieve the position of the ellipsoid center in parent frame
% (thorax_offset)
scapThor = model_preScaled.getJointSet().get(2);
position_in_parent = scapThor.get_frames(0).get_translation().getAsMat()';

% retrieve the current parameters for the thoracoscapular joint
radiiProp = scapThor.getPropertyByName('thoracic_ellipsoid_radii_x_y_z');
r_x = PropertyHelper.getValueVec3(radiiProp, 0);
r_y = PropertyHelper.getValueVec3(radiiProp, 1);

oInPvec = PhysicalOffsetFrame.safeDownCast(scapThor.getParentFrame()).get_orientation;
rot_z = oInPvec.get(2);

% build the initial guess for the overall optimization vector
s0 = [position_in_parent, r_x, r_y, rot_z];

%% Optimizer options 
%options = optimset;
if(print_flag)
    options = optimset('PlotFcns',@optimplotfval);
end
options.TolX = 5e-3;
options.TolFun = 1e-3;
options.MaxIter=300;

tic

%% Optimizer 
f = @(p) computeParamFitToData_Ellipsoid(p, markersSamp, model_preScaled, ikTool, print_flag);
[sFinal, fValFinal_Scaling, ScaleOptExitFlag, outputScaleOpt] = fminsearch(f, s0, options);


%% Apply final scale factors and save model if optimized
if(exist('sFinal', 'var')) 
    x = sFinal;

    tOptScale = toc;
    disp(['Found optimal set of scale factors in ',num2str(tOptScale/60),' minutes.'])
    
    model_preScaled.initSystem();

    % update ellipsoid center position
    scapThor = model_preScaled.getJointSet().get(2);
    scapThor.get_frames(0).set_translation(Vec3(x(1), x(2), x(3)));

    % update ellipsoid radii with optimal values
    radiiProp = scapThor.getPropertyByName('thoracic_ellipsoid_radii_x_y_z');
    PropertyHelper.setValueVec3(x(4), radiiProp, 0);
    PropertyHelper.setValueVec3(x(5), radiiProp, 1);
    PropertyHelper.setValueVec3(x(4), radiiProp, 2);
    
    % update ellipsoid orientation with optimal values
    ellipsoidFrame = PhysicalOffsetFrame.safeDownCast(scapThor.getParentFrame());
    ellipsoidOrientation = ellipsoidFrame.get_orientation;
    orientationXAxis = ellipsoidOrientation.get(0);
    orientationYAxis = ellipsoidOrientation.get(1);

    ellipsoidFrame.set_orientation(Vec3(orientationXAxis, orientationYAxis, x(6)));

    % Update parameters of the wrapping object representing the ellipsoid
    ellipsoid = model_preScaled.get_BodySet().get('thorax').get_WrapObjectSet().get('EllipsoidSurface');

    ellipsoid.set_translation(Vec3(x(1), x(2), x(3)));

    ellipsoidDim = ellipsoid.getPropertyByName('dimensions');
    PropertyHelper.setValueVec3(x(4), ellipsoidDim, 0);
    PropertyHelper.setValueVec3(x(5), ellipsoidDim, 1);
    PropertyHelper.setValueVec3(x(4), ellipsoidDim, 2);
    
    ellipsoid.set_xyz_body_rotation(Vec3(orientationXAxis, orientationYAxis, x(6)));

    % update the AC constraint

   model_preScaled = UpdateAC(model_preScaled);

    %% Save scaled model
    modelName = [model_preScaled.getName().toCharArray()', '_scaled'];
    model_preScaled.setName(modelName);
    model_preScaled.print([modelName, '.osim']);
    
    % Note: inspecting the model (.osim file) it is possible to retrieve which
    % are the values of the overall scaling applied
end


%% Evaluate the RMSE of the fitted model
% retrieve ellipsoid radii values
radius_x = PropertyHelper.getValueVec3(radiiProp, 0);
radius_y = PropertyHelper.getValueVec3(radiiProp, 1);

ikMarkersFileName = '_ik_model_marker_locations.sto';
markersIK_data = dlmread(ikMarkersFileName, '', 7, 0);
%timeIK = markersIK_data(:,1);
modelData = markersIK_data(:,2:end);
clear markersIK_data;
error = markersSamp - modelData;

distances = computeDistancesFromDisplacementsXYZ(error);
RMSE_final_model_markers = sqrt(mean((distances.^2),'omitnan'));
RMSE_final_avg_mm = 1000 * mean(RMSE_final_model_markers);

error_per_marker = mean(distances, 1)
labelsExp
beep
