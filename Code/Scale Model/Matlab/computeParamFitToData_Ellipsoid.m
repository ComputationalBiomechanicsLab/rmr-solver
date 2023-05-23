function [sumSqrError, model, RMSEdistance] = computeParamFitToData_Ellipsoid(x, expData, fixedModel, ikTool, print)
% Compute the error of a model's fit to experimental data according to the
% input parameters, x. Parameters are the scale factors of the bodies in the model.
% USAGE error = computeParamFitToData_Ellipsoid.m(parameters, expData, fixedModel, ikTool, print_flag)
% Ajay Seth  2015
%
% modified: Italo Belli (i.belli@tudelft.nl) 2022, adding print option:
% * print == True -> print information while optimizing
% * print == False -> do not print information while optimizing
%
% x is composed of 7 elements: 
% 3D position of the center in parent frame 
% radius_x, radius_y, radius_z(=radius_x)
% rotation_y, rotation_z
% All the properties describing the joint are optimized in a single loop.

import org.opensim.modeling.*

model = fixedModel.clone();

% Initialize the model
model.initSystem();

%% Edit thoracoscapular parameters

% edit the ellipsoid center position
scapThor = model.getJointSet().get(2);
scapThor.get_frames(0).set_translation(Vec3(x(1), x(2), x(3)));

% edit the ellipsoid radii
radiiProp = scapThor.getPropertyByName('thoracic_ellipsoid_radii_x_y_z');
PropertyHelper.setValueVec3(x(4), radiiProp, 0);
PropertyHelper.setValueVec3(x(5), radiiProp, 1);
PropertyHelper.setValueVec3(x(4), radiiProp, 2);

% edit the ellipsoid orientation (along Y and Z axis of the thorax_frame)
ellipsoidFrame = PhysicalOffsetFrame.safeDownCast(scapThor.getParentFrame());
ellipsoidOrientation = ellipsoidFrame.get_orientation;
orientationXAxis = ellipsoidOrientation.get(0);
orientationYAxis = ellipsoidOrientation.get(1);
ellipsoidFrame.set_orientation(Vec3(orientationXAxis, orientationYAxis, x(6)));

if(print)
    disp('ThoracoScapular Joint parameters:')
    disp(['center position: [',num2str(x(1)),' , ',num2str(x(2)),' , ',num2str(x(3)),']'])
    disp(['radii: [',num2str(x(4)),' , ',num2str(x(5)),' , ',num2str(x(4)),']'])
    disp(['rotation: [Z= ', num2str(x(6)), ' ]'])
    disp('')
end

% Update points in AC constraint
model = UpdateAC(model);

%% 3. Run Inverse Kinematics
ikTool.setModel(model);
ikTool.run();

%% 4. Determine SSE and RMSE of IK solution
% Load the IK marker results
ikMarkersFileName = '_ik_model_marker_locations.sto';
markersIK_data = dlmread(ikMarkersFileName, '', 7, 0);
%timeIK = markersIK_data(:,1);
modelData = markersIK_data(:,2:end);
clear markersIK_data;
error = expData - modelData;

sumSqrError = sum(sum(error.^2,'omitnan'),'omitnan');
if (print)
    fprintf('sumSqrError: %g m\n', sumSqrError);
end

distances = computeDistancesFromDisplacementsXYZ(error);
RMSEdistance = sqrt(mean((distances.^2),'omitnan'));
if (print)
    fprintf('Mean RMSE: %g mm\n\n\n', 1000*mean(RMSEdistance));
end
