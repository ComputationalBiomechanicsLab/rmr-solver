function [sumSqrError, model, RMSEdistance] = computeScaleFitToData_ThoraxScapula(x, expData, fixedModel, ikTool, scaleSet, print)
% Compute the error of a model's fit to experimental data according to the
% input parameters, x. Parameters are the scale factors of the bodies in the model.
% USAGE error = computeModelFitToData.m(parameters, expData, fixedModel, ikTool, scaleSet)
% Ajay Seth  2015
% modified: Italo Belli (i.belli@tudelft.nl) 2022, adding print option:
% * print == True -> print information while optimizing
% * print == False -> do not print information while optimizing

import org.opensim.modeling.*

model = fixedModel.clone();

%% 1. Edit scale factors

% thorax
factors = Vec3(x(1), x(2), x(3));
scaleSet.get(0).setScaleFactors(factors);
segName = scaleSet.get(0).getSegmentName();
segFactors = scaleSet.get(0).getScaleFactors().toString();
if (print)
    fprintf('\n body %s factor: %s\n', segName.toChar(), segFactors.toChar());
end

% scapula
factors = Vec3(x(4), x(5), x(6));
scaleSet.get(2).setScaleFactors(factors);
segName = scaleSet.get(2).getSegmentName();
segFactors = scaleSet.get(2).getScaleFactors().toString();
if (print)
    fprintf('body %s factor: %s\n', segName.toChar(), segFactors.toChar());
end

% Update model
state = model.initSystem();
model.scale(state, scaleSet, true);


%% 2. Update ScapuloThoracicJoint and AC constraint

% Apply scale factors to ScapuloThoracicJoint properties
% EllipsoidSurface = model.getBodySet().get('thorax').getWrapObjectSet().get('EllipsoidSurface');
% EllipsDim = EllipsoidSurface.getPropertyByName('dimensions');
% scapThor = model.getJointSet().get(2);
% pf = PhysicalOffsetFrame.safeDownCast(scapThor.getParentFrame());
% pf.set_translation(EllipsoidSurface.get_translation());
% 
% radiiProp = scapThor.getPropertyByName('thoracic_ellipsoid_radii_x_y_z');
% PropertyHelper.setValueVec3( PropertyHelper.getValueVec3(EllipsDim,0) , radiiProp, 0);
% PropertyHelper.setValueVec3( PropertyHelper.getValueVec3(EllipsDim,1) , radiiProp, 1);
% PropertyHelper.setValueVec3( PropertyHelper.getValueVec3(EllipsDim,2) , radiiProp, 2);

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
