function [sumSqrError, model, RMSEdistance] = computeScaleFitToData(x, expData, fixedModel, ikTool, scaleSet, print)
% Compute the error of a model's fit to experimental data according to the
% input parameters, x. Parameters are the scale factors of the bodies in the model.
% USAGE error = computeModelFitToData.m(parameters, expData, fixedModel,
% ikTool, scaleSet, print)
% Ajay Seth  2015
%
% modified Italo Belli (i.belli@tudelft.nl) 2022
% added print option:
% * print == true -> print information while optimizing
% * print == false -> do not print information while optimizing

import org.opensim.modeling.*

model = fixedModel.clone();

%% 1. Edit scale factors

% thorax
factors = Vec3(x(1), x(2), x(3));
scaleSet.get(0).setScaleFactors(factors);
segName = scaleSet.get(0).getSegmentName();
segFactors = scaleSet.get(0).getScaleFactors().toString();
if print
    fprintf('\n body %s factor: %s\n', segName.toChar(), segFactors.toChar());
end

% clavicle
factors = Vec3(x(4), x(4), x(4));
scaleSet.get(1).setScaleFactors(factors);
segName = scaleSet.get(1).getSegmentName();
segFactors = scaleSet.get(1).getScaleFactors().toString();
if print
    fprintf('body %s factor: %s\n', segName.toChar(), segFactors.toChar());
end

% scapula
factors = Vec3(x(5), x(6), x(7));
scaleSet.get(2).setScaleFactors(factors);
segName = scaleSet.get(2).getSegmentName();
segFactors = scaleSet.get(2).getScaleFactors().toString();
if print
    fprintf('body %s factor: %s\n', segName.toChar(), segFactors.toChar());
end

% humerus
factors = Vec3(x(8), x(8), x(8));
scaleSet.get(3).setScaleFactors(factors);
segName = scaleSet.get(3).getSegmentName();
segFactors = scaleSet.get(3).getScaleFactors().toString();
if print
    fprintf('body %s factor: %s\n', segName.toChar(), segFactors.toChar());
end

% radius
factors = Vec3(x(9), x(9), x(9));
scaleSet.get(4).setScaleFactors(factors);
segName = scaleSet.get(4).getSegmentName();
segFactors = scaleSet.get(4).getScaleFactors().toString();
if print
    fprintf('body %s factor: %s\n', segName.toChar(), segFactors.toChar());
end

% ulna
factors = Vec3(x(9), x(9), x(9));
scaleSet.get(5).setScaleFactors(factors);
segName = scaleSet.get(5).getSegmentName();
segFactors = scaleSet.get(5).getScaleFactors().toString();
if print
    fprintf('body %s factor: %s\n', segName.toChar(), segFactors.toChar());
end

% Update model
state = model.initSystem();
model.scale(state, scaleSet, true);


%% 2. Update ScapuloThoracicJoint and AC constraint

% Apply scale factors to ScapuloThoracicJoint properties
EllipsoidSurface = model.getBodySet().get('thorax').getWrapObjectSet().get('EllipsoidSurface');
EllipsDim = EllipsoidSurface.getPropertyByName('dimensions');
scapThor = model.getJointSet().get(2);
pf = PhysicalOffsetFrame.safeDownCast(scapThor.getParentFrame());
% cf = PhysicalOffsetFrame.safeDownCast(scapThor.getChildFrame());
pf.set_translation(EllipsoidSurface.get_translation());
% pf.set_orientation();
% cf.set_translation();
% cf.set_orientation();
% 
radiiProp = scapThor.getPropertyByName('thoracic_ellipsoid_radii_x_y_z');
PropertyHelper.setValueVec3( PropertyHelper.getValueVec3(EllipsDim,0) , radiiProp, 0);
PropertyHelper.setValueVec3( PropertyHelper.getValueVec3(EllipsDim,1) , radiiProp, 1);
PropertyHelper.setValueVec3( PropertyHelper.getValueVec3(EllipsDim,2) , radiiProp, 2);

% Update points in AC constraint
model = UpdateAC(model);
% % point constraint representing acromioclavicular joint
% ac = PointConstraint.safeDownCast(model.getComponent('/constraintset/AC'));
% 
% % constraint point in the clavicle and scapula
% cpInClav = model.getComponent('/bodyset/clavicle/cp_in_clavicle');
% cpInScap = model.getComponent('/bodyset/scapula/cp_in_scapula');
% cpInClav = PhysicalOffsetFrame.safeDownCast(cpInClav);
% cpInScap = PhysicalOffsetFrame.safeDownCast(cpInScap);
% 
% % update constraint locations based on the scaled offset frame locations
% ac.set_location_body_1(cpInClav.get_translation());
% ac.set_location_body_2(cpInScap.get_translation());


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
if print
    fprintf('sumSqrError: %g m\n', sumSqrError);
end

distances = computeDistancesFromDisplacementsXYZ(error);
RMSEdistance = sqrt(mean((distances.^2),'omitnan'));

if print
    fprintf('Mean RMSE: %g mm\n\n\n', 1000*mean(RMSEdistance));
end
