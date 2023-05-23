function [sumSqrError, model, RMSEdistance] = ...
    computeJointFitToMotion_reduced(x, time, markerData, markerNames, coordData, coordNames, fixedModel, print)
%function [RMSEdistance, model, sumSqrError] = ...
    %computeJointFitToMotion_reduced(x, time, markerData, markerNames, coordData, coordNames, fixedModel)
    
% Compute the error of a model's fit to experimental data according to the
% input parameters, p. Parameters could be scale factors and other model properties.
% USAGE error = computeJointFitToMotion(parameters, markerData, fixedModel, motion)
% Ajay Seth  2015
%
% modified Italo Belli (i.belli@tudelft.nl) 2022
% added print option:
% * print == True -> print information while optimizing
% * print == False -> do not print information while optimizing

import org.opensim.modeling.*

model = fixedModel.clone();
model.initSystem();

if(print)
    disp('ThoracoScapular Joint parameters:')
    disp(['ellipsoid radii: [',num2str(x(1)),' , ',num2str(x(2)),' , ',num2str(x(1)),']'])
    disp(['ellipsoid rotation: [ Y= ', num2str(x(3)),' , Z= ', num2str(x(4)), ' ]'])
    disp('')
end

% get scapulothoracic joint
scapThor = model.getJointSet().get(2);

% update the ellipsoid radii
radiiProp = scapThor.getPropertyByName('thoracic_ellipsoid_radii_x_y_z');
PropertyHelper.setValueVec3(x(1), radiiProp, 0);
PropertyHelper.setValueVec3(x(2), radiiProp, 1);
PropertyHelper.setValueVec3(x(1), radiiProp, 2);

% update the ellipsoid orientation (along Y and Z axis of the thorax_frame)
ellipsoidFrame = PhysicalOffsetFrame.safeDownCast(scapThor.getParentFrame());
ellipsoidOrientation = ellipsoidFrame.get_orientation;
orientationXAxis = ellipsoidOrientation.get(0);

ellipsoidFrame.set_orientation(Vec3(orientationXAxis, x(3), x(4)));


[~, markerErrsLocal] = ...
    computeModelMarkerOffsetsFromMotion(time, markerData, markerNames, coordData, coordNames, model);

sumSqrError = sum(sum(markerErrsLocal.^2,'omitnan'),'omitnan');
if(print)
    fprintf('sumSqrError: %g m\n', sumSqrError);
end

distances = computeDistancesFromDisplacementsXYZ(markerErrsLocal);
RMSEdistance = sqrt(mean((distances.^2),'omitnan'));
if(print)
    fprintf('Mean RMSE: %g mm\n', 1000*mean(RMSEdistance));
end
