function [sumSqrError, model, RMSEdistance] = ...
    computeJointFitToMotion(p, time, markerData, markerNames, coordData, coordNames, fixedModel, print)
%function [RMSEdistance, model, sumSqrError] = ...
    %computeJointFitToMotion(p, time, markerData, markerNames, coordData, coordNames, fixedModel, print)
    
% Compute the error of a model's fit to experimental data according to the
% input parameters, p. Parameters could be scale factors and other model properties.
% USAGE error = computeJointFitToMotion(parameters, markerData, fixedModel, motion)
% Ajay Seth  2015
%
% modified Italo Belli (i.belli@tudelft.nl) 2022
% added print option:
%
% * print == True -> print information while optimizing
% * print == False -> do not print information while optimizing

import org.opensim.modeling.*

model = fixedModel.clone();
model.initSystem();

x = p(1:15);
if(print)
    disp('ThoracoScapular Joint parameters:')
    disp(['Translation in parent frame: [',num2str(x(1)),', ',num2str(x(2)),', ',num2str(x(3)),']'])
    disp(['Orientation in parent frame: [',num2str(x(4)),', ',num2str(x(5)),', ',num2str(x(6)),']'])
    disp(['Translation in child frame: [',num2str(x(7)),', ',num2str(x(8)),', ',num2str(x(9)),']'])
    disp(['Orientation in child frame: [',num2str(x(10)),', ',num2str(x(11)),', ',num2str(x(12)),']'])
    disp(['Thoracoscapular elipsoid radii: [',num2str(x(13)),', ',num2str(x(14)),', ',num2str(x(15)),']'])
    disp('')
end

% Edit scapulothoracic joint parameters
scapThor = model.getJointSet().get(2);
pf = PhysicalOffsetFrame.safeDownCast(scapThor.getParentFrame());
cf = PhysicalOffsetFrame.safeDownCast(scapThor.getChildFrame());
pf.set_translation(Vec3(x(1), x(2), x(3)));
pf.set_orientation(Vec3(x(4), x(5), x(6)));
cf.set_translation(Vec3(x(7), x(8), x(9)));
cf.set_orientation(Vec3(x(10), x(11), x(12)));

radiiProp = scapThor.getPropertyByName('thoracic_ellipsoid_radii_x_y_z');
PropertyHelper.setValueVec3(x(13), radiiProp, 0);
PropertyHelper.setValueVec3(x(14), radiiProp, 1);
PropertyHelper.setValueVec3(x(15), radiiProp, 2);

% optionally modify the glenohumeral joint axes 
if (length(p) > 16) 
    x = p(16:end);
    glenoHum = model.getJointSet().get(3);
    pf = PhysicalOffsetFrame.safeDownCast(glenoHum.getParentFrame());
    cf = PhysicalOffsetFrame.safeDownCast(glenoHum.getChildFrame());
    pf.set_translation(Vec3(x(1), x(2), x(3)));
    pf.set_orientation(Vec3(x(4), x(5), x(6)));
    cf.set_translation(Vec3(x(7), x(8), x(9)));
    cf.set_orientation(Vec3(x(10), x(11), x(12)));
    
    if(print)
        disp('GlenoHumeral Joint parameters:')
        disp(['Translation in parent frame: [',num2str(x(1)),', ',num2str(x(2)),', ',num2str(x(3)),']'])
        disp(['Orientation in parent frame: [',num2str(x(4)),', ',num2str(x(5)),', ',num2str(x(6)),']'])
        disp(['Translation in child frame: [',num2str(x(7)),', ',num2str(x(8)),', ',num2str(x(9)),']'])
        disp(['Orientation in child frame: [',num2str(x(10)),', ',num2str(x(11)),', ',num2str(x(12)),']'])
        disp('')
    end
end



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
