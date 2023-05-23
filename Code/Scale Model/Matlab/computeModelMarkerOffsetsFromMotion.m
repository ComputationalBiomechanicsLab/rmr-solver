function [offsetsLocal, markerErrsLocal] = ...
    computeModelMarkerOffsetsFromMotion(time, markersExp, names, Qs, coordNames, model)
% Givent the motion (coordinates as a function of time) of the model,
% compute the mean offset between model markers and experimental markers in  
% the local (body) frame for each marker in the model (in the order
% indicated in markerNames).
% USAGE [offsetsLocal, markerNames] = ...
% computeModelMarkerOffsetsFromMotion(time, markers, markerNames, Qs, coordNames, model)
% Ajay Seth and Ricardo Matias 2015

import org.opensim.modeling.*

% number of markers
nm = length(names);

nt = length(time);
nqik = length(coordNames);

state = model.initSystem();
coords = model.getCoordinateSet();
nqm = coords.getSize();

markers = model.getMarkerSet();

ground = model.getGround();

% convert Qs to radians
degToRadMat = ones(1, nqik);
for I = 1:nqik,
    if (coords.get(I-1).getMotionType() == 'Rotational')
        degToRadMat(I) = pi/180;
    end
end

Qs = (ones(nt,1)*degToRadMat) .* Qs;  

% compute the markers errors in body local frame for each marker
markerErrsLocal = zeros(nt, 3*nm);
emInBody = Vec3(); % experimental marker (em) in Body local frame
mmInBody = Vec3(); % model marker (mm) in Body local frame
%markersInGround = zeros(nt, 3*nm);
for ixt = 1:nt,
    % set the pose of the model base on the IK solution
    for ixq = 1:nqik,
        coords.get(coordNames(ixq)).setValue(state, Qs(ixt, ixq), ixq == nqik);
    end
    % verify that cooridnates are updated and changing the pose
    Pcom = model.calcMassCenterPosition(state);
    
    for ixm = 1:nm,
        marker = markers.get(names(ixm));
        % Double check to verify that model marker transformed to ground
        % matches reproted marker locations from inverse kinematics
        % engine.transformPosition(state, marker.getBody(),...
        %                        marker.getOffset(), ground, markerInGround);
        % pInG = [markerInGround.get(0) markerInGround.get(1) markerInGround.get(2)];                   
        % markersInGround(ixt, 3*(ixm-1)+1:3*(ixm-1)+3) = pInG;
        
        % Model marker location in the Body
        mmInBody = marker.get_location();
        mInB = [mmInBody.get(0) mmInBody.get(1) mmInBody.get(2)];
        
        % Transform experimental marker in ground to be expressed in the
        % marker Body's reference frame
        emInBody.set(0, markersExp(ixt, 3*(ixm-1)+1));
        emInBody.set(1, markersExp(ixt, 3*(ixm-1)+2));
        emInBody.set(2, markersExp(ixt, 3*(ixm-1)+3));
        %engine.transformPosition(state, ground, emInBody,...
        %                         marker.getBody(), emInBody);
        mf = marker.getParentFrame();
        emInBody = ground.findStationLocationInAnotherFrame(state, emInBody, mf);
        eInB = [emInBody.get(0) emInBody.get(1) emInBody.get(2)];
        
        % Error expressed in local marker's Body frame
        markerErrsLocal(ixt, 3*(ixm-1)+1:3*(ixm-1)+3) = eInB - mInB;
    end
end

offsetsLocal = mean(markerErrsLocal);





















