function [angle, Vec_H2GC] = get_glenoid_status(model, state)
% This function returns parameters describing the status of the
% glenohumeral joint in the thoracoscapular model
% Input:
% * model: opensim thoracoscapular model, that must be already provided
% with markers on the glenoid center, humerus head and glenoid edge 
% (in this order, and they should be the last ones in the markerset)
% * state: state of the model 
% Outputs:
% * angle: the maximum angle representing the cone in which the reaction forces must
%     be contained is returned
% * Vec_H2GC: 3D vector defined between the humeral head center (origin)
%     and the glenoid center. It is expressed in the ground frame

import org.opensim.modeling.*;

% get the markerset of the model and the number of markers
mkrs = model.getMarkerSet;
nmkrs = mkrs.getSize; 

% manually hardcode the markers that we want (last 3 in the MarkerSet)
G_Cent = mkrs.get(nmkrs-3); 
HH_Cent = mkrs.get(nmkrs-2); 
G_Edge = mkrs.get(nmkrs-1); 

% get the location in ground of the three markers
G_Cent_Loc = G_Cent.getLocationInGround(state).getAsMat()';
HH_Cent_Loc = HH_Cent.getLocationInGround(state).getAsMat()';
G_Edge_Loc = G_Edge.getLocationInGround(state).getAsMat()';

% define the vector from the glenoid center to the humerus head
Vec_H2GC = G_Cent_Loc - HH_Cent_Loc;

% define the vector from the glenoid edge to the humerus head
Vec_H2GE = G_Edge_Loc - HH_Cent_Loc;

% get the cosine fo the angle between the two vectors
CosTheta = max(min(dot(Vec_H2GC,Vec_H2GE)/(norm(Vec_H2GC)*norm(Vec_H2GE)),1),-1);

% find the maximum angle to be returned
angle = acosd(CosTheta);

% get additional informations about the glenohumeral joint
Glenoid_Rad = norm(G_Cent_Loc-G_Edge_Loc);
Head2Glen_dist = norm(HH_Cent_Loc-G_Cent_Loc);