function [c ceq] = jntrxncon_OpenSimJRF(x,params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

import org.opensim.modeling.*;

model = params.model; 
state = params.state;
acts = params.acts;
nActs = params.nActs;
nMuscles = params.nMuscles;
AMuscForce = params.AMuscForce;
PMuscForce = params.PMuscForce;
maxAngle = params.maxAngle;
glen = params.glen;
Vec_H2GC = params.Vec_H2GC;
useControls = params.useControls;

% get the corresponding force for each muscle
Forces = AMuscForce.*x(1:nMuscles)' + PMuscForce;

% Apply the muscle force directly 
for i = 1:nMuscles
    acts{i}.overrideActuation(state, true); 
    acts{i}.setOverrideActuation(state, Forces(i)); 
end

% Inizialize the CoordinateActuators to produce the required effect
if ~useControls  % if we want to set the force
    for k = nMuscles+1:length(acts)
        acts{k}.setOverrideActuation(state, x(k));
    end
else             % if we want to set the controls
    modelControls = params.modelControls;
    for k = nMuscles+1:nActs
        acts{k}.setControls(Vector(1, x(k)), modelControls);
    end
    model.realizeVelocity(state);
    model.setControls(state, modelControls);
end

model.realizeAcceleration(state) ; % moves stage back to acceleration after setting control function

% computing the reaction force direction
v = glen.calcReactionOnParentExpressedInGround(state).get(1).getAsMat();

% evaluate the relative angle between the reaction force and Vec_H2GC
cosTheta = max(min(dot(Vec_H2GC,v)/(norm(Vec_H2GC)*norm(v)),1),-1);
rel_angle = real(acosd(cosTheta));

% value of the constraint violation
c = (rel_angle/maxAngle)^2 - 1;
ceq = 0; 
