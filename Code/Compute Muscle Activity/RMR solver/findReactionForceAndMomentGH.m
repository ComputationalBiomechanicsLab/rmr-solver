function [force, moment]  = findReactionForceAndMomentGH(x,params)
% This function returns the joint reaction force and moment at the gleno-humeral
% (GH) joint, given the model state and the forces exerted by the Muscles (and
% CoorindateActuators). Most of the necessary inputs are passed as a single
% parameter struct, with the followign fields:
% * model: the OpenSim Model considered;
% * state: the State of the model
% * acts: actuators in the model, where the i-th element is retrieved as
%         acts(i+1) = ScalarActuator.safeDownCast(model.getActuators().get(i));
% * muscles: muscles present in the model, muscles = model.getMuscles()
% * AMuscForce: vector containing the muscle multipliers to find the muscle
%               force dependent from the activation level (for each muscle)
% * PMuscForce: vector containing the muscle multipliers to find the
%               passive muscle force (for each muscle)
% * glen: glenohumeral joint (model.getJointSet.get('GlenoHumeral');
% * useControls: flag to indicate whether to use .setControls (1) for the
%                CoordinateActuators, or to overwrite their actuation with
%                .setOverrideActuation (0). In the second case, remeber
%                that the actuation should be overrideable!
% * modelControls: model controls (model.getControls(state)) used if the
%                  corresponding flag is set to 1.
%
% Overall, the input to the function are:
% x: vector containing the muscle activations (for each muscle) and the
%    controls for the CoordinateActuators (for each one of them)
% params: collecting all the parameters presented above
%
% The required outputs are found as a linear combination of the muscle
% activations in the model, so this computation is true in the Static
% Optimization case.
%
% author: Italo Belli (i.belli@tudelft.nl) june 2022

import org.opensim.modeling.*;

model = params.model; 
state = params.state;
acts = params.acts;
muscles = params.muscles;
AMuscForce = params.AMuscForce;
PMuscForce = params.PMuscForce;
glen = params.glen;
useControls = params.useControls;

nMuscles = muscles.getSize();

% Inizialize the muscles to produce the required forces
MuscleForces = AMuscForce.*x(1:nMuscles)' + PMuscForce; % This calculates for the 33 muscles of the model
for k = 1:nMuscles
    muscle = muscles.get(k-1);
    muscle.setOverrideActuation(state, MuscleForces(k));
end

% Inizialize the CoordinateActuators to produce the required effect
if ~useControls  % if we want to set the force
    for k = nMuscles+1:length(acts)
        acts{k}.setOverrideActuation(state, x(k));
    end
else             % if we want to set the controls
    modelControls = params.modelControls;
    for k = nMuscles+1:length(acts)
        acts{k}.setControls(Vector(1, x(k)), modelControls);
    end
    model.realizeVelocity(state);
    model.setControls(state, modelControls);
end

% Realize the model to the acceleration stage
model.realizeAcceleration(state);

% get moment and force at the GlenoHumeral joint
moment = glen.calcReactionOnParentExpressedInGround(state).get(0).getAsMat();
force = glen.calcReactionOnParentExpressedInGround(state).get(1).getAsMat();
    