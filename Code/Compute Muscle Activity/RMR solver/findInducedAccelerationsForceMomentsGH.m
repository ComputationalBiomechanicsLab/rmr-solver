function [simulatedAccelerations, force, moment] = findInducedAccelerationsForceMomentsGH(x,params)
% This function returns the simulated accelerations for each coordinate of
% a Model, given the model state and the forces exerted by the Muscles (and
% CoorindateActuators). It also returns the joint reaction force and moment 
% at the gleno-humeral (GH) joint. Most of the necessary inputs are passed 
% as a single parameter struct, with the following fields:
% * model: the OpenSim Model considered;
% * state: the State of the model
% * acts: actuators in the model, where the i-th element is retrieved as
%         acts(i+1) = ScalarActuator.safeDownCast(model.getActuators().get(i));
% * muscles: muscles present in the model, muscles = model.getMuscles()
% * numMuscles: number of muscles of the model, numMuscles = muscles.getSize();
% * useMuscles: flag to indicate whether the muscles are used (1) or they
%               are ignored in the computation of the accelerations (0)
%               If ignored, the next 2 fields can be omitted.
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
coords = params.coords;
coordNames = params.coordNames;
muscles = params.muscles;
numMuscles = params.numMuscles;
glen = params.glen;
useMuscles = params.useMuscles;
useControls = params.useControls;


% Inizialize the muscles to produce the required forces
if useMuscles
    AcMuscForce= params.AMuscForce;
    PasMuscForce = params.PMuscForce;
    MuscleForces = AcMuscForce.*x(1:numMuscles)' + PasMuscForce; % This calculates for the 33 muscles of the model
    for k = 1:numMuscles
        muscle = muscles.get(k-1);
        muscle.setOverrideActuation(state, MuscleForces(k));
    end
end

% Inizialize the CoordinateActuators to produce the required effect
if ~useControls  % if we want to set the force
    for k = numMuscles+1:length(acts)
        acts{k}.setOverrideActuation(state, x(k));
    end
else             % if we want to set the controls
    modelControls = params.modelControls;
    for k = numMuscles+1:length(acts)
        acts{k}.setControls(Vector(1, x(k)), modelControls);
    end
    model.realizeVelocity(state);
    model.setControls(state, modelControls);
end

% Realize the model to the acceleration stage
model.realizeAcceleration(state);

% retrieve the simulated accelerations for each coordinate
simulatedAccelerations = zeros(length(coordNames),1);
for j = 1:length(coordNames)
    coord = coords.get(coordNames{j});
    simulatedAccelerations(j) = coord.getAccelerationValue(state);
end

% get moment and force at the GlenoHumeral joint
moment = glen.calcReactionOnParentExpressedInGround(state).get(0).getAsMat();
force = glen.calcReactionOnParentExpressedInGround(state).get(1).getAsMat();