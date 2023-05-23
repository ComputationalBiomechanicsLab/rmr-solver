function simulatedAccelerations = findInducedAccelerations(x,params)
% This function returns the simulated accelerations for each coordinate of
% a Model, given the model state and the forces exerted by the Muscles (and
% CoorindateActuators). Most of the necessary inputs are passed as a single
% parameter struct, with the followign fields:
% * model: the OpenSim Model considered;
% * state: the State of the model
% * coords: the CoordinateSet of the model (retrieved as model.getCoordinateSet())
% * coordNames: a cell object of dimension numCoords x 1, filled with the
%               coord names as a char datatype
% * acts: actuators in the model, where the i-th element is retrieved as
%         acts(i+1) = ScalarActuator.safeDownCast(model.getActuators().get(i));
% * muscles: muscles present in the model, muscles = model.getMuscles()
% * useMuscles: flag to indicate whether the muscles are used (1) or they
%               are ignored in the computation of the accelerations (0)
% * AMuscForce: vector containing the muscle multipliers to find the muscle
%               force dependent from the activation level (for each muscle)
% * PMuscForce: vector containing the muscle multipliers to find the
%               passive muscle force (for each muscle)
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
% The accelerations for each of the coordinates are returned.
%
% author: Sagar Joshi (s.d.joshi@tudelft.nl) 2022

import org.opensim.modeling.*;

model = params.model;
state = params.state;
coords = params.coords;
coordNames = params.coordNames;
acts = params.acts;
muscles = params.muscles;
nMuscles = muscles.getSize();
useMuscles = params.useMuscles;
useControls = params.useControls;

% Inizialize the muscles to produce the required forces
if useMuscles
    AcMuscForce= params.AMuscForce;
    PasMuscForce = params.PMuscForce;
    MuscleForces = AcMuscForce.*x(1:nMuscles)' + PasMuscForce; % This calculates for the 33 muscles of the model
    for k = 1:nMuscles
        muscle = muscles.get(k-1);
        muscle.setOverrideActuation(state, MuscleForces(k));
    end
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

% retreive the simulated accelerations for each coordinate
simulatedAccelerations = zeros(length(coordNames),1);
for j = 1:length(coordNames)
    coord = coords.get(coordNames{j});
    simulatedAccelerations(j) = coord.getAccelerationValue(state);
end
