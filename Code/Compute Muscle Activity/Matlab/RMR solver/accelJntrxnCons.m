function [c,ceq] = accelJntrxnCons(x,params, has_coord_acts)
% This function compares the accelerations obtained by a forward dynamics
% simulation of the model with the desired acceleration provided, for the
% timestep considered. The inputs are:
% - x: the vector containing muscle activations and coordinateActuators
%      controls (if those are present on the model)
% - params: struct containing data required (see below) like the model
%           considered, the coodinate set, the acceleration ground truths...
% - has_coord_acts : flag that discriminates whether coordActs are present
%                    or not in the model
%
% The acceleration data (simulated and ground truth) are compared, and the
% value of the joint reaction constraint between humerus head and glenoid
% is also returned.
%
% author: Italo Belli (i.belli@tudelft.nl) June 2022


import org.opensim.modeling.*;

model = params.model;                        % model (assumed to be a thoracoscapular shoulder model with markers already)
state = params.state;                        % state of the model
coords = params.coords;                      % coordinates of the model
coordNames = params.coordNames;              % coordinate names
desiredAccelerations = params.AccelData;     % accelerations ground truth
muscles = params.muscles;                    % muscles in the model
nMuscles = params.nMuscles;                  % number of muscles in the model
AcMuscForce = params.AMuscForce;             % active muscle force vector (for all the muscles)
PasMuscForce = params.PMuscForce;            % passive muscle force vector (for all the muscles)
maxAngle = params.maxAngle;                  % maximum angle allowed between for the reaction force on the glenoid
glen = params.glen;                          % glenoid joint

if has_coord_acts    % check presence of coord actuators
    nActs = params.nActs;
    acts = params.acts;
    model_controls = params.model_controls;
else
    nActs = 0;
end

% set the coordinate actuator controls
for k = 1:nActs
    acts(nMuscles+k).setControls(Vector(1, x(nMuscles+k)), model_controls);
end

MuscleForces = AcMuscForce.*x(1:nMuscles)'+ PasMuscForce;                 % This calculates the total force/tension in the muscle, for all the muscles in the model

% override the muscle forces
for k = 1:nMuscles
    if has_coord_acts
        acts(k).setOverrideActuation(state, MuscleForces(k));
    else
        muscle = muscles.get(k-1);
        muscle.setOverrideActuation(state, MuscleForces(k));
    end
end

model.realizeVelocity(state);

if has_coord_acts
    model.setControls(state, model_controls);
end

model.realizeAcceleration(state);

simulatedAccelerations = zeros(1, length(coordNames));

for j = 1:length(coordNames)
    coord = coords.get(coordNames{j});
    simulatedAccelerations(j) = coord.getAccelerationValue(state);
end

ceq = desiredAccelerations - simulatedAccelerations;
% ceq = (desiredAccelerations - simulatedAccelerations).^2;    % performs much worse...


% computing the reaction force direction
v = glen.calcReactionOnParentExpressedInGround(state).get(1).getAsMat();

% get the vector Vec_H2GC between humeral head and the glenoid center
[~, Vec_H2GC] = get_glenoid_status(model, state);

% evaluate the relative angle between the reaction force and Vec_H2GC
cosTheta = max(min(dot(Vec_H2GC,v)/(norm(Vec_H2GC)*norm(v)),1),-1);
rel_angle = real(acosd(cosTheta));

% value of the constraint violation
c = (rel_angle/maxAngle)^2 - 1;
