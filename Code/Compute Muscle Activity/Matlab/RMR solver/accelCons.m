function [c,ceq] = accelCons(x,params, has_coord_acts)
% This function compares the accelerations obtained by a forward dynamics
% simulation of the model with the desired acceleration provided (that
% comes from the trial).

import org.opensim.modeling.*;

state = params.state;
model = params.model;
coords = params.coords;
coordNames = params.coordNames;
desiredAccelerations = params.AccelData;
muscles = params.muscles;
nMuscles = params.nMuscles;

if has_coord_acts    % check presence of coord actuators
    nActs = params.nActs;
    acts = params.acts;
    model_controls = params.model_controls;
else
    nActs = 0;
end

AcMuscForce= params.AMuscForce;
PasMuscForce = params.PMuscForce;

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

c = [];
ceq = desiredAccelerations - simulatedAccelerations;
% ceq = (desiredAccelerations - simulatedAccelerations).^2;    % performs much worse...
