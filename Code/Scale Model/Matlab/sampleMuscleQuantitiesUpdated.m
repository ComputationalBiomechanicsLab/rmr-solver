%-------------------------------------------------------------------------%
% Copyright (c) 2015 Modenese L., Ceseracciu, E., Reggiani M., Lloyd, D.G.%
%                                                                         %
% Licensed under the Apache License, Version 2.0 (the "License");         %
% you may not use this file except in compliance with the License.        %
% You may obtain a copy of the License at                                 %
% http://www.apache.org/licenses/LICENSE-2.0.                             %
%                                                                         %
% Unless required by applicable law or agreed to in writing, software     %
% distributed under the License is distributed on an "AS IS" BASIS,       %
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or         %
% implied. See the License for the specific language governing            %
% permissions and limitations under the License.                          %
%                                                                         %
%    Author: Luca Modenese, March 2013                                    %
%                        revised July 2014 (improved management of joints)%
%                        modified May 2015 (management of max resolution) %
%    email:    l.modenese@sheffield.ac.uk                                 %
%-------------------------------------------------------------------------%
%
% Given as INPUT an OpenSim muscle OSModel, OSMuscle, a muscle variable and a nr
% of evaluation points this function returns as
% musOutput a vector of the muscle variable of interest
% obtained by sampling the ROM of the joint spanned by the muscle in
% N_EvalPoints evaluation points.
% For multidof joint the combinations of ROMs are considered.
% For multiarticular muscles the combination of ROM are considered.
% The script is totally general because based on generating strings of code
% correspondent to the encountered code. The strings are evaluated at the end.
%
% IMPORTANT 1
% The function can decrease the N_EvalPoints if there are too many dof
% involved (ASSUMPTION is that a better sampling will be vanified by the
% huge amount of data generated). This is an option that can be controlled
% by the user by deciding if to use it or not (setting limit_discr = 0/1)
% and, in case the sampling is limited, by setting a lower limit to the
% discretization.
%
% IMPORTANT 2
% Another check is done on the dofs: only INDEPENDENT coordinate are
% considered. This is fundamental for patellofemoral joint that both in
% LLLM and Arnold's model are constrained dof, dependent on the knee
% flexion angle. This function assumes to be used with Arnold's model.
% currentState is initialState;
%
% IMPORTANT 3
% At the purposes of the muscle optimizer it is important that here the
% model is initialize every time. So there are no risk of working with an
% old state (important for Schutte muscles for instance, where we observed
% that it was necessary to re-initialize the muscle after updating the
% muscle parameters).

%% NOTE: Changes in this file:
% - A setting was added to limitting the amount of evaluated model poses 
% per muscle to limit runaway computation times for muscles actuating many
% degrees of freedom whilst ensuring enough evaluations are made in
% optimization.
%
% - The original getJointsSpannedByMuscle.m function has been replaced by
% the new JointsSpannedByMuscle.m function. This function serves the same
% purpose and functions with OpenSim model 4.0+   
%
% Alterations by Joran de Vet

function musOutput = sampleMuscleQuantitiesUpdated(osimModel,OSMuscle,muscleQuant, N_EvalPoints,nPosesLimit)  %#ok<STOUT>

%======= SETTINGS ======
% limit (1) or not (0) the discretization of the joint space sampling
limit_discr = 0;
% minimum angular discretization
min_increm_in_deg = 2.5;
% limit of number of evaluated poses per evaluated muscle
% nPosesLimit = 1000;
%=======================

% initialize the model
currentState = osimModel.initSystem();


%% NOTE: this is a new matlab function written for the same purpose
% getting the joint crossed by a muscle
[muscleCrossedJointSet,muscleCrossedCoordSet,CoordsJointsList] = getJointsSpannedByMuscleUpdated(osimModel, OSMuscle);


%% NOTE: code herafter is added to limit runaway of the number of poses evaluated for muscles that actuated multiple degrees of freedom
% For example: the latissimus dorsi in the thoracoscapular shoulder mondel
% crosses 6 degrees of freedom.

N_Coords = size(muscleCrossedCoordSet,1);

while (N_EvalPoints^N_Coords) > nPosesLimit
    N_EvalPoints = N_EvalPoints - 1;
end

disp([char(OSMuscle.getName()),' in model: ',char(osimModel),' influences ',num2str(N_Coords),' DOFs, ',num2str(N_EvalPoints),' EvalPoints used.'])

% index for effective dofs
n_dof = 1;
DOF_Index = [];
for n_joint = 1:size(muscleCrossedJointSet,1)

    % current joint
    curr_joint = muscleCrossedJointSet{n_joint};

    %%%%%%%%%%%%%%%%%%%%%% NOTE: joint.getCoordinateSet() no longer
    %%%%%%%%%%%%%%%%%%%%%% functions from OpenSim 4.0 and onwards. Original
    %%%%%%%%%%%%%%%%%%%%%% code is commented, workaround written below
    % get CoordinateSet for the crossed joint
    % curr_joint_CoordinateSet = osimModel.getJointSet().get(curr_joint).getCoordinateSet();
    curr_joint_CoordinateSet = muscleCrossedCoordSet((strcmp(char(curr_joint),CoordsJointsList)==1));

    % Initial estimation of the nr of Dof of the CoordinateSet for that
    % joint before checking for locked and constraint dofs.
    %nDOF = osimModel.getJointSet().get(curr_joint).getCoordinateSet().getSize();
    nDOF = size(curr_joint_CoordinateSet,1);

    % skip welded joint and removes welded joint from muscleCrossedJointSet
    if nDOF == 0;
        continue;
    end

    % calculating effective dof for that joint
    effect_DOF = nDOF;
    for n_coord = 0:nDOF-1

        % get coordinate
        %curr_coord = curr_joint_CoordinateSet.get(n_coord);
        %curr_coord_name = char(curr_coord.getName());
        curr_coord_name = curr_joint_CoordinateSet{n_coord+1};
        curr_coord = osimModel.getCoordinateSet().get(curr_coord_name);

        % skip dof if locked
        if curr_coord.getLocked(currentState) == 1
            continue;
        end

        % if coordinate is constrained then the independent coordinate and
        % associated joint will be listed in the sampling "map"
        if curr_coord.isConstrained(currentState) && ~curr_coord.getLocked(currentState)
            constraint_coord_name = curr_coord_name;
            % finding the independent coordinate
            [ind_coord_name, ind_coord_joint_name] = getIndipCoordAndJoint(osimModel, constraint_coord_name); %#ok<NASGU>
            % updating the coordinate name to be saved in the list
            curr_coord_name = ind_coord_name;
            effect_DOF = effect_DOF-1;
            % ignoring constrained dof if they point to an independent
            % coordinate that has already been stored
            if sum(ismember(DOF_Index, osimModel.getCoordinateSet.getIndex(curr_coord_name)))>0
                continue
            end
            % skip dof if independent coordinate locked (the coord
            % correspondent to the name needs to be extracted)
            if osimModel.getCoordinateSet.get(curr_coord_name).getLocked(currentState)
                continue;
            end
        end

        % NB: DOF_Index is used later in the string generated code.
        % CRUCIAL: the index of dof now is model based ("global") and
        % different from the joint based used until now.
        DOF_Index(n_dof) = osimModel.getCoordinateSet.getIndex(curr_coord_name);

        % necessary update/reload the curr_coord to avoid problems with
        % dependent coordinates
        curr_coord = osimModel.getCoordinateSet.get(DOF_Index(n_dof));

        % Getting the values defining the range
        jointRange(1) = curr_coord.getRangeMin();
        jointRange(2) = curr_coord.getRangeMax();

        % Storing range of motion conveniently
        CoordinateBoundaries(n_dof) = {jointRange}; %#ok<NASGU>

        % increments in the variables when sampling the mtl space.
        % Increments are different for each dof and based on N_eval.
        %Defining the increments
        degIncrem(n_dof) = (jointRange(2)-jointRange(1))/(N_EvalPoints-1);

        % limit or not the discretization of the joint space sampling
        if limit_discr==1
            % a limit to the increase can be set though
            if degIncrem(n_dof)<(min_increm_in_deg/180*pi)
                degIncrem(n_dof)=(min_increm_in_deg/180*pi); %#ok<*AGROW>
            end
        end

        % updating index of list of dof
        n_dof = n_dof+1;
    end
end


% initializes the counter to save the results
iter = 1;  %#ok<NASGU>

% assigns an initial and a final value variable for each dof X
% calling them setAngleStartDofX and setLimitDof respectively.
for n_instr = 1:size(DOF_Index,2);
    % setting up the variables
    eval(['setAngleStartDof',num2str(n_instr),' = CoordinateBoundaries{',num2str(n_instr),'}(1);'])
    eval(['setLimitDof',num2str(n_instr),' = CoordinateBoundaries{',num2str(n_instr),'}(2);'])
end


% setting up for loops in order to explore all the possible combination of
% joint angles (looping on all the dofs of each joint for all the joint
% crossed by the muscle).
% The model pose is updated via: " coordToUpd.setValue(currentState,setAngleDof)"
% The right dof to update is chosen via: "coordToUpd = osimModel.getCoordinateSet.get(n_instr)"
stringToExcute1 = '';
for n_instr = 1:size(DOF_Index,2);
    stringToExcute1 = [stringToExcute1,[' for setAngleDof',num2str(n_instr),'=setAngleStartDof',num2str(n_instr),':degIncrem(',num2str(n_instr),'):setLimitDof',num2str(n_instr),';   coordToUpd = osimModel.getCoordinateSet.get(DOF_Index(',num2str(n_instr),'));    coordToUpd.setValue(currentState,setAngleDof',num2str(n_instr),');']]; %#ok<AGROW>
end

% calculating muscle length for the muscle
switch muscleQuant
    case 'MTL'
            stringToExcute2 = 'musOutput(iter) = OSMuscle.getGeometryPath.getLength(currentState);';
    case 'LfibNorm'
            stringToExcute2 =  ['OSMuscle.setActivation(currentState,1.0);',...
                                ' osimModel.equilibrateMuscles(currentState);',...
                                ' musOutput(iter) = OSMuscle.getNormalizedFiberLength(currentState);'];
    case 'Lten'
            stringToExcute2 =  ['OSMuscle.setActivation(currentState,1.0);',...
                                ' osimModel.equilibrateMuscles(currentState);',...
                                ' musOutput(iter) = OSMuscle.getTendonLength(currentState);'];
    case 'Ffib'
            stringToExcute2 =  ['OSMuscle.setActivation(currentState,1.0);',...
                                ' osimModel.equilibrateMuscles(currentState);',...
                                ' musOutput(iter) = OSMuscle.getActiveFiberForce(currentState);'];
    case 'all'
            stringToExcute2 = [ 'OSMuscle.setActivation(currentState,1.0);',...
                                ' osimModel.equilibrateMuscles(currentState);',...
                                'musOutput(iter,1) = OSMuscle.getGeometryPath.getLength(currentState);',...
                                ' musOutput(iter,2) = OSMuscle.getNormalizedFiberLength(currentState);',...
                                ' musOutput(iter,3) = OSMuscle.getTendonLength(currentState);',...
                                ' musOutput(iter,4) = OSMuscle.getActiveFiberForce(currentState);',...
                                ' musOutput(iter,5) = OSMuscle.getPennationAngle(currentState);'];
end

% updating iteration
stringToExcute3 =['iter = iter+1;',...
%     'angles.colheaders(iter-1) = {char(osimModel.getCoordinateSet.get(DOF_Index(',num2str(n_instr),')).getName)};'...
%     'angles.data(iter-1,1) = setAngleDof',num2str(n_instr),';'...
    ''];

% closing the loops started in stringToExcute1 with appropriate nr of end
stringToExcute4 = '';
for n_instr = 1:size(DOF_Index,2);
    stringToExcute4 =[stringToExcute4,'end;']; %#ok<AGROW>
end

% assembling the code string
Code = [stringToExcute1,stringToExcute2,stringToExcute3,stringToExcute4];

% evaluate the code
eval(Code)

end
