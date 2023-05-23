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
%    Author:   Luca Modenese, July 2013                                   %
%    email:    l.modenese@sheffield.ac.uk                                 %
% ----------------------------------------------------------------------- %
%
% Given as INPUT a muscle OSMuscleName from an OpenSim model, this function
% returns the OUTPUT structure jointNameSet containing the OpenSim jointNames
% crossed by the OSMuscle.
%
% It works through the following steps:
%   1) extracts the GeometryPath
%   2) loops through the single points, determining the body they belong to
%   3) stores the bodies to which the muscle points are attached to
%   4) determines the nr of joints based on body indexes
%   5) stores the crossed OpenSim joints in the output structure named jointNameSet
%
% NB this function return the crossed joints independently on the
% constraints applied to the coordinates. Eg patello-femoral is considered as a
% joint, although in Arnold's model it does not have independent
% coordinates, but it is moved in dependency of the knee flexion angle.


function [jointNameSet,coordNameSet,coordjointSet,varargout] = getJointsSpannedByMuscleUpdated(osimModel, OSMuscleName)

% just in case the OSMuscleName is given as java string
OSMuscleName = char(OSMuscleName);

%useful initializations
JointSet= osimModel.getJointSet();
muscle  = osimModel.getMuscles().get(OSMuscleName);


%% Get list of bodies muscle attaches to
% Extracting the PathPointSet via GeometryPath
musclePath = muscle.getGeometryPath();
musclePathPointSet = musclePath.getPathPointSet();

% for loops to get the attachment bodies
n_body = 1;
muscleAttachBodies = '';
%muscleAttachIndex = [];
for n_point = 0:musclePathPointSet.getSize()-1

    % get the current muscle point
    currentAttachBody = char(musclePathPointSet.get(n_point).getBodyName());

    %Initialize
    if n_point ==0;
        previousAttachBody = currentAttachBody;
        muscleAttachBodies{n_body} = currentAttachBody;
        %muscleAttachIndex(n_body) = BodySet.getIndex(currentAttachBody);
        n_body = n_body+1;
    end

    % building a vectors of the bodies attached to the muscles
    if ~strncmp(currentAttachBody,previousAttachBody, size(char(currentAttachBody),2))
        muscleAttachBodies{n_body} = currentAttachBody;
        %muscleAttachIndex(n_body) = BodySet.getIndex(currentAttachBody);
        previousAttachBody = currentAttachBody;
        n_body = n_body+1;
    end
end

%% Get list of joints crossed by muscle
% Compare Distal- and Proximal muscle attachment bodies to Parent- and
% Child bodies of joints to determine if joints are crossed.

DistalBodyName = muscleAttachBodies{end};
ProximalBodyName= muscleAttachBodies{1};

% initialize output structures
NoDofjointNameSet = {};
jointNameSet = [];
coordNameSet = [];
coordjointSet = [];

JointName = [];
ParentBodyName = [];
ChildBodyName = [];

    % Loop through joints, get array with attachment bodies
    for i = 1:JointSet.getSize()
        Joint = JointSet.get(i-1);
        
        JointName = [JointName;{char(Joint.getName())}];
        
        ParentBodyName = [ParentBodyName;{char(Joint.getParentFrame().findBaseFrame().getName())}];
        ChildBodyName  = [ChildBodyName;{char(Joint.getChildFrame().findBaseFrame().getName())}];
        
    end
    
    
    %% Loop through Constraints, temporarily add joint without degrees of freedom for PointConstraint
    ConstraintNames = [];
    ConstraintSet = osimModel.getConstraintSet();
    nConstraints = ConstraintSet.getSize();
    
    if nConstraints > 0
        % I there are constraints, loop throught constraints
        for i = 0:nConstraints-1
        Constraint = ConstraintSet.get(i);
        % Get classname constraint
        ConstraintClass = Constraint.getConcreteClassName().toString();
        
        % only do stuff if it is a PointConstraint
            if strcmp(ConstraintClass,'PointConstraint')
            % Get bodies in constraint
            %Downcasted = PointConstraint().safeDownCast(Constraint);
            
            ConstraintName = char(Constraint.getName());
            
            Body1 = char(osimModel.getComponent(Constraint.getPropertyByName('socket_body_1').toString()).getName());
            Body2 = char(osimModel.getComponent(Constraint.getPropertyByName('socket_body_2').toString()).getName());
            
            
            % Introduce this constraint as a joint without degrees of
            % freedom
            JointName = [JointName;{ConstraintName}];
            ConstraintNames = [ConstraintNames; {ConstraintName}];
            
            ParentBodyName = [ParentBodyName;{Body1}];
            ChildBodyName = [ChildBodyName;{Body2}];
            
            else
                % other constraints not supported
                continue
            end
        
        end
    end
    
    %% Find joints crossed by muscle
    
    %% Single-jointed
     % Array filled with all the parent and child bodies per joint
    JointBodies = [ParentBodyName,ChildBodyName];
    
    % find distal and proximal bodies in JointBodies structure
    DistalPos=strcmp(JointBodies,DistalBodyName);
    ProximalPos=strcmp(JointBodies,ProximalBodyName);
    
    % check if muscle spans a single joint; its proximal and distal body are
    % identical to the joints parent and child bodies
    % THIS PART ADDS FUNCTIONALITY FOR SERRATUS ANTERIOR MUSCLES DUE TO
    % PROXIMAL NATURE SITES OF MUSCLE INSERTION SITE
    singleJoint = find(sum(ProximalPos+DistalPos,2)==2);
    
    % fill output structure if singleJoint is found
    jointNameSet = JointName(singleJoint);
    
    
    %% double-jointed
    % If no single-spanned joint is found for muscle, check for
    % double-jointedness by finding intermediate body. 
    % THIS PART ADDS FUNCTIONALITY FOR CONSTRAINT-JOINTS
    if isempty(singleJoint)
        % check joints where either the Distal- or Proximal body is found
        check = find(sum(ProximalPos+DistalPos,2)==1);
        
        % Remove Distal- and Proximal bodies from array to form array of
        % possible intermediate bodies
        checkParent = erase(ParentBodyName(check),[{DistalBodyName};{ProximalBodyName}]);
        checkChild = erase(ChildBodyName(check),[{DistalBodyName};{ProximalBodyName}]);
        
        % Find intermediate body and both spanned joints
        IntermediateBody = intersect(checkParent,checkChild);
        % Remove empty cells from resulting vector
        IntermediateBody = IntermediateBody(~cellfun('isempty',IntermediateBody));
        
        if ~isempty(IntermediateBody)
        jointNameSet = [JointName(check(find(strcmp(ChildBodyName(check),IntermediateBody))));JointName(check(find(strcmp(ParentBodyName(check),IntermediateBody))))];
        % disp(['Double-jointedness detected for muscle: ',OSMuscleName]);
        end
       
    end
    
    %% multi-jointed
    % THIS PART ADDS FUNCTIONALITY FOR JOINTS CROSSING MORE THAN TWO JOINTS
    if isempty(jointNameSet)
    CurrentBodyName = DistalBodyName;
    
    while ~strcmp(CurrentBodyName,ProximalBodyName)
        
        % Break loop if CurrentBodyName is "ground" and empty out
        % jointNameSet
        if strcmp(CurrentBodyName,'ground')
            jointNameSet = [];
            break
        end
        
        % Find joint with CurrentBodyName as ChildBodyName
        CurrentJoint = find(strcmp(CurrentBodyName,ChildBodyName));
        
        % In case multiple joints are found, only use first. Later ones are
        % constraints at the end of the evaluated vector.
        if length(CurrentJoint)>1
            warning(['Multiple valid crossed joints found for muscle ',OSMuscleName,'.\n Optimization might be invalid if ',char(JointName(CurrentJoint(end))),' is not a constraint-joint.'])
            CurrentJoint = CurrentJoint(1);
        end
        
        % Add JointName of CurrentJoint to jointNameSet
        jointNameSet = [jointNameSet; JointName(CurrentJoint)];
        
        % Update CurrentBodyName to ParentBodyName
        CurrentBodyName = ParentBodyName(CurrentJoint);
    end
    
    end
    

    % THIS HAS THE SAME FUNCTIONALITY AS CODE ABOVE BUT SWAPS PROXIMAL-
    % AND DISTAL BODY NAMES, SHOULD BE REDUNDANT
    % If no joints are found, start loop from other body 
    if isempty(jointNameSet)
        
        CurrentBodyName = ProximalBodyName;
        
  	while ~strcmp(CurrentBodyName,DistalBodyName)
        
        % Break loop if CurrentBodyName is "ground" and empty out
        % jointNameSet
        if strcmp(CurrentBodyName,'ground')
            jointNameSet = [];
            break
        end
        
        % Find joint with CurrentBodyName as ChildBodyName
        CurrentJoint = find(strcmp(CurrentBodyName,ChildBodyName));
        
        % In case multiple joints are found, only use first. Later ones are
        % constraints at the end of the evaluated vector.
        if length(CurrentJoint)>1
            warning(['Multiple valid crossed joints found for muscle ',OSMuscleName,'. Optimization might be invalid if ',char(JointName(CurrentJoint(end))),' is not a constraint-joint.'])
            CurrentJoint = CurrentJoint(1);
        end
        
        % Add JointName of CurrentJoint to jointNameSet
        jointNameSet = [jointNameSet; JointName(CurrentJoint)];
        
        % Update CurrentBodyName to ParentBodyName
        CurrentBodyName = ParentBodyName(CurrentJoint);
       
    end
   
    end
    
       
    
   %% Create output arrays
   
   for j = 1:length(jointNameSet)
       
       % Check if constraint-joints are crossed
       if ismember(jointNameSet(j),ConstraintNames)
           
           % disp([OSMuscleName,' crosses constraint-joint ', jointNameSet(j)])
           NoDofjointNameSet = [NoDofjointNameSet;{jointNameSet(j)}];
           % delete constraint-joint from jointNameSet
           jointNameSet = erase(jointNameSet,jointNameSet(j));
               
       else
        nCoordsJoint = JointSet.get(jointNameSet(j)).numCoordinates();
            if nCoordsJoint == 0
                 NoDofjointNameSet = [NoDofjointNameSet;{jointNameSet(j)}];
                else
                    for n = 0:nCoordsJoint-1
                    coordNameSet = [coordNameSet;{char(JointSet.get(jointNameSet(j)).get_coordinates(n))}];
                    coordjointSet = [coordjointSet; jointNameSet(j)];
                    end      
            end
       end
   end
    
   % remove all empty entries due to removed constraints from jointNameSet
   jointNameSet = jointNameSet(~cellfun('isempty',jointNameSet));
   
%%%%%%%%%%%

if isempty(jointNameSet)
    error(['No joint detected for muscle ',OSMuscleName]);
end
if  ~isempty(NoDofjointNameSet)
    for n_v = 1:length(NoDofjointNameSet)
        % display(['Joint ',NoDofjointNameSet{n_v},' has no dof.'])
    end
end

varargout = NoDofjointNameSet;
end