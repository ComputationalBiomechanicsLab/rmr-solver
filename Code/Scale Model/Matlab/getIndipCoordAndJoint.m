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
%    Author:   Luca Modenese, August 2014                                 %
%    email:    l.modenese@sheffield.ac.uk                                 %
% ----------------------------------------------------------------------- %
%
% Function that given a dependent coordinate finds the independent
% coordinate and the associated joint. The function assumes that the
% constraint is a CoordinateCoupleConstraint as used by Arnold, Delp and
% LLLM. The function can be useful to manage the patellar joint for instance.

function [ind_coord_name, ind_coord_joint_name] = getIndipCoordAndJoint(osimModel, constraint_coord_name)

% importing OpenSim libraries
import org.opensim.modeling.*
% Load scapulothoracic joint library
Model.LoadOpenSimLibrary('C:\Users\joran\Documents\OpenSim\4.1\Models\ShoulderDatabase ThoracoscapularShoulderModel\ScapulothoracicJointPlugin40\WinX64\\ScapulothoracicJointPlugin40_WinX64');
 
% in case constrained coordinate is java string. Quality of life change
% added
constraint_coord_name = char(constraint_coord_name);

% init state
s = osimModel.initSystem;

% get coordinate
constraint_coord = osimModel.getCoordinateSet.get(constraint_coord_name);

% double check: if not constrained then function returns
if ~(constraint_coord.isConstrained(s))
    display([char(constraint_coord.getName),' is not a constrained coordinate.'])
    return
end

% otherwise search through the constraints
for n = 0:osimModel.getConstraintSet.getSize-1

    % get current constraint
    curr_constr = osimModel.getConstraintSet.get(n);
    
    % this function assumes that the constraint will be a coordinate
    % coupler contraint ( Arnold's model and LLLM uses this)
    % cast down constraint
    curr_constr_casted = CoordinateCouplerConstraint.safeDownCast(curr_constr);

    % get dep coordinate and check if it is the coord of interest
    %dep_coord_name = char(curr_constr_casted.getDependentCoordinateName());

    if strncmp(constraint_coord_name,dep_coord_name,length(constraint_coord_name))
        %     if curr_constr_casted.getIndependentCoordinateNames().getSize
        ind_coord_name_set = curr_constr_casted.getIndependentCoordinateNames();

        % extract independent coordinate and independent joint to which the
        % coordinate refers
        if ind_coord_name_set.getSize==1
            ind_coord_name = curr_constr_casted.getIndependentCoordinateNames().get(0);
            ind_coord_joint_name = osimModel.getCoordinateSet.get(ind_coord_name).getJoint.getName;
            return
        elseif ind_coord_name_set.getSize>1
            error('getIndipCoordAndJoint.m. The CoordinateCouplerConstraint has more than one indipendent coordinate and this is not managed by this function yet.')
        end
    end
end