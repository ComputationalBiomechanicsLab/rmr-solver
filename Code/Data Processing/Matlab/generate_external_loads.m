function external_loads = generate_external_loads(force_magnitude, force_direction_in_ground, application_point_trajectory_in_ground, body_of_application, frequency, loads_name, file_name)
% this function allows to create and save an ExternalLoads object given:
% * force_magnitude: a (constant) value of force in N to be applied
% * force_direction_in_ground: direction of the force expressed in the
%     ground reference frame, with sign (string like '+X', '-Z', ...)
% * application_point_trajectory_in_ground: array of Nx3, describing 3D
%     trajectory (with N points) of the point at which the force is applied
% * body_of_applicatio: string with the name of the body to which the force
%     will be applied
% * frequency: frequency of sampling of the force value (in Hz)
% * loads_name: name of the ExternalLoads object created
% * file_name: name of the file to be created
% 
% The output is:
% * external_loads: the ExternalLoads object created
% Also, an XML file (saved in the working directory, as file_name.xml) is
% created and can be just loaded into an Inverse Dynamic tool (or similar tools) 
% as idtool.setExternalLoadsFileName(file_name.xml)
%
% author: Italo Belli (i.belli@tudelft.nl) 2022

import org.opensim.modeling.*;

%% Create the .mot file where forces and moments are stored

% get the length of the trajectory
len_trj = size(application_point_trajectory_in_ground, 1);

% create the timestamps for the experiment
time_stamps = (1:len_trj)*1/frequency;

% get the force direction in ground
force_sign = force_direction_in_ground(1);
force_axis = force_direction_in_ground(2);

if strcmp(force_sign,'-')
    force_magnitude = - force_magnitude;
end

force_vector = ones(len_trj,1)*force_magnitude;

force_moments_matrix = zeros(len_trj, 6);

if strcmp(force_axis, 'X')
    force_moments_matrix(:,1) = force_vector;
elseif strcmp(force_axis, 'Y')
    force_moments_matrix(:,2) = force_vector;
elseif strcmp(force_axis, 'Z')
    force_moments_matrix(:,3) = force_vector;
end

data = [time_stamps', force_moments_matrix, application_point_trajectory_in_ground];

% create  the motion file to be used in the ExternalLoads .xml file
columnNames = ["F_x" "F_y" "F_z" "M_x" "M_y" "M_z" "p_x" "p_y" "p_z"];

writeMotStoData(data, 1, columnNames, file_name);

mot_file_name = append(file_name, '.mot');

% creating storage object
force_storage = Storage(mot_file_name, false);
force_storage.setName(mot_file_name(1:end-4));
force_storage.print(mot_file_name);

%% Create the .xml file
external_loads = ExternalLoads();
external_loads.setName(loads_name)
external_loads.setDataFileName(mot_file_name);

% create the external force object
external_force = ExternalForce();
external_force.setName('Force');
external_force.set_applied_to_body(body_of_application);
external_force.set_force_expressed_in_body('ground');
external_force.set_point_expressed_in_body('ground');
external_force.set_force_identifier('F_');
external_force.set_point_identifier('p_');
external_force.set_torque_identifier('M_');
external_force.set_data_source_name(mot_file_name);
external_force.setDataSource(force_storage);

% add the external force to the external loads file
external_loads.cloneAndAppend(external_force);

xml_file_name = append(file_name, '.xml');
external_loads.print(xml_file_name);