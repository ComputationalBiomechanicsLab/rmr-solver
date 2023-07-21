function external_force = generate_external_force_3D(params, time)
% This function is used to generate a 3D external force to be applied to an
% OpenSim model. The function produces a .mot file containing the
% information of the required force, that can be later applied to the
% model.
% Inputs are:
% * params: struct containing information on the external force. It
%           should contain the following fields:
%           * EF_point_of_application: point on which the force is applied 
%                                     (the body is specified in this
%                                     function)
%           * EF_filename: name of the .mot file where the force
%                           is saved
%           * EF.x: x component of the force
%           * EF.y: y component of the force
%           * EF.z: z component of the force (all must have the same length)
%
% * time: time steps corresponing to the x,y,z components of the force
%        (dimensions must be consistent, at each time isntant must be
%         associated one triplet x,y,z for the force value)
%
% Author: Irene Beck 2023

import org.opensim.modeling.*;

point_of_application = params.EF_point_of_application;
file_name = params.EF_filename;

% Set point of application. Define body and relative position in body-fixed
% frame of force application
if strcmp(point_of_application, 'elbow')
    body = 'humerus';
    marker_location = [0.00109386 -0.240965 -0.00999386];
end

if strcmp(point_of_application, 'thorax')
    body = 'thorax';
    marker_location = [0.0 0.0 0.0];
end

% TODO: add here the identifier that you are willing to use, and specify
% the 
if strcmp(point_of_application, 'your_identifier')
    body = 'hand';                    % I guess you would apply the force to the hand/wrist?
    marker_location = [0.0 0.0 0.0];  % To be modified with the correct position
                                      % of the marker in the body frame
end

len_traj = size(time, 1);

force_vector_x = repelem(params.EF.x, floor(len_traj/length(params.EF.x)));
force_vector_y = repelem(params.EF.y, floor(len_traj/length(params.EF.y)));
force_vector_z = repelem(params.EF.z, floor(len_traj/length(params.EF.z)));

% here we are setting the moments to be all 0. In case this is not the
% case, they should be set here too!
force_moments_matrix = zeros(len_traj, 6);

force_moments_matrix(1:numel(force_vector_x),1) = force_vector_x;       % force along x axis    
force_moments_matrix(1:numel(force_vector_y),2) = force_vector_y;       % force along y axis
force_moments_matrix(1:numel(force_vector_z),3) = force_vector_z;       % force along z axis

% the moments should be added here, as the last 3 columns of the matrix

% The marker position is kept fixed, as it is defined in the body frame
force_point = repmat(marker_location, [len_traj,1]);

data = [time, force_moments_matrix, force_point];

% create  the motion file to be used in the ExternalLoads .xml file
columnNames = ["F_x" "F_y" "F_z" "M_x" "M_y" "M_z" "p_x" "p_y" "p_z"];

writeMotStoData(data, 1, columnNames, file_name);

mot_file_name = append(file_name, '.mot');

% creating storage object
force_storage = Storage(mot_file_name, false);
force_storage.setName(mot_file_name(1:end-4));
force_storage.print(mot_file_name);

% Create external force object
external_force = ExternalForce();
external_force.setName('Force');
external_force.set_applied_to_body(body);
external_force.set_force_expressed_in_body("ground");       % TODO: here it depends from the coordinate frame in which the force is expressed
external_force.set_point_expressed_in_body(body);
external_force.set_force_identifier('F_');
external_force.set_point_identifier('p_');
external_force.set_torque_identifier('M_');
external_force.set_data_source_name(mot_file_name);         % this line and the next "link" the storage file to the external force
external_force.setDataSource(force_storage);    

xml_file_name = append(file_name, '.xml');
external_force.print(xml_file_name);                        % print the xml file where the ExternalForce is defined