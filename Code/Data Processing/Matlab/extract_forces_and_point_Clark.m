function extract_forces_and_point_Clark(mot_output_name, trc_file_Clark, csv_file_Clark)
% This function is intended for extracting the forces (N and Nm) from
% Clark's dataset, together with the point of application, and save a .mot
% file that can be used to create an ExternalForce object with OpenSim.
%
% NOTE: it is still missing the correct conversion between the raw singals
% contained in the database and the true force values!!
%
% example of usage:
% extract_forces_and_point_Clark('3_exp.mot', 'U:\PTbot\DataShare\dataset_Clark\MoCap\TRC files\ACC\Processed\3_processed.trc' , 'U:\PTbot\DataShare\dataset_Clark\SMAP\SMAP 4 Raw\ACC\ACC\ACC\3.CSV')
%
% author: Italo Belli (i.belli@tudelft.nl) 2022

%% Load and process marker data
% first, we load the data from the trc file. We retain just the information
% about the last two markers (MCP5 and MCP2), and from those we build the 
% virtual marker representing the position of application of the force and
% moments.
[markersExp, timesExp, labelsExp, unitsExp] = readTRC(trc_file_Clark);

% check if the markers that we need are present in the trc file
if (max(strcmp(labelsExp,'MCP5')) * max(strcmp(labelsExp, 'MCP2')) == 0)
    error("No hand markers (MCP5 or MCP2) present in .trc file")
end

% convert the marker data to meters if needed
if unitsExp == "mm"
     markersExp = markersExp./1000;
     unitsExp = "m";
end

sample_time_trc = timesExp(2)-timesExp(1);
rate_trc = 1/sample_time_trc;

% extract the marker positions that we want (TODO: which reference frame is
% used here????)
mcp5 = markersExp(:, 49:51);
mcp2 = markersExp(:, 52:54);
clear markersExp

% find the position of the application point in the middle of the hand
application_point = (mcp2+mcp5)/2;

%% Load and process force and moment data
% load the csv file containing the force data
csv_content = table2array(readtable(csv_file_Clark, 'VariableNamingRule', 'preserve'));
csv_content = csv_content(2:end, :);

% extract the frame and subframe
frame_id = csv_content(:,1);
subframe_id = csv_content(:,2);

% get the rate at which the csv data
rate_csv = rate_trc*(max(subframe_id+1));

% extract signals related to forces (v_F) and moments (v_M) from the table
% They are most likely voltages and should be converted to values in N and
% Nm respectively
v_F = csv_content(:, 17:19);
v_M = csv_content(:, 20:22);

v_F_norm = sqrt(v_F(:,1).^2 + v_F(:,2).^2 + v_F(:,3).^2);
v_M_norm = sqrt(v_M(:,1).^2 + v_M(:,2).^2 + v_M(:,3).^2);

v_F_significant = v_F_norm(5*rate_csv:6*rate_csv);   % considering just the time window 5-6s as in Meszaros and Dickenson's paper

% at this point, we would need to find a way to transform the raw voltages
% into forces. This could be done by comparing the values of
% v_F_significant to the expected force, and finding a conversion
% coefficient

clear csv_content

% since the rate of the marker data and that of the force data is
% different (and higher) we need to resample the former.
padding_initial = [application_point(1,1)*ones(rate_trc,1), application_point(1,2)*ones(rate_trc,1), application_point(1,3)*ones(rate_trc,1)];
padding_end = [application_point(end,1)*ones(rate_trc,1), application_point(end,2)*ones(rate_trc,1), application_point(end,3)*ones(rate_trc,1)];
application_point_padded = [padding_initial; application_point; padding_end];

application_point_padded_resampled = resample(application_point_padded, round(rate_csv), round(rate_trc));
application_point_padded_resampled(1:rate_csv,:) = [];
application_point_padded_resampled(end-rate_csv+1:end,:) = [];
application_point_resampled = application_point_padded_resampled;

% also the timestamps need to be recreated
timestamps_resampled = linspace(timesExp(1), timesExp(end), size(application_point_resampled,1))';

% writing the .mot file
data = [timestamps_resampled, v_F, v_M, application_point_resampled];
columnNames = ["F_x" "F_y" "F_z" "M_x" "M_y" "M_y" "p_x" "p_y" "p_z"];
writeMotStoData(data, 1, columnNames, mot_output_name);