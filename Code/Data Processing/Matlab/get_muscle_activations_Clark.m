% This script processes the EMG data coming from the experiments recorded
% in Clark's dataset, and saves the corresponding muscle activation
% histories in a user-chosen folder. It processes all of the 120 CSV files
% that are associated to the subject, and saves a corresponding .mat file
% storing the muscle activations.
% The processing of the EMG data is done according to what described at
% section 2.4 of the paper "The effects of hand force variation on shoulder
% muscle activation during submaximal exertions" (Meszaros et al., 2017).
% Please check the USER INPUTS section of the code to specify the required
% data.
% Note that we first check if the required folder (containing the results)
% already exists. In case, this script does not run and indicates that the
% results are present already.
%
% REMARK : the naming of the muscles is fixed in the function
% get_MVC_values_Clark.m
%
% The muscles that are present in the dataset are (in order):
% 'Anterior Deltoid', 'Middle Deltoid', 'Posterior Deltoid', 'Biceps', 
% 'Triceps', 'Infra, 'Supra', 'PecMajorC', 'PecMajorS', 'Lats', 'Sert', 
% 'TrapLow', 'TrapMid', 'TrapUp'.

clear; clc;

%% USER INPUTS
path_to_subject_folder = 'U:\PTbot\DataShare\dataset_Clark\SMAP\SMAP 4 Raw\ACC\ACC\ACC';
path_to_results = 'U:\PTbot\DataShare\dataset_Clark\SMAP\Muscle Activations';
name_result_folder = 'ACC';

% cut-off frequency Hz of the highpass 4th order Butterworth filter applied
cutoff_high = 30;

% cut-off frequency Hz of the lowpass 4th order Butterworth filter applied
cutoff_low = 4;

% acquisition frequency for the data (always 1500 Hz for Clark's)
fs = 1500;

% here we select all of the 120 files to be processed
list_of_file_inds = [1:120];

%% PERFORM INITIAL CHECKS
% check if the required folder already exists, and create it if needed
fullfileName = fullfile(path_to_results, name_result_folder);
[status, msg, ~] = mkdir(fullfileName);

if size(msg,1)>0
    disp(msg)
    error("It appears that the folder that you want to create exists already")
end

%% initialize Matlab path
% set the path current folder to be the one where this script is contained
mfile_name = mfilename('fullpath');
[pathstr,name,ext] = fileparts(mfile_name);
cd(pathstr);

%% GET MVC VALUES FROM SUBJECT FOLDER
disp("Computing MVC values")
array_MVC = get_MVC_values_Clark(path_to_subject_folder);

%% LOOP THROUGH ALL OF THE CSV FILES AND PROCESS THEM
disp('Processing CSV files')

for file_ind = list_of_file_inds
    % define the name of the muscle activation file to be printed
    name_file_muscle_activ = string(file_ind);

    % open the correct CSV file in the workspace (checking if it is there)
    csv_file = append(string(file_ind), ".CSV");
    search_string = fullfile(path_to_subject_folder, csv_file);
    file_info = dir(search_string);

    csv_content = table2array(readtable(fullfile(file_info.folder, file_info.name), 'VariableNamingRule', 'preserve'));
    emg_data = csv_content(2:end, 3:16);

    % NOTE: from now on, we assume that the muscle names and order is the 
    % same as indicated at the beginning of this script
    
    % find and remove the DC bias from the signal
    bias_DC = mean(emg_data,1);
    emg_data = emg_data - bias_DC;

    % highpass Butterworth filter of 4th order (cutoff_high) is applied 
    nyquist_rate = fs/2;
    [b_h,a_h]=butter(4, cutoff_high/nyquist_rate, 'high');
    emg_data = filter(b_h, a_h, emg_data);

    % rectification of the emg data
    emg_data = abs(emg_data);

    % lowpass Butterworth filter of 4th order (cutoff_low) is applied
    [b_l,a_l]=butter(4, cutoff_low/nyquist_rate, 'low');
    emg_data = filter(b_l, a_l, emg_data);

    % normalize the data wrt the MVC values (get the muscle activations)
    % looping over all the muscles
    muscle_activations = zeros(size(emg_data));

    for i=1:size(emg_data,2)
        muscle_activations(:,i) = emg_data(:,i)/array_MVC(i);
    end

    % save the muscle activations to a .mat file (in the right folder)
    current_dir = pwd;
    cd(fullfileName);
    save(name_file_muscle_activ, 'muscle_activations')
    cd(pathstr);
end

% print also an additional file to explain the order of the muscles
cd(fullfileName);
fileID = fopen('muscle_order.txt','w');
fprintf(fileID,'The muscle activations are reported in the following order (column-wise) \n');
fprintf(fileID,'1. Anterior Deltoid \n');
fprintf(fileID,'2. Middle Deltoid \n');
fprintf(fileID,'3. Posterior Deltoid \n');
fprintf(fileID,'4. Biceps \n');
fprintf(fileID,'5. Triceps \n');
fprintf(fileID,'6. Infra \n');
fprintf(fileID,'7. Supra \n');
fprintf(fileID,'8. PecMajorC \n');
fprintf(fileID,'9. PecMajorS \n');
fprintf(fileID,'10. Lats \n');
fprintf(fileID,'11. Sert \n');
fprintf(fileID,'12. TrapLow \n');
fprintf(fileID,'13. TrapMid \n');
fprintf(fileID,'14. TrapUp \n');
fclose(fileID);