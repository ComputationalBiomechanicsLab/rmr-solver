function augment_motion_data(number_new_files)
% this function is used to generate more data starting from the motion data
% that we already have

% set the path current folder to be the one where this script is contained
mfile_name          = mfilename('fullpath');
[pathstr,~,~]  = fileparts(mfile_name);
cd(pathstr);

% getting path to other folders in this repo
addpath(pathstr)
cd ..\..\..\
path_to_repo = pwd;
addpath(path_to_repo)
addpath(fullfile(path_to_repo, 'Code\Data Processing\Matlab\'))

% select the files over which to operate the "data augmentation"
[files,path_mot] = uigetfile('*.mot', 'Select the .mot files to analyse', path_to_repo, 'MultiSelect','on');

if iscell(files)
    num_files = size(files, 2);
else
    num_files = 1;
    files = {files};
end

% divide the number of the required files by the available ones, and vith
% integrer division and modulus
num_file_per_available_data = fix(number_new_files/num_files);
num_missing_files = mod(number_new_files, num_files);

% for reproducibility
rng('default')

num_files_to_generate = num_file_per_available_data + num_missing_files;
for input_file=1:num_files
    % Load the trc file to be considered
    [~, experiment_name] = fileparts(files(input_file));
    q = read_motionFile(fullfile(path_mot, files{input_file}));
    labels = q.labels;
    labels_coord = labels(2:end);
    data = q.data;
    time = data(:,1);
    coordinate_data = data(:, 2:end);
    num_timeInstants = q.nr;
    num_coords = q.nc-1;
    mean_values = mean(coordinate_data,1);
    max_values = max(coordinate_data,[], 1);
    for j=1:num_files_to_generate
        new_coordinate_values = coordinate_data;
        for coord_index=1:num_coords
            change_coordValue = abs(abs(max_values(coord_index))-abs(mean_values(coord_index)));
            if change_coordValue>5
                noise_sd = change_coordValue/50;
                noise = normrnd(0,noise_sd, num_timeInstants, 1);
                % low pass filter : Butterworth filter to filter the noise
                Fs = 1/(time(2)-time(1));
                NN = 2; % filter order
                fc = 2; % cut off frequency (Hz)
                Wn = fc/(Fs/2); % normalized cut off frequency 
                [B,A] = butter(NN,Wn);
                noise_filtered = filter(B,A,noise); % noise (after LP filtering)
                new_coordinate_values(:,coord_index) = new_coordinate_values(:,coord_index) + noise_filtered; 
            end
        end
        new_file_name = [experiment_name, '_', num2str(j)];
        writeMotStoData([time,new_coordinate_values], 1, labels_coord, new_file_name);      
    end
    num_files_to_generate = num_file_per_available_data;
end
