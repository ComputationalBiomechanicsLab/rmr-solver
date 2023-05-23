function augment_marker_data(number_new_files)
% this function is used to generate more data starting from the marker data
% that we already have

noise_sd = 0.005;    % standard deviation of the white noise to add to the marker trajectories, in meters
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
[files,path_trc] = uigetfile('*.trc', 'Select the .trc files to analyse', path_to_repo, 'MultiSelect','on');

if iscell(files)
    num_files = size(files, 2);
else
    num_files = 1;
    files = {files};
end

% divide the number of the required files by the available ones, and with
% integrer division and modulus
num_file_per_available_data = fix(number_new_files/num_files);
num_missing_files = mod(number_new_files, num_files);

% for reproducibility
rng('default')

num_files_to_generate = num_file_per_available_data + num_missing_files;
for input_file=1:num_files
    % Load the trc file to be considered
    [~, experiment_name] = fileparts(files(input_file));
    [markersExp, timesExp, labels, unitsExp] = readTRC(fullfile(path_trc, files{input_file}));
    % convert everything to SI units
    if strcmp(unitsExp, 'mm')
        markersExp = markersExp/1000;
        unitsExp = 'm';
    end

    % get information about the data
    rate = 1/(timesExp(2)-timesExp(1));
    [num_timeInstants, num_traj] = size(markersExp);
    frames = [1:num_timeInstants]';
    mean_values = mean(markersExp,1);
    max_values = max(markersExp,[], 1);

%     noise_sd = abs(max_values-mean_values)/15;

    % generate the augmented data, by adding noise to the marker
    % trajectories (every 1D trajectory is treated independently)
    for j=1:num_files_to_generate
        new_marker_traj = markersExp;
        for traj_index=1:num_traj
            if (abs(abs(max_values(traj_index))-abs(mean_values(traj_index)))/abs(mean_values(traj_index)))>0.1
                noise = normrnd(0,noise_sd, num_timeInstants, 1);
                % low pass filter : Butterworth filter to filter the noise
                NN = 2; % filter order
                fc = 1; % cut off frequency (Hz)
                Wn = fc/(rate/2); % normalized cut off frequency 
                [B,A] = butter(NN,Wn);
                noise_filtered = filter(B,A,noise); % noise (after LP filtering)
                new_marker_traj(:,traj_index) = new_marker_traj(:,traj_index) + noise_filtered; 
            end
        end
        new_file_name = [experiment_name, num2str(j), '.trc'];
        writeMarkersToTRC(new_file_name, new_marker_traj, labels, rate, frames, timesExp, unitsExp);      
    end
    num_files_to_generate = num_file_per_available_data;
end
