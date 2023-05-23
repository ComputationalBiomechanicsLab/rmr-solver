function [dataAverage,dataStdup,dataStdlow, target_data] = evaluate_mean_std_RMR_multipleFiles(path_results, index_muscle, task, files, save_data_struct) 
% This function is used to align the results of the RMR along
% different trials. It leverages scripts from 2019 to compute the standard
% deviations across repetitions of the same task
% INPUTS
% - path_results : path from which to read the data
% - index_muscle: index of the muscle to be considered
% - task: name of the task whose activations need to be aligned
% - files: files to be considered
% - save_data_struct: flag indicating whether aligned processed data should
%                     be saved. Can be a string to specify additional info
%
% author: Italo Belli (i.belli@tudelft.nl) 2023

% get the number of files we are dealing with
if iscell(files)
    num_files = size(files, 2);
else
    num_files = 1;
    files = {files};
end

% pre-allocate cell array to load the results from the files
experiment = cell(num_files,1);

% load results from the files
for index_file=1:num_files
    experiment{index_file} = load(fullfile(path_results, files{index_file}));
end

% load the matrix containing the relevant time instants for each of the
% experiments
load('Results\timings_matrix_experiments.mat')

% preallocate matrix to contain the important indexes defining the phases
% of the motion
index_motion_phase = zeros(num_files, 3);

% get the various important indexes (from the precomputed matrix), and
% store them in the matrix
commom_initial_file_name = 'muscle_activations_Seth2019_experiment_';
for index_file = 1:num_files
    code_name_experiment = extractBefore(extractAfter(files{index_file}, commom_initial_file_name), '_');
    switch code_name_experiment
        case 'ABD01'
            index_motion_phase(index_file,1) = timings_matrix(1,1); index_motion_phase(index_file,2) =timings_matrix(1,2); index_motion_phase(index_file,3) =timings_matrix(1,3);
        case 'ABD02'
            index_motion_phase(index_file,1) = timings_matrix(2,1); index_motion_phase(index_file,2) =timings_matrix(2,2); index_motion_phase(index_file,3) =timings_matrix(2,3);
        case 'ABD03'
            index_motion_phase(index_file,1) = timings_matrix(3,1); index_motion_phase(index_file,2) =timings_matrix(3,2); index_motion_phase(index_file,3) =timings_matrix(3,3);
        case 'ABD21'
            index_motion_phase(index_file,1) = timings_matrix(4,1); index_motion_phase(index_file,2) =timings_matrix(4,2); index_motion_phase(index_file,3) =timings_matrix(4,3);
        case 'ABD22'
            index_motion_phase(index_file,1) = timings_matrix(5,1); index_motion_phase(index_file,2) =timings_matrix(5,2); index_motion_phase(index_file,3) =timings_matrix(5,3);
        case 'ABD23'
            index_motion_phase(index_file,1) = timings_matrix(6,1); index_motion_phase(index_file,2) =timings_matrix(6,2); index_motion_phase(index_file,3) =timings_matrix(6,3);
        case 'FLX01'
            index_motion_phase(index_file,1) = timings_matrix(7,1); index_motion_phase(index_file,2) =timings_matrix(7,2); index_motion_phase(index_file,3) =timings_matrix(7,3);
        case 'FLX02'
            index_motion_phase(index_file,1) = timings_matrix(8,1); index_motion_phase(index_file,2) =timings_matrix(8,2); index_motion_phase(index_file,3) =timings_matrix(8,3);
        case 'FLX03'
            index_motion_phase(index_file,1) = timings_matrix(9,1); index_motion_phase(index_file,2) =timings_matrix(9,2); index_motion_phase(index_file,3) =timings_matrix(9,3);
        case 'FLX21'
            index_motion_phase(index_file,1) = timings_matrix(10,1); index_motion_phase(index_file,2) =timings_matrix(10,2); index_motion_phase(index_file,3) =timings_matrix(10,3);
        case 'FLX22'
            index_motion_phase(index_file,1) = timings_matrix(11,1); index_motion_phase(index_file,2) =timings_matrix(11,2); index_motion_phase(index_file,3) =timings_matrix(11,3);
        case 'FLX23'
            index_motion_phase(index_file,1) = timings_matrix(12,1); index_motion_phase(index_file,2) =timings_matrix(12,2); index_motion_phase(index_file,3) =timings_matrix(12,3);
        case 'SHRUG01'
            index_motion_phase(index_file,1) = timings_matrix(13,1); index_motion_phase(index_file,2) =timings_matrix(13,2); index_motion_phase(index_file,3) =timings_matrix(13,3);
        case 'SHRUG02'
            index_motion_phase(index_file,1) = timings_matrix(14,1); index_motion_phase(index_file,2) =timings_matrix(14,2); index_motion_phase(index_file,3) =timings_matrix(14,3);
        case 'SHRUG03'
            index_motion_phase(index_file,1) = timings_matrix(15,1); index_motion_phase(index_file,2) =timings_matrix(15,2); index_motion_phase(index_file,3) =timings_matrix(15,3);
        case 'SHRUG21'
            index_motion_phase(index_file,1) = timings_matrix(16,1); index_motion_phase(index_file,2) =timings_matrix(16,2); index_motion_phase(index_file,3) =timings_matrix(16,3);
        case 'SHRUG22'
            index_motion_phase(index_file,1) = timings_matrix(17,1); index_motion_phase(index_file,2) =timings_matrix(17,2); index_motion_phase(index_file,3) =timings_matrix(17,3);
        case 'SHRUG23'
            index_motion_phase(index_file,1) = timings_matrix(18,1); index_motion_phase(index_file,2) =timings_matrix(18,2); index_motion_phase(index_file,3) =timings_matrix(18,3);
        otherwise
            warning("The initial, maximum and final indexes defining the movement were not found. Select the correct files (expected name is: 'muscle_activations_Seth2019_experiment_CODENAME_i.mat'")
            warning("You might want to select other files, or change the hardcoded stuff at the beginning of this swith statement to work with your data")
            error("could not proceed") 
    end
end

% transform time into indexes
frequency_data = zeros(num_files,1);
for index_file=1:num_files
    frequency_data(index_file) = experiment{index_file}.frequency_solution;
end

if min(frequency_data) ~= max(frequency_data)
    warning('frequencies of the data are not the same')
end

% frequency_data = frequency_data(1);

startIndex_exp = zeros(num_files,1);
endIndex_exp = zeros(num_files,1);
maxIndex_exp = zeros(num_files,1);

for index_file = 1:num_files
    startIndex_exp(index_file) = max(round(index_motion_phase(index_file,1) * frequency_data(index_file)),1);
    endIndex_exp(index_file) = round(index_motion_phase(index_file,2) * frequency_data(index_file));
    maxIndex_exp(index_file) = round(index_motion_phase(index_file,3) * frequency_data(index_file));
end

data = cell(num_files,1);
dataPercMot = cell(num_files,1);
sampleSize = zeros(num_files,1);
% align data to up-down movements
for index_file = 1:num_files
    [data{index_file},dataPercMot{index_file}]=dataNorm(experiment{index_file}.xsol,startIndex_exp(index_file),endIndex_exp(index_file),maxIndex_exp(index_file),index_muscle);
    sampleSize(index_file) = length(data{index_file});
end

% get the average and standard deviations across the experiments, for
% the muscle considered
[~, index_shortestExp]=min(sampleSize);

target_data = dataPercMot{index_shortestExp};
indexes_data = cell(num_files,1);

for index_file=1:num_files
    [dataPercMot{index_file}, indexes_data{index_file}] = unique(dataPercMot{index_file});
end

% interpolate the data
data_aln = cell(num_files,1);

for index_file=1:num_files
    data_aln{index_file}=interp1(dataPercMot{index_file},data{index_file}(indexes_data{index_file}),target_data,'linear');
end

if save_data_struct
    name_saved_file = append(task, '_', save_data_struct, '_muscle_', char(string(index_muscle)));
    aligned_data = [];
        for index_file=1:num_files
            aligned_data = [aligned_data; data_aln{index_file}];
        end
    save(name_saved_file, 'aligned_data')
end

dataAverage = mean(aligned_data);stdev = std(aligned_data);
dataStdup = dataAverage+stdev;
dataStdlow = dataAverage-stdev;
