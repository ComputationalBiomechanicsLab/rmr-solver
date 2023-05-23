function struct_aligned_activations = align_activations(path_to_folder, nMusc)
% This function allows to align the muscle activations for different repetitions
% of the same task, stored in a given folder. The activations, computed
% with the RMR solver, are aligned using suitable time indexes, precomputed
% on the basis of the experimental study (saved earlier and available in
% timings_matrix_experiments.mat)
% The output is a struct with nMusc fields (where nMusc is the number of muscles in
% the model), each one containing the activations of a specific muscle in
% the model during the various repetitions of the task.

% get all the relevant files in the folder
file_list = dir(fullfile(path_to_folder, '*.mat'));
file_list = {file_list.name}';

num_files = size(file_list,1);

% get info from file name
folder_names_cell = regexp(path_to_folder, '\', 'split');
is_GH_enforced = folder_names_cell{end};

% load the matrix containing the relevant time instants for each of the
% experiments
load('Results\timings_matrix_experiments.mat');

% pre-allocate cell array to load the results from the files
experiment = cell(num_files,1);

% load results from the files
for index_file=1:num_files
    experiment{index_file} = load(fullfile(path_to_folder, file_list{index_file}));
end

% preallocate matrix to contain the important indexes defining the phases
% of the motion
index_motion_phase = zeros(num_files, 3);

% get the various important indexes (from the precomputed matrix), and
% store them in index_motion_phase
for index_file = 1:num_files
    fields_name_experiment = regexp(file_list{index_file}, append('_', is_GH_enforced, '_'), 'split');
    code_name_experiment = extractBefore(fields_name_experiment{end}, '.mat');
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
            warning("You might want to select other files, or change the hardcoded stuff at the beginning of this switch statement to work with your data")
            error("could not proceed") 
    end
end

% transform time into indexes
frequency_data = zeros(num_files,1);
for index_file=1:num_files
    frequency_data(index_file) = experiment{index_file}.frequency_solution;
end

if min(frequency_data) ~= max(frequency_data)
    % issue warning if the frequencies of the solutions considered are not
    % the same
    warning('frequencies of the data are not the same')
end

startIndex_exp = zeros(num_files,1);
endIndex_exp = zeros(num_files,1);
maxIndex_exp = zeros(num_files,1);

for index_file = 1:num_files
    startIndex_exp(index_file) = max(round(index_motion_phase(index_file,1) * frequency_data(index_file)),1);
    endIndex_exp(index_file) = round(index_motion_phase(index_file,2) * frequency_data(index_file));
    maxIndex_exp(index_file) = round(index_motion_phase(index_file,3) * frequency_data(index_file));
end

% initialize struct to store aligned activations muscle by muscle
struct_aligned_activations = cell(nMusc,1);

% consider each muscle separately and save its aligned activations in a
% struct field
for index_muscle=1:nMusc

    % initialize variables
    data = cell(num_files,1);
    dataPercMot = cell(num_files,1);
    sampleSize = zeros(num_files,1);
    
    % align data to up-down movements
    for index_file = 1:num_files
        [data{index_file},dataPercMot{index_file}]=dataNorm(experiment{index_file}.xsol,startIndex_exp(index_file),endIndex_exp(index_file),maxIndex_exp(index_file),index_muscle);
        sampleSize(index_file) = length(data{index_file});
    end
    
    % select the shortest experiment
    [~, index_shortestExp]=min(sampleSize);
    
    target_data = dataPercMot{index_shortestExp};
    indexes_data = cell(num_files,1);
    
    for index_file=1:num_files
        [dataPercMot{index_file}, indexes_data{index_file}] = unique(dataPercMot{index_file});
    end
    
    % interpolate the muscle activations for each experiment
    data_aln = zeros(num_files, min(sampleSize));
    
    for index_file=1:num_files
        data_aln(index_file, :)=interp1(dataPercMot{index_file},data{index_file}(indexes_data{index_file}),target_data,'linear');
    end

    struct_aligned_activations{index_muscle} = data_aln;
end
