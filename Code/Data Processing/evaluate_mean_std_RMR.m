function [dataAverage,dataStdup,dataStdlow, target_data] = evaluate_mean_std_RMR(path_results, index_muscle, task, save_data_struct) 
% This function is used to align the results of the RMR solver along
% different trials. It leverages scripts from 2019 to compute the standard
% deviations across repetitions of the same task
% INPUTS
% - path_results : path from which to read the data
% - index_muscle: index of the muscle to be considered
% - task: task considered (ABD_0kg, ABD_2kg, FLX_0kg, etc ..)
% - save_data_struct: flag indicating whether aligned processed data should
%                     be saved. Can be a string to specify additional info
%
% author: Italo Belli (i.belli@tudelft.nl) 2023


% load RMR results from the 3 experiments
experiment1 = load(fullfile(path_results, (strcat('muscle_activations_Seth2019_experiment_', task(1:end-4), task(end-2), '1'))));
experiment2 = load(fullfile(path_results, (strcat('muscle_activations_Seth2019_experiment_', task(1:end-4), task(end-2), '2'))));
experiment3 = load(fullfile(path_results, (strcat('muscle_activations_Seth2019_experiment_', task(1:end-4), task(end-2), '3'))));

% load the matrix containing the relevant time instants for each of the
% experiments
load('Results\timings_matrix_experiments.mat')

% get the various important indexes (from the precomputed 
if strcmpi(task, 'ABD_0kg')
    start_time_1 = timings_matrix(1,1); end_time_1 =timings_matrix(1,2); max_time_1 =timings_matrix(1,3);
    start_time_2 = timings_matrix(2,1); end_time_2 =timings_matrix(2,2); max_time_2 =timings_matrix(2,3);
    start_time_3 = timings_matrix(3,1); end_time_3 =timings_matrix(3,2); max_time_3 =timings_matrix(3,3);

elseif strcmpi(task, 'ABD_2kg')
    start_time_1 = timings_matrix(4,1); end_time_1 =timings_matrix(4,2); max_time_1 =timings_matrix(4,3);
    start_time_2 = timings_matrix(5,1); end_time_2 =timings_matrix(5,2); max_time_2 =timings_matrix(5,3);
    start_time_3 = timings_matrix(6,1); end_time_3 =timings_matrix(6,2); max_time_3 =timings_matrix(6,3);

elseif strcmpi(task, 'FLX_0kg')
    start_time_1 = timings_matrix(7,1); end_time_1 =timings_matrix(7,2); max_time_1 =timings_matrix(7,3);
    start_time_2 = timings_matrix(8,1); end_time_2 =timings_matrix(8,2); max_time_2 =timings_matrix(8,3);
    start_time_3 = timings_matrix(9,1); end_time_3 =timings_matrix(9,2); max_time_3 =timings_matrix(9,3);

elseif strcmpi(task, 'FLX_2kg')
    start_time_1 = timings_matrix(10,1); end_time_1 =timings_matrix(10,2); max_time_1 =timings_matrix(10,3);
    start_time_2 = timings_matrix(11,1); end_time_2 =timings_matrix(11,2); max_time_2 =timings_matrix(11,3);
    start_time_3 = timings_matrix(12,1); end_time_3 =timings_matrix(12,2); max_time_3 =timings_matrix(12,3);

elseif strcmpi(task, 'SHRUG_0kg')
    start_time_1 = timings_matrix(13,1); end_time_1 =timings_matrix(13,2); max_time_1 =timings_matrix(13,3);
    start_time_2 = timings_matrix(14,1); end_time_2 =timings_matrix(14,2); max_time_2 =timings_matrix(14,3);
    start_time_3 = timings_matrix(15,1); end_time_3 =timings_matrix(15,2); max_time_3 =timings_matrix(15,3);

elseif strcmpi(task, 'SHRUG_2kg')
    start_time_1 = timings_matrix(16,1); end_time_1 =timings_matrix(16,2); max_time_1 =timings_matrix(16,3);
    start_time_2 = timings_matrix(17,1); end_time_2 =timings_matrix(17,2); max_time_2 =timings_matrix(17,3);
    start_time_3 = timings_matrix(18,1); end_time_3 =timings_matrix(18,2); max_time_3 =timings_matrix(18,3);

end

% transform time into indexes
frequency_data_1 = experiment1.frequency_solution;
frequency_data_2 = experiment2.frequency_solution;
frequency_data_3 = experiment3.frequency_solution;

if (frequency_data_1 ~= frequency_data_2) OR (frequency_data_1 ~= frequency_data_3)
    error('frequencies of the data are not the same')
end
frequency_data = frequency_data_1;

startIndex_exp1 = max(round(start_time_1 * frequency_data),1);
endIndex_exp1 = round(end_time_1 * frequency_data);
maxIndex_exp1 = round(max_time_1 * frequency_data);

startIndex_exp2 = max(round(start_time_2 * frequency_data),1);
endIndex_exp2 = round(end_time_2 * frequency_data);
maxIndex_exp2 = round(max_time_2 * frequency_data);

startIndex_exp3 = max(round(start_time_3 * frequency_data),1);
endIndex_exp3 = round(end_time_3 * frequency_data);
maxIndex_exp3 = round(max_time_3 * frequency_data);

% align data to up-down movements
[data1,dataPercMot1]=dataNorm(experiment1.xsol,startIndex_exp1,endIndex_exp1,maxIndex_exp1,index_muscle);

[data2,dataPercMot2]=dataNorm(experiment2.xsol,startIndex_exp2,endIndex_exp2,maxIndex_exp2,index_muscle);

[data3,dataPercMot3]=dataNorm(experiment3.xsol,startIndex_exp3,endIndex_exp3,maxIndex_exp3,index_muscle);

% get the average and standard deviations across the 3 experiments, for
% the muscle considered
sampleSize=min([length(data1),length(data2),length(data3)]);

if sampleSize == length(data1)
    target_data=dataPercMot1;
elseif sampleSize == length(data2)
    target_data=dataPercMot2;
else
    target_data=dataPercMot3;
end

[dataPercMot1, indexes_data1] = unique(dataPercMot1);
[dataPercMot2, indexes_data2] = unique(dataPercMot2);
[dataPercMot3, indexes_data3] = unique(dataPercMot3);

data1_aln=interp1(dataPercMot1,data1(indexes_data1),target_data,'linear');
data2_aln=interp1(dataPercMot2,data2(indexes_data2),target_data,'linear');
data3_aln=interp1(dataPercMot3, data3(indexes_data3),target_data,'linear');

if save_data_struct
    name_saved_file = append(task, '_', save_data_struct, '_muscle_', char(string(index_muscle)));
    aligned_data = [data1_aln; data2_aln; data3_aln];
    save(name_saved_file, 'aligned_data')
end

[dataAverage,dataStdup,dataStdlow] = meanstd(data1_aln,data2_aln,data3_aln);
