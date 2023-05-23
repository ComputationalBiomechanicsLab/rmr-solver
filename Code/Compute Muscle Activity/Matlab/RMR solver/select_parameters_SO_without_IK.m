% Script used to analyze the sensitivity of the solution of the static
% optimization (SO) to different conditions.

close all; clear all; clc; beep off;

% Import the OpenSim libraries.
import org.opensim.modeling.*;

% Select model
modelFile = 'C:\Users\italobelli\Desktop\GitHub\PTbot_officialCBL\Personal_Results\TSM_Ajay2019_2kgWeight.osim';
model = Model(modelFile);

% Select the experimental data to be considered
subject_considered = 'Ajay2019';

experiment1 = 'C:\Users\italobelli\Desktop\Ajay old studies\FrontiersPubMaterials-latest\ThoracoscapularShoulderPaperMaterials\ThoracoscapularShoulderPaperMaterials\Results\IK\flx21_IK.mot';
experiment2 = 'C:\Users\italobelli\Desktop\Ajay old studies\FrontiersPubMaterials-latest\ThoracoscapularShoulderPaperMaterials\ThoracoscapularShoulderPaperMaterials\Results\IK\abd21_IK.mot';
experiment3 = 'C:\Users\italobelli\Desktop\Ajay old studies\FrontiersPubMaterials-latest\ThoracoscapularShoulderPaperMaterials\ThoracoscapularShoulderPaperMaterials\Results\IK\shrug21_IK.mot';

% Downsampling
time_interval = 10;

% Flags 
dynamic_bounds = true;
enforce_GH_constraint = true;

%% Solving SO
fprintf('Performing SO 1st experiment \n');
[optimizationStatus1, unfeasibility_flags1, tOptim1, result_file_SO1] = CSO_analysis_without_IK(subject_considered, model, experiment1, time_interval, dynamic_bounds, enforce_GH_constraint);
fprintf('\n Solved with %i unfeasible solutions \n \n \n', sum(unfeasibility_flags1));

fprintf('Performing SO 2nd experiment \n');
[optimizationStatus2, unfeasibility_flags2, tOptim2, result_file_SO2] = CSO_analysis_without_IK(subject_considered, model, experiment2, time_interval, dynamic_bounds, enforce_GH_constraint);
fprintf('\n Solved with %i unfeasible solutions \n \n \n', sum(unfeasibility_flags2));

fprintf('Performing SO 3rd experiment \n');
[optimizationStatus3, unfeasibility_flags3, tOptim3, result_file_SO3] = CSO_analysis_without_IK(subject_considered, model, experiment3, time_interval, dynamic_bounds, enforce_GH_constraint);
fprintf('\n Solved with %i unfeasible solutions \n \n \n', sum(unfeasibility_flags3));

%% Visualizing the results and comparing them to the EMG data
compare_results_CMC2019vsSO(result_file_SO1, 'C:\Users\italobelli\Desktop\GitHub\PTbot_officialCBL\Personal_Results\muscle activations 2019 study Ajay\CMC2019 results\flx2x', 'flx21');
visualize_max_activation(result_file_SO1);

compare_results_CMC2019vsSO(result_file_SO2, 'C:\Users\italobelli\Desktop\GitHub\PTbot_officialCBL\Personal_Results\muscle activations 2019 study Ajay\CMC2019 results\abd2x', 'abd21');
visualize_max_activation(result_file_SO2);

compare_results_CMC2019vsSO(result_file_SO3, 'C:\Users\italobelli\Desktop\GitHub\PTbot_officialCBL\Personal_Results\muscle activations 2019 study Ajay\CMC2019 results\shrug2x', 'shrug21');
visualize_max_activation(result_file_SO3);


