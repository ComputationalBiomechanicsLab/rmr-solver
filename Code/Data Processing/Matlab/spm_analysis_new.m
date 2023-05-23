% This script performs a statistical analysis of the muscle activations
% found with the RMR solver, allowing a comparison between activations for
% the same muscles found under different conditions. 
%
% In its current structure, it is employed to analyze the difference in
% activations induced by the glenohumeral (GH) stability condition in the
% solver, when simulating shoulder movements.
% The results of the simulations with the GH constraint are assumed to be 
% stored into a single folder (named withGH), and the results for the 
% simulations without in another one (named noGH).
% The naming convention for the results of the RMR is assumed to be:
% "muscle_activations_Seth2019_2kgWeight1_noGH_ABD21.mat"
% where the last part of the name specifies the experimental data (ABD21)
% used as input to the RMR analysis.
% Resulting activations for the same muscle are aligned across different
% experiments, and compared using a statistical parametric mapping (SPM)
% paired t-test.
%
% The external package SPM1D is required, check detailed instructions on
% how to get it and configure it at https://spm1d.org/install/InstallationMatlab.html

clear; close all; clc

%% Parameters
% number of muscles whose activation was evaluated
nMusc = 33;

% Flag to decide whether the activations of different bundles composing the 
% rotator cuff muscles should be averaged before the SPM analysis
average_RC_activation = true;

% p-value for significance in the SPM paired t-test
p_value = 0.01;

%% Script
% set the path current folder to be the one where this script is contained
mfile_name          = mfilename('fullpath');
[pathstr,~,~]  = fileparts(mfile_name);
cd(pathstr);

% change it to path to repo
addpath(pathstr)
cd '../../..'
path_to_repo = pwd;
addpath(path_to_repo)
addpath(append(path_to_repo, '/Code/Compute Muscle Activity/Matlab'))

% select folders where the RMR results are stored
% select folder of the results obtained considering GH stability (withGH)
path_withGH = uigetdir(path_to_repo, 'Select foder of results withGH');

% select folder of the results obtained without GH stability (noGH)
path_noGH = uigetdir(fullfile(path_withGH, '..'), 'Select foder of results noGH');

% get aligned activations for the two cases
struct_aligned_activations_withGH = align_activations(path_withGH, nMusc);
struct_aligned_activations_noGH = align_activations(path_noGH, nMusc);

% index and corresponding names of the muscles in the RMR solution
if average_RC_activation
    muscle_RMR = [1, 2, 3, 12, 13, 14, 18, 6, 16, 21, 25, 10, 11, 24, 34, 35, 36];   % last three don't exist in the model, are averages of the rotator-cuff elements
    muscle_Name_list={'Trapezius, scapula middle';
        'Trapezius, scapula superior';
        'Trapezius, scapula inferior';
        'Deltoid anterior';
        'Deltoid posterior';
        'Deltoid middle';
        'Pectoralis major clavical';
        'Serratus anterior (Medialis)';
        'Latissimus Dorsi (Medialis)';
        'Teres major';
        'Teres minor';
        'Levator Scapulae';
        'Coracobrachialis';
        'Pectoralis minor';
        'Infraspinatus';
        'Subscapularis';
        'Supraspinatus'};
else
    muscle_RMR = [1, 2, 3, 12, 13, 14, 18, 6, 23, 16, 21, 22, 25, 26, 27, 28, 29, 30, 10, 11, 24];
    muscle_Name_list={'Trapezius, scapula middle';
        'Trapezius, scapula superior';
        'Trapezius, scapula inferior';
        'Deltoid anterior';
        'Deltoid posterior';
        'Deltoid middle';
        'Pectoralis major clavical';
        'Serratus anterior (Medialis)';
        'Infraspinatus (Superior)';
        'Latissimus Dorsi (Medialis)';
        'Teres major';
        'Infraspinatus (Inferior)';
        'Teres minor';
        'Subscapularis (Superior)';
        'Subscapularis (Medialis)';
        'Subscapularis (Inferior)';
        'Supraspinatus (Posterior)';
        'Supraspinatus (Anterior)';
        'Levator Scapulae';
        'Coracobrachialis';
        'Pectoralis minor';};
end

% Number of muscles of interest
nMusc_to_analyze = max(size(muscle_RMR));

if average_RC_activation
    % preallocate new structs, to host also the averaged activations
    aux_struct_aligned_activations_withGH = cell(nMusc+3,1);
    aux_struct_aligned_activations_noGH = cell(nMusc+3,1);

    for index_muscle=1:nMusc
        aux_struct_aligned_activations_withGH{index_muscle} = struct_aligned_activations_withGH{index_muscle};
        aux_struct_aligned_activations_noGH{index_muscle} = struct_aligned_activations_noGH{index_muscle};
    end
    
    % re-assign the structs to the new augmented variables, that will now
    % be filled with averaged activations for the rotator cuff in last 3
    % fields
    struct_aligned_activations_noGH = aux_struct_aligned_activations_noGH;
    struct_aligned_activations_withGH = aux_struct_aligned_activations_withGH;

    clear aux_struct_aligned_activations_withGH aux_struct_aligned_activations_noGH;

    % evaluate the averages across the rotator-cuffs
    % 1. Infraspinatus (divided in 2 bundles) - index:22, 23
    struct_aligned_activations_withGH{nMusc+1} = (struct_aligned_activations_withGH{22}+struct_aligned_activations_withGH{23})/2;
    struct_aligned_activations_noGH{nMusc+1} = (struct_aligned_activations_noGH{22}+struct_aligned_activations_noGH{23})/2;

    % 2. Subscapularis (divided in 3 bundles) - index: 26, 27, 28
    struct_aligned_activations_withGH{nMusc+2} = (struct_aligned_activations_withGH{26}+struct_aligned_activations_withGH{27}+struct_aligned_activations_withGH{28})/3;
    struct_aligned_activations_noGH{nMusc+2} = (struct_aligned_activations_noGH{26}+struct_aligned_activations_noGH{27}+struct_aligned_activations_noGH{28})/3;

    % 3. Supraspinatus (divided in 2 bundles) - index: 29, 30
    struct_aligned_activations_withGH{nMusc+3} = (struct_aligned_activations_withGH{29}+struct_aligned_activations_withGH{30})/2;
    struct_aligned_activations_noGH{nMusc+3} = (struct_aligned_activations_noGH{29}+struct_aligned_activations_noGH{30})/2;

    % augment nMusc to take into account the newly introduced variables
    nMusc = nMusc+3;
end

% performing SPM paired t-test
ratio_significance = zeros(nMusc_to_analyze,1);
significantly_different_muscles_SPM = cell(nMusc_to_analyze,1);

exp_size =  size(struct_aligned_activations_withGH{1}, 2);
mae_raw = zeros(nMusc_to_analyze, size(struct_aligned_activations_withGH{1}, 1));
effect_size_L1 = zeros(nMusc_to_analyze, 1);	
peak_difference = zeros(nMusc_to_analyze, 1);

for index_name_muscle = 1:nMusc_to_analyze
    if index_name_muscle == 2
        doSmthg =0;
    end
    % get the muscle considered
    index_muscle = muscle_RMR(index_name_muscle);
    muscle_name = muscle_Name_list{index_name_muscle};

    % define more convenient names for analysis
    muscle_activations_withGH = struct_aligned_activations_withGH{index_muscle};
    muscle_activations_noGH = struct_aligned_activations_noGH{index_muscle};
    spm = spm1d.stats.ttest_paired(muscle_activations_withGH, muscle_activations_noGH);   % performing a paired t-test
    spmi = spm.inference(p_value);

    % evaluate the percentage of the motion in which the SPM paired t-test
    % detects significant difference between the two sets of results
    significance_time_instants = find(abs(spmi.z)>spmi.zstar);
    ratio_significance(index_name_muscle) = size(significance_time_instants,2)/size(muscle_activations_withGH,2);

    % plot mean and standard deviation
    fig = figure;
    subplot(121)
    spm1d.plot.plot_meanSD(muscle_activations_withGH, 'color', [0.3 0.3 1]);
    hold on
    spm1d.plot.plot_meanSD(muscle_activations_noGH, 'color','g');
    ylim([0, 0.5])
    xticks([0, round(size(muscle_activations_withGH,2)/4), round(size(muscle_activations_withGH,2)/2), round(3*size(muscle_activations_withGH,2)/4), size(muscle_activations_withGH,2)])
    xticklabels({'0', '25', '50', '75', '100'})
    xlim([0, size(muscle_activations_withGH,2)])
    xlabel('% of motion')
    title(append('Mean and SD  (', muscle_name, ')'))

    % plot SPM results
    subplot(122)
    spmi.plot();
    spmi.plot_threshold_label();
    spmi.plot_p_values();
    xticks([0, round(size(muscle_activations_withGH,2)/4), round(size(muscle_activations_withGH,2)/2), round(3*size(muscle_activations_withGH,2)/4), size(muscle_activations_withGH,2)])
    xticklabels({'0', '25', '50', '75', '100'})
    xlim([0, size(muscle_activations_withGH,2)])
    xlabel('% of motion')
    title('Hypothesis test')

    % plot effect size as the difference between the means (withGH-noGH)
    mean_withGH = mean(muscle_activations_withGH);
    mean_noGH = mean(muscle_activations_noGH);
    effect_size = mean_withGH-mean_noGH;
    effect_size_abs = abs(effect_size);
    fig.WindowState = 'maximized';

    fig2 = figure;
    subplot(321)
    plot(effect_size_abs, 'LineWidth', 2, 'Color', 'black')
    xticks([0, round(size(muscle_activations_withGH,2)/4), round(size(muscle_activations_withGH,2)/2), round(3*size(muscle_activations_withGH,2)/4), size(muscle_activations_withGH,2)])
    xticklabels({'0', '25', '50', '75', '100'})
    xlim([0, size(muscle_activations_withGH,2)])
    xlabel('% of motion')
    title(append('Effect size  ', muscle_name))
    ylim([0, 0.3])
    fig2.WindowState = 'maximized';


    effect_size_L1(index_name_muscle) = mean(abs(effect_size));
    peak_difference(index_name_muscle) = max(abs(effect_size));
end

% display information about the significant differences detected by SPM in
% the activations
name_task = split(path_noGH, '\');
fprintf(append("Analising: ", name_task{end-1}, '\n \n'))
fprintf(append('Muscles and related effect size (L1 norm): \n \n'))

for index_name_muscle = 1:nMusc_to_analyze
        fprintf("%s  (%.2f) \n", string(muscle_Name_list{index_name_muscle}), effect_size_L1(index_name_muscle))
end

fprintf(append('\n \n Muscles and peak difference: \n \n'))

for index_name_muscle = 1:nMusc_to_analyze
        fprintf("%s  (%.2f) \n", string(muscle_Name_list{index_name_muscle}), peak_difference(index_name_muscle))
end