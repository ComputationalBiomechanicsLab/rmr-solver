% Statistic Parametric Mapping (SPM) analysis to compare the muscle
% activations predicted by RaMReS in the following two cases:
% - glenohumeral constraint is enforced
% - glenohumeral constraint is not enforced
%
% The SPM packaged is provided for free at https://spm1d.org/ and can be
% installed from there. Note that this script does not run "as is", but should
% be copied in the folder where the rest of the SPM utilities are located
% (as explained in the spm1d.org installation steps). 
% It is licensed with GPL-3.0 license
%
% Author: Italo Belli (i.belli@tudelft.nl) 2023

clear; clc;

automatic_analysis = true;
task_for_automatic_analysis = 'SHRUG_2kg';

aggregate_rotator_cuffs = true;

% set the path current folder to be the one where this script is contained
mfile_name          = mfilename('fullpath');
[pathstr,~,~]  = fileparts(mfile_name);
cd(pathstr);

path_GH = 'C:\Users\italobelli\Desktop\GitHub\PTbot_officialCBL\Personal_Results\noise_on_elbow\SHRUG_2kg\aligned_activations\3samplesPerExp\withGH';
path_noGH = 'C:\Users\italobelli\Desktop\GitHub\PTbot_officialCBL\Personal_Results\noise_on_elbow\SHRUG_2kg\aligned_activations\3samplesPerExp\noGH';

% index and corresponding names of the muscles in the solution
if aggregate_rotator_cuffs
    muscle_RaMReS = [1, 2, 3, 12, 13, 14, 18, 6, 16, 21, 25, 10, 11, 24, 34, 35, 36];   % last three don't exist in the model, are averages of the rotator-cuff elements
    muscle_Name={'Trapezius, scapula middle';
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
    muscle_RaMReS = [1, 2, 3, 12, 13, 14, 18, 6, 23, 16, 21, 22, 25, 26, 27, 28, 29, 30, 10, 11, 24];
    muscle_Name={'Trapezius, scapula middle';
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

if ~automatic_analysis
    % select first file
    [file_GH,~] = uigetfile('*.mat', 'Select the file to analyse (with GH)', path_GH, 'MultiSelect','off');
    
    % select second file
    [file_noGH, ~] = uigetfile('*.mat', 'Select the file to analyse (no GH)', path_noGH, 'MultiSelect','off');

    num_files = 1;
    file_noGH = {file_noGH};
    file_GH = {file_GH};
else
    num_files = max(size(muscle_RaMReS));
    file_GH = cell(num_files,1);
    file_noGH = cell(num_files,1);

    if aggregate_rotator_cuffs
        % evaluate the averages across the rotator-cuffs, if not present
        % 1. infraspinatus (divided in 2 bundles)
        if ~isfile(fullfile(path_GH,append(task_for_automatic_analysis, '_withGH_muscle_34.mat')))
            data1 = load(fullfile(path_GH,append(task_for_automatic_analysis, '_withGH_muscle_22.mat')));
            data1 = data1.aligned_data; 
            data2 = load(fullfile(path_GH,append(task_for_automatic_analysis, '_withGH_muscle_23.mat')));
            data2 = data2.aligned_data; 
            aligned_data = (data1+data2)/2;
            name_file_avg_infraspinatus = fullfile(path_GH, append(task_for_automatic_analysis, '_withGH_muscle_34.mat'));
            save(name_file_avg_infraspinatus, "aligned_data");
        end

        if ~isfile(fullfile(path_noGH,append(task_for_automatic_analysis, '_noGH_muscle_34.mat')))
            data1 = load(fullfile(path_noGH,append(task_for_automatic_analysis, '_noGH_muscle_22.mat')));
            data1 = data1.aligned_data; 
            data2 = load(fullfile(path_noGH,append(task_for_automatic_analysis, '_noGH_muscle_23.mat')));
            data2 = data2.aligned_data; 
            aligned_data = (data1+data2)/2;
            name_file_avg_infraspinatus = fullfile(path_noGH, append(task_for_automatic_analysis, '_noGH_muscle_34.mat'));
            save(name_file_avg_infraspinatus, "aligned_data");
        end

        % 2. Subscapularis (divided in 3 bundles)
        if ~isfile(fullfile(path_GH,append(task_for_automatic_analysis, '_withGH_muscle_35.mat')))
            data1 = load(fullfile(path_GH,append(task_for_automatic_analysis, '_withGH_muscle_26.mat')));
            data1 = data1.aligned_data;
            data2 = load(fullfile(path_GH,append(task_for_automatic_analysis, '_withGH_muscle_27.mat')));
            data2 = data2.aligned_data;
            data3 = load(fullfile(path_GH,append(task_for_automatic_analysis, '_withGH_muscle_28.mat')));
            data3 = data3.aligned_data;
            aligned_data = (data1+data2+data3)/3;
            name_file_avg_subscap = fullfile(path_GH, append(task_for_automatic_analysis, '_withGH_muscle_35.mat'));
            save(name_file_avg_subscap, "aligned_data");
        end
        if ~isfile(fullfile(path_noGH,append(task_for_automatic_analysis, '_noGH_muscle_35.mat')))
            data1 = load(fullfile(path_noGH,append(task_for_automatic_analysis, '_noGH_muscle_26.mat')));
            data1 = data1.aligned_data;
            data2 = load(fullfile(path_noGH,append(task_for_automatic_analysis, '_noGH_muscle_27.mat')));
            data2 = data2.aligned_data;
            data3 = load(fullfile(path_noGH,append(task_for_automatic_analysis, '_noGH_muscle_28.mat')));
            data3 = data3.aligned_data;
            aligned_data = (data1+data2+data3)/3;
            name_file_avg_subscap = fullfile(path_noGH, append(task_for_automatic_analysis, '_noGH_muscle_35.mat'));
            save(name_file_avg_subscap, "aligned_data");
        end

        % 2. Supraspinatus (divided in 2 bundles)
        if ~isfile(fullfile(path_GH,append(task_for_automatic_analysis, '_withGH_muscle_36.mat')))
            data1 = load(fullfile(path_GH,append(task_for_automatic_analysis, '_withGH_muscle_29.mat')));
            data1 = data1.aligned_data;
            data2 = load(fullfile(path_GH,append(task_for_automatic_analysis, '_withGH_muscle_30.mat')));
            data2 = data2.aligned_data;
            aligned_data = (data1+data2)/2;
            name_file_avg_supraspin = fullfile(path_GH, append(task_for_automatic_analysis, '_withGH_muscle_36.mat'));
            save(name_file_avg_supraspin, "aligned_data");
        end
        if ~isfile(fullfile(path_noGH,append(task_for_automatic_analysis, '_noGH_muscle_36.mat')))
            data1 = load(fullfile(path_noGH,append(task_for_automatic_analysis, '_noGH_muscle_29.mat')));
            data1 = data1.aligned_data;
            data2 = load(fullfile(path_noGH,append(task_for_automatic_analysis, '_noGH_muscle_30.mat')));
            data2 = data2.aligned_data;
            aligned_data = (data1+data2)/2;
            name_file_avg_supraspin = fullfile(path_noGH, append(task_for_automatic_analysis, '_noGH_muscle_36.mat'));
            save(name_file_avg_supraspin, "aligned_data");
        end
    end

    for index_file = 1:num_files
        file_GH{index_file} = append(task_for_automatic_analysis, '_withGH_muscle_', num2str(muscle_RaMReS(index_file)), '.mat');
        file_noGH{index_file} = append(task_for_automatic_analysis, '_noGH_muscle_', num2str(muscle_RaMReS(index_file)), '.mat');
    end
end

ratio_significance = [];
threshold_significance = 0.1;
significantly_different_muscles_SPM = [];

for index_file = 1:num_files

    withGH = load(fullfile(path_GH, file_GH{index_file}));
    withGH = withGH.aligned_data;
    
    noGH = load(fullfile(path_noGH, file_noGH{index_file}));
    noGH = noGH.aligned_data;

    % get the muscle considered
    muscle_1 = extractBefore(file_GH{index_file}, '.');
    muscle_1 = split(muscle_1, '_');
    muscle_1 = str2num(muscle_1{end});
    muscle_2 = extractBefore(file_noGH{index_file}, '.');
    muscle_2 = split(muscle_2, '_');
    muscle_2 = str2num(muscle_2{end});
    if muscle_1==muscle_2
        muscle_index = muscle_1;
        muscle_name = muscle_Name{find(muscle_RaMReS==muscle_index)};
    else
        error("different muscles compared")
    end
    
    spm = spm1d.stats.ttest_paired(withGH, noGH);   % performing a paired t-test
    spmi = spm.inference(0.01);

    % evaluate the percentage of the motion in which the SPM paired t-test
    % detects significant difference between the two sets of results
    significance_time_instants = find(abs( spmi.z)>spmi.zstar);
    ratio_significance(index_file) = size(significance_time_instants,2)/size(withGH,2);
    
    if ratio_significance(index_file)>threshold_significance
        significantly_different_muscles_SPM = [significantly_different_muscles_SPM, {muscle_name}];
    end

    %%% plot mean and SD:
    fig = figure;
    subplot(121)
    spm1d.plot.plot_meanSD(withGH, 'color', [0.3 0.3 1]);
    hold on
    spm1d.plot.plot_meanSD(noGH, 'color','g');
    ylim([0, 0.5])
    xticks([0, round(size(withGH,2)/4), round(size(withGH,2)/2), round(3*size(withGH,2)/4), size(withGH,2)])
    xticklabels({'0', '25', '50', '75', '100'})
    xlim([0, size(withGH,2)])
    xlabel('% of motion')
    title(append('Mean and SD  (', muscle_name, ')'))
    %%% plot SPM results:
    subplot(122)
    spmi.plot();
    spmi.plot_threshold_label();
    spmi.plot_p_values();
    xticks([0, round(size(withGH,2)/4), round(size(withGH,2)/2), round(3*size(withGH,2)/4), size(withGH,2)])
    xticklabels({'0', '25', '50', '75', '100'})
    xlim([0, size(withGH,2)])
    xlabel('% of motion')
    title('Hypothesis test')
    fig.WindowState = 'maximized';
%     name_fig = append(muscle_name, '_', task_for_automatic_analysis, '_SPM_pairedTTEst.svg');
%     saveas(fig, name_fig)
end

% display information about the significant differences detected by SPM in
% the activations
fprintf(append('The muscles that the test detected as different for more than ', num2str(threshold_significance*100), '%% of the motion are: \n \n'))
significantly_different_muscles_SPM'