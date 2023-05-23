function compare_results_CMC2019vsSO(results_SO, folder_results_CMC2019, experiment)
% This function allows to compare the results of our Static Optimization
% procedure with the ones that were obtained in teh study
% "Muscle contributions to upper-extremity movement and work from a
% musculoskeletal model of the human shoulder", 2019.
% It takes as input the results of CSO_accelerationConstLin_TSM.m
% (results_SO) and the results of CMC2019 (assumed to be saved in the
% correct form).
% "experiment" discriminates which experiment is considered
%
% example usage:
% compare_results_CMC2109vsSO('C:\Users\italobelli\Desktop\GitHub\PTbot_officialCBL\Personal_Results\muscle activations 2019 study Ajay\SO results\shrug2x\muscle_activations_Ajay2019__experiment_SHRUG23.mat', 'C:\Users\italobelli\Desktop\GitHub\PTbot_officialCBL\Personal_Results\muscle activations 2019 study Ajay\CMC2019 results\shrug2x', 'shrug23')

load(results_SO);

if strcmpi(experiment, 'flx01')
    start_time = 0.634;
    end_time = 5.926;
    label_x = '% flexion';
elseif strcmpi(experiment, 'flx02')
    start_time = 0.938;
    end_time = 6.720;
    label_x = '% flexion';
elseif strcmpi(experiment, 'flx03')
    start_time = 1.077;
    end_time = 6.737;
    label_x = '% flexion';
elseif strcmpi(experiment, 'flx21')
    start_time = 0.050;
    end_time = 5.230;
    label_x = '% flexion';
elseif strcmpi(experiment, 'flx22')
    start_time = 1.061;
    end_time = 7.297;
    label_x = '% flexion';
elseif strcmpi(experiment, 'flx23')
    start_time = 1.519;
    end_time = 7.114;
    label_x = '% flexion';
elseif strcmpi(experiment, 'abd01')
    start_time = 2.925;
    end_time = 9.787;
    label_x = '% abduction';
elseif strcmpi(experiment, 'abd02')
    start_time = 1.315;
    end_time = 6.875;
    label_x = '% abduction';
elseif strcmpi(experiment, 'abd03')
    start_time = 0.1;
    end_time = 6.666;
    label_x = '% abduction';
elseif strcmpi(experiment, 'abd21')
    start_time = 0.552;
    end_time = 7.005;
    label_x = '% abduction';
elseif strcmpi(experiment, 'abd22')
    start_time = 0.81;
    end_time = 7.43;
    label_x = '% abduction';
elseif strcmpi(experiment, 'abd23')
    start_time = 0.717;
    end_time = 7.223;
    label_x = '% abduction';
elseif strcmpi(experiment, 'shrug01')
    start_time = 0.894;
    end_time = 3.12;
    label_x = '% shrug';
elseif strcmpi(experiment, 'shrug02')
    start_time = 0.084;
    end_time = 2.857;
    label_x = '% shrug';
elseif strcmpi(experiment, 'shrug03')
    start_time = 0.71;
    end_time = 2.94;
    label_x = '% shrug';
elseif strcmpi(experiment, 'shrug21')
    start_time = 0.83;
    end_time = 2.828;
    label_x = '% shrug';
elseif strcmpi(experiment, 'shrug22')
    start_time = 0.7;
    end_time = 3.191;
    label_x = '% shrug';
elseif strcmpi(experiment, 'shrug23')
    start_time = 1.15;
    end_time = 3.38;
    label_x = '% shrug';
end

start_index_SO = max(round(start_time * frequency_solution),1);
end_index_SO = round(end_time * frequency_solution);
step = 1/(end_index_SO-start_index_SO);

muscles_SO = xsol(start_index_SO:end_index_SO, :);
length_SO = size(muscles_SO,1);

% allocate array to contain MAE, MAPE, correlations
MAE_SO = zeros(1,11);
MAPE_SO = zeros(1,11);
correlation_SO = zeros(2,11);
MAE_CMC = zeros(1,11);
MAPE_CMC = zeros(1,11);
correlation_CMC = zeros(2,11);

%% PROCESS and PLOT the RESULTS -------------------------------------------------------
max_act = 1.1;
lim_axis = 1;       % flag to limit the axis between 0-200 (X) and 0-max_act (Y)

scatter_downsampled_point = false;

%% (1) - Trapezius Scapula Middle
muscle_SO = muscles_SO(:,1);
name_CMC2019_results = append(folder_results_CMC2019, '\2019_trapezius_scapula_middle_EMG.mat');
load(name_CMC2019_results);

% retrieve the length of the experimental EMG data and CMC2019 solution
length_EMG = size(betweenx_emg,2)/2;
length_CMC = size(CMCaverage,2);

% get the average of EMG values
average_EMG = (betweeny_emg(1:length_EMG)+flip(betweeny_emg(length_EMG+1:end)))/2;

downFactor_EMG =  round(size(average_EMG,2)/size(muscle_SO,1));
downFactor_CMC =  round(size(CMCaverage,2)/size(muscle_SO,1));

emg_low_freq = downsample(average_EMG, downFactor_EMG);
cmc_low_freq = downsample(CMCaverage, downFactor_CMC);

if ~((size(emg_low_freq,2) == size (cmc_low_freq,2)) && (size(emg_low_freq,2)==size(muscle_SO,1)))
    warning("Downsampling leads to un-equal arrays.. proceeding anyway, some samples will be discared")
end
data_length = min([size(emg_low_freq,2), size(muscle_SO,1), size(cmc_low_freq,2)]);

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(emg_low_freq(1:data_length), muscle_SO(1:data_length), 'coeff');
correlation_SO(1,1) = max(corrs);
correlation_SO(2,1) = lags(find(ismember(corrs, max(corrs))));
[corrs, lags] = xcorr(emg_low_freq(1:data_length), cmc_low_freq(1:data_length), 'coeff');
correlation_CMC(1,1) = max(corrs);
correlation_CMC(2,1) = lags(find(ismember(corrs, max(corrs))));

aux1 = 0;
aux2 = 0;
aux3 = 0;
aux4 = 0;
for i=1:data_length
    % building the MAE elements
    aux1 = aux1 + abs(emg_low_freq(i)-muscle_SO(i));
    aux2 = aux2 + abs(emg_low_freq(i)-cmc_low_freq(i));
    % building the MAPE elements
    aux3 = aux3 + abs((emg_low_freq(i)-muscle_SO(i))/emg_low_freq(i));
    aux4 = aux4 + abs((emg_low_freq(i)-cmc_low_freq(i))/emg_low_freq(i));
end

MAE_SO(1) = aux1/data_length;
MAE_CMC(1) = aux2/data_length;
MAPE_SO(1) = 100/data_length*aux3;
MAPE_CMC(1) = 100/data_length*aux4;

% PLOT the comparison
time_SO = (0:step:1)*target_CMC(end);

f_muscle_activ = figure;
subplot(3,4,1)
hold on
set(gca,'fontsize',15,'LineWidth',2)
h_mae=fill(betweenx_emg,betweeny_emg,[0.8 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_EMG, EMGstdup,'Color',[0.8 0.8 0.8])
h_mae=fill(betweenx,betweeny,[1 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_CMC, CMCaverage,'r','LineWidth',2)
plot(time_SO, muscle_SO,'Color' , '#7E2F8E', 'LineWidth',2)
if lim_axis
    axis([0 200 0 max_act])
end
% legend('EMG data', '', '', 'CMC 2019', 'SO result')
xlabel(label_x)
ylabel('activation')
title('Trap Scapula Mid')

% scatter the downsampled points
if scatter_downsampled_point
    inds_emg_downsampled = find(ismember(average_EMG, emg_low_freq));
    inds_cmc_downsampled = find(ismember(CMCaverage, cmc_low_freq));
    
    scatter(time_SO(1:data_length), muscle_SO(1:data_length), 'magenta', 'filled')
    scatter(time_SO(1:data_length), emg_low_freq(1:data_length), 'black', 'filled')
    scatter(time_SO(1:data_length), cmc_low_freq(1:data_length), 'red', 'filled')
end
hold off

%% (2) - Trapezius Scapula Superior
muscle_SO = muscles_SO(:,2);
name_CMC2019_results = append(folder_results_CMC2019, '\2019_trapezius_scapula_superior_EMG.mat');
load(name_CMC2019_results);
average_EMG = (betweeny_emg(1:length_EMG)+flip(betweeny_emg(length_EMG+1:end)))/2;

% RESAMPLE the EMG data and CMC solution to same frequency as SO solution
emg_low_freq = downsample(average_EMG, downFactor_EMG);
cmc_low_freq = downsample(CMCaverage, downFactor_CMC);

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(emg_low_freq(1:data_length), muscle_SO(1:data_length), 'coeff');
correlation_SO(1,2) = max(corrs);
correlation_SO(2,2) = lags(find(ismember(corrs, max(corrs))));
[corrs, lags] = xcorr(emg_low_freq(1:data_length), cmc_low_freq(1:data_length), 'coeff');
correlation_CMC(1,2) = max(corrs);
correlation_CMC(2,2) = lags(find(ismember(corrs, max(corrs))));

aux1 = 0;
aux2 = 0;
aux3 = 0;
aux4 = 0;
for i=1:data_length
    % building the MAE elements
    aux1 = aux1 + abs(emg_low_freq(i)-muscle_SO(i));
    aux2 = aux2 + abs(emg_low_freq(i)-cmc_low_freq(i));
    % building the MAPE elements
    aux3 = aux3 + abs((emg_low_freq(i)-muscle_SO(i))/emg_low_freq(i));
    aux4 = aux4 + abs((emg_low_freq(i)-cmc_low_freq(i))/emg_low_freq(i));
end

MAE_SO(2) = aux1/data_length;
MAE_CMC(2) = aux2/data_length;
MAPE_SO(2) = 100/data_length*aux3;
MAPE_CMC(2) = 100/data_length*aux4;

% PLOT the comparison
time_SO = (0:step:1)*target_CMC(end);

subplot(3,4,2)
hold on
set(gca,'fontsize',15,'LineWidth',2)
h_mae=fill(betweenx_emg,betweeny_emg,[0.8 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_EMG, EMGstdup,'Color',[0.8 0.8 0.8])
h_mae=fill(betweenx,betweeny,[1 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_CMC, CMCaverage,'r','LineWidth',2)
plot(time_SO, muscle_SO,'Color' , '#7E2F8E', 'LineWidth',2)
if lim_axis
    axis([0 200 0 max_act])
end
% legend('EMG data', '', '', 'CMC 2019', 'SO result')
xlabel(label_x)
ylabel('activation')
title('Trap Scapula Sup')

% scatter the downsampled points
if scatter_downsampled_point
    scatter(time_SO(1:data_length), muscle_SO(1:data_length), 'magenta', 'filled')
    scatter(time_SO(1:data_length), emg_low_freq(1:data_length), 'black', 'filled')
    scatter(time_SO(1:data_length), cmc_low_freq(1:data_length), 'red', 'filled')
end
hold off

%% (3) - Trapezius Scapula Inferior
muscle_SO = muscles_SO(:,3);
name_CMC2019_results = append(folder_results_CMC2019, '\2019_trapezius_scapula_inferior_EMG.mat');
load(name_CMC2019_results);
average_EMG = (betweeny_emg(1:length_EMG)+flip(betweeny_emg(length_EMG+1:end)))/2;

% RESAMPLE the EMG data and CMC solution to same frequency as SO solution
emg_low_freq = downsample(average_EMG, downFactor_EMG);
cmc_low_freq = downsample(CMCaverage, downFactor_CMC);

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(emg_low_freq(1:data_length), muscle_SO(1:data_length), 'coeff');
correlation_SO(1,3) = max(corrs);
correlation_SO(2,3) = lags(find(ismember(corrs, max(corrs))));
[corrs, lags] = xcorr(emg_low_freq(1:data_length), cmc_low_freq(1:data_length), 'coeff');
correlation_CMC(1,3) = max(corrs);
correlation_CMC(2,3) = lags(find(ismember(corrs, max(corrs))));

aux1 = 0;
aux2 = 0;
aux3 = 0;
aux4 = 0;
for i=1:data_length
    % building the MAE elements
    aux1 = aux1 + abs(emg_low_freq(i)-muscle_SO(i));
    aux2 = aux2 + abs(emg_low_freq(i)-cmc_low_freq(i));
    % building the MAPE elements
    aux3 = aux3 + abs((emg_low_freq(i)-muscle_SO(i))/emg_low_freq(i));
    aux4 = aux4 + abs((emg_low_freq(i)-cmc_low_freq(i))/emg_low_freq(i));
end

MAE_SO(3) = aux1/data_length;
MAE_CMC(3) = aux2/data_length;
MAPE_SO(3) = 100/data_length*aux3;
MAPE_CMC(3) = 100/data_length*aux4;

% PLOT the comparison
time_SO = (0:step:1)*target_CMC(end);

subplot(3,4,3)
hold on
set(gca,'fontsize',15,'LineWidth',2)
h_mae=fill(betweenx_emg,betweeny_emg,[0.8 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_EMG, EMGstdup,'Color',[0.8 0.8 0.8])
h_mae=fill(betweenx,betweeny,[1 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_CMC, CMCaverage,'r','LineWidth',2)
plot(time_SO, muscle_SO,'Color' , '#7E2F8E', 'LineWidth',2)
if lim_axis
    axis([0 200 0 max_act])
end
% legend('EMG data', '', '', 'CMC 2019', 'SO result')
xlabel(label_x)
ylabel('activation')
title('Trap Scapula Inf')

% scatter the downsampled points
if scatter_downsampled_point
    scatter(time_SO(1:data_length), muscle_SO(1:data_length), 'magenta', 'filled')
    scatter(time_SO(1:data_length), emg_low_freq(1:data_length), 'black', 'filled')
    scatter(time_SO(1:data_length), cmc_low_freq(1:data_length), 'red', 'filled')
end
hold off

%% (4) - Serratus Anterior
muscle_SO = muscles_SO(:,6);
name_CMC2019_results = append(folder_results_CMC2019, '\2019_serratus_anterior_EMG.mat');
load(name_CMC2019_results);
average_EMG = (betweeny_emg(1:length_EMG)+flip(betweeny_emg(length_EMG+1:end)))/2;

% RESAMPLE the EMG data and CMC solution to same frequency as SO solution
emg_low_freq = downsample(average_EMG, downFactor_EMG);
cmc_low_freq = downsample(CMCaverage, downFactor_CMC);

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(emg_low_freq(1:data_length), muscle_SO(1:data_length), 'coeff');
correlation_SO(1,4) = max(corrs);
correlation_SO(2,4) = lags(find(ismember(corrs, max(corrs))));
[corrs, lags] = xcorr(emg_low_freq(1:data_length), cmc_low_freq(1:data_length), 'coeff');
correlation_CMC(1,4) = max(corrs);
correlation_CMC(2,4) = lags(find(ismember(corrs, max(corrs))));

aux1 = 0;
aux2 = 0;
aux3 = 0;
aux4 = 0;
for i=1:data_length
    % building the MAE elements
    aux1 = aux1 + abs(emg_low_freq(i)-muscle_SO(i));
    aux2 = aux2 + abs(emg_low_freq(i)-cmc_low_freq(i));
    % building the MAPE elements
    aux3 = aux3 + abs((emg_low_freq(i)-muscle_SO(i))/emg_low_freq(i));
    aux4 = aux4 + abs((emg_low_freq(i)-cmc_low_freq(i))/emg_low_freq(i));
end

MAE_SO(4) = aux1/data_length;
MAE_CMC(4) = aux2/data_length;
MAPE_SO(4) = 100/data_length*aux3;
MAPE_CMC(4) = 100/data_length*aux4;

% PLOT the comparison
time_SO = (0:step:1)*target_CMC(end);

subplot(3,4,4)
hold on
set(gca,'fontsize',15,'LineWidth',2)
h_mae=fill(betweenx_emg,betweeny_emg,[0.8 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_EMG, EMGstdup,'Color',[0.8 0.8 0.8])
h_mae=fill(betweenx,betweeny,[1 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_CMC, CMCaverage,'r','LineWidth',2)
plot(time_SO, muscle_SO,'Color' , '#7E2F8E', 'LineWidth',2)
if lim_axis
    axis([0 200 0 max_act])
end
% legend('EMG data', '', '', 'CMC 2019', 'SO result')
xlabel(label_x)
ylabel('activation')
title('Serratus Ant')

% scatter the downsampled points
if scatter_downsampled_point
    scatter(time_SO(1:data_length), muscle_SO(1:data_length), 'magenta', 'filled')
    scatter(time_SO(1:data_length), emg_low_freq(1:data_length), 'black', 'filled')
    scatter(time_SO(1:data_length), cmc_low_freq(1:data_length), 'red', 'filled')
end
hold off

%% (5) - Deltoid Anterior
muscle_SO = muscles_SO(:,12);
name_CMC2019_results = append(folder_results_CMC2019, '\2019_deltoid_anterior_EMG.mat');
load(name_CMC2019_results);
average_EMG = (betweeny_emg(1:length_EMG)+flip(betweeny_emg(length_EMG+1:end)))/2;

% RESAMPLE the EMG data and CMC solution to same frequency as SO solution
emg_low_freq = downsample(average_EMG, downFactor_EMG);
cmc_low_freq = downsample(CMCaverage, downFactor_CMC);

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(emg_low_freq(1:data_length), muscle_SO(1:data_length), 'coeff');
correlation_SO(1,5) = max(corrs);
correlation_SO(2,5) = lags(find(ismember(corrs, max(corrs))));
[corrs, lags] = xcorr(emg_low_freq(1:data_length), cmc_low_freq(1:data_length), 'coeff');
correlation_CMC(1,5) = max(corrs);
correlation_CMC(2,5) = lags(find(ismember(corrs, max(corrs))));

aux1 = 0;
aux2 = 0;
aux3 = 0;
aux4 = 0;
for i=1:data_length
    % building the MAE elements
    aux1 = aux1 + abs(emg_low_freq(i)-muscle_SO(i));
    aux2 = aux2 + abs(emg_low_freq(i)-cmc_low_freq(i));
    % building the MAPE elements
    aux3 = aux3 + abs((emg_low_freq(i)-muscle_SO(i))/emg_low_freq(i));
    aux4 = aux4 + abs((emg_low_freq(i)-cmc_low_freq(i))/emg_low_freq(i));
end

MAE_SO(5) = aux1/data_length;
MAE_CMC(5) = aux2/data_length;
MAPE_SO(5) = 100/data_length*aux3;
MAPE_CMC(5) = 100/data_length*aux4;

% PLOT the comparison
time_SO = (0:step:1)*target_CMC(end);

subplot(3,4,5)
hold on
set(gca,'fontsize',15,'LineWidth',2)
h_mae=fill(betweenx_emg,betweeny_emg,[0.8 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_EMG, EMGstdup,'Color',[0.8 0.8 0.8])
h_mae=fill(betweenx,betweeny,[1 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_CMC, CMCaverage,'r','LineWidth',2)
plot(time_SO, muscle_SO,'Color' , '#7E2F8E', 'LineWidth',2)
if lim_axis
    axis([0 200 0 max_act])
end
% legend('EMG data', '', '', 'CMC 2019', 'SO result')
xlabel(label_x)
ylabel('activation')
title('Deltoid Ant')

% scatter the downsampled points
if scatter_downsampled_point
    scatter(time_SO(1:data_length), muscle_SO(1:data_length), 'magenta', 'filled')
    scatter(time_SO(1:data_length), emg_low_freq(1:data_length), 'black', 'filled')
    scatter(time_SO(1:data_length), cmc_low_freq(1:data_length), 'red', 'filled')
end
hold off

%% (6) - Deltoid Posterior
% muscle_SO = muscles_SO(:,13); % correct muscle 
muscle_SO =  muscles_SO(:,14); % incorrect, to compare with CMC 2019
warning('COMPARING DELTOID POSTERIOR TO DeltoideusScapula_M, as in 2019 CMC')
name_CMC2019_results = append(folder_results_CMC2019, '\2019_deltoid_posterior_EMG.mat');
load(name_CMC2019_results);
average_EMG = (betweeny_emg(1:length_EMG)+flip(betweeny_emg(length_EMG+1:end)))/2;

% RESAMPLE the EMG data and CMC solution to same frequency as SO solution
emg_low_freq = downsample(average_EMG, downFactor_EMG);
cmc_low_freq = downsample(CMCaverage, downFactor_CMC);

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(emg_low_freq(1:data_length), muscle_SO(1:data_length), 'coeff');
correlation_SO(1,6) = max(corrs);
correlation_SO(2,6) = lags(find(ismember(corrs, max(corrs))));
[corrs, lags] = xcorr(emg_low_freq(1:data_length), cmc_low_freq(1:data_length), 'coeff');
correlation_CMC(1,6) = max(corrs);
correlation_CMC(2,6) = lags(find(ismember(corrs, max(corrs))));

aux1 = 0;
aux2 = 0;
aux3 = 0;
aux4 = 0;
for i=1:data_length
    % building the MAE elements
    aux1 = aux1 + abs(emg_low_freq(i)-muscle_SO(i));
    aux2 = aux2 + abs(emg_low_freq(i)-cmc_low_freq(i));
    % building the MAPE elements
    aux3 = aux3 + abs((emg_low_freq(i)-muscle_SO(i))/emg_low_freq(i));
    aux4 = aux4 + abs((emg_low_freq(i)-cmc_low_freq(i))/emg_low_freq(i));
end

MAE_SO(6) = aux1/data_length;
MAE_CMC(6) = aux2/data_length;
MAPE_SO(6) = 100/data_length*aux3;
MAPE_CMC(6) = 100/data_length*aux4;

% PLOT the comparison
time_SO = (0:step:1)*target_CMC(end);

subplot(3,4,6)
hold on
set(gca,'fontsize',15,'LineWidth',2)
h_mae=fill(betweenx_emg,betweeny_emg,[0.8 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_EMG, EMGstdup,'Color',[0.8 0.8 0.8])
h_mae=fill(betweenx,betweeny,[1 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_CMC, CMCaverage,'r','LineWidth',2)
plot(time_SO, muscle_SO,'Color' , '#7E2F8E', 'LineWidth',2)
if lim_axis
    axis([0 200 0 max_act])
end
% legend('EMG data', '', '', 'CMC 2019', 'SO result')
xlabel(label_x)
ylabel('activation')
title('Deltoid Post')

% scatter the downsampled points
if scatter_downsampled_point
    scatter(time_SO(1:data_length), muscle_SO(1:data_length), 'magenta', 'filled')
    scatter(time_SO(1:data_length), emg_low_freq(1:data_length), 'black', 'filled')
    scatter(time_SO(1:data_length), cmc_low_freq(1:data_length), 'red', 'filled')
end
hold off

%% (7) - Deltoid Middle
muscle_SO = muscles_SO(:,14);
name_CMC2019_results = append(folder_results_CMC2019, '\2019_deltoid_middle_EMG.mat');
load(name_CMC2019_results);
average_EMG = (betweeny_emg(1:length_EMG)+flip(betweeny_emg(length_EMG+1:end)))/2;

% RESAMPLE the EMG data and CMC solution to same frequency as SO solution
emg_low_freq = downsample(average_EMG, downFactor_EMG);
cmc_low_freq = downsample(CMCaverage, downFactor_CMC);

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(emg_low_freq(1:data_length), muscle_SO(1:data_length), 'coeff');
correlation_SO(1,7) = max(corrs);
correlation_SO(2,7) = lags(find(ismember(corrs, max(corrs))));
[corrs, lags] = xcorr(emg_low_freq(1:data_length), cmc_low_freq(1:data_length), 'coeff');
correlation_CMC(1,7) = max(corrs);
correlation_CMC(2,7) = lags(find(ismember(corrs, max(corrs))));

aux1 = 0;
aux2 = 0;
aux3 = 0;
aux4 = 0;
for i=1:data_length
    % building the MAE elements
    aux1 = aux1 + abs(emg_low_freq(i)-muscle_SO(i));
    aux2 = aux2 + abs(emg_low_freq(i)-cmc_low_freq(i));
    % building the MAPE elements
    aux3 = aux3 + abs((emg_low_freq(i)-muscle_SO(i))/emg_low_freq(i));
    aux4 = aux4 + abs((emg_low_freq(i)-cmc_low_freq(i))/emg_low_freq(i));
end

MAE_SO(7) = aux1/data_length;
MAE_CMC(7) = aux2/data_length;
MAPE_SO(7) = 100/data_length*aux3;
MAPE_CMC(7) = 100/data_length*aux4;

% PLOT the comparison
time_SO = (0:step:1)*target_CMC(end);

subplot(3,4,7)
hold on
set(gca,'fontsize',15,'LineWidth',2)
h_mae=fill(betweenx_emg,betweeny_emg,[0.8 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_EMG, EMGstdup,'Color',[0.8 0.8 0.8])
h_mae=fill(betweenx,betweeny,[1 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_CMC, CMCaverage,'r','LineWidth',2)
plot(time_SO, muscle_SO,'Color' , '#7E2F8E', 'LineWidth',2)
if lim_axis
    axis([0 200 0 max_act])
end
% legend('EMG data', '', '', 'CMC 2019', 'SO result')
xlabel(label_x)
ylabel('activation')
title('Deltoid Mid')

% scatter the downsampled points
if scatter_downsampled_point
    scatter(time_SO(1:data_length), muscle_SO(1:data_length), 'magenta', 'filled')
    scatter(time_SO(1:data_length), emg_low_freq(1:data_length), 'black', 'filled')
    scatter(time_SO(1:data_length), cmc_low_freq(1:data_length), 'red', 'filled')
end
hold off

%% (8) - Latissimus Dorsi
muscle_SO1 = muscles_SO(:,15);
muscle_SO2 = muscles_SO(:,16);
muscle_SO3 = muscles_SO(:,17);
name_CMC2019_results = append(folder_results_CMC2019, '\2019_latissimus_dorsi_EMG.mat');
load(name_CMC2019_results);
average_EMG = (betweeny_emg(1:length_EMG)+flip(betweeny_emg(length_EMG+1:end)))/2;

% RESAMPLE the EMG data and CMC solution to same frequency as SO solution
emg_low_freq = downsample(average_EMG, downFactor_EMG);
cmc_low_freq = downsample(CMCaverage, downFactor_CMC);

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(emg_low_freq(1:data_length), muscle_SO1(1:data_length), 'coeff');
correlation_SO(1,8) = max(corrs);
correlation_SO(2,8) = lags(find(ismember(corrs, max(corrs))));
[corrs, lags] = xcorr(emg_low_freq(1:data_length), cmc_low_freq(1:data_length), 'coeff');
correlation_CMC(1,8) = max(corrs);
correlation_CMC(2,8) = lags(find(ismember(corrs, max(corrs))));

aux1 = 0;
aux2 = 0;
aux3 = 0;
aux4 = 0;
for i=1:data_length
    % building the MAE elements
    aux1 = aux1 + abs(emg_low_freq(i)-muscle_SO1(i));
    aux2 = aux2 + abs(emg_low_freq(i)-cmc_low_freq(i));
    % building the MAPE elements
    aux3 = aux3 + abs((emg_low_freq(i)-muscle_SO1(i))/emg_low_freq(i));
    aux4 = aux4 + abs((emg_low_freq(i)-cmc_low_freq(i))/emg_low_freq(i));
end

MAE_SO(8) = aux1/data_length;
MAE_CMC(8) = aux2/data_length;
MAPE_SO(8) = 100/data_length*aux3;
MAPE_CMC(8) = 100/data_length*aux4;

% PLOT the comparison
time_SO = (0:step:1)*target_CMC(end);

subplot(3,4,8)
hold on
set(gca,'fontsize',15,'LineWidth',2)
h_mae=fill(betweenx_emg,betweeny_emg,[0.8 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_EMG, EMGstdup,'Color',[0.8 0.8 0.8])
h_mae=fill(betweenx,betweeny,[1 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_CMC, CMCaverage,'r','LineWidth',2)
plot(time_SO, muscle_SO1,'Color' , '#7E2F8E', 'LineWidth',2)
plot(time_SO, muscle_SO2, 'LineWidth',2)
plot(time_SO, muscle_SO3, 'LineWidth',2)
if lim_axis
    axis([0 200 0 max_act])
end
% legend('EMG data', '', '', 'CMC 2019', 'SO result')
xlabel(label_x)
ylabel('activation')
title('Lats Dorsi')

% scatter the downsampled points
if scatter_downsampled_point
    scatter(time_SO(1:data_length), muscle_SO1(1:data_length), 'magenta', 'filled')
    scatter(time_SO(1:data_length), emg_low_freq(1:data_length), 'black', 'filled')
    scatter(time_SO(1:data_length), cmc_low_freq(1:data_length), 'red', 'filled')
end
hold off

%% (9) - Pectoralis Major Clavicle
muscle_SO = muscles_SO(:,18);
name_CMC2019_results = append(folder_results_CMC2019, '\2019_pectoralis_major_clavical_EMG.mat');
load(name_CMC2019_results);
average_EMG = (betweeny_emg(1:length_EMG)+flip(betweeny_emg(length_EMG+1:end)))/2;

% RESAMPLE the EMG data and CMC solution to same frequency as SO solution
emg_low_freq = downsample(average_EMG, downFactor_EMG);
cmc_low_freq = downsample(CMCaverage, downFactor_CMC);

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(emg_low_freq(1:data_length), muscle_SO(1:data_length), 'coeff');
correlation_SO(1,9) = max(corrs);
correlation_SO(2,9) = lags(find(ismember(corrs, max(corrs))));
[corrs, lags] = xcorr(emg_low_freq(1:data_length), cmc_low_freq(1:data_length), 'coeff');
correlation_CMC(1,9) = max(corrs);
correlation_CMC(2,9) = lags(find(ismember(corrs, max(corrs))));

aux1 = 0;
aux2 = 0;
aux3 = 0;
aux4 = 0;
for i=1:data_length
    % building the MAE elements
    aux1 = aux1 + abs(emg_low_freq(i)-muscle_SO(i));
    aux2 = aux2 + abs(emg_low_freq(i)-cmc_low_freq(i));
    % building the MAPE elements
    aux3 = aux3 + abs((emg_low_freq(i)-muscle_SO(i))/emg_low_freq(i));
    aux4 = aux4 + abs((emg_low_freq(i)-cmc_low_freq(i))/emg_low_freq(i));
end

MAE_SO(9) = aux1/data_length;
MAE_CMC(9) = aux2/data_length;
MAPE_SO(9) = 100/data_length*aux3;
MAPE_CMC(9) = 100/data_length*aux4;

% PLOT the comparison
time_SO = (0:step:1)*target_CMC(end);

subplot(3,4,9)
hold on
set(gca,'fontsize',15,'LineWidth',2)
h_mae=fill(betweenx_emg,betweeny_emg,[0.8 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_EMG, EMGstdup,'Color',[0.8 0.8 0.8])
h_mae=fill(betweenx,betweeny,[1 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_CMC, CMCaverage,'r','LineWidth',2)
plot(time_SO, muscle_SO,'Color' , '#7E2F8E', 'LineWidth',2)
if lim_axis
    axis([0 200 0 max_act])
end
% legend('EMG data', '', '', 'CMC 2019', 'SO result')
xlabel(label_x)
ylabel('activation')
title('Pect Major Clav')

% scatter the downsampled points
if scatter_downsampled_point
    scatter(time_SO(1:data_length), muscle_SO(1:data_length), 'magenta', 'filled')
    scatter(time_SO(1:data_length), emg_low_freq(1:data_length), 'black', 'filled')
    scatter(time_SO(1:data_length), cmc_low_freq(1:data_length), 'red', 'filled')
end
hold off

%% (10) - Teres Major
muscle_SO =muscles_SO(:,21);
name_CMC2019_results = append(folder_results_CMC2019, '\2019_teres_major_EMG.mat');
load(name_CMC2019_results);
average_EMG = (betweeny_emg(1:length_EMG)+flip(betweeny_emg(length_EMG+1:end)))/2;

% RESAMPLE the EMG data and CMC solution to same frequency as SO solution
emg_low_freq = downsample(average_EMG, downFactor_EMG);
cmc_low_freq = downsample(CMCaverage, downFactor_CMC);

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(emg_low_freq(1:data_length), muscle_SO(1:data_length), 'coeff');
correlation_SO(1,10) = max(corrs);
correlation_SO(2,10) = lags(find(ismember(corrs, max(corrs))));
[corrs, lags] = xcorr(emg_low_freq(1:data_length), cmc_low_freq(1:data_length), 'coeff');
correlation_CMC(1,10) = max(corrs);
correlation_CMC(2,10) = lags(find(ismember(corrs, max(corrs))));

aux1 = 0;
aux2 = 0;
aux3 = 0;
aux4 = 0;
for i=1:data_length
    % building the MAE elements
    aux1 = aux1 + abs(emg_low_freq(i)-muscle_SO(i));
    aux2 = aux2 + abs(emg_low_freq(i)-cmc_low_freq(i));
    % building the MAPE elements
    aux3 = aux3 + abs((emg_low_freq(i)-muscle_SO(i))/emg_low_freq(i));
    aux4 = aux4 + abs((emg_low_freq(i)-cmc_low_freq(i))/emg_low_freq(i));
end

MAE_SO(10) = aux1/data_length;
MAE_CMC(10) = aux2/data_length;
MAPE_SO(10) = 100/data_length*aux3;
MAPE_CMC(10) = 100/data_length*aux4;

% PLOT the comparison
time_SO = (0:step:1)*target_CMC(end);

subplot(3,4,10)
hold on
set(gca,'fontsize',15,'LineWidth',2)
h_mae=fill(betweenx_emg,betweeny_emg,[0.8 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_EMG, EMGstdup,'Color',[0.8 0.8 0.8])
h_mae=fill(betweenx,betweeny,[1 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_CMC, CMCaverage,'r','LineWidth',2)
plot(time_SO, muscle_SO,'Color' , '#7E2F8E', 'LineWidth',2)
if lim_axis
    axis([0 200 0 max_act])
end
% legend('EMG data', '', '', 'CMC 2019', 'SO result')
xlabel(label_x)
ylabel('activation')
title('Teres Major')

% scatter the downsampled points
if scatter_downsampled_point
    scatter(time_SO(1:data_length), muscle_SO(1:data_length), 'magenta', 'filled')
    scatter(time_SO(1:data_length), emg_low_freq(1:data_length), 'black', 'filled')
    scatter(time_SO(1:data_length), cmc_low_freq(1:data_length), 'red', 'filled')
end
hold off

%% (11) - Infraspinatus
muscle_SO1 = muscles_SO(:,22);
muscle_SO2 = muscles_SO(:,23);
name_CMC2019_results = append(folder_results_CMC2019, '\2019_infraspinatus_EMG.mat');
load(name_CMC2019_results);
average_EMG = (betweeny_emg(1:length_EMG)+flip(betweeny_emg(length_EMG+1:end)))/2;

% RESAMPLE the EMG data and CMC solution to same frequency as SO solution
emg_low_freq = downsample(average_EMG, downFactor_EMG);
cmc_low_freq = downsample(CMCaverage, downFactor_CMC);

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(emg_low_freq(1:data_length), muscle_SO1(1:data_length), 'coeff');
correlation_SO(1,11) = max(corrs);
correlation_SO(2,11) = lags(find(ismember(corrs, max(corrs))));
[corrs, lags] = xcorr(emg_low_freq(1:data_length), cmc_low_freq(1:data_length), 'coeff');
correlation_CMC(1,11) = max(corrs);
correlation_CMC(2,11) = lags(find(ismember(corrs, max(corrs))));

aux1 = 0;
aux2 = 0;
aux3 = 0;
aux4 = 0;
for i=1:data_length
    % building the MAE elements
    aux1 = aux1 + abs(emg_low_freq(i)-muscle_SO1(i));
    aux2 = aux2 + abs(emg_low_freq(i)-cmc_low_freq(i));
    % building the MAPE elements
    aux3 = aux3 + abs((emg_low_freq(i)-muscle_SO1(i))/emg_low_freq(i));
    aux4 = aux4 + abs((emg_low_freq(i)-cmc_low_freq(i))/emg_low_freq(i));
end

MAE_SO(11) = aux1/data_length;
MAE_CMC(11) = aux2/data_length;
MAPE_SO(11) = 100/data_length*aux3;
MAPE_CMC(11) = 100/data_length*aux4;

% PLOT the comparison
time_SO = (0:step:1)*target_CMC(end);

subplot(3,4,11)
hold on
set(gca,'fontsize',15,'LineWidth',2)
h_mae=fill(betweenx_emg,betweeny_emg,[0.8 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_EMG, EMGstdup,'Color',[0.8 0.8 0.8])
h_mae=fill(betweenx,betweeny,[1 0.8 0.8]);
set(h_mae,'EdgeColor','none')
plot(target_CMC, CMCaverage,'r','LineWidth',2)
plot(time_SO, muscle_SO1,'Color' , '#7E2F8E', 'LineWidth',2)
plot(time_SO, muscle_SO2, 'LineWidth',2)
if lim_axis
    axis([0 200 0 max_act])
end

% scatter the downsampled points
if scatter_downsampled_point
    scatter(time_SO(1:data_length), muscle_SO1(1:data_length), 'magenta', 'filled')
    scatter(time_SO(1:data_length), emg_low_freq(1:data_length), 'black', 'filled')
    scatter(time_SO(1:data_length), cmc_low_freq(1:data_length), 'red', 'filled')
end
hold off

xlabel(label_x)
ylabel('activation')
title('Infra')

%% Create Legend
subplot(3,4,12)
hold on
set(gca,'XColor', 'none','YColor','none')
set(gca, 'color', 'none');
set(gca,'fontsize',15,'LineWidth',2)
plot(1, 1,'Color',[0.8 0.8 0.8],'LineWidth',3)
plot(1,1, 'r', 'LineWidth',3)
plot(1, 1,'Color' , '#7E2F8E', 'LineWidth',3)

legend('EMG data','CMC 2019', 'SO result');
legend('Location','east')

f_muscle_activ.WindowState = 'maximized';
name_fig_muscle_activ = append(experiment, '_SO_EMG_CMC.png');
saveas(f_muscle_activ, name_fig_muscle_activ)

%% VISUALIZE MAE and CORR for ACTIVATIONS
muscle_name_list = ["Trap Scap Mid", "Trap Scap Sup", "Trap Scap Inf", "Serratus Ant", "Deltoid Ant", "Deltoid Post", "Deltoid Mid", "Lats Dorsi", "Pect Major Clav", "Teres Major", "Infra"];
figure
h_mae = heatmap(MAE_SO', 'ColorLimits', [0,0.4], 'Colormap', jet(4));
colorbar
ax = gca;
ax.YDisplayLabels = muscle_name_list;
set(gcf,'position',[0,0,300,650])
title(experiment)

nameFig =  append(experiment, '_MAE',  '.png');
saveas(h_mae, nameFig)

figure
h_corr = heatmap(correlation_SO(1,:)', 'ColorLimits', [-1,1], 'Colormap', hot);
colorbar
ax = gca;
ax.YDisplayLabels = muscle_name_list;
set(gcf,'position',[300,0,300,650])
title(experiment)

nameFig =  append(experiment, '_CORR',  '.png');
saveas(h_corr, nameFig)

%% RETURN AGGREGATED METRICS ---------------------------------------------------------
MAE_CMC_all = mean(MAE_CMC);
MAPE_CMC_all = mean(MAPE_CMC);
correlation_CMC_all = mean(correlation_CMC(1,:));

MAE_SO_all = mean(MAE_SO);
MAPE_SO_all = mean(MAPE_SO);
correlation_SO_all = mean(correlation_SO(1,:));

% visualize on screen
string_MAE = sprintf('MAE                %.2f                   %.2f \n', MAE_CMC_all, MAE_SO_all);
string_MAPE = sprintf('MAPE               %.2f                 %.2f \n', MAPE_CMC_all, MAPE_SO_all);
string_CORR = sprintf('CORRELATION        %.2f                   %.2f \n', correlation_CMC_all, correlation_SO_all);
disp('                    CMC2019                SO         ')
fprintf(string_MAE);
fprintf(string_MAPE);
fprintf(string_CORR);

%% Save the metrics related to every muscle
metrics.muscleNames = muscle_name_list;
metrics.mae_SO = MAE_SO;
metrics.mae_CMC = MAE_CMC;
metrics.corr_SO = correlation_SO;
metrics.corr_CMC = correlation_CMC;

name_file = append(experiment, '_metrics_correctMuscles');
save(name_file, 'metrics');
