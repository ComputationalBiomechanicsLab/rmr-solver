function [MAE, correlation] = compare_muscle_activations(file_results_opensim, file_filtered_EMG, starting_sample)
% This function allows to compare the results of the muscle activations
% coming from its two inputs:
% * file_results_opensim = file (fullpath+name) containing the results of
%   the static optimization as returned in CustomStaticOptimization_adapted.m
% * file_filtered_EMG = file (fullpath+name) containing the muscle
%   activations coming from filtered EMG data
% * starting sample = first sample of the data to be considered
%
% The function plots the comparisons for all the muscles. Note that the
% frequency of the data in the first file is retrieved from the file
% itself, while for the EMG data we assume it to be recorded at 1500 Hz (as
% this is the standard frequency in Clark's dataset).
% The function also computes and returns the value of the Mean Absolute
% Error (MAE) between the predicted muscle activations and the filtered
% EMG, as well as the correlation coefficient between the two signals.

if nargin>2
    initial_sample = starting_sample;
else
    initial_sample = 1;
end

fit_to_1 = 1;   % flag to plot muscle activations between 0 and 1

freq_Clark = 1500; % Hz

% load the .mat file containing the muscle activations, names of muscles
% and frequency of the data ('xsol', 'muscle_order', 'frequency_solution')
load(file_results_opensim, 'xsol', 'muscle_order', 'frequency_solution');
ma_opensim = xsol(:, 1:33);    % keeep the values coming from real muscles

% load the .mat file containing the filtered EMG data
% ('muscle_activations')
load(file_filtered_EMG);
ma_EMG = muscle_activations;

% define the final time of the muscle activations coming from opensim
end_time = size(xsol,1)/frequency_solution;

% crop accordingly the filtered EMG data
end_sample_EMG = freq_Clark * end_time;
ma_EMG = ma_EMG(1:end_sample_EMG, :);

% resample the EMG data to match 'frequency_solution'
step = round(freq_Clark/frequency_solution);
ma_EMG = ma_EMG(1:step:end, :);

% make sure the two sets of data are comparable
new_size = size(ma_opensim, 1);
ma_EMG = ma_EMG(1:new_size, :);

% allocate array to contain MAE, and correlation coefficients
% (assuming 14 muscles, Clark's dataset)
MAE = zeros(1,14);
correlation = zeros(2,14); % second row indicates the lag value at which the corresponding correlation is found

%% PROCESS and PLOT the RESULTS -------------------------------------------

side = ceil(sqrt(size(ma_EMG,2)));

pgc = linspace(0, 100, size(ma_EMG,1)-initial_sample+1);

figure
%% 1. Anterior Deltoid
% comparing 'DeltoideusClavicle_A' to 'Anterior Deltoid'

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(ma_EMG(initial_sample:end, 1), ma_opensim(initial_sample:end, 12), 'coeff');
correlation(1,1) = max(corrs);
correlation(2,1) = lags(find(ismember(corrs, max(corrs))));

aux = 0;
for i=initial_sample:new_size
    % building the MAE elements
    aux = aux + abs(ma_EMG(i, 1)-ma_opensim(i,12));
end
MAE(1) = aux/(new_size-initial_sample+1);

subplot(side,side,1)
title('Anterior Deltoid')
hold on
plot(pgc, ma_EMG(initial_sample:end,1));
plot(pgc, ma_opensim(initial_sample:end,12));
legend('Anterior Deltoid\_EMG', 'DeltoideusClavicle\_A\_OS')
xlabel("% of motion")
ylabel("muscle activation")
if fit_to_1
    ylim([0 1])
end
grid on 
hold off

%% 2. Middle Deltoid
% comparing 'DeltoideusScapula_M' to 'Middle Deltoid'

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(ma_EMG(initial_sample:end, 2), ma_opensim(initial_sample:end, 14), 'coeff');
correlation(1,2) = max(corrs);
correlation(2,2) = lags(find(ismember(corrs, max(corrs))));

aux = 0;
for i=initial_sample:new_size
    % building the MAE elements
    aux = aux + abs(ma_EMG(i, 2)-ma_opensim(i,14));
end
MAE(2) = aux/(new_size-initial_sample+1);

subplot(side,side,2)
title('Middle Deltoid')
hold on
plot(pgc, ma_EMG(initial_sample:end,2));
plot(pgc, ma_opensim(initial_sample:end,14));
legend('Anterior Deltoid\_EMG', 'DeltoideusClavicle\_M\_OS')
if fit_to_1
    ylim([0 1])
end
grid on 
hold off

%% 3. Posterior Deltoid
% comparing 'DeltoideusScapula_P' to 'Posterior Deltoid'

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(ma_EMG(initial_sample:end, 3), ma_opensim(initial_sample:end, 13), 'coeff');
correlation(1,3) = max(corrs);
correlation(2,3) = lags(find(ismember(corrs, max(corrs))));

aux = 0;
for i=initial_sample:new_size
    % building the MAE elements
    aux = aux + abs(ma_EMG(i, 3)-ma_opensim(i,13));
end
MAE(3) = aux/(new_size-initial_sample+1);

subplot(side,side,3)
title('Posterior Deltoid')
hold on
plot(pgc, ma_EMG(initial_sample:end,3));
plot(pgc, ma_opensim(initial_sample:end,13));
legend('Posterior Deltoid\_EMG', 'DeltoideusClavicle\_P\_OS')
if fit_to_1
    ylim([0 1])
end
grid on 
hold off

%% 4. Biceps
% comparing 'BIC_long' and 'BIC_brevis' to 'Biceps' 
% # with our current model this might not be accurate since the elbow is
% actuated individually by a CoordinateActuator

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(ma_EMG(initial_sample:end, 4), ma_opensim(initial_sample:end, 33), 'coeff');
correlation(1,4) = max(corrs);
correlation(2,4) = lags(find(ismember(corrs, max(corrs))));

aux = 0;
for i=initial_sample:new_size
    % building the MAE elements
    aux = aux + abs(ma_EMG(i, 4)-ma_opensim(i,33));
end
MAE(4) = aux/(new_size-initial_sample+1);

subplot(side,side,4)
title("Biceps")
hold on
plot(pgc, ma_EMG(initial_sample:end,4))
plot(pgc, ma_opensim(initial_sample:end,32))
plot(pgc, ma_opensim(initial_sample:end,33))
legend('biceps\_EMG', 'BIC\_long\_OS', 'BIC\_brevis\_OS')
if fit_to_1
    ylim([0 1])
end
grid on
hold off

%% 5. Triceps
% comparing 'TRIlong' to 'Triceps'

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(ma_EMG(initial_sample:end, 5), ma_opensim(initial_sample:end, 31), 'coeff');
correlation(1,5) = max(corrs);
correlation(2,5) = lags(find(ismember(corrs, max(corrs))));

aux = 0;
for i=initial_sample:new_size
    % building the MAE elements
    aux = aux + abs(ma_EMG(i, 5)-ma_opensim(i,31));
end
MAE(5) = aux/(new_size-initial_sample+1);

subplot(side,side,5)
title("Triceps")
hold on
plot(pgc, ma_EMG(initial_sample:end,5))
plot(pgc, ma_opensim(initial_sample:end,31))
legend('Triceps\_EMG', 'TRIlong\_OS')
if fit_to_1
    ylim([0 1])
end
grid on
hold off

%% 6. Infra
% comparing 'Infraspinatus_I' and 'Infraspinatus_S' to 'Infra'

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(ma_EMG(initial_sample:end, 6), ma_opensim(initial_sample:end, 22), 'coeff');
correlation(1,6) = max(corrs);
correlation(2,6) = lags(find(ismember(corrs, max(corrs))));

aux = 0;
for i=initial_sample:new_size
    % building the MAE elements
    aux = aux + abs(ma_EMG(i, 6)-ma_opensim(i,22));
end
MAE(6) = aux/(new_size-initial_sample+1);

subplot(side,side,6)
title("Infra")
hold on
plot(pgc, ma_EMG(initial_sample:end,6))
plot(pgc, ma_opensim(initial_sample:end,22))
plot(pgc, ma_opensim(initial_sample:end,23))
legend('Infra\_EMG', 'Infra\_I\_OS', 'Infra\_S\_OS')
if fit_to_1
    ylim([0 1])
end
grid on
hold off

%% 7. Supra
% comparing 'Supraspinatus_P' and 'Supraspinatus_A' to 'Supra'

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(ma_EMG(initial_sample:end, 7), ma_opensim(initial_sample:end, 29), 'coeff');
correlation(1,7) = max(corrs);
correlation(2,7) = lags(find(ismember(corrs, max(corrs))));

aux = 0;
for i=initial_sample:new_size
    % building the MAE elements
    aux = aux + abs(ma_EMG(i, 7)-ma_opensim(i,29));
end
MAE(7) = aux/(new_size-initial_sample+1);

subplot(side,side,7)
title("Supfra")
hold on
plot(pgc, ma_EMG(initial_sample:end,7))
plot(pgc, ma_opensim(initial_sample:end,29))
plot(pgc, ma_opensim(initial_sample:end,30))
legend('Supra\_EMG', 'Supraspinatus\_P\_OS', 'Supraspinatus\_A\_OS')
if fit_to_1
    ylim([0 1])
end
grid on
hold off

%% 8. PecMajorC
% comparing 'PectoralisMajorClavicle_S' and 'PecMajorC'

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(ma_EMG(initial_sample:end, 8), ma_opensim(initial_sample:end, 18), 'coeff');
correlation(1,8) = max(corrs);
correlation(2,8) = lags(find(ismember(corrs, max(corrs))));

aux = 0;
for i=initial_sample:new_size
    % building the MAE elements
    aux = aux + abs(ma_EMG(i, 8)-ma_opensim(i,18));
end
MAE(8) = aux/(new_size-initial_sample+1);

subplot(side,side,8)
title("PecMajorC (clavicular insertion of PecMajor)")
hold on
plot(pgc, ma_EMG(initial_sample:end,8))
plot(pgc, ma_opensim(initial_sample:end,18))
legend('PecMajorC\_EMG', 'PectoralisMajorClavicle\_S')
if fit_to_1
    ylim([0 1])
end
grid on
hold off

%% 9. PecMajorS
% comparing "PectoralisMajorThorax_I" and "PectoralisMajorThorax_M" to 'PecMajorS'

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(ma_EMG(initial_sample:end, 9), ma_opensim(initial_sample:end, 19), 'coeff');
correlation(1,9) = max(corrs);
correlation(2,9) = lags(find(ismember(corrs, max(corrs))));

aux = 0;
for i=initial_sample:new_size
    % building the MAE elements
    aux = aux + abs(ma_EMG(i, 9)-ma_opensim(i,19));
end
MAE(9) = aux/(new_size-initial_sample+1);

subplot(side,side,9)
title("PecMajorS (sternal insertion of PecMajor)")
hold on
plot(pgc, ma_EMG(initial_sample:end,9))
plot(pgc, ma_opensim(initial_sample:end,19))
plot(pgc, ma_opensim(initial_sample:end,20))
legend('PecMajorS\_EMG', 'PectoralisMajorThorax\_I\_OS', 'PectoralisMajorThorax\_M\_OS')
if fit_to_1
    ylim([0 1])
end
grid on
hold off

%% 10. Lats
% comparing "LatissimusDorsi_S", "LatissimusDorsi_M", "LatissimusDorsi_I"
% to 'Lats'

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(ma_EMG(initial_sample:end, 10), ma_opensim(initial_sample:end, 15), 'coeff');
correlation(1,10) = max(corrs);
correlation(2,10) = lags(find(ismember(corrs, max(corrs))));

aux = 0;
for i=initial_sample:new_size
    % building the MAE elements
    aux = aux + abs(ma_EMG(i, 10)-ma_opensim(i,15));
end
MAE(10) = aux/(new_size-initial_sample+1);

subplot(side,side,10)
title("Lats (latissimus dorsi)")
hold on
plot(pgc, ma_EMG(initial_sample:end,10))
plot(pgc, ma_opensim(initial_sample:end,15))
plot(pgc, ma_opensim(initial_sample:end,16))
plot(pgc, ma_opensim(initial_sample:end,17))
legend('Lats\_EMG', 'LatissimusDorsi\_S\_OS', 'LatissimusDorsi\_M\_OS', 'LatissimusDorsi\_I\_OS')
if fit_to_1
    ylim([0 1])
end
grid on
hold off

%% 11. Sert
% comparing "SerratusAnterior_I", "SerratusAnterior_M", "SerratusAnterior_S"    

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(ma_EMG(initial_sample:end, 11), ma_opensim(initial_sample:end, 5), 'coeff');
correlation(1,11) = max(corrs);
correlation(2,11) = lags(find(ismember(corrs, max(corrs))));

aux = 0;
for i=initial_sample:new_size
    % building the MAE elements
    aux = aux + abs(ma_EMG(i, 11)-ma_opensim(i,5));
end
MAE(11) = aux/(new_size-initial_sample+1);

subplot(side,side,11)
title("Sert (serratus anterior")
hold on
plot(pgc, ma_EMG(initial_sample:end,11))
plot(pgc, ma_opensim(initial_sample:end,5))
plot(pgc, ma_opensim(initial_sample:end,6))
plot(pgc, ma_opensim(initial_sample:end,7))
legend('Sert\_EMG', 'SerratusAnterior\_I\_OS' , 'SerratusAnterior\_M\_OS', 'SerratusAnterior\_S\_OS')
if fit_to_1
    ylim([0 1])
end
grid on
hold off

%% 12. TrapLow
% comparing 'TrapeziusScapula_I' to 'TrapLow'    

% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(ma_EMG(initial_sample:end, 12), ma_opensim(initial_sample:end, 3), 'coeff');
correlation(1,12) = max(corrs);
correlation(2,12) = lags(find(ismember(corrs, max(corrs))));

aux = 0;
for i=initial_sample:new_size
    % building the MAE elements
    aux = aux + abs(ma_EMG(i, 12)-ma_opensim(i,3));
end
MAE(12) = aux/(new_size-initial_sample+1);

subplot(side,side,12)
title("TrapLow")
hold on
plot(pgc, ma_EMG(initial_sample:end,12))
plot(pgc, ma_opensim(initial_sample:end,3))
legend('TrapLow\_EMG', 'TrapeziusScapula\_I\_OS')
if fit_to_1
    ylim([0 1])
end
grid on
hold off

%% 13. TrapMid
% comparing 'TrapeziusScapula_M' with 'TrapMid'
    
% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(ma_EMG(initial_sample:end, 13), ma_opensim(initial_sample:end, 1), 'coeff');
correlation(1,13) = max(corrs);
correlation(2,13) = lags(find(ismember(corrs, max(corrs))));

aux = 0;
for i=initial_sample:new_size
    % building the MAE elements
    aux = aux + abs(ma_EMG(i, 13)-ma_opensim(i,1));
end
MAE(13) = aux/(new_size-initial_sample+1);

subplot(side,side,13)
title("TrapMid")
hold on
plot(pgc, ma_EMG(initial_sample:end,13))
plot(pgc, ma_opensim(initial_sample:end,1))
legend('TrapMid\_EMG', 'TrapeziusScapula\_M\_OS')
if fit_to_1
    ylim([0 1])
end
grid on
hold off

%% 14. TrapUp
% comparing 'TrapeziusScapula_s', 'TrapeziusClavicle_S' to 'TrapUp'
    
% getting the correlation at different lags - retain max correlation and
% corresponding lag
[corrs, lags] = xcorr(ma_EMG(initial_sample:end, 14), ma_opensim(initial_sample:end, 2), 'coeff');
correlation(1,14) = max(corrs);
correlation(2,14) = lags(find(ismember(corrs, max(corrs))));

aux = 0;
for i=initial_sample:new_size
    % building the MAE elements
    aux = aux + abs(ma_EMG(i, 14)-ma_opensim(i,2));
end
MAE(14) = aux/(new_size-initial_sample+1);

subplot(side,side,14)
title("TrapUp")
hold on
plot(pgc, ma_EMG(initial_sample:end,14))
plot(pgc, ma_opensim(initial_sample:end,2))
plot(pgc, ma_opensim(initial_sample:end,4))
legend('TrapUp\_EMG', 'TrapeziusScapula\_S\_OS', 'TrapeziusClavicle\_S\_OS')
if fit_to_1
    ylim([0 1])
end
xlabel("% of movement")
ylabel("muscle activation")
grid on
hold off

%% RETURN METRICS ---------------------------------------------------------
MAE_all = mean(MAE);
correlation_all = mean(correlation(1,:));

% visualize on screen
string_MAE = sprintf('MAE                %.2f', MAE_all);
string_CORR = sprintf('CORRELATION        %.2f', correlation_all);
disp(string_MAE)
disp(string_CORR)
