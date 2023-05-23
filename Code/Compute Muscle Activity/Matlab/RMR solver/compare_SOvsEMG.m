activation_file_SO = 'C:\Users\italobelli\Desktop\GitHub\PTbot_officialCBL\Personal_Results\muscle_activations_ACC_experiment_6.mat';
activation_file_EMG = 'U:\PTbot\DataShare\dataset_Clark\SMAP\Muscle Activations\ACC\107.mat';

load(activation_file_EMG);
ma_EMG = muscle_activations;

load(activation_file_SO, 'xsol', 'muscle_order', 'frequency_solution');

% define the time in which the CMC was run
time_start = 0;
time_end = 7;

% transform time bounds into sample index for SO solution

initial_sample = time_start * frequency_solution+1;
final_sample = time_end * frequency_solution;

% get the correct interval of data from SO solution
ma_SO = xsol(initial_sample:final_sample, 1:33);
reserv_acts_SO = xsol(initial_sample:final_sample, 34:end);

% get the correct samples from EMG solution
% crop accordingly the filtered EMG data
start_time_EMG = 1500 * time_start+1;
end_sample_EMG = 1500 * time_end;
ma_EMG = ma_EMG(start_time_EMG:end_sample_EMG, :);

% resample the EMG data to match 'frequency_solution'
step = round(1500/frequency_solution);
ma_EMG = ma_EMG(1:step:end, :);

% Plot the comparison
side = ceil(sqrt(size(ma_EMG,2)));
figure

% 1. Anterior Deltoid
% comparing 'DeltoideusClavicle_A' to 'Anterior Deltoid'
subplot(side,side,1)
title('Anterior Deltoid')
hold on
plot(ma_EMG(:,1));
plot(xsol(:,12));
legend('Anterior Deltoid\_EMG','DeltoideusClavicle\_A\_SO')
xlabel("samples @ 5 Hz at steady state")
ylabel("muscle activation")
ylim([0, 1])
grid on 
hold off

% 2. Middle Deltoid
% comparing 'DeltoideusScapula_M' to 'Middle Deltoid'
subplot(side,side,2)
title('Middle Deltoid')
hold on
plot(ma_EMG(:,2));
plot(xsol(:,14));
legend('Anterior Deltoid\_EMG', 'DeltoideusScapula\_M\_SO')
ylim([0, 1])
grid on 
hold off

% 3. Posterior Deltoid
% comparing 'DeltoideusScapula_P' to 'Posterior Deltoid'
subplot(side,side,3)
title('Posterior Deltoid')
hold on
plot(ma_EMG(:,3));
plot(xsol(:,13));
legend('Posterior Deltoid\_EMG', 'DeltoideusScapula\_P\_SO')
ylim([0, 1])
grid on 
hold off

% 4. Biceps
% comparing 'BIC_long' and 'BIC_brevis' to 'Biceps'
subplot(side,side,4)
title("Biceps")
hold on
plot(ma_EMG(:,4))
plot(xsol(:,32))
plot(xsol(:,33))
legend('biceps\_EMG', 'BIC\_long\_SO', 'BIC\_brevis\_SO')
ylim([0, 1])
grid on
hold off

% 5. Triceps
% comparing 'TRIlong' to 'Triceps'
subplot(side,side,5)
title("Triceps")
hold on
plot(ma_EMG(:,5))
plot(xsol(:,31))
legend('Triceps\_EMG', 'TRIlong\_SO')
ylim([0, 1])
grid on
hold off

% 6. Infra
% comparing 'Infraspinatus_I' and 'Infraspinatus_S' to 'Infra'
subplot(side,side,6)
title("Infra")
hold on
plot(ma_EMG(:,6))
plot(xsol(:,23))
legend('Infra\_EMG', 'Infra\_S\_SO')
ylim([0, 1])
grid on
hold off

% 7. Supra
% comparing 'Supraspinatus_P' and 'Supraspinatus_A' to 'Supra'
subplot(side,side,7)
title("Supfra")
hold on
plot(ma_EMG(:,7))
plot(xsol(:,29))
plot(xsol(:,30))
legend('Supra\_EMG', 'Supraspinatus\_P\_SO', 'Supraspinatus\_A\_SO')
ylim([0, 1])
grid on
hold off

% 8. PecMajorC
% comparing 'PectoralisMajorClavicle_S' and 'PecMajorC'
subplot(side,side,8)
title("PecMajorC (clavicular insertion of PecMajor)")
hold on
plot(ma_EMG(:,8))
plot(xsol(:,18))
legend('PecMajorC\_EMG', 'PectoralisMajorClavicle\_S\_SO')
ylim([0, 1])
grid on
hold off

% 9. PecMajorS
% comparing "PectoralisMajorThorax_I" and "PectoralisMajorThorax_M" to 'PecMajorS'
subplot(side,side,9)
title("PecMajorS (sternal insertion of PecMajor)")
hold on
plot(ma_EMG(initial_sample:end,9))
plot(xsol(:,19))
plot(xsol(:,20))
legend('PecMajorS\_EMG', 'PectoralisMajorThorax\_I\_SO', 'PectoralisMajorThorax\_M\_SO')
ylim([0, 1])
grid on
hold off

% 10. Lats
% comparing "LatissimusDorsi_S", "LatissimusDorsi_M", "LatissimusDorsi_I"
% to 'Lats'
subplot(side,side,10)
title("Lats (latissimus dorsi)")
hold on
plot(ma_EMG(:,10))
plot(xsol(:,15))
plot(xsol(:,16))
plot(xsol(:,17))
legend('Lats\_EMG', 'LatissimusDorsi\_S\_CMC', 'LatissimusDorsi\_M\_SO', 'LatissimusDorsi\_I\_SO')
ylim([0, 1])
grid on
hold off

% 11. Sert
% comparing "SerratusAnterior_I", "SerratusAnterior_M", "SerratusAnterior_S"    
subplot(side,side,11)
title("Sert (serratus anterior")
hold on
plot(ma_EMG(:,11))
plot(xsol(:,5))
legend('Sert\_EMG', 'SerratusAnterior\_I\_SO')
ylim([0, 1])
grid on
hold off

% 12. TrapLow
% comparing 'TrapeziusScapula_I' to 'TrapLow'
subplot(side,side,12)
title("TrapLow")
hold on
plot(ma_EMG(:,12))
plot(xsol(:,3))
legend('TrapLow\_EMG', 'TrapeziusScapula\_I\_SO')
ylim([0, 1])
grid on
hold off

% 13. TrapMid
% comparing 'TrapeziusScapula_M' with 'TrapMid'
subplot(side,side,13)
title("TrapMid")
hold on
plot(ma_EMG(:,13))
plot(xsol(:,1))
legend()
ylim([0, 1])
grid on
hold off

% 14. TrapUp
% comparing 'TrapeziusScapula_s', 'TrapeziusClavicle_S' to 'TrapUp'
subplot(side,side,14)
title("TrapUp")
hold on
plot(ma_EMG(:,14))
plot(xsol(:,2))
plot(xsol(:,4))
legend('TrapUp\_EMG', 'TrapeziusScapula\_S\_SO', 'TrapeziusClavicle\_S\_SO')
ylim([0, 1])
xlabel("samples @ 5 Hz at steady state")
ylabel("muscle activation")
grid on
hold off