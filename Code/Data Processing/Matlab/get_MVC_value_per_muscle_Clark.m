function value_MVC = get_MVC_value_per_muscle_Clark(csv_file_Clark, muscle)
% this function can be used to get the MVC value of a particular muscle
% from Clark's dataset. The inputs are
% * csv_file: fullpath+name of the csv file in which the EMG values are
%             stored;
% * muscle: string containing the name of the muscle to focus on
%           (acceptable strings are specified below)
% The output of the function is a single (positive) float number, indicating
% the maximum EMG value [V] recorded in the file for that particular
% muscle.
%
% Note that the procedure for finding the MVC involves removing the DC bias
% from the data, using a high pass Butterworth filter, rectifying the data, 
% using a lowpass Butterworth filter and then filtering with a moving
% averge (following the filtering done by Clark in section 2.4 of his
% paper). The parameters adopted for the filters are highlighted at the
% beginning of the function - if needed they can be changed.
%
% The muscles that are present in the dataset are:
% 'Anterior Deltoid', 'Middle Deltoid', 'Posterior Deltoid', 'Biceps', 
% 'Triceps', 'Infra, 'Supra', 'PecMajorC', 'PecMajorS', 'Lats', 'Sert', 
% 'TrapLow', 'TrapMid', 'TrapUp'.
% The csv files can be found in the DatShare folder of the PTbot shared drive
% under 'dataset_Clark/SMAP/SMAP 4 RAW/(name of the subject)'


%% Stating the frequency values adopted here [Hz]
% cut-off frequency Hz of the highpass 4th order Butterworth filter applied
cutoff_high = 30;

% cut-off frequency Hz of the lowpass 4th order Butterworth filter applied
cutoff_low = 4;

% acquisition frequency for the data (always 1500 Hz for Clark's)
fs = 1500;

%% Checking the inputs
% if needed, we can automatize the selection of the column corresponding to
% a certain muscle
muscles_in_order = ["Anterior Deltoid", "Middle Deltoid", "Posterior Deltoid", ...
                     "Biceps", "Triceps", "Infra", "Supra", "PecMajorC", "PecMajorS", ...
                     "Lats", "Sert", "TrapLow", "TrapMid", "TrapUp"];

% check if the muscle is present in the list of expected ones
aux = strcmp(muscle, muscles_in_order);
if max(aux)==0
    error("The specified muscle is not present in this file")
end

%% Find MVCs
index_muscle = find(aux==1);

% read as a table the data from .csv file, and retain the EMG part
csv_content = table2array(readtable(csv_file_Clark, 'VariableNamingRule', 'preserve', 'VariableNamesLine', 4));
emg_data = csv_content(2:end, 3:16);

% select the colum corresponding to the requested muscle, and return its
% maximum (absolute) value as MVC
required_muscle = emg_data(:, index_muscle);

%% 
% y = fft(required_muscle);
% n = length(required_muscle);          % number of samples
% f = (0:n-1)*(fs/n);     % frequency range
% power = abs(y).^2/n;    % power of the DFT
% 
% figure
% subplot(2,3,1)
% plot(f,power)
% xlabel('Frequency')
% ylabel('Power')
% grid on
% title("raw EMG frequency component")
%%

% find and remove the DC bias from the signal
bias_DC = mean(required_muscle,1);
required_muscle = required_muscle - bias_DC;

%%

% y = fft(required_muscle);
% n = length(required_muscle);          % number of samples
% f = (0:n-1)*(fs/n);     % frequency range
% power = abs(y).^2/n;    % power of the DFT
% 
% subplot(2,3,2)
% plot(f,power)
% xlabel('Frequency')
% ylabel('Power')
% grid on
% title("debiased EMG frequency component")

%% 

% highpass Butterworth filter of 4th order (cutoff_high) is applied 
nyquist_rate = fs/2;
[b_h,a_h]=butter(4, cutoff_high/nyquist_rate, 'high');
required_muscle = filter(b_h, a_h, required_muscle);

%%
% y = fft(required_muscle);
% n = length(required_muscle);          % number of samples
% f = (0:n-1)*(fs/n);     % frequency range
% power = abs(y).^2/n;    % power of the DFT
% 
% subplot(2,3,3)
% plot(f,power)
% xlabel('Frequency')
% ylabel('Power')
% grid on
% title("highpass EMG frequency component")

%%

% rectification of the emg data
required_muscle = abs(required_muscle);

%%
% y = fft(required_muscle);
% n = length(required_muscle);          % number of samples
% f = (0:n-1)*(fs/n);     % frequency range
% power = abs(y).^2/n;    % power of the DFT
% 
% subplot(2,3,4)
% plot(f,power)
% xlabel('Frequency')
% ylabel('Power')
% grid on
% title("rectified EMG frequency component")
% 
%%

% lowpass Butterworth filter of 4th order (cutoff_low) is applied
[b_l,a_l]=butter(4, cutoff_low/nyquist_rate, 'low');
required_muscle = filter(b_l, a_l, required_muscle);

%%
% y = fft(required_muscle);
% n = length(required_muscle);          % number of samples
% f = (0:n-1)*(fs/n);     % frequency range
% power = abs(y).^2/n;    % power of the DFT
% 
% subplot(2,3,5)
% plot(f,power)
% xlabel('Frequency')
% ylabel('Power')
% grid on
% title("lowpass EMG frequency component")
%%

% the signal is averaged with a moving average of 500 ms. Since the signal
% is acquired at 1500 Hz, this means averaging with a window of 750 samples
window_samples = 750;
required_muscle_filt = movmean(required_muscle, window_samples);

%%
% y = fft(required_muscle_filt);
% n = length(required_muscle_filt);          % number of samples
% f = (0:n-1)*(fs/n);     % frequency range
% power = abs(y).^2/n;    % power of the DFT
% 
% subplot(2,3,6)
% plot(f,power)
% xlabel('Frequency')
% ylabel('Power')
% grid on
% title("movmean EMG frequency component")
%%

value_MVC = max(required_muscle_filt);
