% author: Irene Beck (I.L.Y.Beck@student.tudelft.nl)

clear all
close all
clc

% Note: EMG signals measured in uV. 

% Define directory path with data (.mat files) and muscles measured during
% the experiment
dir_path = 'U:\staff-umbrella\PTbot\Experiments\2021-10-12_MicahTest2\Data\AnalogData\';

% Middle Deltoid, Anterior Deltoid, Posterior Deltoid, Pectoralis Major
% Clavicular, Pectoralis Major Sternocostal, Serratus Anterior, LongHead of
% Biceps, Upper Trapezius, Middle Trapezius, Lower Trapezius,
% Infraspinatus, LongHead of Triceps, Latissimus Dorsi

%!! NB: They have to be stated in the same order as the applied EMG sensors
muscles = ["MD", "AD","PD","PM CL","PM ST","SA","LHB","UT","MT","LT","ISP","LHT","LAT"];

% Processing parameters
samplerate = 2000;
window = 300; %ms. Used for moving average window filtering
window_frames = window*samplerate/1000; % Amount of frames needed for moving average window filtering
%% Get MVC values
% Load all files in directory path starting with 'MVC' and with a .mat
% format
MVC_search_string = fullfile(dir_path, 'MVC*.mat'); % Search query
MVC_Files = dir(MVC_search_string); % List files in directory

% Create arrays with field names to be later used in structure
mvc_string = 'mvc_' + string(linspace(1,length(MVC_Files),length(MVC_Files))); % string array with field names for MVC
emg_string = 'emg_channel_' + string(linspace(1,length(muscles),length(muscles))); % string array with field names for channels per MVC
mat_string = 'MVC_' + string(linspace(1, length(MVC_Files), length(MVC_Files))) + '_Standing'; % string array with .mat names
time_string = 'Time_' + string(linspace(1,length(MVC_Files),length(MVC_Files))); % string array with field names for time vectors

% Extract EMG data from analog data and determine maximum value per muscle
[MVC_EMG_Struct, MVC_filtered_Struct, MVC_movavg_Struct] = extract_data(MVC_Files, mvc_string, emg_string, mat_string, time_string, samplerate, window_frames);
[MVC, MVC_movavg] = extract_mvc(MVC_filtered_Struct, MVC_movavg_Struct, MVC_Files, mvc_string, emg_string);

%% Get EMG data for motions and normalize
% Load all files in directory path starting with 'motion' and with a .mat
% format
Motion_search_string = fullfile(dir_path, 'motion*.mat');
Motion_Files = dir(Motion_search_string);

% Create arrays with field names to be later used in structure
motion_string = 'Motion_' + string(linspace(1,length(Motion_Files),length(Motion_Files))); % string array with field names for motion = analog_string
emg_string = 'emg_channel_' + string(linspace(1,length(muscles),length(muscles))); % string array with field names for channels per motion
mat_string = 'motion_' + string(linspace(1, length(Motion_Files), length(Motion_Files))); % string array with .mat names
time_string = 'Time_' + string(linspace(1,length(Motion_Files),length(Motion_Files))); % string array with field names for time vectors

% Extract EMG data from analog data and filter + normalize
[Motion_EMG_Struct, Motion_filtered_Struct, Motion_movavg_Struct] = extract_data(Motion_Files, motion_string, emg_string, mat_string, time_string, samplerate, window_frames);
[Motion_max, Motion_movavg_max] = extract_mvc(Motion_filtered_Struct, Motion_movavg_Struct, Motion_Files, motion_string, emg_string);


% Determine maximum muscle activation based on both the MVC and motion
% experiments. Use the combined maximum value to normalize motion data
MVC_Motion_movavg_max = [MVC_movavg, Motion_movavg_max];
Total_movavg_max = max(MVC_Motion_movavg_max, [], 2);
Motion_normalized_Struct = normalize_emg(Total_movavg_max, Motion_movavg_Struct, motion_string);
%% Plotting figures HARDCODED
for i = 1:13
%     % MVC Raw
%     figure(1)
%     hold on
%     subplot(7,2,i)
%     plot(MVC_filtered_Struct.Time_1, MVC_EMG_Struct.mvc_1.Data.(emg_string(i)))  
%     title(string(muscles(i)))
%     sgtitle('MVC: Raw EMG signal [uV]')
%     xlabel('time [s]')    
%     
%     % MVC filtered (2x Butterworth) and with average window
%     figure(2)
%     subplot(7,2,i)
%     plot(MVC_filtered_Struct.Time_1, MVC_filtered_Struct.mvc_1(i,:))
%     hold on
%     plot(MVC_filtered_Struct.Time_1, MVC_movavg_Struct.mvc_1(i,:)) 
%     title(string(muscles(i)))
%     sgtitle('MVC: Filtered EMG signal [uV]')
%     xlabel('time [s]')
% 
%     % Normalized EMG of motions
%     figure(3)
%     hold on
%     subplot(7,2,i)
%     plot(MVC_filtered_Struct.Time_1, MVC_normalized_Struct.mvc_1(i,:))
%     title(string(muscles(i)))
%     sgtitle('normalized EMG data [-]')
%     xlabel('time [s]')

    % Motion RAW
    figure(4)
    hold on
    subplot(7,2,i)
    plot(Motion_filtered_Struct.Time_3, Motion_EMG_Struct.Motion_3.Data.(emg_string(i)))
    title(string(muscles(i))) 
    sgtitle('Raw EMG signal [uV]')
    xlabel('time [s]')    

    % Motion filtered (2x Butterworth) and with average window
    figure(5)
    subplot(7,2,i)
    plot(Motion_filtered_Struct.Time_3, Motion_filtered_Struct.Motion_3(i,:))
    hold on
    plot(Motion_movavg_Struct.Time_3, Motion_movavg_Struct.Motion_3(i,:)) 
    title(string(muscles(i)))
    sgtitle('Filtered EMG signal [uV]')
    xlabel('time [s]')
    
    % Normalized EMG of motions
    figure(6)
    hold on
    subplot(7,2,i)
    plot(Motion_filtered_Struct.Time_3, Motion_normalized_Struct.Motion_3(i,:))
    ylim([0, 1])
    title(string(muscles(i)))
    sgtitle('normalized EMG data [-]')
    xlabel('time [s]')
end

%% FUNCTIONS %%

%% extract_data(): Extract data + filter to only have smooth EMG signal
function [EMG_Struct, EMG_filtered_Struct, EMG_movavg_Struct] = extract_data(list_of_files, analog_string, emg_string, mat_string, time_string, samplerate, window_frames)
    time = zeros(length(list_of_files),1);
    for k = 1: length(list_of_files)
        % Find all files matching search query in specified directory [1]
        baseFileName = list_of_files(k).name;
        fullFileName = fullfile(list_of_files(k).folder, baseFileName);        

        % Save .mat data in structure for all MVC data
        All_Analog_Struct.(analog_string(k)) = load(fullFileName);

        % Filter data to only have EMG channels (1, 11, 21, etc)
        for i = 1:(length(emg_string))
            EMG_Struct.(analog_string(k)).Data.(emg_string(i)) = All_Analog_Struct.(analog_string(k)).(mat_string(k)).Analog.Data(10*(i-1) + 1,:);  
        end
        
        % Combine EMG Fields in to single matrix. Filter Matrix using
        % Butterworth filters
        EMG_Matrix_Struct.(analog_string(k)) = cell2mat(struct2cell(EMG_Struct.(analog_string(k)).Data));
        EMG_filtered_Struct.(analog_string(k)) = butterworth(EMG_Matrix_Struct.(analog_string(k)), size(EMG_Matrix_Struct.(analog_string(k)),1), samplerate);   

        % Moving average window
        EMG_movavg_Struct.(analog_string(k)) = movmean(EMG_filtered_Struct.(analog_string(k)), window_frames, 2);
        
        % Create time vector
        time(k) = length(EMG_filtered_Struct.(analog_string(k)))/samplerate;
        EMG_filtered_Struct.(time_string(k)) = linspace(0, time(k), length(EMG_filtered_Struct.(analog_string(k))));
        EMG_movavg_Struct.(time_string(k)) = linspace(0, time(k), length(EMG_movavg_Struct.(analog_string(k))));
    end

end

%% extract_mvc(): Extract MVC values
function [MVC, MVC_movavg] = extract_mvc(EMG_filtered_Struct, EMG_movavg_Struct, list_of_files, analog_string, emg_string)
    MVC_filt_max = zeros(length(emg_string), length(analog_string));
    MVC_movavg_max = zeros(length(emg_string), length(analog_string));
    for k = 1: length(list_of_files)

        % Find maximum for each MVC measurement, results in matrix of size
        % numMVCs x numEMGs
        MVC_filt_max(:,k) = max(EMG_filtered_Struct.(analog_string(k)), [], 2);
        MVC_movavg_max(:,k) = max(EMG_movavg_Struct.(analog_string(k)), [], 2);

        % Maximum of maximum. Find resulting MVC for each muscle over all MVC
        % measurements
        MVC = max(MVC_filt_max, [], 2);
        MVC_movavg = max(MVC_movavg_max, [], 2);
    end
end

%% normalize_emg(): Normalize EMG data
function [EGM_normalized_Struct] = normalize_emg(MVC_movavg, EMG_filtered_Struct, analog_string)
    for m = 1:length(analog_string)
        EGM_normalized_Struct.(analog_string(m)) = (EMG_filtered_Struct.(analog_string(m))) ./ MVC_movavg;
    end
end

%% butterworth(): Create butterworth filters as used in Thoracoscapular paper by Ajay
function[EMG_filt] = butterworth(EMG_dat, numEMGs, samplerate)
    % EMG_dat = zeros(numEMGs,length(EMG_dat));
    data_filtered = zeros(numEMGs,length(EMG_dat));
    data_rect = zeros(numEMGs,length(EMG_dat));
    EMG_filt = zeros(numEMGs,length(EMG_dat));

    %4th order bandpass butterworth filter 
    fc = [20 400]; %cut-off frequency
    fs = samplerate; %sample rate
    [b_h,a_h]=butter(4,fc/(fs/2),'bandpass');

    %2nd order lowpass buttworth filter
    fc = 10; %cut-off frequency is 2.5?
    fs = samplerate; %sample rate
    [b_l,a_l]=butter(2,fc/(fs/2));

    %Filter and rectify each of the EMG channels
        for i = 0:numEMGs-1
            data_raw = EMG_dat(i+1,:);
            data_filtered(i+1,:)=filter(b_h,a_h,data_raw);        %highpass filtering
            data_rect(i+1,:)=abs(data_filtered(i+1,:));           %rectify signal
            EMG_filt(i+1,:) = filtfilt(b_l,a_l,data_rect(i+1,:)); %zero phase filtering with lowpass
        end
end