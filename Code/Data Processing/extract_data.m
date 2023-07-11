% author: Irene Beck I.L.Y.Beck@student.tudelft.nl
%% FUNCTION: Extract data from Qualysis .mat outputs. Filters to only take the EMG channels.
% Filters EMG data using function 'butterworth' to only have smooth EMG signal
% WARNING: Function is specific for Qualysis output in combination with Delsys (10 channels per EMG sensor)
% May not work if used on different formats (CSV or C3D) or different amount of channels on EMG sensor
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
        EMG_filtered_Struct.(analog_string(k)) = butterworth_filter(EMG_Matrix_Struct.(analog_string(k)), size(EMG_Matrix_Struct.(analog_string(k)),1), samplerate);   

        % Moving average window
        EMG_movavg_Struct.(analog_string(k)) = movmean(EMG_filtered_Struct.(analog_string(k)), window_frames, 2);
        
        % Create time vector
        time(k) = length(EMG_filtered_Struct.(analog_string(k)))/samplerate;
        EMG_filtered_Struct.(time_string(k)) = linspace(0, time(k), length(EMG_filtered_Struct.(analog_string(k))));
        EMG_movavg_Struct.(time_string(k)) = linspace(0, time(k), length(EMG_movavg_Struct.(analog_string(k))));
    end

end