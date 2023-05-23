% author: Irene Beck I.L.Y.Beck@student.tudelft.nl
%% FUNCTION: Extract MVC values
% Takes the maximum of each EMG sensor for all trials that were used as
% input. Then takes maximum per EMG sensor again to have the maximum
% activation per muscle as measured by the EMG sensors
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
