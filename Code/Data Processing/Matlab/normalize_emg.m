% author: Irene Beck I.L.Y.Beck@student.tudelft.nl
%% FUNCTION: Normalize EMG data
% Normalizes EMG data by taking it as a part of the MVC (maximum voluntary
% contraction) -> scales between 0 and 1.
function [EGM_normalized_Struct] = normalize_emg(MVC_movavg, EMG_filtered_Struct, analog_string)
for m = 1:length(analog_string)
    EGM_normalized_Struct.(analog_string(m)) = (EMG_filtered_Struct.(analog_string(m))) ./ MVC_movavg;
end
end