% author: Irene Beck I.L.Y.Beck@student.tudelft.nl
%% FUNCTION: Create butterworth filters as used in Thoracoscapular paper by Ajay
% Raw EMG data is sampled first trough a bandpass filter, rectified and
% then through a lowpass filter. Resulting in filtered EMG data for further
% data processing.
function[EMG_filt] = butterworth_filter(EMG_dat, numEMGs, samplerate)
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
        data_rect(i+1,:)=abs(data_filtered(i+1,:));          %rectify signal
        EMG_filt(i+1,:) = filtfilt(b_l,a_l,data_rect(i+1,:)); %zero phase filtering with lowpass
    end
%     min_filt = min(min(data_rect))
%     min_emg = min(min(EMG_filt))
end