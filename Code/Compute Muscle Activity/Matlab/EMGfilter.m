function datafit=EMGfilter(data_raw,a_h,b_h,a_l,b_l,muscle_CMC)
data_filtered=filter(b_h,a_h,data_raw);
data_rect=abs(data_filtered);
datafit = filtfilt(b_l,a_l,data_rect);
if muscle_CMC==67||muscle_CMC==57
 datafit = datafit-min(datafit);
end