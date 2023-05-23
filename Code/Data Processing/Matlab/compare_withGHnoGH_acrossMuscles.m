% This script is used to compute the root mean square error between the	
% muscle activations obtained when enforcing the GH constraint and those	
% obtained without.
clear	

% set the path current folder to be the one where this script is contained	
mfile_name          = mfilename('fullpath');	
[pathstr,~,~]  = fileparts(mfile_name);	
cd(pathstr);	

addpath(pathstr)	
cd '..\..\..'	
path_to_repo = pwd;	
addpath(path_to_repo)	

folder_to_consider = 'C:\Users\italobelli\Desktop\GitHub\PTbot_officialCBL\Personal_Results\noise_on_elbow\ABD_2kg\aligned_activations\3samplesPerExp';	

withGH_folder = fullfile(folder_to_consider, 'withGH');	
noGH_folder = fullfile(folder_to_consider, 'noGH');	

% select and load files obtained enforcing GH constraint	
[files_GH,path_GH] = uigetfile('*.mat', 'Select the .mat files to analyse (withGH)', withGH_folder, 'MultiSelect','on');	

if iscell(files_GH)	
    num_files_GH = size(files_GH, 2);	
else	
    num_files_GH = 1;	
end	

for i=1:num_files_GH	
    aux_load = load(append(path_GH, files_GH{i}));	
    withGH_data(:,:,i) = aux_load.aligned_data;	
end	

% select and load files obtained without including GH constraint	
[files_noGH,path_noGH] = uigetfile('*.mat', 'Select the .mat files to analyse (noGH)', noGH_folder, 'MultiSelect','on');	

if iscell(files_noGH)	
    num_files_noGH = size(files_noGH, 2);	
else	
    num_files_noGH = 1;	
end	

for i=1:num_files_noGH	
    aux_load = load(append(path_noGH, files_noGH{i}));	
    noGH_data(:,:,i) = aux_load.aligned_data;	
end	

% now check if the dimensions are the same	
if ~(size(withGH_data, 3)==size(noGH_data, 3))	
    error("Sizes must agree! Different number of muscles considered")	
end	
if ~(size(withGH_data, 1)==size(noGH_data, 1))	
    error("Sizes must agree! Different number of repetitions considered")	
end	

% now that we are sure that the dimensions match we can just use a more
% general term
num_file = num_files_noGH;
num_muscles = size(noGH_data, 3);	
exp_size = size(noGH_data,2);
num_reps = size(noGH_data, 1);


% compute the individual RMSE and then their average, per muscle	
rmse_total = zeros(num_muscles, 1);	
rmse_raw = zeros(num_muscles, num_reps);

for index_muscle=1:num_muscles	
    for index_reps = 1:num_reps
        rmse_raw(index_muscle, index_reps) = sqrt(sum((noGH_data(index_reps,:,index_muscle) - withGH_data(index_reps,:,index_muscle)).^2)/exp_size);	
    end
end	

for index_muscle = 1:num_muscles
    rmse_total(index_muscle) = mean(rmse_raw(index_muscle, :));
end	

% now compute the MAE value corresponding to each muscle
mae_raw = zeros(num_muscles, num_reps);
mae_total = zeros(num_muscles, 1);	


for index_muscle=1:num_muscles	
    for index_reps = 1:num_reps
        mae_raw(index_muscle, index_reps) = sum(abs(noGH_data(index_reps,:,index_muscle) - withGH_data(index_reps,:,index_muscle)))/exp_size;	
    end
end	

for index_muscle = 1:num_muscles
    mae_total(index_muscle) = mean(mae_raw(index_muscle, :));
end	

% now compute the RMS value corresponding to each muscle	
rms_withGH = zeros(num_muscles,num_reps);	
rms_noGH = zeros(num_muscles,num_reps);	
	

for index_muscle=1:num_muscles	
    for index_reps = 1:num_reps
        rms_withGH(index_muscle, index_reps) = rms(withGH_data(index_reps,:,index_muscle));	
        rms_noGH(index_muscle, index_reps) = rms(noGH_data(index_reps,:,index_muscle));
    end
end	

[t_test_h_results, t_test_p_values]=ttest(rms_withGH, rms_noGH, 'Alpha', 0.05, 'Tail', 'both', 'Dim', 2);	

% % now compute the cross correlation between the two conditions	
% cross_corr_total = zeros(num_muscles, 1);	
% 
% for index_muscle=1:num_muscles	
%     cross_corr(index_muscle,1) = xcorr(withGH_data(1,:,index_muscle), noGH_data(1,:, index_muscle), 0, "normalized");	
%     cross_corr(index_muscle, 2) = xcorr(withGH_data(2,:,index_muscle), noGH_data(2,:, index_muscle), 0, "normalized");	
%     cross_corr(index_muscle, 3) = xcorr(withGH_data(3,:,index_muscle), noGH_data(3,:, index_muscle), 0, "normalized");	
%     cross_corr_total(index_muscle) = (cross_corr(index_muscle, 1)+cross_corr(index_muscle, 2)+cross_corr(index_muscle, 3))/3;	
% end	
