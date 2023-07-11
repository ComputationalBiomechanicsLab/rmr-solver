% This script allows to augment the motion data which we have by
% considering only the elbow markers as starting point for injecting noise
% in the movement and generating motions which are different yet similar to
% the ones in our dataset.

original_trc_file = 'SHRUG23.trc';
original_path = 'C:\Users\italobelli\Desktop\GitHub\PTbot_officialCBL\ExperimentalData\Markers';
marker_names = ["centelbow", "EpL", "centpxt8", "centijc7", "ij"];
saving_path = 'C:\Users\italobelli\Desktop\GitHub\PTbot_officialCBL\Personal_Results';

% extract only the required markers from the original file
% new_trc = extract_markers_by_name(original_trc_file, original_path, marker_names, saving_path);
new_trc = original_trc_file;

% downsample them to a reasonable frequency
downsampling_freq = 50;
[~, ~, ~, downsampled_file_name] = get_downsampled_trc_data(new_trc, ...
                                                            saving_path, ...
                                                            downsampling_freq, ...
                                                            true);

% add white noise to the markers on the elbow
cutoff_freq = 0.5; 
markers_to_be_modified = ["centelbow", "EpL"];
mean_wn = 0;
sigma_wn = 0.01;

% loop to generate multiple instances of the noisy files
num_files = 2;
noisy_file_name = cell(num_files, 1);

for index_file = 1:num_files
    file_name = [extractBefore(original_trc_file, '.trc'), '_noisy', num2str(index_file), '.trc'];
    [markers_noisy, noisy_file_name{index_file}, max_deviation]  = add_equal_white_noise_to_marker(downsampled_file_name, ...
                                                                                                   saving_path, ...
                                                                                                   markers_to_be_modified, ...
                                                                                                   mean_wn, ...
                                                                                                   sigma_wn, ...
                                                                                                   cutoff_freq, ...
                                                                                                   true, ...
                                                                                                   file_name);
    % if the level of injected noise is too small or too big, consider the
    % same file again and augment it with a different level of noise
    if max_deviation<0.01 || max_deviation > 0.06
        index_file = index_file - 1;
    end
end

%%
% up-sample the noisy markers to original frequency
upsampling_freq = downsampling_freq;
interpolating_method = 2; 

for index_file=1:num_files
    get_upsampled_trc_data(noisy_file_name{index_file}, ...
                           saving_path, ...
                           upsampling_freq, ...
                           interpolating_method, ...
                           true, ...
                           noisy_file_name{index_file});
end

% TODO: process the multiple noisy files directly with an ad-hoc IK routine (where the
% scapula movement needs to be heavily penanlized). After that, run RMR
% solver on the resulting motion!