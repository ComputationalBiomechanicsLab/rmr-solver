% quick script to average the activations for serratus anterior in RMR
% results

[files, path_folder] = uigetfile('*.mat', 'MultiSelect','on');

if iscell(files)
    num_files = size(files, 2);
else
    num_files = 1;
    files = {files};
end

for index_file = 1:num_files
    file_name = string(files{index_file});
    curr_file = load(fullfile(path_folder, files{index_file}));
    if ~isfield(curr_file, 'xsol')
        warning("%s has no field named xsol -> will be ignored", name_files{index_file})
    elseif ~isfield(curr_file, 'muscle_order')
        warning("%s has no field named muscle_order -> will be ignored", name_files{index_file})
    elseif ~isfield(curr_file, 'frequency_solution')
        warning("%s has no field named frequency_solution -> will be ignored", name_files{index_file})
    else
        xsol_old = curr_file.xsol;
        serratus_avg = (xsol_old(:,5)+xsol_old(:,6)+xsol_old(:,7))/3; 
        xsol_new = [xsol_old, serratus_avg];
        xsol = xsol_new;

        muscle_order = curr_file.muscle_order;
        frequency_solution = curr_file.frequency_solution;
        save(fullfile(path_folder, file_name), 'xsol', 'muscle_order', 'frequency_solution');
    end
end
