% just a script checking that the results of the RMR analysis contained in
% the selected folder make sense.
% The amount of times in which the RMR solver did not converge is returned.
% It is also possible to plot the muscle activations so that we can perform
% a visual check if they make sense.
clear; clc; close all
%% Parameters
check_RMR_unfeasibilities = true;  % if true, the unfeasibility_flags field of each .mat file is checked
                                   % if any unfeasible solution is found,
                                   % the name of the corresponding file is
                                   % returned

plot_activations = true;           % if true, the activations are visualized

num_muscles = 33;                  % number of muscles in the results

save_plots = false;                 % whether the plots will be saved or not
                                   % If false, they are displayed to the
                                   % user instead.

num_plots = 50;                    % just num_plots results are shown at a time
                                   % the user can inspect them, close them
                                   % and then the remaining will be plotted
                                   % too. Do not increase the number too
                                   % much to avoid overloading the RAM.
                                   % Meaningful only if save_lots=false.


%% Script
% set the path current folder to be the one where this script is contained
mfile_name          = mfilename('fullpath');
[pathstr,~,~]  = fileparts(mfile_name);
cd(pathstr)

% getting path to other folders in this repo
addpath(pathstr)
cd ..\..\..\
path_to_repo = pwd;
addpath(path_to_repo)

% Select files where RMR results have been stored
[name_files, results_path] = uigetfile('*.mat', 'Select the RMR result files to consider', path_to_repo, 'MultiSelect','on');

if save_plots
    % get the directory to save images to, if required
    saving_path=uigetdir(results_path, 'Select folder to save muscle activation plots to');
end

if iscell(name_files)
    num_files = size(name_files, 2);
else
    num_files = 1;
    name_files = {name_files};
end

% check that all the activations contained in the result_files have been
% found properly, without any unfeasibility reported by RMR
num_unfeasible = 0;
unreliable_exp = [];

if check_RMR_unfeasibilities
    for index_file=1:num_files
        curr_file = load(fullfile(results_path, name_files{index_file}));
        if ~isfield(curr_file, 'unfeasibility_flags')
            warning("%s has no field named unfeasibility_flags -> discarded from analysis", name_files{index_file})
        else
            unfeasibilities = sum(curr_file.unfeasibility_flags);
            num_unfeasible = num_unfeasible+unfeasibilities;
            if unfeasibilities
                unreliable_exp = [unreliable_exp, string(name_files{index_file})];
            end
        end
    end
end

fprintf("The total number of unreliable experiment results among the ones you selected is: %i \n\n", num_unfeasible)

% prints the file names of the unreliable experiments that have been
% selected, so the user can doublecheck those more accurately
for index_unreliable_exp = 1:size(unreliable_exp,2)
    fprintf("%s \n", unreliable_exp(index_unreliable_exp))
end

% now the resulting muscle activations are plotted, num_plots at a time.
% After that the user can close all the plots themselves to get the next
% ones to be visualized
if plot_activations

    if ~save_plots
        % initialize fields of dialog box for printing figures
        group = "Updates";
        pref = "Conversion";
        title_box = "Plotting is paused";
        quest = ["Inspect the plots that have been generated.", ...
                 "Once you are done, close all the plots and press 'Continue' to display more results"];
    end

    % loop through the files containing RMR results
    for index_file = 1:num_files
        curr_file = load(fullfile(results_path, name_files{index_file}));

        % set visibility of the figure (to be saved or to be shown now)
        if save_plots
            f = figure('visible', 'off', 'Position', get(0, 'Screensize'));
        else 
            f = figure();
            f.WindowState = 'maximized';
        end

        muscle_names = curr_file.muscle_order;
        pgc = linspace(0, 100, size(curr_file.xsol,1));
        for index_muscle = 1:num_muscles
           subplot(5,8,index_muscle)
           hold on
           plot(pgc,curr_file.xsol(:,index_muscle),'b-')
           ylim([0 1])
           title(muscle_names(index_muscle))
           hold off
        end
        sgtitle("Muscle Activations")
        
        % save the plots in the saving_path, or pause displaying too many
        % plots at a time
        if save_plots
            saveas(f, append(fullfile(saving_path, name_files{index_file}(1:end-4)), '.png'))
        else
            if mod(index_file, num_plots)==0 && index_file~=0
                uigetpref(group,pref,title_box,quest,"Continue");
            end
        end
    end
end
