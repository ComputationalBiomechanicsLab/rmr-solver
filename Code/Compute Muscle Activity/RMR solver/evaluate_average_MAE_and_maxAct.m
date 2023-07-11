function evaluate_average_MAE_and_maxAct(folderPath)
% This function is used to process results from the SO. In particular,
% taking as an inpute a path, it returns the average values of the maximum
% activation of each muscle over the tasks present as subfolders in the
% folder. Also, the average MAE is returned (taken from comparison with
% ground truth EMG data)

tasks = ["abd01", "abd02", "abd03", "abd21", "abd22", "abd23", ...
         "flx01", "flx02", "flx03", "flx21", "flx22", "flx23", ...
         "shrug01", "shrug02", "shrug03", "shrug21", "shrug22", "shrug23"];

% inizialize the number of tasks that have been effectively evaluated
num_tasks = 0;

MAE = zeros(1,11);
correlation = zeros(1,11);
maxActivations = zeros(1,33);

% switching to the right folder to analyze
cd(folderPath)

for i=1:size(tasks,2)
    task = tasks(i);
    if isfolder(task)
        subfolder = append(folderPath, '\', task);
        cd(subfolder)

        % load and add the max activations 
        file_maxAct = append(task, '_maxActivation');
        maxActFromFile = load(file_maxAct);
        maxActivations = maxActivations + maxActFromFile.metrics.maxActivations;

        % load add the MAE
        file_mae = append(task, '_metrics_correctMuscles');
        maeFromFile = load(file_mae);
        MAE = MAE + maeFromFile.metrics.mae_SO;

        % load and add the correlation
        correlation = correlation + maeFromFile.metrics.corr_SO;
        
        num_tasks = num_tasks + 1;
        cd(folderPath)
    end
end

avg_maxActivations = maxActivations/num_tasks;
avg_mae = MAE/num_tasks;
avg_corr = correlation/num_tasks;

% visualize the averaged max activations
figure
h_max = heatmap(avg_maxActivations', 'ColorLimits', [0,1], 'Colormap', jet);
colorbar
set(gcf,'position',[0,0,350,650])
title('averaged max activations')

ax = gca;
ax.YDisplayLabels = maxActFromFile.metrics.muscleNames(1:33);


% visualize the averaged MAE
figure
h_mae = heatmap(avg_mae', 'ColorLimits', [0,0.4], 'Colormap', jet(4));
colorbar
ax = gca;
ax.YDisplayLabels = maeFromFile.metrics.muscleNames;
set(gcf,'position',[0,0,300,650])
title('averaged MAE')

metrics.avg_max_act = avg_maxActivations;
metrics.avg_mae = avg_mae;
save('metrics', 'metrics')

nameFig =  append(task{1}(1:end-1),  '_avg_max_act.png');
saveas(h_max, nameFig)

nameFig =  append(task{1}(1:end-1),  '_avg_mae.png');
saveas(h_mae, nameFig)
