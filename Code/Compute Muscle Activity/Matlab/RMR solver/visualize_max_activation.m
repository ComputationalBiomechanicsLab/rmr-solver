function visualize_max_activation(file_SO_results)
% This function plots the activation levels resulting from Static Optimization (SO)
% in a single table, using a colour map to quickly capture the level of (maximum) activation 
% that each muscle is reaching during the movement which was used to
% generate the results.
% It requires as input 
% * file_SO_results: filepath+filename of the results of a run of SO, as
%                    returned by the CSO_accelerationConst*.m scripts
% It prints a summarizing table for the data, to help their analysis
%
% author: Italo Belli (i.belli@tudelft.nl), 2022

% set names for figure and results
nameExp = extractBefore(file_SO_results, '.');
nameExp = extractAfter(nameExp, 'experiment_');
nameFig =  append(nameExp, '_maxActivation', '.png');
nameFile = append(nameExp, '_maxActivation');

% load data and save the required fields
data = load(file_SO_results);
activations = data.xsol;
muscle_order = data.muscle_order;

% if Thoracoscapular Model is considered
muscleNames = ["Trap Scap Mid", "Trap Scap Sup", "Trap Scap Inf", "Trap Clav Sup", "Serr Ant Inf", "Serr Ant Mid", "Serr Ant Sup", "Rhomb Sup", "Rhomb Inf", ...
               "Levat Scap", "Coracobrac", "Deltoid Ant", "Deltoid Post", "Deltoid Mid", "Lat Dorsi Sup", "Lat DOrsi Mid", "Lat Dorsi Inf", "Pect Major Clav", ...
               "Pect Major Thor Inf", "Pect Major Thor Mid", "Teres Major", "Infra Inf", "Infra Sup", "Pect Minor", "Teres Minor", "Subscap Sup", "Subscap Mid",...
               "Subscap Inf", "Supraspin Post", "Supraspin Ant", "Triceps", "Bic Long", "Bic Brevis"];

% keep just data regarding real muscles, not coordinateActuators
activations_muscles = activations(:,1:33);
level_coordAct = activations(:,34:end);

% get the max activation per muscle
max_activations = max(activations_muscles, [], 1);

figure
h = heatmap(max_activations', 'ColorLimits', [0,1], 'Colormap', jet); %, 'YLabel', muscle_order(1:33))
colorbar
set(gcf,'position',[0,0,350,650])
title(nameExp)

ax = gca;
ax.YDisplayLabels = muscleNames;

saveas(h, nameFig)

metrics.muscleNames = muscle_order;
metrics.maxActivations = max_activations;
metrics.maxCoordActs = level_coordAct;
save(nameFile, 'metrics');