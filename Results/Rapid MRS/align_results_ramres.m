function align_results_ramres()
% This function is used to group and align the predicted activations for
% every muscle in the model considered by RaMReS, and save the result in
% various files (one per muscle, according to the convention
% TASK_NAME_index_muscle.m
%
% author Italo Belli (i.belli@tudelft.nl) 2023

% set the path current folder to be the one where this script is contained
mfile_name          = mfilename('fullpath');
[pathstr,~,~]  = fileparts(mfile_name);
cd(pathstr);

% change it to path to repo
addpath(pathstr)
cd '..\..'
path_to_repo = pwd;
addpath(path_to_repo)
addpath(append(path_to_repo, '\Code\Compute Muscle Activity\Matlab'))

% select the files
[files,path] = uigetfile('*.mat', 'Select the files to analyse (3 files per task)', pathstr, 'MultiSelect','on');

% indicate the task considered
task = 'ABD_0kg';

% select if the glenohumeral constraint was considered
answer = questdlg('Is the GH constrait included?', ...
	'Select an option', ...
	'Yes','No', 'No');

switch answer
    case 'Yes'
        disp('Gh constraint is considered')
        GH_enforced = 'withGH';
    case 'No'
        disp('GH constraint not considered')
        GH_enforced = 'noGH';
end

for index_muscle=1:33
    evaluate_mean_std_RaMReS(path, index_muscle, task, GH_enforced);
end
