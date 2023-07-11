% This script performs a batch CMC analysis of the various experiments in
% the dataset (specified with taskNames). A more synthetic analysis can be
% achieved by specifying only the tasks of interest in the list.
%
% Adapted from original script released with the paper:
%
%     "Muscle Contributions to Upper-Extremity Movement and Work From a Musculoskeletal Modelof the Human Shoulder, Seth et al. 2019"
%
% (available at https://simtk.org/projects/thoracoscapular)
%
% Modified 2023 by Italo Belli (i.belli@tudelft.nl)

import org.opensim.modeling.*;

% uncomment for older versions of OpenSim
% Model.LoadOpenSimLibrary('..\..\ScapulothoracicJointPlugins40\WinX64\ScapulothoracicJointPlugin40_WinX64');

% set the path current folder to be the one where this script is contained
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

% getting path to other folders in this repo
addpath(pathstr)
cd ..\..\..\
path_to_repo = pwd;
addpath(path_to_repo)
addpath(fullfile(path_to_repo, 'Code\Data Processing\'))
cd(pathstr);

% where to save the results of CMC analysis
saving_path = fullfile(path_to_repo, '\Personal_Results\');

% specify task names to run CMC for (can be a subset of the following)
taskNames={'abd01';'abd02';'abd03';'flx01';'flx02';'flx03';'shrug01'; 'shrug02';'shrug03'; ...
           'abd21';'abd22';'abd23';'flx21';'flx22';'flx23';'shrug21';'shrug22';'shrug23'};

% load the two models, where the hand-mass and irnetia have been modified
% The models are specific for CMC as some coordinates have been locked for this analysis.
modelBase = Model(fullfile(path_to_repo, 'OpenSim Models/for CMC/TSM_subject_CMC_noWeight.osim'));
model2kgBase = Model(fullfile(path_to_repo, 'OpenSim Models/for CMC/TSM_subject_CMC_2kgWeight.osim'));

cmcBase = CMCTool('CMC_setup.xml', false);

% create task set for shrugging cases where we'll ignore plane_elv and
% axial_rot coordinates that are poorly defined when humerus is vertical.
ts = CMC_TaskSet(cmcBase.getTaskSetFileName());
tpe = ts.get('plane_elv');
ts.remove(tpe);
tar = ts.get('axial_rot');
ts.remove(tar);
tempShrugTaskFileName = 'temp_shrug_CMC_Tracking_Tasks.xml';
ts.print(tempShrugTaskFileName);

currentDir = pwd();

for i=1:length(taskNames)
    % start with a fresh Tool for each trial
    cmc = CMCTool('CMC_setup.xml', false);
    % give each Tool its own copy (clone) of the necessary Model
    if i<10
        model = modelBase.clone();
    else
        model = model2kgBase.clone();
    end
    % CMC that requires reserve actuators to be appended first
    cmc.updateModelForces(model, 'CMC_setup.xml');
    cmc.setModel(model);

    % do not track plane_elv and axial_rot during shrugging, locking them
    if(strmatch('shrug', taskNames{i}))
        m = cmc.getModel();
        m.getCoordinateSet().get('plane_elv').set_default_value(-0.433725);
        m.getCoordinateSet().get('plane_elv').set_locked(true);
        m.getCoordinateSet().get('axial_rot').set_default_value(0.8125346);
        m.getCoordinateSet().get('axial_rot').set_locked(true);
        cmc.setTaskSetFileName(tempShrugTaskFileName);
    end 
    
    motion=[fullfile(path_to_repo, 'Results\IK solutions\',taskNames{i}),'.mot'];
    cmc.setDesiredKinematicsFileName(motion);
    cmc.setStartTime(0);
    cmc.setFinalTime(15);
    output=[saving_path, taskNames{i}];
    cmc.setResultsDir(output);
    setupfile=['CMC_setup_', taskNames{i}, '.xml'];
    cmc.print(setupfile);
    tic;
    success(i)=cmc.run();       % store whether the CMC analysis concluded successfully (1) or not (0)
    time(i)=toc;                % store the time required for each analysis
    
    [collabels, data] = readStoFile(motion);
    duration = data(end,1);
    
    disp(['CMC for ' taskNames{i} ' took ' num2str(time(i)) 's for ' num2str(duration) 's.']);
    disp(['A compution to real-time ratio of: ' num2str(time(i)/duration)]);
    
    if success(i)==0
        continue;
    else
        cmc.delete();
        model.delete();
    end
end    
