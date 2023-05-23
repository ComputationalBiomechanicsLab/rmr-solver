import org.opensim.modeling.*;
taskNames={'abd01';'abd02';'abd03';'flx01';'flx02';'flx03';'shrug01';'shrug02';'shrug03';...
    'abd21';'abd22';'abd23';'flx21';'flx22';'flx23';'shrug21';'shrug22';'shrug23'};

Model.LoadOpenSimLibrary('..\..\ScapulothoracicJointPlugins40\WinX64\ScapulothoracicJointPlugin40_WinX64');

modelBase = Model('../New_Scaled_FinalTest.osim');
model2kgBase = modelBase.clone();
hand = model2kgBase.getBodySet().get('hand');
% add 2kg to hand body
hand.setMass(hand.getMass()+2.0);

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

for i=1:length(taskNames),
    % start with a fresh Tool for each trial
    cmc = CMCTool('CMC_setup.xml', false);
    % give each Tool its own copy (clone) of the necessary Model
    if i<10
        model = modelBase.clone();
    else
        model = model2kgBase.clone();
    end
    % Bug/odd behavior in CMC that requires reserve actuators to be
    % appended first
    cmc.updateModelForces(model, 'CMC_setup.xml');
    cmc.setModel(model);
    % do not track plane_elv and axial_rot during shrugging. Lock these.
    if(strmatch('shrug', taskNames{i}))
        m = cmc.getModel();
        m.getCoordinateSet().get('plane_elv').set_default_value(-0.433725);
        m.getCoordinateSet().get('plane_elv').set_locked(true);
        m.getCoordinateSet().get('axial_rot').set_default_value(0.8125346);
        m.getCoordinateSet().get('axial_rot').set_locked(true);
        cmc.setTaskSetFileName(tempShrugTaskFileName);
    end 
        
    motion=['..\IK\output\',taskNames{i},'_IK.mot'];
    cmc.setDesiredKinematicsFileName(motion);
    cmc.setStartTime(0);
    cmc.setFinalTime(15);
    output=[currentDir, '/', taskNames{i}];
    cmc.setResultsDir(output);
    setupfile=['CMC_setup_', taskNames{i}, '.xml'];
    cmc.print(setupfile);
    tic;
    success=cmc.run();
    time(i)=toc;
    
    [collabels, data] = readStoFile(motion);
    duration = data(end,1);
    
    disp(['CMC for ' taskNames{i} ' took ' num2str(time(i)) 's for ' num2str(duration) 's.']);
    disp(['A compution to real-time ratio of: ' num2str(time(i)/duration)]);
    
    if success==0
        continue;
    else
        cmc.delete();
        model.delete();
    end
end    
