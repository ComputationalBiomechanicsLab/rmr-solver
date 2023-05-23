%-------------------------------------------------------------------------%
% Copyright (c) 2015 Modenese L., Ceseracciu, E., Reggiani M., Lloyd, D.G.%
%                                                                         %
% Licensed under the Apache License, Version 2.0 (the "License");         %
% you may not use this file except in compliance with the License.        %
% You may obtain a copy of the License at                                 %
% http://www.apache.org/licenses/LICENSE-2.0.                             %
%                                                                         %
% Unless required by applicable law or agreed to in writing, software     %
% distributed under the License is distributed on an "AS IS" BASIS,       %
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or         %
% implied. See the License for the specific language governing            %
% permissions and limitations under the License.                          %
%                                                                         %
%    Author:   Luca Modenese, January 2015                                %
%    email:    l.modenese@sheffield.ac.uk                                 %
% ----------------------------------------------------------------------- %
%
% This function optimizes the muscle parameters as described in Modenese L,
% Ceseracciu E, Reggiani M, Lloyd DG (2015). Estimation of
% musculotendon parameters for scaled and subject specific musculoskeletal
% models using an optimization technique. Journal of Biomechanics (submitted)
% and prints the results to command window.
% Also it stores information about the optimization in the structure SimInfo
%
%% Changes in this file:
%
% - Added a line to load plugin library for ThoracoScapular Shoulder Model
%
% - Added a workaround to properly function in case the reference muscle
% has 0 tendon slack length; for example, serratus muscle-elements in the
% Thoracoscapular Shoulder Model
%
% - Added an input that allows for selection of a single solver as well as
% the native method: 'modenese' (native solver), 'lsq', 'lsqnonneg', and
% 'twostep'.
%
% - Added an input to limit the amount of evaluated poses per muscle, in
% case a high number of degrees of freedom are actuated by the muscle,
% causing an exponential in evaluated poses and computational time.
%
% - Added Optimization time per muscle to SimInfo output struct
%
% - Added automatic saving of optimized model and SimInfo struct
%
% Alterations by Joran de Vet 

function [osimModel_opt, SimInfo] = optimMuscleParamsUpdated(osimModel_ref_filepath, osimModel_targ_filepath, N_eval,nPosesLimit, log_folder,solver)

% import opensim libraries
import org.opensim.modeling.*
% Load scapulothoracic joint library
Model.LoadOpenSimLibrary('C:\OpenSim 4.1\plugins\ScapulothoracicJointPlugin40_WinX64.dll');

% results file identifier
res_file_id_exp = ['_N',num2str(N_eval),'_nLimit',num2str(nPosesLimit)];

% import models
osimModel_ref   = Model(osimModel_ref_filepath);
osimModel_targ  = Model(osimModel_targ_filepath);

% models details
[~, name, ext]   = fileparts(osimModel_targ_filepath);

% assigning new name to the model
osimModel_opt_name  = [name,'_opt',res_file_id_exp,ext];
osimModel_targ.setName(osimModel_opt_name);

% set current directory to log_folder
cd(log_folder)

% initializing log file
log_filepath = fullfile(log_folder,[name,'_opt',res_file_id_exp,'.log']);
% cleaning file (otherwise it appends)
fopen(log_filepath,'w+'); fclose all;
% set up diary, without writing at the moment
diary(log_filepath); diary off

% get muscles
muscles = osimModel_ref.getMuscles;
muscles_scaled = osimModel_targ.getMuscles;

% init model
% si = osimModel_ref.initSystem;
% si_scaled = osimModel_targ.initSystem;

% initialize with recognizable values
LmOptLts_opt = ones(muscles.getSize,2)*(-1000);

for n_mus = 0:muscles.getSize-1
    tic

    % current muscle name (here so that it is possible to choose a single
    % muscle when developing.
    curr_mus_name = char(muscles.get(n_mus).getName);
    display(['processing mus ',num2str(n_mus+1),': ',char(curr_mus_name)]);

    % import muscles
    curr_mus = muscles.get(n_mus);
    curr_mus_scaled = muscles_scaled.get(curr_mus_name);

    % extracting the muscle parameters from reference model
    LmOptLts = [curr_mus.getOptimalFiberLength, curr_mus.getTendonSlackLength];
    PenAngleOpt = curr_mus.getPennationAngleAtOptimalFiberLength();
    Mus_ref = sampleMuscleQuantitiesUpdated(osimModel_ref,curr_mus,'all',N_eval,nPosesLimit);

    % calculating minimum fiber length before having pennation 90 deg
    % acos(0.1) = 1.47 red = 84 degrees, chosen as in OpenSim
    limitPenAngle = acos(0.1);
    % this is the minimum length the fiber can be for geometrical reasons.
    LfibNorm_min = sin(PenAngleOpt)/sin(limitPenAngle);
    % LfibNorm as calculated above can be shorter than the minimum length
    % at which the fiber can generate force (taken to be 0.5 Zajac 1989)
    if (LfibNorm_min<0.5)==1
        LfibNorm_min = 0.5;
    end

    % checking the muscle configuration that do not respect the condition.
    LfibNorm_ref = Mus_ref(:,2);
    okList = (LfibNorm_ref>LfibNorm_min);

    % keeping only acceptable values
    LfibNorm_ref        = LfibNorm_ref(okList);
    LtenNorm_ref        = Mus_ref(okList,3)/LmOptLts(2);
    %Ffib               = Mus_ref(okList,4)/curr_mus.getMaxIsometricForce;
    MTL_ref             = Mus_ref(okList,1);
    penAngle_ref        = Mus_ref(okList,5);
    LfibNormOnTen_ref   = LfibNorm_ref.*cos(penAngle_ref);

    % in the target only MTL is needed for all muscles
    MTL_targ = sampleMuscleQuantitiesUpdated(osimModel_targ,curr_mus_scaled,'MTL',N_eval,nPosesLimit);
    evalTotPoints = length(MTL_targ);
    MTL_targ = MTL_targ(okList)';
    evalOkPoints  = length(MTL_targ);


    
    %% NOTE: Incorporated functionality for cases where Lts = 0 in the reference model
    % This is the case for some muscle-elements in the thoracoscapular
    % shoulder model representing the serratus anterior.
    % Lts = 0 should be maintained on the optimization
    
    if LmOptLts(2) < 1e-6   % check if tendon slack length is zero
        
        % In case of zero tendon slack length, the problem to be solved is:
        % [Lnorm*cos(penAngle)]*[Lmopt] = MTL as tendon length equals zero
        % for every pose. Calculated normalized tendon lengths are NaN in
        % this case, and thus this column is removed from the problem.
        
        A = LfibNormOnTen_ref;
        b = MTL_targ;
        % calculate least-square solution for LmOpt, maintaining zero Lts
        x = [A\b ; 0];
        LmOptLts_opt(n_mus+1,:) = x;
        
        
    else
        
    %% NOTE: Setting added to also use one of the three solver methods exclusively    
    
    % The problem to be solved is:
    % [LmNorm*cos(penAngle) LtNorm]*[Lmopt Lts]' = MTL;
    % written as Ax = b
    A = [LfibNormOnTen_ref LtenNorm_ref];
    b = MTL_targ;
     
    switch solver
        case 'modenese'
    % ===== LINSOL =======
    % solving the problem to calculate the muscle param
    x = A\b;
    LmOptLts_opt(n_mus+1,:) = x;

    % checking the results
    if min(x)<=0
        diary on
        % informing the user
        display(['Negative value estimated for muscle parameter of muscle ',curr_mus_name]);
        display( '                         Lm Opt        Lts'   );
        display(['Template model       : ',num2str(LmOptLts)]);
        display(['Optimized param      : ',num2str(LmOptLts_opt(n_mus+1,:))]);

        % ===== IMPLEMENTING CORRECTIONS IF ESTIMATION IS NOT CORRECT =======
        % first try lsqnonlin
        x = lsqnonneg(A,b);
        LmOptLts_opt(n_mus+1,:) = x;
        display(['Opt params (lsqnonneg): ',num2str(LmOptLts_opt(n_mus+1,:))]);
        % In our tests, if something goes wrong is generally tendon slack
        % length becoming negative or zero because tendon length doesn't change
        % throughout the range of motion, so lowering the rank of A.
        % if x(2)==0
        if min(x)<=0
            if (max(Mus_ref(okList,3))-min(Mus_ref(okList,3)))<0.0001
                display('Tendon length not changing throughout range of motion')
            end
            % calculating proportion of tendon and fiber
%            Lfib_fraction = LfibNormOnTen*LmOptLts(1)./MTL;
            Lten_fraction = Mus_ref(okList,3)./MTL_ref;
            Lten_targ = (Lten_fraction.*MTL_targ);

            % first round: optimizing Lopt maintaing the proportion of
            % tendon as in the reference model
            A_1 = LfibNormOnTen_ref;
            b_1 = (MTL_targ-Lten_targ);
            x(1) = A_1\b_1;

            % second round: using the optimized Lopt to recalculate Lts
            b_2 = MTL_targ-A_1*x(1);
            A_2 = LtenNorm_ref;
            x(2) = A_2\b_2;
            LmOptLts_opt(n_mus+1,:) = x;
        end
        diary off
    end
    
         case 'lsq'
             % ===== LINSOL =======
    % solving the problem to calculate the muscle param
    x = A\b;
    LmOptLts_opt(n_mus+1,:) = x;
                 
         case 'lsqnonneg'
             % ===== LSGNONNEG =======
    % solving the problem to calculate the muscle param
             x = lsqnonneg(A,b);
        LmOptLts_opt(n_mus+1,:) = x;     
               
         case 'twostep'
             % ===== TWO-STEP Solver =======
             
    % calculating proportion of tendon and fiber
%            Lfib_fraction = LfibNormOnTen*LmOptLts(1)./MTL;
            Lten_fraction = Mus_ref(okList,3)./MTL_ref;
            Lten_targ = (Lten_fraction.*MTL_targ);
            
            % initiatie output vector x, as it has not been defined yet
            % before when using this setting.
            x = zeros(2,1);

            % first round: optimizing Lopt maintaing the proportion of
            % tendon as in the reference model
            A_1 = LfibNormOnTen_ref;
            b_1 = (MTL_targ-Lten_targ);
            x(1) = A_1\b_1;

            % second round: using the optimized Lopt to recalculate Lts
            b_2 = MTL_targ-A_1*x(1);
            A_2 = LtenNorm_ref;
            x(2) = A_2\b_2;
            LmOptLts_opt(n_mus+1,:) = x;
    end
    
    end


    % Here tests about '\' against optimizers were implemented.

    %% NOTE: In case of zero tendon (slack) length, fval check is evaluated differently
    if LmOptLts(2) < 1e-6   % check if tendon slack length is zero 
        
    % In this case, A is a vector of normalized fiber length as all
    % normalized tendon lengths calculated are NaN  
        
    % calculating the error (sum of squared errors)
    fval = norm(A*x(1)-b).^2.0;
    
    else
    % calculating the error (sum of squared errors)
    fval = norm(A*x-b).^2.0;
    end

    % update muscles from scaled model
    curr_mus_scaled.setOptimalFiberLength(LmOptLts_opt(n_mus+1,1));
    curr_mus_scaled.setTendonSlackLength(LmOptLts_opt(n_mus+1,2));

    
    % Record duration of optimization
    OptTime = toc;
    % PRINT LOGS
    diary on
    display('  ');
    display(['Calculated optimized muscle parameters for ', char(curr_mus),' in ',num2str(OptTime),' seconds.'])
    display( '                         Lm Opt        Lts'   );
    display(['Template model       : ',num2str(LmOptLts)]);
    display(['Optimized param      : ',num2str(LmOptLts_opt(n_mus+1,:))]);
    display(['Nr of eval points    : ',num2str(evalOkPoints), '/',num2str(evalTotPoints),' used'])
    display(['fval                 : ',num2str(fval)]);
    display(['var from template [%]: ',num2str(100*(abs(LmOptLts-LmOptLts_opt(n_mus+1,:)))./LmOptLts),'%'])
    display('  ');
    diary off

    % SIMULATION INFO AND RESULTS
    SimInfo.colheader(n_mus+1)               = {char(curr_mus)};
    SimInfo.LmOptLts_ref(1:2,n_mus+1)        = LmOptLts;
    SimInfo.LmOptLts_opt(1:2,n_mus+1)        = LmOptLts_opt(n_mus+1,:);
    SimInfo.varPercLmOptLts(1:2,n_mus+1)     = 100*(abs(LmOptLts-LmOptLts_opt(n_mus+1,:)))./LmOptLts;
    SimInfo.sampledEvalPoints(n_mus+1)       = evalOkPoints;
    SimInfo.usedEvalPoints(n_mus+1)          = evalTotPoints;
    SimInfo.fval(n_mus+1)                    = fval;
    % Added tracking of optimization time
    SimInfo.OptTime(n_mus+1)                 = OptTime;
end

% assigning optimized model as output
osimModel_opt = osimModel_targ;


%% Quality of life changes

% saving optimized model
fileoutpath = [name,'_Optimized.osim'];
osimModel_opt.print(fileoutpath);

% saving SimInfo as .mat structure
InfoStructName = ['OptInfo_',name,res_file_id_exp,'.mat'];
save(InfoStructName,'SimInfo');

end
