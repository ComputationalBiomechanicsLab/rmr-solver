% ----------------------------------------------------------------------- %
% The OpenSim API is a toolkit for musculoskeletal modeling and           %
% simulation. See http://opensim.stanford.edu and the NOTICE file         %
% for more information. OpenSim is developed at Stanford University       %
% and supported by the US National Institutes of Health (U54 GM072970,    %
% R24 HD065690) and by DARPA through the Warrior Web program.             %
%                                                                         %
% Copyright (c) 2005-2018 Stanford University and the Authors             %
% Author(s): James Dunne                                                  %
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
% ----------------------------------------------------------------------- %

%% Custom version of c3dExport.m
% It takes as input a .c3d file obtained with a motion capture system such
% as Qualisys, where just marker positions where recorder (no force data).
% It saves then to a .trc file whose name is specified by the user.
%
% author: Italo Belli (i.belli@tudelft.nl), on 7/12/2021
clear 
clc

% Load OpenSim libs
import org.opensim.modeling.*

% Get the path to a C3D file
[filename, path] = uigetfile('*.c3d');
c3dpath = fullfile(path,filename);

% Construct an opensimC3D object with input c3d path
% Constructor takes full path to c3d file and an integer for forceplate
% representation (1 = COP). The integer must be input for the function to work,
% but it does not really matter since the force data are not present.
c3d = osimC3D(c3dpath,1);

% Get the number of marker trajectories
nTrajectories = c3d.getNumTrajectories();

% Get the marker data rate
rMakers = c3d.getRate_marker();

% Get Start and end time
t0 = c3d.getStartTime();
tn = c3d.getEndTime();

% Rotate the data (to comply with OpenSim convention for the frames)
c3d.rotateData('x',-90)

% Write marker data to trc file, with some examples before
%
% c3d.writeTRC()                       Write to dir of input c3d.
% c3d.writeTRC('Walking.trc')          Write to dir of input c3d with defined file name.
% c3d.writeTRC('C:/data/Walking.trc')  Write to defined path input path.

dlgtitle = 'Input';
dims = [1 80];
filename_out = strsplit( filename , '.' );
filename_out = append(filename_out{1} ,'.trc');
definput = {filename_out};
filename_out = inputdlg({'Choose a name for your file (rember to include the .trc extension!):'},dlgtitle,dims,definput);
c3d.writeTRC(char(filename_out));