% ----------------------------------------------------------------------- %
% The OpenSim API is a toolkit for musculoskeletal modeling and           %
% simulation. See http://opensim.stanford.edu and the NOTICE file         %
% for more information. OpenSim is developed at Stanford University       %
% and supported by the US National Institutes of Health (U54 GM072970,    %
% R24 HD065690) and by DARPA through the Warrior Web program.             %
%                                                                         %   
% Copyright (c) 2005-2012 Stanford University and the Authors             %
% Author(s): Edith Arnold                                                 %  
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

% setupAndRunIKBatchExample.m                                                 
% Author: Edith Arnold

% This example script runs multiple inverse kinematics trials for the model Subject01. 
% All input files are in the folder ../Matlab/testData/Subject01
% To see the results load the model and ik output in the GUI.

dirorig = cd;

trialsForIK = dir('*.trc');

results_folder = [dirorig];

nTrials = size(trialsForIK);

% Loop through the trials
parfor trial= 1:nTrials(1)
    
    % Pull in the modeling classes straight from the OpenSim distribution
    import org.opensim.modeling.*
    
    % Initialize IKtool object
    ikTool = InverseKinematicsTool([dirorig '\GenericIKSetupFile.xml']);

    % Load the model and initialize
    model = Model([dirorig '\Patient4_optModel.osim']);
    model.initSystem();

    % Tell Tool to use the loaded model
    ikTool.setModel(model);
    
    % Get the name of the file for this trial
    markerFile = trialsForIK(trial).name;
    
    % Create name of trial from .trc file name
    name = regexprep(markerFile,'.trc','');
    fullpath = ([dirorig '\' markerFile]);
    
    % Get trc data to determine time range
    markerData = MarkerData(fullpath);
    
    % Get initial and intial time 
    initial_time = markerData.getStartFrameTime();
    final_time = markerData.getLastFrameTime();

    % Setup the ikTool for this trial
    ikTool.setName(name);
    ikTool.setMarkerDataFileName(fullpath);
    ikTool.setStartTime(initial_time);
    ikTool.setEndTime(final_time);
    ikTool.setOutputMotionFileName([results_folder '\' name '_ik.mot']);

%     % Save the settings in a setup file
%     outfile = ['Setup_IK_' name '.xml'];
%     ikTool.print([genericSetupPath '\' outfile]);
    
    fprintf(['Performing IK on cycle # ' num2str(trial) '\n']);
    % Run IK
    ikTool.run();

end
