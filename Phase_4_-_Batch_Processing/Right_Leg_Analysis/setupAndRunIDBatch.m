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

IKtrials = dir('*ik.mot');
GRFtrials = dir('*GRF.mot');

results_folder = [dirorig];

nTrials = size(IKtrials);

% Loop through the trials
for trial= 600:nTrials(1)
    
    activeExternalForceFile = 'ExternalLoadsFile.xml';
    
    % Pull in the modeling classes straight from the OpenSim distribution
    import org.opensim.modeling.*
    
    % Load the model and initialize
    model = Model([dirorig '\Patient4_optModel.osim']);
    model.initSystem();
    
    % Create an external loads object
    extLoadsObject = ExternalLoads(model,activeExternalForceFile);
    
    idTool = InverseDynamicsTool([dirorig '\GenericIDSetupFile.xml']);
    idTool.setResultsDir(results_folder);
    
    % Tell Tool to use the loaded model
    idTool.setModel(model);
    
    % Get the name of the file for this trial
    kinematicsFile = IKtrials(trial).name;
    GRFFile = GRFtrials(trial).name;
    
    % Create name of trial from .trc file name
    name = kinematicsFile;
    fullpath = ([dirorig '\' kinematicsFile]);
    fullpathGRF = ([dirorig '\' GRFFile]);
    
    % Set the grf filename
    extLoadsObject.setDataFileName(fullpathGRF);
    
    extLoadsObject.print([dirorig '\' activeExternalForceFile]);
    
    % Print updated external load xml file
    
    % Get trc data to determine time range
    motData = Storage(fullpath);
    GRFData = Storage(fullpathGRF);
    
    % Get initial and intial time
    initial_time = GRFData.getFirstTime();
    final_time = GRFData.getLastTime();

    % Setup the ikTool for this trial
    idTool.setName(name);
    idTool.setLowpassCutoffFrequency(7/((final_time-initial_time)*101/141));
    idTool.setCoordinatesFileName(fullpath);
    idTool.setStartTime(initial_time);
    idTool.setEndTime(final_time);
    idTool.setOutputGenForceFileName(strrep(name, 'ik.mot', 'id'));
    %     idTool.getExternalLoads.setDataFileName(fullpathGRF);
    %     idTool.updExternalLoads;
    %     % Save the settings in a setup file
    %     outfile = ['Setup_IK_' name '.xml'];
    %     idTool.print([genericSetupPath '\' outfile]);
    
    fprintf(['\n==== Performing ID on cycle # ' num2str(trial) '====\n\n'])
    % Run IK
    idTool.run();
    
end
