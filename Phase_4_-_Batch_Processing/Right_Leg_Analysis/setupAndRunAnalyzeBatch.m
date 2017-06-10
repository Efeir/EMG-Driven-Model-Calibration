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

trialsForAn = dir('*ik.mot');
grfTrials = dir('*GRF.mot');

results_folder = [dirorig];

nTrials = size(trialsForAn);
% matlabpool open 10

TrialDivision = ceil(nTrials(1)/10);

% Loop through the trials
%spmd
% for trial= 1:nTrials(1)
spmd  
    % Pull in the modeling classes straight from the OpenSim distribution
    import org.opensim.modeling.*
    
    anTool = AnalyzeTool([dirorig '\GenericAnalyzeTool.xml']);

    % Load the model and initialize
    model = Model([dirorig '\Patient4_optModel.osim']);
    model.initSystem();

    % Tell Tool to use the loaded model
    anTool.setModel(model);
    
    %for trial= (labindex-1)*TrialDivision+1:min([labindex*TrialDivision nTrials(1)])
    for trial= 1:nTrials(1)
        % Get the name of the file for this trial
        motionFile = trialsForAn(trial).name;

        % Create name of trial from .trc file name
        name = regexprep(motionFile,'_ik.mot','');
        fullpath = ([dirorig '\' motionFile]);

        % Set results folder
        results_folder = [dirorig '\' name '_Analyses'];

        % Get trc data to determine time range
        coordData = Storage(fullpath);
        % Get GRF data for the time vector, since ik sometimes omits trials, an
        % extra frame was added that should be ignored
        grfData = Storage([dirorig '\' grfTrials(trial).name]);

        % Get initial and intial time 
        initial_time = grfData.getFirstTime();
        final_time = grfData.getLastTime();

        % Setup the analyzeTool for this trial
        anTool.setName(name);
        anTool.setLowpassCutoffFrequency(7/((final_time-initial_time)*101/141));
        anTool.setResultsDir(results_folder);
        anTool.setInitialTime(initial_time);
        anTool.setFinalTime(final_time);
        anTool.setCoordinatesFileName(fullpath);

    %     % Save the settings in a setup file
    %     outfile = ['Setup_IK_' name '.xml'];
    %     analyzeTool.print([genericSetupPath '\' outfile]);

        fprintf(['\n=====Performing muscle analysis on cycle # ' num2str(trial) '=====\n\n'])
        % Run IK
        anTool.run();

    end
end

%parpool close

