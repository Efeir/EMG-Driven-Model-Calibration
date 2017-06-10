function compileMuscleData_l_v3
% For the left leg only
% v2 compiles data for duplicate muscles that have multiple heads or for
% similar muscles such as semiten and semimem

dirorig = cd;
Files = dir('*l.xls');
nMusc = 16;

MovementSpeedOld = [];

load muscle_analysis_l_results.mat

for i = 1:length(Files)
    
    muscleAnalFileName = Files(i).name;
    
    TrialName = strrep(muscleAnalFileName,'_muscle_analysis_l.xls','')
    FolderName = [TrialName '_Analyses'];
    motionFileName = strrep(muscleAnalFileName,'_muscle_analysis_l.xls','_ik.mot');
    statesFileName = strrep(muscleAnalFileName,'_muscle_analysis_l.xls','_StatesReporter_states.sto');
    IDFileName = strrep(muscleAnalFileName,'_muscle_analysis_l.xls','_id.sto');
    emgFileName = strrep(muscleAnalFileName,'_muscle_analysis_l.xls','_EMG.mot');
    
    % Extract Joint Coordinates
    [coords, coordinateLabels, Time] = readMotFile(motionFileName);
    TimeLongest = Time(1:141);

    % Extract Model States (includes joint velocities and positions
    cd(FolderName);
    [states, statesLabels, TimeBad] = readMotFile(statesFileName);
    TimeLongest = Time(1:141);
    cd(dirorig);
    
    % Extract Moment Arms
    [HipMA] = muscAnalysisData(i).MADataHipFlex;
    HipAAMA = muscAnalysisData(i).MADataHipAA;
    KneeMA = muscAnalysisData(i).MADataKnee;
    AnkleMA = muscAnalysisData(i).MADataAnkle;
    SubtalarMA = muscAnalysisData(i).MADataSubtalar;
    
    % Extract Fiber Lengths and velocities
    %     Lvf = xlsread(FileName, 'Fiber Velocity');
    %     Lmvtilda = xlsread(FileName, 'Normalized Fiber Velocity');
    Lmt = muscAnalysisData(i).MTL;
    
    % Load inverse dynamics torques
    IDnames = {'Hip FE Moment' 'Hip AA Moment' 'Knee FE Moment' 'Ankle FE Moment' 'Subtalar IE Moment'};
    [IDloads, IDLabels, TimeID] = readMotFile(IDFileName);
    IDloads = IDloads(:,[7 8 10 12 13]+8);
%     IDloads = filterData(IDloads, mean(diff(TimeID)),6);
%     IDloads(:,1) = SplineSmooth(TimeLongest, IDloads(:,1), .999999);
%     IDloads(:,2) = SplineSmooth(TimeLongest, IDloads(:,2), .999999);
%     IDloads(:,3) = SplineSmooth(TimeLongest, IDloads(:,3), .999999);
%     IDloads(:,4) = SplineSmooth(TimeLongest, IDloads(:,4), .999999);
%     IDloads(:,5) = SplineSmooth(TimeLongest, IDloads(:,5), .999999);
%     IDloads = IDloads(:,:);
    
    % Load emg data
    [emgData, emgLabels, TimeEMG] = readMotFile(emgFileName);
    % Reorder EMG data to match order found in open sim file
    % Original order
    % 17RAddLong	18RGlutMax	19RGlutMed	20RIliopsoas	21RSemiMemb	22RBicFemLong	...
    %23RRecFem	24RVastMed	25RVastLat	26RGasMed	27RTibAnt	28RTibPost	...
    %29RPeroneusLong	30RSol	31RExtDigLong	32RFlexDigLong
    emgLabels = emgLabels(2:17);
    emgData = emgData(:,1:16);
    
    % Take desired DOFs from coordinates (left leg only, sagittal plane)
    coords = coords(1:141,[7 8 10 12 13 9]+8);
%     jointVels = states(1:141,[54 55 57 59 60 56]);
    coordNames = {'Hip FE' 'Hip AA' 'Knee FE' 'Ankle FE' 'Subtalar IE' 'Hip IE'};
    jointVels = coords;
    [coords(:,1), jointVels(:,1)] = SplineSmooth(TimeLongest(1:length(coords)), coords(:,1), .9999);
    [coords(:,2), jointVels(:,2)] = SplineSmooth(TimeLongest(1:length(coords)), coords(:,2), .9999);
    [coords(:,3), jointVels(:,3)] = SplineSmooth(TimeLongest(1:length(coords)), coords(:,3), .9999);
    [coords(:,4), jointVels(:,4)] = SplineSmooth(TimeLongest(1:length(coords)), coords(:,4), .9999);
    [coords(:,5), jointVels(:,5)] = SplineSmooth(TimeLongest(1:length(coords)), coords(:,5), .9999);
    [coords(:,6), jointVels(:,6)] = SplineSmooth(TimeLongest(1:length(coords)), coords(:,6), .9999);
    
    % Take desired muscles values (right leg, 16 muscles)
    HipMA = HipMA(:,[1:35]+35);
    HipAAMA = HipAAMA(:,[1:35]+35);
    KneeMA = KneeMA(:,[1:35]+35);
    AnkleMA = AnkleMA(:,[1:35]+35);
    SubtalarMA = SubtalarMA(:,[1:35]+35);
    MuscNames = muscAnalysisData(i).ColumnLabels([2:36]+35);
    Lmt = Lmt(:,[1:35]+35);
    nMusc = length([1:35]);
    
    Vmt = zeros(size(Lmt));
    for j = 1:nMusc
        [junk, Vmt(:,j)] = SplineSmooth(TimeLongest(1:length(coords)), Lmt(:,j), .99999);
    end
    
    MovementSpeed = TrialName(10:20);
    if strcmp(MovementSpeed, MovementSpeedOld)
        count = count+1;
    else
        count = 1;
    end
    
    MovementSpeedOld = MovementSpeed;
    
    eval(['etmData.' MovementSpeed '(' num2str(count) ')' '.TrialName = TrialName;']);
    eval(['etmData.' MovementSpeed '(' num2str(count) ')' '.Time = TimeLongest;']);
    eval(['etmData.' MovementSpeed '(' num2str(count) ')' '.EMG = emgData;']);
    eval(['etmData.' MovementSpeed '(' num2str(count) ')' '.EMGlabels = emgLabels;']);
    eval(['etmData.' MovementSpeed '(' num2str(count) ')' '.OpensimMuscleLabels = MuscNames;']);
    eval(['etmData.' MovementSpeed '(' num2str(count) ')' '.HipFEMomentArms = HipMA;']);
    eval(['etmData.' MovementSpeed '(' num2str(count) ')' '.HipAAMomentArms = HipAAMA;']);
    eval(['etmData.' MovementSpeed '(' num2str(count) ')' '.KneeMomentArms = KneeMA;']);
    eval(['etmData.' MovementSpeed '(' num2str(count) ')' '.AnkleMomentArms = AnkleMA;']);
    eval(['etmData.' MovementSpeed '(' num2str(count) ')' '.SubtalarMomentArms = SubtalarMA;']);
    eval(['etmData.' MovementSpeed '(' num2str(count) ')' '.MuscleTendonLengths = Lmt;']);
    eval(['etmData.' MovementSpeed '(' num2str(count) ')' '.MuscleTendonVelocities = Vmt;']);
    eval(['etmData.' MovementSpeed '(' num2str(count) ')' '.CoordinateLabels = coordNames;']);
    eval(['etmData.' MovementSpeed '(' num2str(count) ')' '.JointAngles = coords;']);
    eval(['etmData.' MovementSpeed '(' num2str(count) ')' '.JointVelocities = jointVels;']);
    eval(['etmData.' MovementSpeed '(' num2str(count) ')' '.InverseDynamicsLoadLabels = IDnames;']);
    eval(['etmData.' MovementSpeed '(' num2str(count) ')' '.InverseDynamicsLoads = IDloads;']);
    
    
end

save Patient4_etmData_left_reprocessed2.mat etmData

function [Data, ColumnLabels, Time] = readMotFile(inFile)

%open file
fid = fopen(inFile);

%skip initial lines
line = 'asdf';
while ~strcmp(upper(line),'ENDHEADER')
    line = fgetl(fid);
end

% Find the labels on the columns
ColumnLabels = fgetl(fid);
ColumnLabels = textscan(ColumnLabels, '%t');
ColumnLabels = ColumnLabels{1}';

% Read numeric data
Data = fscanf(fid, '%f');

% Find length of data
nColumns = length(ColumnLabels);
nRows = length(Data)/nColumns;

% Reshape data into matrix
Data = reshape(Data, nColumns, nRows)';

% Extract time vector
Time = Data(:,1);
Data(:,1) = [];

fclose(fid);