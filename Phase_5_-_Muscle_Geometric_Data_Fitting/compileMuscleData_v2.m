function compileMuscleData_v2
% For the left leg only, but since bilateral symmetry is maintained it
% works for both legs
% v2 compiles data for duplicate muscles that have multiple heads or for
% similar muscles such as semiten and semimem

dirorig = cd;
Files = dir('*.xls');
nMusc = 16;

MovementSpeedOld = [];

for i = 1:length(Files)
    
    muscleAnalFileName = Files(i).name;
    
    TrialName = 'muscle_analysis_test1.xls'
    FolderName = [TrialName ' Muscle Analysis'];
    motionFileName = 'MuscGeometrySampling.mot';

    % Extract Joint Coordinates
    [coords, coordinateLabels, Time] = readMotFile(motionFileName);
    TimeLongest = Time(1:1000);
    
    % Extract Moment Arms
    [HipMA, txt] = xlsread(muscleAnalFileName, 'Hip Flexion Moment Arms', '', 'basic');
    HipAAMA = xlsread(muscleAnalFileName, 'Hip Adduction Moment Arms', '', 'basic');
    KneeMA = xlsread(muscleAnalFileName, 'Knee Moment Arms', '', 'basic');
    AnkleMA = xlsread(muscleAnalFileName, 'Ankle Moment Arms', '', 'basic');
    SubtalarMA = xlsread(muscleAnalFileName, 'Subtalar Moment Arms', '', 'basic');
    
    % Extract Fiber Lengths and velocities
    Lmt = xlsread(muscleAnalFileName, 'Muscle Tendon Length', '', 'basic');
    
    % Take desired DOFs from coordinates (left leg only, sagittal plane)
    coords = coords(1:1000,[15 16 18 20 21 17]);
    coordNames = {'Hip FE' 'Hip AA' 'Knee FE' 'Ankle FE' 'Subtalar IE' 'Hip IE'};
  
    % Take desired muscles values (right leg, 16 muscles)
    HipMA = HipMA(:,[2:36]);
    HipAAMA = HipAAMA(:,[2:36]);
    KneeMA = KneeMA(:,[2:36]);
    AnkleMA = AnkleMA(:,[2:36]);
    SubtalarMA = SubtalarMA(:,[2:36]);

    MuscNames = txt([2:36]);
    Lmt = Lmt(:,[2:36]);
    nMusc = length([2:36]);
    
    eval(['MuscleGeoFitData.TrialName = TrialName;']);
    eval(['MuscleGeoFitData.Time = TimeLongest;']);
    eval(['MuscleGeoFitData.OpensimMuscleLabels = MuscNames;']);
    eval(['MuscleGeoFitData.HipFEMomentArms = HipMA;']);
    eval(['MuscleGeoFitData.HipAAMomentArms = HipAAMA;']);
    eval(['MuscleGeoFitData.KneeMomentArms = KneeMA;']);
    eval(['MuscleGeoFitData.AnkleMomentArms = AnkleMA;']);
    eval(['MuscleGeoFitData.SubtalarMomentArms = SubtalarMA;']);
    eval(['MuscleGeoFitData.MuscleTendonLengths = Lmt;']);
    eval(['MuscleGeoFitData.CoordinateLabels = coordNames;']);
    eval(['MuscleGeoFitData.JointAngles = coords;']);
    
end

save Patient4_MuscleGeoFitData.mat MuscleGeoFitData

function [Data, ColumnLabels, Time] = readMotFile(inFile)

%open file
fid = fopen(inFile);

%skip initial lines
line = 'asdfsadfsadfsadf';
while ~strcmp(upper(line(1:9)),'ENDHEADER')
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