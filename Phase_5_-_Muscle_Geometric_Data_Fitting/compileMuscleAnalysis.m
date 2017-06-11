function compileMuscleAnalysis

dirOrig = cd;
 
    % Move to muscle analysis folder
    cd([dirOrig '\Muscle Analysis']);
    
    name = 'Patient4_optModel'
    
    % Create moment arm file names of interest
    MAFileHipFlex = [name '_MuscleAnalysis_MomentArm_hip_flexion_r.sto'];
    MAFileHipAA = [name '_MuscleAnalysis_MomentArm_hip_adduction_r.sto'];
    MAFileKnee = [name '_MuscleAnalysis_MomentArm_knee_angle_r.sto'];
    MAFileAnkle = [name '_MuscleAnalysis_MomentArm_ankle_angle_r.sto'];
    MAFileSubtalar = [name '_MuscleAnalysis_MomentArm_subtalar_angle_r.sto'];
    
    % Create names for other quantities of interest
    FLFile = [name '_MuscleAnalysis_FiberLength.sto'];
    FVFile = [name '_MuscleAnalysis_FiberVelocity.sto'];
    FLtildaFile = [name '_MuscleAnalysis_NormalizedFiberLength.sto'];
    FVtildaFile = [name '_MuscleAnalysis_NormFiberVelocity.sto'];
    alphaFile = [name '_MuscleAnalysis_PennationAngle.sto'];
    MTLengthFile = [name '_MuscleAnalysis_Length.sto'];

    % Read data files
    [MADataHipFlex, ColumnLabels, Time] = readMotFile(MAFileHipFlex);
    [MADataHipAA] = readMotFile(MAFileHipAA);
    [MADataKnee] = readMotFile(MAFileKnee);
    [MADataAnkle] = readMotFile(MAFileAnkle);
    [MADataSubtalar] = readMotFile(MAFileSubtalar);
    
    [FL] = readMotFile(FLFile);
    [FV] = readMotFile(FVFile);
    [FLtilda] = readMotFile(FLtildaFile);
    [FVtilda] = readMotFile(FVtildaFile);
    [alpha] = readMotFile(alphaFile);
    [MTL] = readMotFile(MTLengthFile);
    
    cd(dirOrig);
    outFile = ['muscle_analysis_test1.xls'];
    
    xlswrite(outFile,[ColumnLabels; num2cell([Time MADataHipFlex])],'Hip Flexion Moment Arms');
    
    xlswrite(outFile,[ColumnLabels; num2cell([Time MADataHipAA])],'Hip Adduction Moment Arms');
    
    xlswrite(outFile,[ColumnLabels; num2cell([Time MADataKnee])],'Knee Moment Arms');
    
    xlswrite(outFile,[ColumnLabels; num2cell([Time MADataAnkle])],'Ankle Moment Arms');
    
    xlswrite(outFile,[ColumnLabels; num2cell([Time MADataSubtalar])],'Subtalar Moment Arms');
    
    xlswrite(outFile,[ColumnLabels; num2cell([Time FL])],'Fiber Length');
    
    xlswrite(outFile,[ColumnLabels; num2cell([Time FV])],'Fiber Velocity');
    
    xlswrite(outFile,[ColumnLabels; num2cell([Time FLtilda])],'Normalized Fiber Length');
    
    xlswrite(outFile,[ColumnLabels; num2cell([Time FVtilda])],'Normalized Fiber Velocity');
    
    xlswrite(outFile,[ColumnLabels; num2cell([Time alpha])],'Pennation Angle');
    
    xlswrite(outFile,[ColumnLabels; num2cell([Time MTL])],'Muscle Tendon Length');


function [Data, ColumnLabels, Time] = readMotFile(inFile)

%open file
fid = fopen(inFile);

%skip initial lines
for i = 1:11
   junk = fgetl(fid); 
   if i == 3
       nRows = str2double(strrep(junk,'nRows=',''));
   elseif i == 4
       nColumns = str2double(strrep(junk,'nColumns=',''));
   end
end

% Find the labels on the columns
ColumnLabels = fgetl(fid);
ColumnLabels = textscan(ColumnLabels, '%t');
ColumnLabels = ColumnLabels{1}';

Data = [];
% Read numeric data
% count = 1;
while length(Data)<nRows*nColumns
    DataPiece = fscanf(fid, '%f');
    junk = fscanf(fid, '%s', 1);
    Data = [Data; DataPiece];
%     if count>1e12
%         error('Infinite Loop')
%     end
%     count = count+1;    
end

Data(Data==-1) = NaN;

% Reshape data into matrix
Data = reshape(Data, nColumns, nRows)';

% Extract time vector
Time = Data(:,1);
Data(:,1) = [];

fclose(fid);