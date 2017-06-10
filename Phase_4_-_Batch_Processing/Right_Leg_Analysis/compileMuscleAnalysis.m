function compileMuscleAnalysis

dirOrig = cd;

Trials = dir('*ik.mot');
nTrials = size(Trials);

muscAnalysisData(nTrials(1)).initialize = 0;

for i = 1:nTrials(1)
    
    name = strrep(Trials(i).name, '_ik.mot', '')
    
    % Move to muscle analysis folder
    cd([dirOrig '\' name '_Analyses']);
    
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
    [muscAnalysisData(i).MADataHipFlex, muscAnalysisData(i).ColumnLabels, muscAnalysisData(i).Time] = readMotFile(MAFileHipFlex);
    [muscAnalysisData(i).MADataHipAA, muscAnalysisData(i).ColumnLabels, muscAnalysisData(i).Time] = readMotFile(MAFileHipAA);
    [muscAnalysisData(i).MADataKnee, muscAnalysisData(i).ColumnLabels, muscAnalysisData(i).Time] = readMotFile(MAFileKnee);
    [muscAnalysisData(i).MADataAnkle, muscAnalysisData(i).ColumnLabels, muscAnalysisData(i).Time] = readMotFile(MAFileAnkle);
    [muscAnalysisData(i).MADataSubtalar, muscAnalysisData(i).ColumnLabels, muscAnalysisData(i).Time] = readMotFile(MAFileSubtalar);
    
    [muscAnalysisData(i).FL, muscAnalysisData(i).ColumnLabels, muscAnalysisData(i).Time] = readMotFile(FLFile);
    [muscAnalysisData(i).FV, muscAnalysisData(i).ColumnLabels, muscAnalysisData(i).Time] = readMotFile(FVFile);
    [muscAnalysisData(i).FLtilda, muscAnalysisData(i).ColumnLabels, muscAnalysisData(i).Time] = readMotFile(FLtildaFile);
    [muscAnalysisData(i).FVtilda, muscAnalysisData(i).ColumnLabels, muscAnalysisData(i).Time] = readMotFile(FVtildaFile);
    [muscAnalysisData(i).alpha, muscAnalysisData(i).ColumnLabels, muscAnalysisData(i).Time] = readMotFile(alphaFile);
    [muscAnalysisData(i).MTL, muscAnalysisData(i).ColumnLabels, muscAnalysisData(i).Time] = readMotFile(MTLengthFile);
    
     %cd(dirOrig);
%outFile = [name '_muscle_analysis_r.xls'];
%xlswrite(outFile,[muscAnalysisData(i).ColumnLabels; num2cell([muscAnalysisData(i).Time muscAnalysisData(i).MADataHipFlex])],'Hip Flexion Moment Arms');
%xlswrite(outFile,[muscAnalysisData(i).ColumnLabels; num2cell([muscAnalysisData(i).Time muscAnalysisData(i).MADataHipAA])],'Hip Adduction Moment Arms');
%xlswrite(outFile,[muscAnalysisData(i).ColumnLabels; num2cell([muscAnalysisData(i).Time muscAnalysisData(i).MADataKnee])],'Knee Moment Arms');
%xlswrite(outFile,[muscAnalysisData(i).ColumnLabels; num2cell([muscAnalysisData(i).Time muscAnalysisData(i).MADataAnkle])],'Ankle Moment Arms');
%xlswrite(outFile,[muscAnalysisData(i).ColumnLabels; num2cell([muscAnalysisData(i).Time muscAnalysisData(i).MADataSubtalar])],'Subtalar Moment Arms');
%xlswrite(outFile,[muscAnalysisData(i).ColumnLabels; num2cell([muscAnalysisData(i).Time muscAnalysisData(i).FL])],'Fiber Length');    
%xlswrite(outFile,[muscAnalysisData(i).ColumnLabels; num2cell([muscAnalysisData(i).Time muscAnalysisData(i).FV])],'Fiber Velocity');  
%xlswrite(outFile,[muscAnalysisData(i).ColumnLabels; num2cell([muscAnalysisData(i).Time muscAnalysisData(i).FLtilda])],'Normalized Fiber Length');
%xlswrite(outFile,[muscAnalysisData(i).ColumnLabels; num2cell([muscAnalysisData(i).Time muscAnalysisData(i).FVtilda])],'Normalized Fiber Velocity');
%xlswrite(outFile,[muscAnalysisData(i).ColumnLabels; num2cell([muscAnalysisData(i).Time muscAnalysisData(i).alpha])],'Pennation Angle');
%xlswrite(outFile,[muscAnalysisData(i).ColumnLabels; num2cell([muscAnalysisData(i).Time muscAnalysisData(i).MTL])],'Muscle Tendon Length');
     
end

cd(dirOrig);
save muscle_analysis_r_results.mat muscAnalysisData

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