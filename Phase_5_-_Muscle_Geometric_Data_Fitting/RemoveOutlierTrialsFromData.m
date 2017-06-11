function RemoveOutlierTrialsFromData

load etmData_refitted_both_4.mat

DataSpeeds = {etmData_l.TMGait_0pt2 etmData_l.TMGait_0pt3 etmData_l.TMGait_0pt4 etmData_l.TMGait_0pt5 etmData_l.TMGait_0pt6 etmData_l.TMGait_0pt7 etmData_l.TMGait_0pt8...
    etmData_r.TMGait_0pt2 etmData_r.TMGait_0pt3 etmData_r.TMGait_0pt4 etmData_r.TMGait_0pt5 etmData_r.TMGait_0pt6 etmData_r.TMGait_0pt7 etmData_r.TMGait_0pt8};

RemovedTrials = 0;
nTrialsAll = 0;
for k = 1:14
    
    Data = DataSpeeds{k};
    
    nTrials = length(Data);
    nTrialsAll = nTrialsAll+nTrials;
    nptsLong = 141;
    nptsShort = 101;
    nMusc = 35;
    SampleStep = 1;
    
    TrialsSpeed = zeros(1,nTrials);
    Time = zeros(nptsLong,nTrials);
    JAngles = zeros(nptsLong,6,nTrials);
    JVels = zeros(nptsLong,6,nTrials);
    HipMA = zeros(nptsLong,nMusc,nTrials);
    HipAAMA = zeros(nptsLong,nMusc,nTrials);
    KneeMA = zeros(nptsLong,nMusc,nTrials);
    AnkleMA = zeros(nptsLong,nMusc,nTrials);
    SubtalarMA = zeros(nptsLong,nMusc,nTrials);
    Lmt = zeros(nptsLong,nMusc,nTrials);
    Vmt = zeros(nptsLong,nMusc,nTrials);
    IDloads = zeros(nptsLong,5,nTrials);
    EMG = zeros(nptsLong,16,nTrials);
    
    for i = 1:nTrials
        
        Time(:,i) = Data(i).Time(1:SampleStep:141);
        JAngles(:,:,i) = Data(i).JointAngles(1:SampleStep:141,1:6);
        JVels(:,:,i) = Data(i).JointVelocities(1:SampleStep:141,1:6);
        HipMA(:,:,i) = Data(i).HipFEMomentArms(1:SampleStep:141,:);
        HipAAMA(:,:,i) = Data(i).HipAAMomentArms(1:SampleStep:141,:);
        KneeMA(:,:,i) = Data(i).KneeMomentArms(1:SampleStep:141,:);
        AnkleMA(:,:,i) = Data(i).AnkleMomentArms(1:SampleStep:141,:);
        SubtalarMA(:,:,i) = Data(i).SubtalarMomentArms(1:SampleStep:141,:);
        Lmt(:,:,i) = Data(i).MuscleTendonLengths(1:SampleStep:141,:);
        Vmt(:,:,i) = Data(i).MuscleTendonVelocities(1:SampleStep:141,:);
        IDloads(:,:,i) = Data(i).InverseDynamicsLoads(1:SampleStep:141,:);
        
        tf = (Time(end,i)-Time(1,i))*(101/141);
        lpf = 7/tf;
        
        EMGtrial = processEMG(Data(i).EMG, 1/1000, 40, lpf);
        
        tf = (length(EMGtrial)-1)/1000;
        
        TimeEMGtrial = 0:(1/1000):tf;
        Time(:,i) = linspace(0,tf,141);
        % Normalize to nptsShort time frames

        EMG(:,:,i) = spline(TimeEMGtrial, EMGtrial', Time(:,i))';
        
    end
    
    EMG = EMG(:,[ones(1,6) 2*ones(1,3) 3*ones(1,6) 4*ones(1,2) 5*ones(1,2) 6*ones(1,2) 7:9 9 10*ones(1,2) 11:12 13*ones(1,3) 14:16],:);
    
    DataStruct.Time = Time;
    DataStruct.TrialsSpeed = TrialsSpeed;
    DataStruct.JAngles = JAngles;
    DataStruct.JVels = JVels;
    DataStruct.HipMA = HipMA;
    DataStruct.HipAAMA = HipAAMA;
    DataStruct.KneeMA = KneeMA;
    DataStruct.AnkleMA = AnkleMA;
    DataStruct.SubtalarMA = SubtalarMA;
    DataStruct.Lmt = Lmt;
    DataStruct.Vmt = Vmt;
    DataStruct.IDloads = IDloads;
    DataStruct.EMG = EMG;
    DataStruct.nTrials = nTrials;
    
    [DataStruct, outliers] = removeOutliers(DataStruct);
    
    fprintf(['\nRemoved ' num2str(nTrials-DataStruct.nTrials) ' trials out of ' num2str(nTrials) '.\n'])
    
    RemovedTrials = RemovedTrials+nTrials-DataStruct.nTrials;
    
    Time = DataStruct.Time;
    TrialsSpeed = DataStruct.TrialsSpeed;
    JAngles = DataStruct.JAngles;
    JVels = DataStruct.JVels;
    HipMA = DataStruct.HipMA;
    HipAAMA = DataStruct.HipAAMA;
    KneeMA = DataStruct.KneeMA;
    AnkleMA = DataStruct.AnkleMA;
    SubtalarMA = DataStruct.SubtalarMA;
    Lmt = DataStruct.Lmt;
    Vmt = DataStruct.Vmt;
    IDloads = DataStruct.IDloads;
    EMG = DataStruct.EMG;
    nTrials = DataStruct.nTrials;
    
    Data(outliers) = [];
    DataSpeeds{k} = Data;
    
end

fprintf(['\nRemoved ' num2str(RemovedTrials) ' trials out of ' num2str(nTrialsAll) '.\n'])

etmData_l.TMGait_0pt2 = DataSpeeds{1};
etmData_l.TMGait_0pt3 = DataSpeeds{2};
etmData_l.TMGait_0pt4 = DataSpeeds{3};
etmData_l.TMGait_0pt5 = DataSpeeds{4};
etmData_l.TMGait_0pt6 = DataSpeeds{5};
etmData_l.TMGait_0pt7 = DataSpeeds{6};
etmData_l.TMGait_0pt8 = DataSpeeds{7};
etmData_r.TMGait_0pt2 = DataSpeeds{8};
etmData_r.TMGait_0pt3 = DataSpeeds{9};
etmData_r.TMGait_0pt4 = DataSpeeds{10};
etmData_r.TMGait_0pt5 = DataSpeeds{11};
etmData_r.TMGait_0pt6 = DataSpeeds{12};
etmData_r.TMGait_0pt7 = DataSpeeds{13};
etmData_r.TMGait_0pt8 = DataSpeeds{14};

save etmData_refitted_both_outliersRemoved.mat etmData_l etmData_r

function emgFiltered = processEMG(emgData, dt, HPcutoff, LPcutoff)
% Function that process EMG data according the the process described by
% Lloyd and Besier (2003)
% Inputs:
%       emgData is a M x N matrix of unfiltered EMG data where M is the
%           number of time frames and N is the number of muscles
%       dt is the time step used in the EMG data
%       HPcutoff is the cutoff frequency for the high pass filter
%       LPcutoff is the cutoff frequency for the low pass filter
% Outputs:
%       emgFiltered is the filtered EMG data

SampRate = 1/mean(dt);

% High pass filter the data
degree = 4;
[z,p,k] = butter(degree, HPcutoff/SampRate*2,'high');
[sos,g] = zp2sos(z,p,k);      % Convert to SOS(second-order-sections) form
EMG_HP = filtfilt(sos,g,emgData);

% Demean
EMG_DM = EMG_HP-ones(size(EMG_HP,1),1)*mean(EMG_HP);

% Rectify
EMG_RF = abs(EMG_DM);

% Low pass filter
degree = 4;

[z,p,k] = butter(degree,LPcutoff/SampRate*2);
[sos,g] = zp2sos(z,p,k);      % Convert to SOS(second-order-sections) form
% pad data with final values
numPad = 0;
EMG_Padded = [ones(numPad,1)*EMG_RF(1,:);EMG_RF;ones(numPad,1)*EMG_RF(end,:)];
EMG_LP = filtfilt(sos,g,EMG_Padded);
emgFiltered = EMG_LP(numPad+1:end-numPad,:);

% Remove any negative EMG values that may still exist
emgFiltered(emgFiltered<0) = 0;

function [Data, outliers] = removeOutliers(DataStruct)

Time = DataStruct.Time;
TrialsSpeed = DataStruct.TrialsSpeed;
JAngles = DataStruct.JAngles;
JVels = DataStruct.JVels;
HipMA = DataStruct.HipMA;
HipAAMA = DataStruct.HipAAMA;
KneeMA = DataStruct.KneeMA;
AnkleMA = DataStruct.AnkleMA;
SubtalarMA = DataStruct.SubtalarMA;
Lmt = DataStruct.Lmt;
Vmt = DataStruct.Vmt;
IDloads = DataStruct.IDloads;
EMG = DataStruct.EMG;
nTrials = DataStruct.nTrials;

%Remove any trials with NaN moment arms
outliers = sign(max(max(isnan(HipMA)))+max(max(isnan(HipAAMA)))+max(max(isnan(KneeMA)))...
    +max(max(isnan(AnkleMA)))+max(max(isnan(SubtalarMA))));

outliers = logical(permute(outliers,[3 1 2]));

%Remove any trials with excessively large EMG values pre-normalization
EMGmax = permute(max(EMG,[],1), [2 3 1]);
EMGmaxmean = mean(EMGmax,2);

outliers = outliers+max(EMGmax>EMGmaxmean*ones(1,nTrials)*2)';

%Remove any trials with strage ID load curves
IDloadsmean = mean(IDloads,3);
IDloadsmean = repmat(IDloadsmean, [1 1 nTrials]);

IDloadsdiff = IDloads-IDloadsmean;
IDloadsdiffstd = std(IDloadsdiff,[],3);
IDloadsdiffstd = repmat(IDloadsdiffstd, [1 1 nTrials]);
outliers = outliers+permute(max(max(abs(IDloadsdiff)>5*IDloadsdiffstd,[],1),[],2), [3 1 2]);

%Remove trials with strange joint angles
JAnglesmean = mean(JAngles,3);
JAnglesmean = repmat(JAnglesmean, [1 1 nTrials]);

JAnglesdiff = JAngles-JAnglesmean;
JAnglesdiffstd = std(JAnglesdiff,[],3);
JAnglesdiffstd = repmat(JAnglesdiffstd, [1 1 nTrials]);
outliers = outliers+permute(max(max(abs(JAnglesdiff)>5*JAnglesdiffstd,[],1),[],2), [3 1 2]);

%Remove trials with strange joint velocities
JVelsmean = mean(JVels,3);
JVelsmean = repmat(JVelsmean, [1 1 nTrials]);

JVelsdiff = JVels-JVelsmean;
JVelsdiffstd = std(JVelsdiff,[],3);
JVelsdiffstd = repmat(JVelsdiffstd, [1 1 nTrials]);
outliers = outliers+permute(max(max(abs(JVelsdiff)>5*JVelsdiffstd,[],1),[],2), [3 1 2]);

%Remove trials with strange hip fe moment arms
HipMAmean = mean(HipMA,3);
HipMAmean = repmat(HipMAmean, [1 1 nTrials]);

HipMAdiff = HipMA-HipMAmean;
HipMAdiffstd = std(HipMAdiff,[],3);
HipMAdiffstd = repmat(HipMAdiffstd, [1 1 nTrials]);
outliers = outliers+permute(max(max(abs(HipMAdiff)>5*HipMAdiffstd,[],1),[],2), [3 1 2]);

%Remove trials with strange hip aa moment arms
HipAAMAmean = mean(HipAAMA,3);
HipAAMAmean = repmat(HipAAMAmean, [1 1 nTrials]);

HipAAMAdiff = HipAAMA-HipAAMAmean;
HipAAMAdiffstd = std(HipAAMAdiff,[],3);
HipAAMAdiffstd = repmat(HipAAMAdiffstd, [1 1 nTrials]);
outliers = outliers+permute(max(max(abs(HipAAMAdiff)>5*HipAAMAdiffstd,[],1),[],2), [3 1 2]);

%Remove trials with strange knee fe moment arms
KneeMAmean = mean(KneeMA,3);
KneeMAmean = repmat(KneeMAmean, [1 1 nTrials]);

KneeMAdiff = KneeMA-KneeMAmean;
KneeMAdiffstd = std(KneeMAdiff,[],3);
KneeMAdiffstd = repmat(KneeMAdiffstd, [1 1 nTrials]);
outliers = outliers+permute(max(max(abs(KneeMAdiff)>5*KneeMAdiffstd,[],1),[],2), [3 1 2]);

%Remove trials with strange ankle moment arms
AnkleMAmean = mean(AnkleMA,3);
AnkleMAmean = repmat(AnkleMAmean, [1 1 nTrials]);

AnkleMAdiff = AnkleMA-AnkleMAmean;
AnkleMAdiffstd = std(AnkleMAdiff,[],3);
AnkleMAdiffstd = repmat(AnkleMAdiffstd, [1 1 nTrials]);
outliers = outliers+permute(max(max(abs(AnkleMAdiff)>5*AnkleMAdiffstd,[],1),[],2), [3 1 2]);

%Remove trials with strange subtalar moment arms
SubtalarMAmean = mean(SubtalarMA,3);
SubtalarMAmean = repmat(SubtalarMAmean, [1 1 nTrials]);

SubtalarMAdiff = SubtalarMA-SubtalarMAmean;
SubtalarMAdiffstd = std(SubtalarMAdiff,[],3);
SubtalarMAdiffstd = repmat(SubtalarMAdiffstd, [1 1 nTrials]);
outliers = outliers+permute(max(max(abs(SubtalarMAdiff)>5*SubtalarMAdiffstd,[],1),[],2), [3 1 2]);

outliers = logical(sign(outliers));

Time(:,outliers) = [];
TrialsSpeed(outliers) = [];
JAngles(:,:,outliers) = [];
JVels(:,:,outliers) = [];
HipMA(:,:,outliers) = [];
HipAAMA(:,:,outliers) = [];
KneeMA(:,:,outliers) = [];
AnkleMA(:,:,outliers) = [];
SubtalarMA(:,:,outliers) = [];
Lmt(:,:,outliers) = [];
Vmt(:,:,outliers) = [];
IDloads(:,:,outliers) = [];
EMG(:,:,outliers) = [];
nTrials = nTrials-sum(outliers);

Data.Time = Time;
Data.TrialsSpeed = TrialsSpeed;
Data.JAngles = JAngles;
Data.JVels = JVels;
Data.HipMA = HipMA;
Data.HipAAMA = HipAAMA;
Data.KneeMA = KneeMA;
Data.AnkleMA = AnkleMA;
Data.SubtalarMA = SubtalarMA;
Data.Lmt = Lmt;
Data.Vmt = Vmt;
Data.IDloads = IDloads;
Data.EMG = EMG;
Data.nTrials = nTrials;
