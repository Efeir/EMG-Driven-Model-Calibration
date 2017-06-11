function MuscleModelOptimization_Release_v1

close all

%% Settings
MakePredictions = 0; % Changes input trials to verification set
numTrialsPerSpeed = 10; % number of trials used for each speed
LoadOldResults = 1; % load previous results and use as initial guess or
LoadFileName = 'paramsOpt_WGA_45678.mat'; % load old model parameter values
RunOpt = 1; % run optimization, 0 will show results for initial guess/loaded results
Speeds = [0.4 0.5 0.6 0.7 0.8]; % which speeds will be loaded
GeometricAdjustments = 1; % Turn geometric adjustments on or off
StepWiseOptimizations = 1; % split optimization into smaller optimizations for better convergence
SaveData = 1; % save data !!!Will overwrite old file if name is not changed and optimization is ran!!!
SaveNameParams = 'paramsOpt_WGA_45678.mat'; % only model parameters get saved here
SaveNameResults = 'results_WGA_45678.mat'; % many results from the optimization get saved here
makeParameterChangePlotsSwitch = 0; % activates design variable value plots
makeIndividualMuscMomentPlotsSwitch = 0; % activates individual muscle moment plots
makeMomentMatchingPlotsSwitch = 1; % activates moment matching plots
makeMuscleParameterPlotsSwitch = 0; % activates lmt, vmt, moment arm plots
makeExcitationsPlotsSwitch = 0; % activates excitation plots
hpf = 40; % high pass filter cutoff for the EMG data
lpf = 3.5; % low pass filter cutoff for the EMG data, will be divided by the period of the gait cycle
numIter = 50; % iteration limit for the stepwise optimizations
numIterLarge = 150; % iteration limit for the final overall optimization
algorithmfmincon = 'sqp'; % optimizer type for fmincon
nJoints = 5; % number of joints in model
Leg = 'both'; % which leg is being optimized, right, left, or both
SampleStep = 1; % resamples data to ever SampleStep points using 1:SampleStep:end. Must be 1, 2, 5 or 10.
nptsLong = (141-1)/SampleStep+1; % number of sample points in gait cycle with padded beginnings and ends
nptsShort = (101-1)/SampleStep+1; % number of sample points in gait cycle without padding

%% Load Data

% Load processed trial data, includes raw EMG, joint angles, joint
% velocities, etc
load etmData.mat
% Load initial guesses for the muscle fits
load MuscleFits_lhs.mat
% Load passive muscle forces and associated joint angles from Silder et
% al. 2007 (with Thelen as corresponding author)
load DigitizedPassiveMoments.mat
JAnglesThelen = JAngles;
PassiveMThelen = PassiveM;

%% Initialize Data

% Compile matrices of trials
Data = [];
if MakePredictions
    if strcmpi(Leg,'both')
        if max((Speeds-0.2)==0)
            Data = [Data etmData_l.TMGait_0pt2(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.3)==0)
            Data = [Data etmData_l.TMGait_0pt3(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.4)==0)
            Data = [Data etmData_l.TMGait_0pt4(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.5)==0)
            Data = [Data etmData_l.TMGait_0pt5(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.6)==0)
            Data = [Data etmData_l.TMGait_0pt6(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.7)==0)
            Data = [Data etmData_l.TMGait_0pt7(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.8)==0)
            Data = [Data etmData_l.TMGait_0pt8(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.2)==0)
            Data = [Data etmData_r.TMGait_0pt2(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.3)==0)
            Data = [Data etmData_r.TMGait_0pt3(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.4)==0)
            Data = [Data etmData_r.TMGait_0pt4(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.5)==0)
            Data = [Data etmData_r.TMGait_0pt5(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.6)==0)
            Data = [Data etmData_r.TMGait_0pt6(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.7)==0)
            Data = [Data etmData_r.TMGait_0pt7(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.8)==0)
            Data = [Data etmData_r.TMGait_0pt8(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
    elseif strcmpi(Leg,'left')
        if max((Speeds-0.2)==0)
            Data = [Data etmData_l.TMGait_0pt2(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.3)==0)
            Data = [Data etmData_l.TMGait_0pt3(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.4)==0)
            Data = [Data etmData_l.TMGait_0pt4(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.5)==0)
            Data = [Data etmData_l.TMGait_0pt5(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.6)==0)
            Data = [Data etmData_l.TMGait_0pt6(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.7)==0)
            Data = [Data etmData_l.TMGait_0pt7(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.8)==0)
            Data = [Data etmData_l.TMGait_0pt8(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
    elseif strcmpi(Leg,'right')
        if max((Speeds-0.2)==0)
            Data = [Data etmData_r.TMGait_0pt2(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.3)==0)
            Data = [Data etmData_r.TMGait_0pt3(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.4)==0)
            Data = [Data etmData_r.TMGait_0pt4(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.5)==0)
            Data = [Data etmData_r.TMGait_0pt5(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.6)==0)
            Data = [Data etmData_r.TMGait_0pt6(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.7)==0)
            Data = [Data etmData_r.TMGait_0pt7(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
        if max((Speeds-0.8)==0)
            Data = [Data etmData_r.TMGait_0pt8(numTrialsPerSpeed+1:2*numTrialsPerSpeed)];
        end
    end
else
    if strcmpi(Leg,'both')
        if max((Speeds-0.2)==0)
            Data = [Data etmData_l.TMGait_0pt2(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.3)==0)
            Data = [Data etmData_l.TMGait_0pt3(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.4)==0)
            Data = [Data etmData_l.TMGait_0pt4(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.5)==0)
            Data = [Data etmData_l.TMGait_0pt5(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.6)==0)
            Data = [Data etmData_l.TMGait_0pt6(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.7)==0)
            Data = [Data etmData_l.TMGait_0pt7(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.8)==0)
            Data = [Data etmData_l.TMGait_0pt8(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.2)==0)
            Data = [Data etmData_r.TMGait_0pt2(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.3)==0)
            Data = [Data etmData_r.TMGait_0pt3(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.4)==0)
            Data = [Data etmData_r.TMGait_0pt4(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.5)==0)
            Data = [Data etmData_r.TMGait_0pt5(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.6)==0)
            Data = [Data etmData_r.TMGait_0pt6(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.7)==0)
            Data = [Data etmData_r.TMGait_0pt7(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.8)==0)
            Data = [Data etmData_r.TMGait_0pt8(1:numTrialsPerSpeed)];
        end
    elseif strcmpi(Leg,'left')
        if max((Speeds-0.2)==0)
            Data = [Data etmData_l.TMGait_0pt2(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.3)==0)
            Data = [Data etmData_l.TMGait_0pt3(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.4)==0)
            Data = [Data etmData_l.TMGait_0pt4(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.5)==0)
            Data = [Data etmData_l.TMGait_0pt5(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.6)==0)
            Data = [Data etmData_l.TMGait_0pt6(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.7)==0)
            Data = [Data etmData_l.TMGait_0pt7(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.8)==0)
            Data = [Data etmData_l.TMGait_0pt8(1:numTrialsPerSpeed)];
        end
    elseif strcmpi(Leg,'right')
        if max((Speeds-0.2)==0)
            Data = [Data etmData_r.TMGait_0pt2(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.3)==0)
            Data = [Data etmData_r.TMGait_0pt3(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.4)==0)
            Data = [Data etmData_r.TMGait_0pt4(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.5)==0)
            Data = [Data etmData_r.TMGait_0pt5(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.6)==0)
            Data = [Data etmData_r.TMGait_0pt6(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.7)==0)
            Data = [Data etmData_r.TMGait_0pt7(1:numTrialsPerSpeed)];
        end
        if max((Speeds-0.8)==0)
            Data = [Data etmData_r.TMGait_0pt8(1:numTrialsPerSpeed)];
        end
    end
end
nTrialsSpeed = numTrialsPerSpeed;

% Determine number of gait speeds to be used in later optimizations
if strcmpi(Leg,'both')
    nSpeeds = 2*length(Speeds);
elseif strcmpi(Leg,'left')||strcmpi(Leg,'right')
    nSpeeds = length(Speeds);
end

% Specify proportion of trials to be used from left or right foot
if strcmpi(Leg,'both')
nTrials_l = 0.5*nSpeeds*numTrialsPerSpeed;
nTrials_r = 0.5*nSpeeds*numTrialsPerSpeed;
elseif strcmpi(Leg,'left')
nTrials_l = nSpeeds*numTrialsPerSpeed;
nTrials_r = 0;
elseif strcmpi(Leg,'right')
nTrials_l = 0;
nTrials_r = nSpeeds*numTrialsPerSpeed;
end

% Determine number of trials used in optimization
nTrials = length(Data);
% Specify number of muscles and EMG signals for each muscle
nMuscEMG = 16;
nMusc = 35;

% Append all individual muscle coefficients into a single dimension array form
coefsFit = MuscleFits.Coefs;
coefsOrig = [];
for i = 1:nMusc
    coefsOrig = [coefsOrig; coefsFit{i}];
end

% Number of coefficients used in the muscle geometry fits
ncoefsFit = length(coefsOrig);

% determine new number of trials and time frames
nFrames = nptsShort*nTrials;

% print number of trials being used to command line
fprintf(['\nUsing ' num2str(nTrials) ' trials for model calibration.\n']);

% Preallocate memory
TrialsSpeed = zeros(1,nTrials); % A matrix which indicates the gait speed of each trial
Time = zeros(nptsLong,nTrials); % Time vector for each trial
JAngles = zeros(nptsShort,6,nTrials); % joint angles for each trial in radians
JVels = zeros(nptsShort,6,nTrials); % joint velocities for each trial in radians/s
HipFEMA = zeros(nptsShort,nMusc,nTrials); % Hip FE moment arms in meters
HipAAMA = zeros(nptsShort,nMusc,nTrials); % Hip AA moment arms in meters
KneeMA = zeros(nptsShort,nMusc,nTrials); % Knee FE moment arms in meters
AnkleMA = zeros(nptsShort,nMusc,nTrials); % Ankle moment arms in meters
SubtalarMA = zeros(nptsShort,nMusc,nTrials); % Subtalar moment arms in meters
Lmt = zeros(nptsShort,nMusc,nTrials); % Muscle tendon lengths in meters
Vmt = zeros(nptsShort,nMusc,nTrials); % Muscle tendon velocities in meters/s
IDloads = zeros(nptsShort,5,nTrials); % Inverse dynamics loads in N-m
TrialName = cell(nTrials,1);

% Store data in new format (Time Frames)x(Muscle or DOF)x(Trial)
k = 1;
for j = 1:nTrials
    
    Time(:,j) = Data(j).Time(1:SampleStep:141);
    JAngles(:,:,j) = Data(j).JointAngles(21:SampleStep:121,1:6)*pi/180;
    JVels(:,:,j) = Data(j).JointVelocities(21:SampleStep:121,1:6)*pi/180;
    HipFEMA(:,:,j) = Data(j).HipFEMomentArms(21:SampleStep:121,:);
    HipAAMA(:,:,j) = Data(j).HipAAMomentArms(21:SampleStep:121,:);
    KneeMA(:,:,j) = Data(j).KneeMomentArms(21:SampleStep:121,:);
    AnkleMA(:,:,j) = Data(j).AnkleMomentArms(21:SampleStep:121,:);
    SubtalarMA(:,:,j) = Data(j).SubtalarMomentArms(21:SampleStep:121,:);
    Lmt(:,:,j) = Data(j).MuscleTendonLengths(21:SampleStep:121,:);
    Vmt(:,:,j) = Data(j).MuscleTendonVelocities(21:SampleStep:121,:);
    IDloads(:,:,j) = Data(j).InverseDynamicsLoads(21:SampleStep:121,:);
    TrialName{j} = Data(j).TrialName;
    
end

% store original curves for comparison purposes
Lmtorig = Lmt;
HipMAorig = HipFEMA;
HipAAMAorig = HipAAMA;
KneeMAorig = KneeMA;
AnkleMAorig = AnkleMA;
SubtalarMAorig = SubtalarMA;

%% Set muscle similarity parameters

% Set similarity pairs that keep parameters for similar muscles close
% Each row is a pair of muscles
optParams.ActivationPairs = ... %inclues time delay, activation/deactivation
    {
    1:6 % Adductors
    7:9 % Gluteus maximumus
    10:12 % Gluteus medius
    13:15 % Gluteus minimus
    16:17 % Iliacus and Psoas
    18:19 % Semimem/semiten
    20:21 % Biceps Fem lh/sh
    22:25 % Quads
    26:27 % Gastrocs
    30:32 % Peroneals
    };

optParams.lmtildaPairs = ... % penalizes differences between lmtilda curves
    {
    7:9 % Gluteus maximumus
    10:12 % Gluteus medius
    13:15 % Gluteus minimus
    23:25 % Quads (not rec fem)
    26:27 % Gastrocs
    };

optParams.HipMAPairs = cell(0,2); % penalizes differences between Hip FE MA curves
optParams.HipAAMAPairs = cell(0,2); % penalizes differences between Hip AA MA curves
optParams.KneeMAPairs = { % penalizes differences between Knee MA curves
    18:19 % Semimem/semiten
    20:21 % Biceps Fem lh/sh
    22:25}; % Gastrocs
optParams.AnkleMAPairs = {[26:27 33] % penalizes differences between Ankle MA curves
    30:31}; % peroneals except tertius which has different moment arm
optParams.SubtalarMAPairs = {[26:27 33] % penalizes differences between Subtalar MA curves
    30:32}; % peroneals including tertius

%% Process the EMG

% Preallocate memory
EMG = zeros(nptsLong,nMuscEMG,nTrials);
TimeEMG = Time;

for j = 1:nTrials
    
    % determine time length from heel strike to heel strike (not including
    % padding on beginning and end of cycle)
    tf = (TimeEMG(end,j)-TimeEMG(1,j))*(nptsShort/nptsLong);
    lpfVar = lpf/tf;
    
    % catch a common step parsing error in preprocessing
    if tf >10
        keyboard
    end
    
    % process the EMG data
    EMGtrial = processEMG(Data(j).EMG, 1/1000, hpf, lpfVar);
    
    % find original final time of trial
    tf = (length(EMGtrial)-1)/1000;
    
    % create normalized emg signals time vectors
    TimeEMGtrial = 0:(1/1000):tf;
    TimeEMG(:,j) = linspace(0,tf,nptsLong);
    % Normalize to nptsLong time frames
    EMG(:,:,j) = spline(TimeEMGtrial, EMGtrial', TimeEMG(:,j))';
    
end

% Duplicate EMG signals for similar muscles
EMG = EMG(:,[ones(1,6) 2*ones(1,3) 3*ones(1,6) 4*ones(1,2) 5*ones(1,2) 6*ones(1,2) 7:9 9 10*ones(1,2) 11:12 13*ones(1,3) 14:16],:);

%% Normalize EMG

% Load EMG scale factors determined for all trials
load EMGmaxima_l.mat
load EMGmaxima_r.mat

% split EMG into left and right side
EMG_l = EMG(:,:,1:nTrials_l);
EMG_r = EMG(:,:,nTrials_l+1:end);

% Remove EMG offset for each trial
EMGmins_l = min(EMG_l(21:121,:,:),[],1);
EMG_l = EMG_l-EMGmins_l(ones(141,1),:,:);
EMGmins_r = min(EMG_r(21:121,:,:),[],1);
EMG_r = EMG_r-EMGmins_r(ones(141,1),:,:);

% activate to determine new scale factors for EMG data
if 0
    
    for i = 1:nMusc
        EMGmin_r(i) = min(min(EMG_r(21:121,i,:)));
        EMGmax_r(i) = max(max(EMG_r(:,i,:)-EMGmin_r(i)));
        EMG_r(:,i,:) = (EMG_r(:,i,:)-EMGmin_r(i))/EMGmax_r(i);
        EMGmin_l(i) = min(min(EMG_l(21:121,i,:)));
        EMGmax_l(i) = max(max(EMG_l(:,i,:)-EMGmin_l(i)));
        EMG_l(:,i,:) = (EMG_l(:,i,:)-EMGmin_l(i))/EMGmax_l(i);
    end
    
    save EMGmaxima_r.mat EMGmin_r EMGmax_r
    save EMGmaxima_l.mat EMGmin_l EMGmax_l
    
else
    % Normalize EMG data
    for i = 1:nMusc
        EMG_l(:,i,:) = (EMG_l(:,i,:)-EMGmin_l(i))/EMGmax_l(i);
        EMG_r(:,i,:) = (EMG_r(:,i,:)-EMGmin_r(i))/EMGmax_r(i);
    end
end

% Recombine EMG data
EMG = cat(3, EMG_l, EMG_r);

% remove EMG values less than 0
EMG(EMG<0) = 0;

%% Create EMG splines

optParams.EMG = permute(EMG, [1 3 2]);

EMGsplines = cell(nMusc,nTrials);
for i = 1:nMusc
    for j = 1:nTrials
        EMGsplines{i,j} = spline(Time(1:4:end,j),EMG(1:4:end,i,j)); % Use every other frame to improve speed when evaluating spline
     end
end
optParams.EMGsplines = EMGsplines;

%% Set initial muscle parameter values

% Muscle parameters from Arnold et al. 2010
lmo =   [10.3 10.8 17.7 15.6 13.8 10.6 14.7 15.7 16.7 7.3 7.3 7.3 6.8 5.6 3.8 10.7 11.7 6.9 19.3 9.8 11 7.6 9.7 9.9 9.9 5.9 5.1 6.8 3.8 4.5 5.1 7.9 4.4 6.9 4.5]/100;
lts =   [3.6 13.0 9.0 22.1 4.8 4.3 5.0 7.3 7.0 5.7 6.6 4.6 1.6 2.6 5.1 9.4 9.7 37.8 24.5 32.2 10.4 34.6 11.2 13 10.6 38.2 40.1 24.1 28.2 14.9 33.3 10.0 28.2 36.7 37.8]/100;
alpha = [6.1 7.1 13.8 11.9 14.7 22.2 21.1 21.9 22.8 20.5 20.5 20.5 10.0 0 1.0 14.3 10.7 15.1 12.9 11.6 12.3 13.9 29.6 18.4 4.5 12.0 9.9 9.6 13.7 11.5 14.1 13.0 28.3 10.8 13.6]*pi/180;

%% Setup muscle force estimation equations, equation from Handsfield et al. 2013
sigma = 610e3;
LegVol = (47*1.7*80.5+1285)/100^3;
VolFraction = [1.47 2.26...
    1.97 1.97 1.97 1.97 ...
    3.52 5.02 3.39 ...
    1.82 1.27 1.45...
    0.45 0.48 0.54...
    2.48 3.80 3.46 2.60 2.92 1.40 3.79 6.06 11.66 3.84 2.11 3.62 1.91 1.49...
    0.53 1.14 0.16 6.21...
    0.97 0.43];
VolFraction = VolFraction/sum(VolFraction);
Vmusc = LegVol*VolFraction;

%% Preconstruct moment calculation matrices to increase optimization speed

% define muscle reference matrix that defines which joint(s) a muscle
% actuates
MuscRef = zeros(1,nMusc);
MuscRef(1:17) = 1; % Hip FE and AA
MuscRef([18:20 22]) = 2; % Hip FE and AA and Knee FE
MuscRef([21 23:25]) = 3; % Knee FE
MuscRef([26:27]) = 4; % Knee FE, Ankle and Subtalar
MuscRef([28:35]) = 5; % Ankle and Subtalar

JAnglesReshaped = reshape(permute(JAngles, [1 3 2]), nptsShort*nTrials, 6);
JVelsReshaped = reshape(permute(JVels, [1 3 2]), nptsShort*nTrials, 6);

optParams.Mat = buildMuscleMatrices(JAnglesReshaped,JVelsReshaped);
optParams.MuscRef = MuscRef;

%% Setup matrices for calculating passive forces from thelen data

% Remove high knee flexion angle from hip moment matching
PassiveMThelen{1} = reshape(PassiveMThelen{1},numel(PassiveMThelen{1})/4,4);
PassiveMThelen{1}(:,4) = [];
PassiveMThelen{1} = PassiveMThelen{1}(:);

% Remove high knee flexion from knee moment matching
PassiveMThelen{2} = reshape(PassiveMThelen{2},numel(PassiveMThelen{2})/4,4);
PassiveMThelen{2}(92:121,:) = [];
PassiveMThelen{2} = PassiveMThelen{2}(:);

% Combine the passive moments
optParams.PassiveMThelen = [PassiveMThelen{1} zeros(length(PassiveMThelen{1}),4);...
    zeros(length(PassiveMThelen{2}),2) PassiveMThelen{2} zeros(length(PassiveMThelen{2}),2);...
    zeros(length(PassiveMThelen{3}),3) PassiveMThelen{3} zeros(length(PassiveMThelen{3}),1)];

% Remove high knee flexion angle from hip moment matching
JAnglesThelen{1} = reshape(JAnglesThelen{1},numel(JAnglesThelen{1})/(4*6),4,6);
JAnglesThelen{1}(:,4,:) = [];
JAnglesThelen{1} = reshape(JAnglesThelen{1},numel(JAnglesThelen{1})/(6),6);

% Remove high knee flexion from knee moment matching
JAnglesThelen{2} = reshape(JAnglesThelen{2},numel(JAnglesThelen{2})/(4*6),4,6);
JAnglesThelen{2}(92:121,:,:) = [];
JAnglesThelen{2} = reshape(JAnglesThelen{2},numel(JAnglesThelen{2})/(6),6);

%Combine the joint angles
JAnglesThelen = [JAnglesThelen{1}; JAnglesThelen{2}; JAnglesThelen{3}]*pi/180;

% Add offset to the angles to account for differences in models
JAnglesThelen(:,1) = JAnglesThelen(:,1)-0*pi/180;
JAnglesThelen(:,2) = JAnglesThelen(:,2)-0*pi/180;
JAnglesThelen(:,4) = JAnglesThelen(:,4)+0*pi/180;

JVelsThelen = zeros(size(JAnglesThelen));

nFramesThelen = length(JAnglesThelen);
nFramesAllThelen = nFramesThelen;

optParams.aThelen = zeros(nFramesAllThelen,35);
optParams.onesColThelen = ones(nFramesAllThelen,1);
optParams.nframesAllThelen = nFramesAllThelen;

optParams.MatThelen = buildMuscleMatrices(JAnglesThelen, JVelsThelen);

%% Assign values to optimization parameters structure

% Assign values to optimization structure
optParams.Time = Time;
optParams.JAngles = reshape(permute(JAngles, [1 3 2]), nptsShort*nTrials, 6);
optParams.JVels = reshape(permute(JVels, [1 3 2]), nptsShort*nTrials, 6);
optParams.HipMA = reshape(permute(HipFEMA, [1 3 2]), nptsShort*nTrials, nMusc);
optParams.HipAAMA = reshape(permute(HipAAMA, [1 3 2]), nptsShort*nTrials, nMusc);
optParams.KneeMA = reshape(permute(KneeMA, [1 3 2]), nptsShort*nTrials, nMusc);
optParams.AnkleMA = reshape(permute(AnkleMA, [1 3 2]), nptsShort*nTrials, nMusc);
optParams.SubtalarMA = reshape(permute(SubtalarMA, [1 3 2]), nptsShort*nTrials, nMusc);
optParams.Lmt = reshape(permute(Lmt, [1 3 2]), nptsShort*nTrials, nMusc);
optParams.Vmt = reshape(permute(Vmt, [1 3 2]), nptsShort*nTrials, nMusc);
optParams.IDloads = reshape(permute(IDloads, [1 3 2]), nptsShort*nTrials, 5);
optParams.lmo = lmo;
optParams.lts = lts;
optParams.alpha = alpha;
optParams.nframesAll = nptsShort*nTrials;
optParams.onesCol = ones(nptsShort*nTrials,1);
optParams.SampleStep = SampleStep;
optParams.nptsShort = nptsShort;
optParams.nptsLong = nptsLong;
optParams.coefsOrig = coefsOrig;
optParams.ncoefsFit = ncoefsFit;
optParams.nTrials_l = nTrials_l;
optParams.nTrials_r = nTrials_r;
optParams.numTrialsPerSpeed = numTrialsPerSpeed;
optParams.nSpeeds = nSpeeds;
optParams.MuscLabels = Data(1).OpensimMuscleLabels;
optParams.GeometricAdjustments = GeometricAdjustments;
optParams.Vmusc = Vmusc;
optParams.sigma = sigma;
optParams.nTrials = nTrials;
optParams.nMusc = nMusc;
optParams.nFrames = nFrames;

%% load precalibrated lmo and lts values
load scaledlmolts.mat

optParams.lmo = lmo;
optParams.lts = lts;

%% Run overall model optimization

% set initial guesses for optimization
if GeometricAdjustments
    guess = [0.5*ones(1,2*nMusc) 1.5*ones(1,1*nMusc) zeros(1,1*nMusc) .25*ones(1,2*nMusc) ones(1,2*nMusc) zeros(1,ncoefsFit)];
    upBounds = [1.25*ones(1,2*nMusc) 3.5*ones(1,1*nMusc) 0.35*ones(1,1*nMusc) 1*ones(1,2*nMusc) 1.25*ones(1,2*nMusc) 1000*ones(1,ncoefsFit)];
    lowBounds = [0*ones(1,2*nMusc) .75*ones(1,1*nMusc) zeros(1,1*nMusc) .05*ones(1,2*nMusc) .75*ones(1,2*nMusc) -1000*ones(1,ncoefsFit)];
else
    guess = [0.5*ones(1,2*nMusc) 1.5*ones(1,1*nMusc) zeros(1,1*nMusc) .25*ones(1,2*nMusc) ones(1,2*nMusc)];
    upBounds = [1.25*ones(1,2*nMusc) 3.5*ones(1,1*nMusc) 0.35*ones(1,1*nMusc) 1*ones(1,2*nMusc) 1.25*ones(1,2*nMusc)];
    lowBounds = [0*ones(1,2*nMusc) .75*ones(1,1*nMusc) zeros(1,1*nMusc) .05*ones(1,2*nMusc) .75*ones(1,2*nMusc)];
end

% setup optimization options
options = optimset('display','iter','MaxIter',numIter,'MaxFunEvals',100000,'Algorithm',algorithmfmincon,'UseParallel','always','Hessian','lbfgs');
optionsLarge = optimset('display','iter','MaxIter',numIterLarge,'MaxFunEvals',100000,'Algorithm',algorithmfmincon,'UseParallel','always','Hessian','lbfgs');

if LoadOldResults
    load(LoadFileName)
    guess = params;
end

% if runopt is set to zero, will output results from guess/loaded
% parameters
if RunOpt

    % Splits up optimization into three phases
    % 1: optimize EMG to activation parameters except 
    % 2: optimize lmo, lts
    % 3: optimize the muscle geometric coefficients
    if StepWiseOptimizations
        for j = 1:2
            optParams.otherParams = guess(:,6*nMusc+1:end);
            [paramsTemp, FVAL] = fmincon(@costFuncFMINCON1, guess(:,1:6*nMusc), [],[],[],[],...
                lowBounds(1:6*nMusc), upBounds(1:6*nMusc),[],options,optParams);
            guess(:,1:6*nMusc) = paramsTemp;
            
            if GeometricAdjustments
                optParams.otherParams = guess(:,[1:6*nMusc 8*nMusc+1:end]);
                [paramsTemp, FVAL] = fmincon(@costFuncFMINCON2, guess(:,6*nMusc+1:8*nMusc), [],[],[],[],...
                    lowBounds(6*nMusc+1:8*nMusc), upBounds(6*nMusc+1:8*nMusc),@nonlcon2,options,optParams);
                guess(6*nMusc+1:8*nMusc) = paramsTemp;
                
                optParams.otherParams = guess(:,1:8*nMusc);
                [paramsTemp, FVAL] = fmincon(@costFuncFMINCON3, guess(:,8*nMusc+1:8*nMusc+ncoefsFit), [],[],[],[],...
                    lowBounds(8*nMusc+1:8*nMusc+ncoefsFit), upBounds(8*nMusc+1:8*nMusc+ncoefsFit),@nonlcon3,options,optParams);
                guess(8*nMusc+1:8*nMusc+ncoefsFit) = paramsTemp;
            else
                optParams.otherParams = guess(:,[1:6*nMusc 8*nMusc+1:end]);
                [paramsTemp, FVAL] = fmincon(@costFuncFMINCON2, guess(:,6*nMusc+1:8*nMusc), [],[],[],[],...
                    lowBounds(6*nMusc+1:8*nMusc), upBounds(6*nMusc+1:8*nMusc),@nonlcon2,options,optParams);
                guess(6*nMusc+1:8*nMusc) = paramsTemp;
            end
        end
    end
    
    % Run final optimization that includes all design variables
    [params, FVAL] = fmincon(@costFuncFMINCON, guess, [],[],[],[],...
        lowBounds, upBounds, @nonlcon, optionsLarge, optParams);
end

%% Save Data

if SaveData && RunOpt
    save(SaveNameParams, 'params');
end

%% Extract new parameter values

% extract values from parameters matrix
tdelay = params(1,1:2*nMusc)/10;
tact = params(1,(1:nMusc)+2*nMusc)/100;
tdeact = 4*tact;
AnonlinOpt = params(1,(1:nMusc)+3*nMusc);
EMGScaleOpt = params(1,(1:2*nMusc)+4*nMusc);
VmaxFactor = 10;
lmoscaleOpt = params(1,(1:nMusc)+6*nMusc);
ltsscaleOpt = params(1,(1:nMusc)+7*nMusc);
if GeometricAdjustments
    coefsOffsetOpt = params(1,(1:ncoefsFit)+8*nMusc)/10;
else
    coefsOffsetOpt = 0;
end

% adjust old muscle parameters values with optimized scale factors/offsets
coefsNew = coefsOrig+coefsOffsetOpt';
ltsOpt = optParams.lts.*ltsscaleOpt;
lmoOpt = optParams.lmo.*lmoscaleOpt;

% add new parameters to optimization structure
optParams.VmaxFactor = VmaxFactor;
optParams.lmoscale = lmoscaleOpt;
optParams.ltsscale = ltsscaleOpt;
optParams.EMGScale = EMGScaleOpt;
optParams.tdelay = tdelay;
optParams.tact = tact;
optParams.tdeact = tdeact;
optParams.Anonlin = AnonlinOpt;

%% Make design variable value plots

% make plots of parameter changes/values
if makeParameterChangePlotsSwitch
    makeParameterChangePlots(optParams);
end

%% Calculate new moments and muscle parameter values

[a, Excitations] = EMGtoAct(optParams);

optParams.a = a;
optParams.lts = ltsOpt;
optParams.lmo = lmoOpt;
optParams.coefs = coefsNew;

[moments, lmtilda, vmtilda, HipFEMA, HipAAMA, KneeMA, AnkleMA, SubtalarMA, ...
    muscMoments, passiveF, Lmt, Vmt, muscleForce] = calcMoments(optParams,1);

% Calculate passive forces at each joint during gait
PassiveHipFEM = sum(passiveF.*HipFEMA,2);
PassiveHipAAM = sum(passiveF.*HipAAMA,2);
PassiveKneeFEM = sum(passiveF.*KneeMA,2);
PassiveAnklePFM = sum(passiveF.*AnkleMA,2);
PassiveAnkleIEM = sum(passiveF.*SubtalarMA,2);

PassiveMGait = [PassiveHipFEM,PassiveHipAAM,PassiveKneeFEM,PassiveAnklePFM,PassiveAnkleIEM];
PassiveMGait = reshape(PassiveMGait, nptsShort,nTrials,nJoints);

optParams.PassiveMGait = PassiveMGait;

%% Calculate errors between model and inverse dynamics

% Calculate RMS errors
HipRMS = sqrt(sum((moments(:,1)-optParams.IDloads(:,1)).^2)/(optParams.nframesAll-1))
HipAARMS = sqrt(sum((moments(:,2)-optParams.IDloads(:,2)).^2)/(optParams.nframesAll-1))
KneeRMS = sqrt(sum((moments(:,3)-optParams.IDloads(:,3)).^2)/(optParams.nframesAll-1))
AnkleRMS = sqrt(sum((moments(:,4)-optParams.IDloads(:,4)).^2)/(optParams.nframesAll-1))
SubtalarRMS = sqrt(sum((moments(:,5)-optParams.IDloads(:,5)).^2)/(optParams.nframesAll-1))

% combine errros into single matrix
RMS = [HipRMS HipAARMS KneeRMS AnkleRMS SubtalarRMS];

% Calculate mean absolute errors
HipMAE = sum(abs(moments(:,1)-optParams.IDloads(:,1))/optParams.nframesAll)
HipAAMAE = sum(abs(moments(:,2)-optParams.IDloads(:,2))/optParams.nframesAll)
KneeMAE = sum(abs(moments(:,3)-optParams.IDloads(:,3))/optParams.nframesAll)
AnkleMAE = sum(abs(moments(:,4)-optParams.IDloads(:,4))/optParams.nframesAll)
SubtalarMAE = sum(abs(moments(:,5)-optParams.IDloads(:,5))/optParams.nframesAll)

% combine errros into single matrix
MAEs = [HipMAE HipAAMAE KneeMAE AnkleMAE SubtalarMAE];

% reshape/permute moments and IDloads into form (time
% frames)/(trials)/(joint)
moments = reshape(moments,nptsShort,nTrials,5);
IDloads = permute(IDloads, [1 3 2]);

% Calculate median errors instead of errors over all trials
HipMAEMed = median(sum(abs(moments(:,:,1)-IDloads(:,:,1)))/nptsShort)
HipAAMAEMed = median(sum(abs(moments(:,:,2)-IDloads(:,:,2)))/nptsShort)
KneeMAEMed = median(sum(abs(moments(:,:,3)-IDloads(:,:,3)))/nptsShort)
AnkleMAEMed = median(sum(abs(moments(:,:,4)-IDloads(:,:,4)))/nptsShort)
SubtalarMAEMed = median(sum(abs(moments(:,:,5)-IDloads(:,:,5)))/nptsShort)

% combine errros into single matrix
MAEMed = [HipMAEMed HipAAMAEMed KneeMAEMed AnkleMAEMed SubtalarMAEMed];

% Calculate max and min values for lmtilda and vmtilda for inspection
maxlmtilda = max(lmtilda)
minlmtilda = min(lmtilda)
maxvmtilda = max(vmtilda)
minvmtilda = min(vmtilda)

%% Make moment matching plots

% for plotting individual muscle contributions to moments
if makeIndividualMuscMomentPlotsSwitch
    makeIndividualMuscMomentPlots(moments,muscMoments,optParams);
end

if makeMomentMatchingPlotsSwitch
    makeMomentMatchingPlots(moments,IDloads,optParams);
end

%% Reshape data into form (time frames)/(trials)/(joint or muscle)

muscMoments = reshape(muscMoments, nptsShort,nTrials,nMusc,5);
passiveF = reshape(passiveF,nptsShort,nTrials,nMusc);
Lmt = reshape(Lmt,nptsShort,nTrials,nMusc);
muscleForce = reshape(muscleForce,nptsShort,nTrials,nMusc);

Lmt = reshape(Lmt,nptsShort,nTrials,nMusc);
lmtilda = reshape(lmtilda,nptsShort,nTrials,nMusc);
vmtilda = reshape(vmtilda,nptsShort,nTrials,nMusc);
HipFEMA = reshape(HipFEMA,nptsShort,nTrials,nMusc);
HipAAMA = reshape(HipAAMA,nptsShort,nTrials,nMusc);
KneeMA = reshape(KneeMA,nptsShort,nTrials,nMusc);
AnkleMA = reshape(AnkleMA,nptsShort,nTrials,nMusc);
SubtalarMA = reshape(SubtalarMA,nptsShort,nTrials,nMusc);
OpensimMuscleLabels = optParams.MuscLabels;

%% Make muscle parameter (lmtilda lmt etc...) plots

if makeMuscleParameterPlotsSwitch
    
    makeMuscleParameterPlots(Lmt, lmtilda, vmtilda, HipFEMA, HipAAMA, KneeMA, AnkleMA, SubtalarMA,Lmtorig,HipMAorig,HipAAMAorig,KneeMAorig,AnkleMAorig,SubtalarMAorig, optParams)
    
end

%% Make muscle excitation plots

if makeExcitationsPlotsSwitch
    makeExcitationPlots(Excitations,a,optParams);
end

%% Calculate R values between moment arms and lmt for each muscle

Lmtorigstacked = reshape(permute(Lmtorig,[1 3 2]),nFrames,35);
Lmtstacked = reshape(permute(Lmt,[1 2 3]),nFrames,35);
HipMAorigstacked = reshape(permute(HipMAorig,[1 3 2]),nFrames,35);
HipMAstacked = reshape(permute(HipFEMA,[1 2 3]),nFrames,35);
HipAAMAorigstacked = reshape(permute(HipAAMAorig,[1 3 2]),nFrames,35);
HipAAMAstacked = reshape(permute(HipAAMA,[1 2 3]),nFrames,35);
KneeMAorigstacked = reshape(permute(KneeMAorig,[1 3 2]),nFrames,35);
KneeMAstacked = reshape(permute(KneeMA,[1 2 3]),nFrames,35);
AnkleMAorigstacked = reshape(permute(AnkleMAorig,[1 3 2]),nFrames,35);
AnkleMAstacked = reshape(permute(AnkleMA,[1 2 3]),nFrames,35);
SubtalarMAorigstacked = reshape(permute(SubtalarMAorig,[1 3 2]),nFrames,35);
SubtalarMAstacked = reshape(permute(SubtalarMA,[1 2 3]),nFrames,35);

for i = 1:nMusc
    R = corrcoef(Lmtstacked(:,i),Lmtorigstacked(:,i));
    RLmt(i) = R(2);
    R = corrcoef(HipMAstacked(:,i),HipMAorigstacked(:,i));
    RHip(i) = R(2);
    R = corrcoef(HipAAMAstacked(:,i),HipAAMAorigstacked(:,i));
    RHipAA(i) = R(2);
    R = corrcoef(KneeMAstacked(:,i),KneeMAorigstacked(:,i));
    RKnee(i) = R(2);
    R = corrcoef(AnkleMAstacked(:,i),AnkleMAorigstacked(:,i));
    RAnkle(i) = R(2);
    R = corrcoef(SubtalarMAstacked(:,i),SubtalarMAorigstacked(:,i));
    RSubtalar(i) = R(2);
end

%% Calculate changes in moment arm values

HipMAorig = permute(HipMAorig, [1 3 2]);
HipAAMAorig = permute(HipAAMAorig, [1 3 2]);
KneeMAorig = permute(KneeMAorig, [1 3 2]);
AnkleMAorig = permute(AnkleMAorig, [1 3 2]);
SubtalarMAorig = permute(SubtalarMAorig, [1 3 2]);

[HipMADiff, index] = max(max(max(abs((HipFEMA-HipMAorig)))));
fprintf(['\nThe ' OpensimMuscleLabels{index} ' hip FE moment arm max change was ' num2str(HipMADiff) '.'])
[HipMADiff, index] = max(mean(mean(abs((HipFEMA-HipMAorig)))));
fprintf(['\nThe ' OpensimMuscleLabels{index} ' hip FE moment arm mean change was ' num2str(HipMADiff) '.\n'])

[HipAAMADiff, index] = max(max(max((HipAAMA-HipAAMAorig))));
fprintf(['\nThe ' OpensimMuscleLabels{index} ' hip AA moment arm max change was ' num2str(HipAAMADiff) '.'])
[HipAAMADiff, index] = max(mean(mean(abs((HipAAMA-HipAAMAorig)))));
fprintf(['\nThe ' OpensimMuscleLabels{index} ' hip AA moment arm mean change was ' num2str(HipAAMADiff) '.\n'])

[KneeMADiff, index] = max(max(max(abs((KneeMA-KneeMAorig)))));
fprintf(['\nThe ' OpensimMuscleLabels{index} ' knee FE moment arm max change was ' num2str(KneeMADiff) '.'])
[KneeMADiff, index] = max(mean(mean(abs((KneeMA-KneeMAorig)))));
fprintf(['\nThe ' OpensimMuscleLabels{index} ' knee FE moment arm mean change was ' num2str(KneeMADiff) '.\n'])

[AnkleMADiff, index] = max(max(max(abs((AnkleMA-AnkleMAorig)))));
fprintf(['\nThe ' OpensimMuscleLabels{index} ' ankle PF moment arm max change was ' num2str(AnkleMADiff) '.'])
[AnkleMADiff, index] = max(mean(mean(abs((AnkleMA-AnkleMAorig)))));
fprintf(['\nThe ' OpensimMuscleLabels{index} ' ankle PF moment arm mean change was ' num2str(AnkleMADiff) '.\n'])

[SubtalarMADiff, index] = max(max(max(abs((SubtalarMA-SubtalarMAorig)))));
fprintf(['\nThe ' OpensimMuscleLabels{index} ' subtalar IE moment arm max change was ' num2str(SubtalarMADiff) '.'])
[SubtalarMADiff, index] = max(mean(mean(abs((SubtalarMA-SubtalarMAorig)))));
fprintf(['\nThe ' OpensimMuscleLabels{index} ' subtalar IE moment arm mean change was ' num2str(SubtalarMADiff) '.\n'])

Lmtorig = permute(Lmtorig, [1 3 2]);
HipMAorig = permute(HipMAorig, [1 3 2]);
HipAAMAorig = permute(HipAAMAorig, [1 3 2]);
KneeMAorig = permute(KneeMAorig, [1 3 2]);
AnkleMAorig = permute(AnkleMAorig, [1 3 2]);
SubtalarMAorig = permute(SubtalarMAorig, [1 3 2]);

OpensimMuscleLabels = optParams.MuscLabels

%% Save all results

k = 1;
for i = 1:nMusc
    if MuscRef(i) == 1
        GeometricCoefs{i} = coefsNew(k:k+18);
        k = k+19;
    elseif MuscRef(i) == 2
        GeometricCoefs{i} = coefsNew(k:k+21);
        k = k+22;
    elseif MuscRef(i) == 3
        GeometricCoefs{i} = coefsNew(k:k+3);
        k = k+4;
    elseif MuscRef(i) == 4
        GeometricCoefs{i} = coefsNew(k:k+12);
        k = k+13;
    elseif MuscRef(i) == 5
        GeometricCoefs{i} = coefsNew(k:k+9);
        k = k+10;
    end
end

% a = reshape(a,101,100,35);

save(SaveNameResults, 'tdelay', 'tact', 'AnonlinOpt', 'lmoscaleOpt', 'ltsscaleOpt',...
    'EMGScaleOpt', 'OpensimMuscleLabels', 'Lmtorig', 'HipMAorig', 'HipAAMAorig',...
    'KneeMAorig', 'AnkleMAorig', 'SubtalarMAorig', 'moments', 'IDloads', 'lmtilda',...
    'vmtilda', 'HipFEMA', 'HipAAMA', 'KneeMA', 'AnkleMA', 'SubtalarMA', 'muscMoments',...
    'passiveF', 'Lmt', 'Vmt', 'muscleForce', 'JAngles','a','Excitations','alpha','lmoOpt','ltsOpt','GeometricCoefs');

function cost = costFuncMain(optParams)
% The main body of the cost function which is called from all other cost
% functions

% get needed values from optParams matrix
nMusc = optParams.nMusc;
IDloads = optParams.IDloads;
nframesAll = optParams.nframesAll;
onesCol = optParams.onesCol;
tact = optParams.tact;
Anonlin = optParams.Anonlin;
ltsscale = optParams.ltsscale;
lmoscale = optParams.lmoscale;
EMGScale = optParams.EMGScale;
tdelay = optParams.tdelay;
nFrames = optParams.nFrames;
nTrials = optParams.nTrials;

% calculate activations and excitations
[a,Excitations] = EMGtoAct(optParams);

optParams.a = a;

% calculate model moments and other muscle parameters
[ModelMoments, lmtilda, ~, HipMA, HipAAMA, KneeMA, AnkleMA, SubtalarMA, ~, passiveF,Lmt] = calcMoments(optParams,1);

% calculate passive moments to match with thelen curves (activation is
% zero)
optParams.a = optParams.aThelen;
optParams.Mat = optParams.MatThelen;
optParams.onesCol = optParams.onesColThelen;
optParams.nframesAll = optParams.nframesAllThelen;

[PassiveM] = calcMoments(optParams,1);

% remove joint that thelen has no data for (hip aa and subtalar)
PassiveM = PassiveM(optParams.PassiveMThelen(:,:)~=0);
PassiveMThelen = optParams.PassiveMThelen(optParams.PassiveMThelen~=0);

% calculate moment errors
costID = (ModelMoments-IDloads);%; LowBoundPenalty(:); UpperBoundPenalty(:)];
costPassive = [(PassiveM-PassiveMThelen)];

% extract muscle pairing information
ActivationPairs = optParams.ActivationPairs;
lmtildaPairs = optParams.lmtildaPairs;
HipMAPairs = optParams.HipMAPairs;
HipAAMAPairs = optParams.HipAAMAPairs;
KneeMAPairs = optParams.KneeMAPairs;
AnkleMAPairs = optParams.AnkleMAPairs;
SubtalarMAPairs = optParams.SubtalarMAPairs;

% Costs for violating activation similarity constraints
ActivationPenalty = 0;
for i = 1:length(ActivationPairs)
    ActivationPenalty = ActivationPenalty+std(EMGScale(ActivationPairs{i}))+std(EMGScale(ActivationPairs{i}+nMusc));
end

% Costs for violating time delay similarity constraints
tdelayPenalty = 0;
for i = 1:length(ActivationPairs)
    tdelayPenalty = tdelayPenalty+std(tdelay(ActivationPairs{i}))+std(tdelay(ActivationPairs{i}+nMusc));
end

% Costs for violating lmtilda similarity constraints
lmtildaPenalty = 0;
for i = 1:length(lmtildaPairs)
    lmtildaPenalty = lmtildaPenalty+std(lmtilda(:,lmtildaPairs{i}),0,2)./((nframesAll)^0.5);
end

% Costs for violating moment arm similarity constraints
MAPenalty = 0;
for i = 1:length(HipMAPairs)
    MAPenalty = MAPenalty+std(HipMA(:,HipMAPairs{i}),[],2)./((nframesAll)^0.5);
end

for i = 1:length(HipAAMAPairs)
    MAPenalty = MAPenalty+std(HipAAMA(:,HipAAMAPairs{i}),[],2)./((nframesAll)^0.5);
end

for i = 1:length(KneeMAPairs)
    MAPenalty = MAPenalty+std(KneeMA(:,KneeMAPairs{i}),[],2)./((nframesAll)^0.5);
end

for i = 1:length(AnkleMAPairs)
    MAPenalty = MAPenalty+std(AnkleMA(:,AnkleMAPairs{i}),[],2)./((nframesAll)^0.5);
end

for i = 1:length(SubtalarMAPairs)
    MAPenalty = MAPenalty+std(SubtalarMA(:,SubtalarMAPairs{i}),[],2)./((nframesAll)^0.5);
end

% cost for having moment arms and muscle lengths change too much from
% initial values
LmtMean = mean(Lmt);
LmtMeanOrig = mean(optParams.Lmt);
HipMAMean = mean(HipMA);
HipMAMeanOrig = mean(optParams.HipMA);
HipAAMAMean = mean(HipAAMA);
HipAAMAMeanOrig = mean(optParams.HipAAMA);
KneeMAMean = mean(KneeMA);
KneeMAMeanOrig = mean(optParams.KneeMA);
AnkleMAMean = mean(AnkleMA);
AnkleMAMeanOrig = mean(optParams.AnkleMA);
SubtalarMAMean = mean(SubtalarMA);
SubtalarMAMeanOrig = mean(optParams.SubtalarMA);

% Calculate max absolute value of muscle geometry curves for normalization
LmtRange = range(optParams.Lmt);
HipMARange = range(optParams.HipMA);
HipAAMARange = range(optParams.HipAAMA);
KneeMARange = range(optParams.KneeMA);
AnkleMARange = range(optParams.AnkleMA);
SubtalarMARange = range(optParams.SubtalarMA);

% Penalize differences in mean values of moment arms
LmtPenMean = (LmtMean-LmtMeanOrig)./LmtRange;
HipMAPenMean = (HipMAMean-HipMAMeanOrig)./HipMARange;
HipAAMAPenMean = (HipAAMAMean-HipAAMAMeanOrig)./HipAAMARange;
KneeMAPenMean = (KneeMAMean-KneeMAMeanOrig)./KneeMARange;
AnkleMAPenMean = (AnkleMAMean-AnkleMAMeanOrig)./AnkleMARange;
SubtalarMAPenMean = (SubtalarMAMean-SubtalarMAMeanOrig)./SubtalarMARange;

% Penalize changes in shape of moment arms
LmtPen = (Lmt-optParams.Lmt+onesCol*(LmtMeanOrig-LmtMean))./(onesCol*LmtRange);
HipMAPen = (HipMA-optParams.HipMA+onesCol*(HipMAMeanOrig-HipMAMean))./(onesCol*HipMARange);
HipAAMAPen = (HipAAMA-optParams.HipAAMA+onesCol*(HipAAMAMeanOrig-HipAAMAMean))./(onesCol*HipAAMARange);
KneeMAPen = (KneeMA-optParams.KneeMA+onesCol*(KneeMAMeanOrig-KneeMAMean))./(onesCol*KneeMARange);
AnkleMAPen = (AnkleMA-optParams.AnkleMA+onesCol*(AnkleMAMeanOrig-AnkleMAMean))./(onesCol*AnkleMARange);
SubtalarMAPen = (SubtalarMA-optParams.SubtalarMA+onesCol*(SubtalarMAMeanOrig-SubtalarMAMean))./(onesCol*SubtalarMARange);

% combine all costs into single vector
cost = [costID(:)/(sqrt(1)*(5*nframesAll)^0.5); ... % moment matching cost, time varying
    (tact(:)-.015)/(.015*nMusc.^0.5);... % activation value costs, time invariant
    ActivationPenalty(:);... % activation pair cost, time invariant
    10*tdelayPenalty(:);... % time delay pair cost, time invariant
    10*tdelay(:)/((2*nMusc).^0.5);
    sqrt(10)*Anonlin(:)/(0.25*nMusc^(0.5));... % anonlin cost, time invariant
    (ltsscale(:)-1)/(.25*nMusc^(0.5));... % lts change cost, time invariant
    (lmoscale(:)-1)/(.25*nMusc^(0.5)); ... % lmo change cost, time invariant
    LmtPenMean(:)/(1*(nMusc).^0.5); ... % lmt mean change cost, time invariant
    HipMAPenMean(:)/(0.5*(nMusc).^0.5); ... % Hip MA mean change cost, time invariant
    HipAAMAPenMean(:)/(0.5*(nMusc).^0.5); ... % Hip AA MA mean change cost, time invariant
    KneeMAPenMean(:)/(0.5*(nMusc).^0.5); ... % Knee MA mean change cost, time invariant
    AnkleMAPenMean(:)/(0.5*(nMusc).^0.5); ... % Ankle MA mean change cost, time invariant
    SubtalarMAPenMean(:)/(0.5*(nMusc).^0.5); ... % Subtalar MA mean change cost, time invariant
    LmtPen(:)/(0.25*(nframesAll*nMusc).^0.5); ... % lmt change cost, time varying
    HipMAPen(:)/(0.25*(nframesAll*nMusc).^0.5); ... % Hip MA change cost, time varying
    HipAAMAPen(:)/(0.25*(nframesAll*nMusc).^0.5); ... % Hip AA MA change cost, time varying
    KneeMAPen(:)/(0.125*(nframesAll*nMusc).^0.5); ... % Knee MA change cost, time varying
    AnkleMAPen(:)/(0.125*(nframesAll*nMusc).^0.5); ... % Ankle MA change cost, time varying
    SubtalarMAPen(:)/(0.125*(nframesAll*nMusc).^0.5); ... % Subtalar MA change cost, time varying
    MAPenalty(:)/.02; ...% moment arm similarity penalty costs, time varying
    sqrt(1e0)*costPassive(:)/(numel(PassiveM)^0.5); ... % minimize passive force, time varying
    lmtildaPenalty(:);... % lmtilda similarity cost, time varying
    (EMGScale(:)-1)/(0.5*(2*nMusc).^0.5)]/sqrt(10);

cost(isnan(cost))=0;

function cost = costFuncFMINCON(guess, optParams)

cost = costFunc(guess, optParams);

cost = sum(cost(:).^2);

function [c,ceq] = nonlcon(guess, optParams)

% no equality constraints
ceq = [];

nMusc = optParams.nMusc;
nframesAll = optParams.nframesAll;

% extract values from optimizer guess
VmaxFactor = 10;
lmo = optParams.lmo.*guess(1,(1:nMusc)+6*nMusc);
lts = optParams.lts.*guess(1,(1:nMusc)+7*nMusc);
if optParams.GeometricAdjustments
    coefsOffset = guess(1,(1:optParams.ncoefsFit)+8*nMusc)/10;
else
    coefsOffset = 0;
end

% extract original moment arms curves
HipMA = optParams.HipMA;
HipAAMA = optParams.HipAAMA;
KneeMA = optParams.KneeMA;
AnkleMA = optParams.AnkleMA;
SubtalarMA = optParams.SubtalarMA;

% write new parameter values to optimization parameters matrix
optParams.lts = lts;
optParams.lmo = lmo;
optParams.VmaxFactor = VmaxFactor;
optParams.coefs = coefsOffset'+optParams.coefsOrig;

% set activations to zero, since they are unnecessary for constraints
a = zeros(nframesAll, nMusc);

optParams.a = a;

% calculate muscle parameter values to be constrained
[~, lmtilda, vmtilda, HipMAnew, HipAAMAnew, KneeMAnew, AnkleMAnew, SubtalarMAnew,~,~,Lmtnew] = calcMoments(optParams,0);

% Remove sign constraints on muscles that can act in two directions at a
% joint (mainly occurs at the hip)
HipMultiDirectional = max(sign(HipMA))+min(sign(HipMA));
HipAAMultiDirectional = max(sign(HipAAMA))+min(sign(HipAAMA));
KneeMultiDirectional = max(sign(KneeMA))+min(sign(KneeMA));
AnkleMultiDirectional = max(sign(AnkleMA))+min(sign(AnkleMA));
SubtalarMultiDirectional = max(sign(SubtalarMA))+min(sign(SubtalarMA));

HipMA(:,HipMultiDirectional==0) = 0;
HipAAMA(:,HipAAMultiDirectional==0) = 0;
KneeMA(:,KneeMultiDirectional==0) = 0;
AnkleMA(:,AnkleMultiDirectional==0) = 0;
SubtalarMA(:,SubtalarMultiDirectional==0) = 0;

% calculate constraint that prevents moment arms from switching signs
HipMACon = optParams.onesCol(:)*mode(-sign(HipMA)).*HipMAnew(:,:);
HipAAMACon = optParams.onesCol(:)*mode(-sign(HipAAMA)).*HipAAMAnew(:,:);
KneeMACon = optParams.onesCol(:)*mode(-sign(KneeMA)).*KneeMAnew(:,:);
AnkleMACon = optParams.onesCol(:)*mode(-sign(AnkleMA)).*AnkleMAnew(:,:);
SubtalarMACon = optParams.onesCol(:)*mode(-sign(SubtalarMA)).*SubtalarMAnew(:,:);

HipMACon(:,sum(HipMACon)==0) = [];
HipAAMACon(:,sum(HipAAMACon)==0) = [];
KneeMACon(:,sum(KneeMACon)==0) = [];
AnkleMACon(:,sum(AnkleMACon)==0) = [];
SubtalarMACon(:,sum(SubtalarMACon)==0) = [];

softmaxlmtildaCon = log(sum(exp(1000*(lmtilda-1.3))))/1000;
softminlmtildaCon = log(sum(exp(1000*(0.3-lmtilda))))/1000;
softminlmtildaCon2 = -log(sum(exp(1000*(1-lmtilda))))/1000;
softmaxlmtildaCon2 = -log(sum(exp(1000*(lmtilda-0.8))))/1000;
softmaxvmtilda = log(sum(exp(1000*vmtilda)))/1000;
softminvmtilda = -log(sum(exp(-1000*vmtilda)))/1000;
softmaxHipMACon = log(sum(exp(10000*HipMACon)))/10000;
softmaxHipAAMACon = log(sum(exp(10000*HipAAMACon)))/10000;
softmaxKneeMACon = log(sum(exp(10000*KneeMACon)))/10000;
softmaxAnkleMACon = log(sum(exp(10000*AnkleMACon)))/10000;
softmaxSubtalarMACon = log(sum(exp(10000*SubtalarMACon)))/10000;

% assign values to c matrix
c = [softmaxlmtildaCon, softminlmtildaCon, softmaxvmtilda-1, softminvmtilda-1,...
    softmaxHipMACon, softmaxHipAAMACon, softmaxKneeMACon, softmaxAnkleMACon, softmaxSubtalarMACon,...
    softmaxlmtildaCon2, softminlmtildaCon2];

c(c<-1000) = -1;
c(c>1000) = 1;

function cost = costFunc(guess, optParams)


nMusc = optParams.nMusc;

% extract values from optimizer guess
tdelay = guess(1,1:2*nMusc)/10;
tact = guess(1,(1:nMusc)+2*nMusc)/100;
tdeact = 4*tact;
Anonlin = guess(1,(1:nMusc)+3*nMusc);
EMGScale = guess(1,(1:2*nMusc)+4*nMusc);
VmaxFactor = 10;
lmoscale = guess(1,(1:nMusc)+6*nMusc);
ltsscale = guess(1,(1:nMusc)+7*nMusc);
if optParams.GeometricAdjustments
    coefsOffset = guess(1,(1:optParams.ncoefsFit)+8*nMusc)/10;
else
    coefsOffset = 0;
end

lts = optParams.lts.*ltsscale;
lmo = optParams.lmo.*lmoscale;

% write new parameter values to optimization parameters matrix
optParams.lts = lts;
optParams.lmo = lmo;
optParams.tdelay = tdelay;
optParams.tact = tact;
optParams.tdeact = tdeact;
optParams.Anonlin = Anonlin;
optParams.VmaxFactor = VmaxFactor;
optParams.EMGScale = EMGScale;
optParams.coefs = coefsOffset'+optParams.coefsOrig;
optParams.ltsscale = ltsscale;
optParams.lmoscale = lmoscale;

cost = costFuncMain(optParams);

function cost = costFuncFMINCON1(guess, optParams)

cost = costFunc1(guess, optParams);

cost = sum(cost(:).^2);

function cost = costFunc1(guess, optParams)

nMusc = optParams.nMusc;

% extract values from optimizer guess
tdelay = guess(1,1:2*nMusc)/10;
tact = guess(1,(1:nMusc)+2*nMusc)/100;
tdeact = 4*tact;
Anonlin = guess(1,(1:nMusc)+3*nMusc);
EMGScale = guess(1,(1:2*nMusc)+4*nMusc);
VmaxFactor = 10;
lmoscale = optParams.otherParams(1,(1:nMusc)+0*nMusc);
ltsscale = optParams.otherParams(1,(1:nMusc)+1*nMusc);
if optParams.GeometricAdjustments
    coefsOffset = optParams.otherParams(1,(1:optParams.ncoefsFit)+2*nMusc)/10;
else
    coefsOffset = 0;
end

lts = optParams.lts.*ltsscale;
lmo = optParams.lmo.*lmoscale;

% write new parameter values to optimization parameters matrix
optParams.lts = lts;
optParams.lmo = lmo;
optParams.tdelay = tdelay;
optParams.tact = tact;
optParams.tdeact = tdeact;
optParams.Anonlin = Anonlin;
optParams.VmaxFactor = VmaxFactor;
optParams.EMGScale = EMGScale;
optParams.coefs = coefsOffset'+optParams.coefsOrig;
optParams.ltsscale = ltsscale;
optParams.lmoscale = lmoscale;

cost = costFuncMain(optParams);

function cost = costFuncFMINCON2(guess, optParams)

cost = costFunc2(guess, optParams);

cost = sum(cost(:).^2);

function [c,ceq] = nonlcon2(guess, optParams)

ceq = [];

nMusc = optParams.nMusc;
nframesAll = optParams.nframesAll;

% extract values from optimizer guess
VmaxFactor = 10;
lmo = optParams.lmo.*guess(1,(1:nMusc)+0*nMusc);
lts = optParams.lts.*guess(1,(1:nMusc)+1*nMusc);
if optParams.GeometricAdjustments
    coefsOffset = optParams.otherParams(1,(1:optParams.ncoefsFit)+6*nMusc)/10;
else
    coefsOffset = 0;
end

% write new parameter values to optimization parameters matrix
optParams.lts = lts;
optParams.lmo = lmo;
optParams.VmaxFactor = VmaxFactor;
optParams.coefs = coefsOffset'+optParams.coefsOrig;

% set activations to zero
a = zeros(nframesAll, nMusc);

optParams.a = a;

% calculate lmtilda and vmtilda
[~, lmtilda, vmtilda] = calcMoments(optParams,0);

softmaxlmtildaCon = log(sum(exp(1000*(lmtilda-1.3))))/1000;
softminlmtildaCon = log(sum(exp(1000*(0.3-lmtilda))))/1000;
softminlmtildaCon2 = -log(sum(exp(1000*(1-lmtilda))))/1000;
softmaxlmtildaCon2 = -log(sum(exp(1000*(lmtilda-0.8))))/1000;
softmaxvmtilda = log(sum(exp(1000*vmtilda)))/1000;
softminvmtilda = -log(sum(exp(-1000*vmtilda)))/1000;

% assign values to c matrix
c = [softmaxlmtildaCon, softminlmtildaCon, softmaxvmtilda-1, softminvmtilda-1,...
   softmaxlmtildaCon2, softminlmtildaCon2];

c(c<-1000) = -1;
c(c>1000) = 1;

function cost = costFunc2(guess, optParams)

nMusc = optParams.nMusc;

% extract values from optimizer guess
tdelay = optParams.otherParams(1,1:2*nMusc)/10;
tact = optParams.otherParams(1,(1:nMusc)+2*nMusc)/100;
tdeact = 4*tact;
Anonlin = optParams.otherParams(1,(1:nMusc)+3*nMusc);
EMGScale = optParams.otherParams(1,(1:2*nMusc)+4*nMusc);
VmaxFactor = 10;
lmoscale = guess(1,(1:nMusc)+0*nMusc);
ltsscale = guess(1,(1:nMusc)+1*nMusc);
if optParams.GeometricAdjustments
    coefsOffset = optParams.otherParams(1,(1:optParams.ncoefsFit)+6*nMusc)/10;
else
    coefsOffset = 0;
end

lts = optParams.lts.*ltsscale;
lmo = optParams.lmo.*lmoscale;

% write new parameter values to optimization parameters matrix
optParams.lts = lts;
optParams.lmo = lmo;
optParams.tdelay = tdelay;
optParams.tact = tact;
optParams.tdeact = tdeact;
optParams.Anonlin = Anonlin;
optParams.VmaxFactor = VmaxFactor;
optParams.EMGScale = EMGScale;
optParams.coefs = coefsOffset'+optParams.coefsOrig;
optParams.ltsscale = ltsscale;
optParams.lmoscale = lmoscale;

cost = costFuncMain(optParams);

function cost = costFuncFMINCON3(guess, optParams)

cost = costFunc3(guess, optParams);

cost = sum(cost(:).^2);

function [c,ceq] = nonlcon3(guess, optParams)

ceq = [];

nMusc = optParams.nMusc;
nframesAll = optParams.nframesAll;

% extract values from optimizer guess
VmaxFactor = 10;
lmo = optParams.lmo.*optParams.otherParams(1,(1:nMusc)+6*nMusc);
lts = optParams.lts.*optParams.otherParams(1,(1:nMusc)+7*nMusc);
coefsOffset = guess/10;

% extract original moment arms curves
HipMA = optParams.HipMA;
HipAAMA = optParams.HipAAMA;
KneeMA = optParams.KneeMA;
AnkleMA = optParams.AnkleMA;
SubtalarMA = optParams.SubtalarMA;

% write new parameter values to optimization parameters matrix
optParams.lts = lts;
optParams.lmo = lmo;
optParams.VmaxFactor = VmaxFactor;
optParams.coefs = coefsOffset'+optParams.coefsOrig;

% set activations to zero, since they are unnecessary for constraints
a = zeros(nframesAll, nMusc);

optParams.a = a;

% calculate muscle parameter values to be constrained
[~, lmtilda, vmtilda, HipMAnew, HipAAMAnew, KneeMAnew, AnkleMAnew, SubtalarMAnew,~,~,Lmtnew] = calcMoments(optParams,0);

% Remove sign constraints on muscles that can act in two directions at a
% joint (mainly occurs at the hip)
HipMultiDirectional = max(sign(HipMA))+min(sign(HipMA));
HipAAMultiDirectional = max(sign(HipAAMA))+min(sign(HipAAMA));
KneeMultiDirectional = max(sign(KneeMA))+min(sign(KneeMA));
AnkleMultiDirectional = max(sign(AnkleMA))+min(sign(AnkleMA));
SubtalarMultiDirectional = max(sign(SubtalarMA))+min(sign(SubtalarMA));

HipMA(:,HipMultiDirectional==0) = 0;
HipAAMA(:,HipAAMultiDirectional==0) = 0;
KneeMA(:,KneeMultiDirectional==0) = 0;
AnkleMA(:,AnkleMultiDirectional==0) = 0;
SubtalarMA(:,SubtalarMultiDirectional==0) = 0;

% calculate constraint that prevents moment arms from switching signs
HipMACon = optParams.onesCol(:)*mode(-sign(HipMA)).*HipMAnew(:,:);
HipAAMACon = optParams.onesCol(:)*mode(-sign(HipAAMA)).*HipAAMAnew(:,:);
KneeMACon = optParams.onesCol(:)*mode(-sign(KneeMA)).*KneeMAnew(:,:);
AnkleMACon = optParams.onesCol(:)*mode(-sign(AnkleMA)).*AnkleMAnew(:,:);
SubtalarMACon = optParams.onesCol(:)*mode(-sign(SubtalarMA)).*SubtalarMAnew(:,:);

HipMACon(:,sum(HipMACon)==0) = [];
HipAAMACon(:,sum(HipAAMACon)==0) = [];
KneeMACon(:,sum(KneeMACon)==0) = [];
AnkleMACon(:,sum(AnkleMACon)==0) = [];
SubtalarMACon(:,sum(SubtalarMACon)==0) = [];

softmaxlmtildaCon = log(sum(exp(1000*(lmtilda-1.3))))/1000;
softminlmtildaCon = log(sum(exp(1000*(0.3-lmtilda))))/1000;
softminlmtildaCon2 = -log(sum(exp(1000*(1-lmtilda))))/1000;
softmaxlmtildaCon2 = -log(sum(exp(1000*(lmtilda-0.8))))/1000;
softmaxvmtilda = log(sum(exp(1000*vmtilda)))/1000;
softminvmtilda = -log(sum(exp(-1000*vmtilda)))/1000;
softmaxHipMACon = log(sum(exp(10000*HipMACon)))/10000;
softmaxHipAAMACon = log(sum(exp(10000*HipAAMACon)))/10000;
softmaxKneeMACon = log(sum(exp(10000*KneeMACon)))/10000;
softmaxAnkleMACon = log(sum(exp(10000*AnkleMACon)))/10000;
softmaxSubtalarMACon = log(sum(exp(10000*SubtalarMACon)))/10000;

% assign values to c matrix
c = [softmaxlmtildaCon, softminlmtildaCon, softmaxvmtilda-1, softminvmtilda-1,...
    softmaxHipMACon, softmaxHipAAMACon, softmaxKneeMACon, softmaxAnkleMACon, softmaxSubtalarMACon,...
    softmaxlmtildaCon2, softminlmtildaCon2];

c(c<-1000) = -1;
c(c>1000) = 1;

function cost = costFunc3(guess, optParams)

nMusc = optParams.nMusc;

% extract values from optimizer guess
tdelay = optParams.otherParams(1,1:2*nMusc)/10;
tact = optParams.otherParams(1,(1:nMusc)+2*nMusc)/100;
tdeact = 4*tact;
Anonlin = optParams.otherParams(1,(1:nMusc)+3*nMusc);
VmaxFactor = 10;
EMGScale = optParams.otherParams(1,(1:2*nMusc)+4*nMusc);
lmoscale = optParams.otherParams(1,(1:nMusc)+6*nMusc);
ltsscale = optParams.otherParams(1,(1:nMusc)+7*nMusc);
coefsOffset = guess/10;

lts = optParams.lts.*ltsscale;
lmo = optParams.lmo.*lmoscale;

% write new parameter values to optimization parameters matrix
optParams.lts = lts;
optParams.lmo = lmo;
optParams.tdelay = tdelay;
optParams.tact = tact;
optParams.tdeact = tdeact;
optParams.Anonlin = Anonlin;
optParams.VmaxFactor = VmaxFactor;
optParams.EMGScale = EMGScale;
optParams.coefs = coefsOffset'+optParams.coefsOrig;
optParams.ltsscale = ltsscale;
optParams.lmoscale = lmoscale;

cost = costFuncMain(optParams);

function [moments, lmtilda, vmtilda, HipMA, HipAAMA, KneeMA, AnkleMA, SubtalarMA, muscMoments, passiveF, Lmt, Vmt, MuscForce] = calcMoments(optParams, DoMomentCalcs)

% extract values from optParams structure
nMusc = optParams.nMusc;
nframesAll = optParams.nframesAll;
alpha = optParams.alpha;
onesCol = optParams.onesCol;
Mat = optParams.Mat;
MuscRef = optParams.MuscRef;
coefs = optParams.coefs;

% MACurvatureOffset = optParams.MACurvatureOffset;
lts = optParams.lts;
lmo = optParams.lmo;
a = optParams.a;

% preallocate memory
HipMA = zeros(nframesAll,nMusc);
HipAAMA = zeros(nframesAll,nMusc);
KneeMA = zeros(nframesAll,nMusc);
AnkleMA = zeros(nframesAll,nMusc);
SubtalarMA = zeros(nframesAll,nMusc);
Lmt = zeros(nframesAll,nMusc);
Vmt = zeros(nframesAll,nMusc);

k = 1;
Ref1 = 1:nframesAll;
Ref2 = nframesAll+1:2*nframesAll;
Ref3 = 2*nframesAll+1:3*nframesAll;
Ref4 = 3*nframesAll+1:4*nframesAll;
Ref5 = 4*nframesAll+1:5*nframesAll;
for i = 1:nMusc
    if MuscRef(i) == 1
        Vec = Mat{1}*coefs(k:k+18);
        Lmt(:,i) = Vec(Ref1);
        Vmt(:,i) = Vec(Ref2);
        HipMA(:,i) = Vec(Ref3);
        HipAAMA(:,i) = Vec(Ref4);
        k = k+19;
    elseif MuscRef(i) == 2
        Vec = Mat{2}*coefs(k:k+21);
        Lmt(:,i) = Vec(Ref1);
        Vmt(:,i) = Vec(Ref2);
        HipMA(:,i) = Vec(Ref3);
        HipAAMA(:,i) = Vec(Ref4);
        KneeMA(:,i) = Vec(Ref5);
        k = k+22;
    elseif MuscRef(i) == 3
        Vec = Mat{3}*coefs(k:k+3);
        Lmt(:,i) = Vec(Ref1);
        Vmt(:,i) = Vec(Ref2);
        KneeMA(:,i) = Vec(Ref3);
        k = k+4;
    elseif MuscRef(i) == 4
        Vec = Mat{4}*coefs(k:k+12);
        Lmt(:,i) = Vec(Ref1);
        Vmt(:,i) = Vec(Ref2);
        KneeMA(:,i) = Vec(Ref3);
        AnkleMA(:,i) = Vec(Ref4);
        SubtalarMA(:,i) = Vec(Ref5);
        k = k+13;
    elseif MuscRef(i) == 5
        Vec = Mat{5}*coefs(k:k+9);
        Lmt(:,i) = Vec(Ref1);
        Vmt(:,i) = Vec(Ref2);
        AnkleMA(:,i) = Vec(Ref3);
        SubtalarMA(:,i) = Vec(Ref4);
        k = k+10;
    end
end

lmtilda = (Lmt-onesCol*lts)./(onesCol*(lmo.*cos(alpha)));
vmtilda = Vmt./(optParams.VmaxFactor*onesCol*(lmo.*cos(alpha)));

Fmax = (optParams.Vmusc./lmo)*optParams.sigma;

if DoMomentCalcs
    
    passiveF = onesCol*(Fmax.*cos(alpha)).*passiveForce(lmtilda);
    MuscForce = onesCol*(Fmax.*cos(alpha)).*a.*FLCurve(lmtilda).*FVCurve(vmtilda)+passiveF;
    
    muscMoments(:,:,1) = HipMA.*MuscForce;
    muscMoments(:,:,2) = HipAAMA.*MuscForce;
    muscMoments(:,:,3) = KneeMA.*MuscForce;
    muscMoments(:,:,4) = AnkleMA.*MuscForce;
    muscMoments(:,:,5) = SubtalarMA.*MuscForce;
    
    moments = permute(sum(muscMoments,2),[1 3 2]);
else
    passiveF = 0;
    moments = 0;
    muscMoments = 0;
end

function [a, Excitations] = EMGtoAct(optParams)

nMusc = optParams.nMusc;
nTrials = optParams.nTrials;
Time = optParams.Time;
tdelay = optParams.tdelay;
EMGsplines = optParams.EMGsplines;
tact(1,1,:) = optParams.tact;
tdeact(1,1,:) = optParams.tdeact;
Anonlin = optParams.Anonlin;
EMGScale = optParams.EMGScale;
onesCol = optParams.onesCol;
nptsShort = optParams.nptsShort;
nptsLong = optParams.nptsLong;
SampleStep = optParams.SampleStep;

EMG = zeros(nptsLong,nTrials,nMusc);
for i = 1:nMusc
    for j = 1:optParams.nTrials_l
        
        pp = EMGsplines{i,j};
        
        br = pp.breaks.';
        cf = pp.coefs;
        
        x_int = Time(:,j)-tdelay(i);
        
        [~, inds] = histc(x_int, [-inf; br(2:end-1); +inf]);
        
        x_shf = x_int - br(inds);
        zero  = ones(size(x_shf));
        one   = x_shf;
        two   = one .* x_shf;
        three = two .* x_shf;
        
        EMG(:,j,i) = sum( [three two one zero] .* cf(inds,:), 2);
        
        %         EMGold(:,j,i) = ppval(Time(:,j)-tdelay(i),EMGsplines{i,j});
        
    end
    for j = (optParams.nTrials_l+1):(optParams.nTrials_l+optParams.nTrials_r)
        
        pp = EMGsplines{i,j};
        
        br = pp.breaks.';
        cf = pp.coefs;
        
        x_int = Time(:,j)-tdelay(i+nMusc);
        
        [~, inds] = histc(x_int, [-inf; br(2:end-1); +inf]);
        
        x_shf = x_int - br(inds);
        zero  = ones(size(x_shf));
        one   = x_shf;
        two   = one .* x_shf;
        three = two .* x_shf;
        
        EMG(:,j,i) = sum( [three two one zero] .* cf(inds,:), 2);
        
        %         EMGold(:,j,i) = ppval(Time(:,j)-tdelay(i+nMusc),EMGsplines{i,j});
    end
end

EMGscales = [ones(optParams.nTrials_l,1)*EMGScale(1:nMusc); ones(optParams.nTrials_r,1)*EMGScale(nMusc+1:end);];
EMGscales = permute(EMGscales, [3 1 2]);
EMG = EMG.*EMGscales(ones(nptsLong,1),:,:);

Excitations = EMG((nptsLong-1)/7+1:nptsLong-(nptsLong-1)/7,:,:);

b2(1,1,:) = 1./tdeact;
b1(1,1,:) = 1./tact-b2;
b3 = b1(ones(nptsLong,1), ones(nTrials,1),:).*EMG+b2(ones(nptsLong,1), ones(nTrials,1),:);
dt(1,:,1) = mean(diff(Time));
b4 = 2*dt(ones(nptsLong,1),:,ones(nMusc,1)).*b3;

EMG = b4.*EMG;

a = zeros(size(EMG));
for j = 3:nptsLong
    a(j,:,:) = (EMG(j,:,:)+4.*a(j-1,:,:)-a(j-2,:,:))./(b4(j,:,:)+3);
end

a = a(((21-1)/SampleStep+1):((121-1)/SampleStep+1),:,:);

a = reshape(a, nptsShort*nTrials, nMusc);
a(a<0) = 0;

% From optimization FindNonlinearityLogFunc.m
nonlincoefs = [29.280183270562596   4.107869238218326 1.000004740962477...
    -7.623282868703527 17.227022969058535   0.884220539986325];

a = (1-Anonlin(onesCol,:)).*a+(Anonlin(onesCol,:)).*...
    (nonlincoefs(4)./(nonlincoefs(1)*(a+nonlincoefs(6)).^nonlincoefs(5)+nonlincoefs(2))+nonlincoefs(3));

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
[b,a] = butter(degree, 2*[HPcutoff 450]/SampRate,'bandpass');
% [sos,g] = zp2sos(z,p,k);      % Convert to SOS(second-order-sections) form
EMG_HP = filtfilt(b,a,emgData);

% Demean
EMG_DM = EMG_HP-ones(size(EMG_HP,1),1)*mean(EMG_HP);

% Rectify
EMG_RF = abs(EMG_DM);

% Low pass filter
[b,a] = butter(degree,LPcutoff/SampRate*2);
EMG_LP = filtfilt(b,a,EMG_RF);
emgFiltered = EMG_LP;

% Remove any negative EMG values that may still exist
emgFiltered(emgFiltered<0) = 0;

function FMtilda = FLCurve(lMtilda)

b11 = 0.814483478343008;
b21 = 1.055033428970575;
b31 = 0.162384573599574;
b41 = 0.063303448465465;
b12 = 0.433004984392647;
b22 = 0.716775413397760;
b32 = -0.029947116970696;
b42 = 0.200356847296188;
b13 = 0.1;
b23 = 1;
b33 = 0.5*sqrt(0.5);
b43 = 0;

FMtilda = b11*exp(-0.5*(lMtilda-b21).^2./(b31+b41*lMtilda).^2)...
    +b12*exp(-0.5*(lMtilda-b22).^2./(b32+b42*lMtilda).^2)+...
    b13*exp(-0.5*(lMtilda-b23).^2./(b33+b43*lMtilda).^2);

function FMtilda = FVCurve(vMtilda)
% Coefficients determined manually

d1 = -8.665833008200368;
d2 = 7.623504943764548;
d3 = 3.565407127548163;
d4 = 0.968466940650055;
d5 = -0.394453556647392;
d6 = 8.169548497993912;

FMtilda = d1+d2*atan(d3+d4*atan(d5+d6*vMtilda));

function PFtilda = passiveForce(lMtilda)

c1 = 3.384513491001194;
c2 = -74.090638944399615;
c3 = -2.443789456896774;

PFtilda = c1*exp(c2*(exp(lMtilda)).^c3);

function [Mat] = buildMuscleMatrices(JAngles,JVels)
% This function creates matrices that when multiplied by the muscle
% geometry coefficients will produce the Lmt, Vmt and moment arm curves
% each matrix has the form [Angle Angle^2 Angle^3 ... other joint ...
% interaction terms]

nFrames = size(JAngles,1);
onesCol = ones(nFrames,1);
zerosMat = zeros(nFrames,3);

% matrices for determining lmt
HipFEMat = [JAngles(:,1) JAngles(:,1).^2 JAngles(:,1).^3];
HipAAMat = [JAngles(:,2) JAngles(:,2).^2 JAngles(:,2).^3];
HipIEMat = [JAngles(:,6) JAngles(:,6).^2 JAngles(:,6).^3];
InteractionFEAA = [JAngles(:,1).*JAngles(:,2) JAngles(:,1).^2.*JAngles(:,2) JAngles(:,1).*JAngles(:,2).^2];
InteractionIEAA = [JAngles(:,6).*JAngles(:,2) JAngles(:,6).^2.*JAngles(:,2) JAngles(:,6).*JAngles(:,2).^2];
InteractionFEIE = [JAngles(:,1).*JAngles(:,6) JAngles(:,1).^2.*JAngles(:,6) JAngles(:,1).*JAngles(:,6).^2];
KneeMat = [JAngles(:,3) JAngles(:,3).^2 JAngles(:,3).^3];
AnkleMat = [JAngles(:,4) JAngles(:,4).^2 JAngles(:,4).^3];
SubtalarMat = [JAngles(:,5) JAngles(:,5).^2 JAngles(:,5).^3];
InteractionAnkleSub = [JAngles(:,4).*JAngles(:,5) JAngles(:,4).^2.*JAngles(:,5) JAngles(:,4).*JAngles(:,5).^2];

% matrices for determining vmt
HipFEVelMat = [JVels(:,1) 2*JVels(:,1).*JAngles(:,1) 3*JVels(:,1).*JAngles(:,1).^2];
HipAAVelMat = [JVels(:,2) 2*JVels(:,2).*JAngles(:,2) 3*JVels(:,2).*JAngles(:,2).^2];
HipIEVelMat = [JVels(:,6) 2*JVels(:,6).*JAngles(:,6) 3*JVels(:,6).*JAngles(:,6).^2];
InteractionFEAAVel = [JVels(:,1).*JAngles(:,2)+JAngles(:,1).*JVels(:,2) 2*JAngles(:,1).*JVels(:,1).*JAngles(:,2)+JAngles(:,1).^2.*JVels(:,2) JVels(:,1).*JAngles(:,2).^2+2*JAngles(:,1).*JVels(:,2).*JAngles(:,2)];
InteractionIEAAVel = [JVels(:,6).*JAngles(:,2)+JAngles(:,6).*JVels(:,2) 2*JAngles(:,6).*JVels(:,6).*JAngles(:,2)+JAngles(:,6).^2.*JVels(:,2) JVels(:,6).*JAngles(:,2).^2+2*JAngles(:,6).*JVels(:,2).*JAngles(:,2)];
InteractionFEIEVel = [JVels(:,1).*JAngles(:,6)+JAngles(:,1).*JVels(:,6) 2*JAngles(:,1).*JVels(:,1).*JAngles(:,6)+JAngles(:,1).^2.*JVels(:,6) JVels(:,1).*JAngles(:,6).^2+2*JAngles(:,1).*JVels(:,6).*JAngles(:,6)];
KneeVelMat = [JVels(:,3) 2*JVels(:,3).*JAngles(:,3) 3*JVels(:,3).*JAngles(:,3).^2];
AnkleVelMat = [JVels(:,4) 2*JVels(:,4).*JAngles(:,4) 3*JVels(:,4).*JAngles(:,4).^2];
SubtalarVelMat = [JVels(:,5) 2*JVels(:,5).*JAngles(:,5) 3*JVels(:,5).*JAngles(:,5).^2];
InteractionAnkleSubVel = [JVels(:,4).*JAngles(:,5)+JAngles(:,4).*JVels(:,5) 2*JAngles(:,4).*JVels(:,4).*JAngles(:,5)+JAngles(:,4).^2.*JVels(:,5) JVels(:,4).*JAngles(:,5).^2+2*JAngles(:,4).*JVels(:,5).*JAngles(:,5)];

% matrices for determining moment arms
HipFEMAMat = -[onesCol 2*JAngles(:,1) 3*JAngles(:,1).^2];
HipAAMAMat = -[onesCol 2*JAngles(:,2) 3*JAngles(:,2).^2];
InteractionFEAAMA_FE = -[JAngles(:,2) 2*JAngles(:,1).*JAngles(:,2) JAngles(:,2).^2];
InteractionFEAAMA_AA = -[JAngles(:,1) JAngles(:,1).^2 2*JAngles(:,1).*JAngles(:,2)];
InteractionFEIEMA_FE = -[JAngles(:,6) 2*JAngles(:,1).*JAngles(:,6) JAngles(:,6).^2];
InteractionIEAAMA_AA = -[JAngles(:,6) JAngles(:,6).^2 2*JAngles(:,6).*JAngles(:,2)];
KneeMAMat = -[onesCol 2*JAngles(:,3) 3*JAngles(:,3).^2];
AnkleMAMat = -[onesCol 2*JAngles(:,4) 3*JAngles(:,4).^2];
SubtalarMAMat = -[onesCol 2*JAngles(:,5) 3*JAngles(:,5).^2];
InteractionAnkleMA = -[JAngles(:,5) 2*JAngles(:,4).*JAngles(:,5) JAngles(:,5).^2];
InteractionSubMA = -[JAngles(:,4) JAngles(:,4).^2 2*JAngles(:,4).*JAngles(:,5)];

% combine all matrices into a structure
Mat{1} = [onesCol HipFEMat HipAAMat HipIEMat InteractionFEAA InteractionIEAA InteractionFEIE;...
    0*onesCol HipFEVelMat HipAAVelMat HipIEVelMat InteractionFEAAVel InteractionIEAAVel InteractionFEIEVel;...
    0*onesCol HipFEMAMat zerosMat zerosMat InteractionFEAAMA_FE zerosMat InteractionFEIEMA_FE;...
    0*onesCol zerosMat HipAAMAMat zerosMat InteractionFEAAMA_AA InteractionIEAAMA_AA zerosMat];
Mat{2} = [onesCol HipFEMat HipAAMat HipIEMat KneeMat InteractionFEAA InteractionIEAA InteractionFEIE;...
    0*onesCol HipFEVelMat HipAAVelMat HipIEVelMat KneeVelMat InteractionFEAAVel InteractionIEAAVel InteractionFEIEVel;...
    0*onesCol HipFEMAMat zerosMat zerosMat zerosMat InteractionFEAAMA_FE zerosMat InteractionFEIEMA_FE;...
    0*onesCol zerosMat HipAAMAMat zerosMat zerosMat InteractionFEAAMA_AA InteractionIEAAMA_AA zerosMat;...
    0*onesCol zerosMat zerosMat zerosMat KneeMAMat zerosMat zerosMat zerosMat];
Mat{3} = [onesCol KneeMat;...
    0*onesCol KneeVelMat;...
    0*onesCol KneeMAMat];
Mat{4} = [onesCol KneeMat AnkleMat SubtalarMat InteractionAnkleSub;...
    0*onesCol KneeVelMat AnkleVelMat SubtalarVelMat InteractionAnkleSubVel;...
    0*onesCol KneeMAMat zerosMat zerosMat zerosMat;...
    0*onesCol zerosMat AnkleMAMat zerosMat InteractionAnkleMA;...
    0*onesCol zerosMat zerosMat SubtalarMAMat InteractionSubMA];
Mat{5} = [onesCol AnkleMat SubtalarMat InteractionAnkleSub;...
    0*onesCol AnkleVelMat SubtalarVelMat InteractionAnkleSubVel;...
    0*onesCol AnkleMAMat zerosMat InteractionAnkleMA;...
    0*onesCol zerosMat SubtalarMAMat InteractionSubMA];

function makeParameterChangePlots(optParams)

tact = optParams.tact;
Anonlin = optParams.Anonlin;
ltsscale = optParams.ltsscale;
lmoscale = optParams.lmoscale;
EMGScale = optParams.EMGScale;
nMusc = optParams.nMusc;

% plot design variables
muscLabels = optParams.MuscLabels;
fig = figure('units','normalized','outerposition',[0 0 1 1]);
bar(tact*1000)
title('Muscle Activation Time Constants')
ylabel('Time (ms)')
set(gca, 'XTick', 1:nMusc, 'XTickLabel', muscLabels, 'FontSize', 10);
%     rotateXLabels(gca,90);
%     saveas(fig, 'ActivationTimeConstants.jpg')

fig = figure('units','normalized','outerposition',[0 0 1 1]);
bar(Anonlin)
title('Muscle Nonlinearity Parameter')
ylabel('Nonlinearity Param')
set(gca, 'XTick', 1:nMusc, 'XTickLabel', muscLabels, 'FontSize', 10);
%     rotateXLabels(gca,90);
%     saveas(fig, 'MuscNonlinParams.jpg')

fig = figure('units','normalized','outerposition',[0 0 1 1]);
bar(ltsscale)
title('Muscle Lts Scale Factors')
ylabel('Scale Factor')
set(gca, 'XTick', 1:nMusc, 'XTickLabel', muscLabels, 'FontSize', 10);
%     rotateXLabels(gca,90);
%     saveas(fig, 'MuscLtsScale.jpg')

fig = figure('units','normalized','outerposition',[0 0 1 1]);
bar(lmoscale)
title('Muscle Lmo Scale Factors')
ylabel('Scale Factor')
set(gca, 'XTick', 1:nMusc, 'XTickLabel', muscLabels, 'FontSize', 10);
%     rotateXLabels(gca,90);
%     saveas(fig, 'MuscLmoScale.jpg')

fig = figure('units','normalized','outerposition',[0 0 1 1]);
bar(EMGScale)
title('Muscle Activation Scale Factors')
ylabel('Scale Factor')
set(gca, 'XTick', 1:nMusc, 'XTickLabel', [muscLabels muscLabels], 'FontSize', 10);

function makeIndividualMuscMomentPlots(moments,muscMoments,optParams)

moments = reshape(moments,optParams.nptsShort*optParams.nTrials,5);
figure
plot(-moments(:,1),'r', 'LineWidth',2)
hold on
plot(-optParams.IDloads(:,1),'b', 'LineWidth',2)
plot(-muscMoments(:,:,1))
title('Hip FE Moment Calibration Results')
xlabel('Time Frames')
ylabel('Moment (N-m)')
legend(['Hip FE Moment ID' 'Hip FE Moment Muscles' optParams.MuscLabels])

figure
plot(-moments(:,2),'r', 'LineWidth',2)
hold on
plot(-optParams.IDloads(:,2),'b', 'LineWidth',2)
plot(-muscMoments(:,:,2))
title('Hip AA Moment Calibration Results')
xlabel('Time Frames')
ylabel('Moment (N-m)')
legend(['Hip AA Moment ID' 'Hip AA Moment Muscles' optParams.MuscLabels])

figure
plot(-moments(:,3),'r', 'LineWidth',2)
hold on
plot(-optParams.IDloads(:,3),'b', 'LineWidth',2)
plot(-muscMoments(:,:,3))
title('Knee FE Moment Calibration Results')
xlabel('Time Frames')
ylabel('Moment (N-m)')
legend(['Knee FE Moment ID' 'Knee FE Moment Muscles' optParams.MuscLabels])

figure
plot(-moments(:,4),'r', 'LineWidth',2)
hold on
plot(-optParams.IDloads(:,4),'b', 'LineWidth',2)
plot(-muscMoments(:,:,4))
title('Ankle Moment Calibration Results')
xlabel('Time Frames')
ylabel('Moment (N-m)')
legend(['Ankle Moment ID' 'Ankle Moment Muscles' optParams.MuscLabels])

figure
plot(-moments(:,5),'r', 'LineWidth',2)
hold on
plot(-optParams.IDloads(:,5),'b', 'LineWidth',2)
plot(-muscMoments(:,:,5))
title('Subtalar Moment Calibration Results')
xlabel('Time Frames')
ylabel('Moment (N-m)')
legend(['Subtalar Moment ID' 'Subtalar Moment Muscles' optParams.MuscLabels])

function makeMomentMatchingPlots(moments,IDloads,optParams)

nTrials_l = optParams.nTrials_l;
nTrials_r = optParams.nTrials_r;
nTrialsPerSpeed = optParams.numTrialsPerSpeed;
PassiveMoments = optParams.PassiveMGait;

for i = 1:optParams.nSpeeds/2
    
    momentsMean_l = permute(mean(moments(:,(i-1)*nTrialsPerSpeed+1:i*nTrialsPerSpeed,:),2),[1 3 2]);
    momentsMin_l = permute(max(moments(:,(i-1)*nTrialsPerSpeed+1:i*nTrialsPerSpeed,:),[],2),[1 3 2]);
    momentsMax_l = permute(min(moments(:,(i-1)*nTrialsPerSpeed+1:i*nTrialsPerSpeed,:),[],2),[1 3 2]);
    PassiveMomentsMean_l = permute(min(PassiveMoments(:,(i-1)*nTrialsPerSpeed+1:i*nTrialsPerSpeed,:),[],2),[1 3 2]);
    
    IDloadsMean_l = permute(mean(IDloads(:,(i-1)*nTrialsPerSpeed+1:i*nTrialsPerSpeed,:),2),[1 3 2]);
    IDloadsMin_l = permute(max(IDloads(:,(i-1)*nTrialsPerSpeed+1:i*nTrialsPerSpeed,:),[],2),[1 3 2]);
    IDloadsMax_l = permute(min(IDloads(:,(i-1)*nTrialsPerSpeed+1:i*nTrialsPerSpeed,:),[],2),[1 3 2]);
    
    momentsMean_r = permute(mean(moments(:,((i-1)*nTrialsPerSpeed+1:i*nTrialsPerSpeed)+nTrials_l,:),2),[1 3 2]);
    momentsMin_r = permute(max(moments(:,((i-1)*nTrialsPerSpeed+1:i*nTrialsPerSpeed)+nTrials_l,:),[],2),[1 3 2]);
    momentsMax_r = permute(min(moments(:,((i-1)*nTrialsPerSpeed+1:i*nTrialsPerSpeed)+nTrials_l,:),[],2),[1 3 2]);
    PassiveMomentsMean_r = permute(min(PassiveMoments(:,((i-1)*nTrialsPerSpeed+1:i*nTrialsPerSpeed)+nTrials_l,:),[],2),[1 3 2]);
    
    IDloadsMean_r = permute(mean(IDloads(:,((i-1)*nTrialsPerSpeed+1:i*nTrialsPerSpeed)+nTrials_l,:),2),[1 3 2]);
    IDloadsMin_r = permute(max(IDloads(:,((i-1)*nTrialsPerSpeed+1:i*nTrialsPerSpeed)+nTrials_l,:),[],2),[1 3 2]);
    IDloadsMax_r = permute(min(IDloads(:,((i-1)*nTrialsPerSpeed+1:i*nTrialsPerSpeed)+nTrials_l,:),[],2),[1 3 2]);
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,5,1), plot(-momentsMean_l(:,1),'r', 'LineWidth',2)
    hold on
    subplot(2,5,1), plot(-IDloadsMean_l(:,1),'k', 'LineWidth',2)
    subplot(2,5,1), plot(-PassiveMomentsMean_l(:,1),'r--', 'LineWidth',2)
    title('Hip FE Moment Calibration Results')
    xlabel('Time Frames')
    ylabel('Moment (N-m)')
    axis([0 100 -50 100])
    set(gca, 'FontSize', 10)
    
    subplot(2,5,2), plot(-momentsMean_l(:,2),'r', 'LineWidth',2)
    hold on
    subplot(2,5,2), plot(-IDloadsMean_l(:,2),'k', 'LineWidth',2)
    subplot(2,5,2), plot(-PassiveMomentsMean_l(:,2),'r--', 'LineWidth',2)
    title('Hip AA Moment Calibration Results')
    xlabel('Time Frames')
    ylabel('Moment (N-m)')
    axis([0 100 -50 100])
    set(gca, 'FontSize', 10)
    
    subplot(2,5,3), plot(-momentsMean_l(:,3),'r', 'LineWidth',2)
    hold on
    subplot(2,5,3), plot(-IDloadsMean_l(:,3),'k', 'LineWidth',2)
    subplot(2,5,3), plot(-PassiveMomentsMean_l(:,3),'r--', 'LineWidth',2)
    title('Knee FE Moment Calibration Results')
    xlabel('Time Frames')
    ylabel('Moment (N-m)')
    axis([0 100 -50 100])
    set(gca, 'FontSize', 10)
    
    subplot(2,5,4), plot(-momentsMean_l(:,4),'r', 'LineWidth',2)
    hold on
    subplot(2,5,4), plot(-IDloadsMean_l(:,4),'k', 'LineWidth',2)
    subplot(2,5,4), plot(-PassiveMomentsMean_l(:,4),'r--', 'LineWidth',2)
    title('Ankle FE Moment Calibration Results')
    xlabel('Time Frames')
    ylabel('Moment (N-m)')
    axis([0 100 -50 100])
    set(gca, 'FontSize', 10)
    
    subplot(2,5,5), plot(-momentsMean_l(:,5),'r', 'LineWidth',2)
    hold on
    subplot(2,5,5), plot(-IDloadsMean_l(:,5),'k', 'LineWidth',2)
    subplot(2,5,5), plot(-PassiveMomentsMean_l(:,5),'r--', 'LineWidth',2)
    title('Subtalar IE Moment Calibration Results')
    xlabel('Time Frames')
    ylabel('Moment (N-m)')
    axis([0 100 -50 100])
    legend({'Average Model Moment' 'Average ID Moment'});% 'Model Min/Max Moment' 'ID Min/Max Moment'})
    set(gca, 'FontSize', 10)
    
    subplot(2,5,6), plot(-momentsMean_r(:,1),'r', 'LineWidth',2)
    hold on
    subplot(2,5,6), plot(-IDloadsMean_r(:,1),'k', 'LineWidth',2)
    subplot(2,5,6), plot(-PassiveMomentsMean_r(:,1),'r--', 'LineWidth',2)
    title('Hip FE Moment Calibration Results')
    xlabel('Time Frames')
    ylabel('Moment (N-m)')
    axis([0 100 -50 100])
    set(gca, 'FontSize', 10)
    
    subplot(2,5,7), plot(-momentsMean_r(:,2),'r', 'LineWidth',2)
    hold on
    subplot(2,5,7), plot(-IDloadsMean_r(:,2),'k', 'LineWidth',2)
    subplot(2,5,7), plot(-PassiveMomentsMean_r(:,2),'r--', 'LineWidth',2)
    title('Hip AA Moment Calibration Results')
    xlabel('Time Frames')
    ylabel('Moment (N-m)')
    axis([0 100 -50 100])
    set(gca, 'FontSize', 10)
    
    subplot(2,5,8), plot(-momentsMean_r(:,3),'r', 'LineWidth',2)
    hold on
    subplot(2,5,8), plot(-IDloadsMean_r(:,3),'k', 'LineWidth',2)
    subplot(2,5,8), plot(-PassiveMomentsMean_r(:,3),'r--', 'LineWidth',2)
    title('Knee FE Moment Calibration Results')
    xlabel('Time Frames')
    ylabel('Moment (N-m)')
    axis([0 100 -50 100])
    set(gca, 'FontSize', 10)
    
    subplot(2,5,9), plot(-momentsMean_r(:,4),'r', 'LineWidth',2)
    hold on
    subplot(2,5,9), plot(-IDloadsMean_r(:,4),'k', 'LineWidth',2)
    subplot(2,5,9), plot(-PassiveMomentsMean_r(:,4),'r--', 'LineWidth',2)
    title('Ankle FE Moment Calibration Results')
    xlabel('Time Frames')
    ylabel('Moment (N-m)')
    axis([0 100 -50 100])
    set(gca, 'FontSize', 10)
    
    subplot(2,5,10), plot(-momentsMean_r(:,5),'r', 'LineWidth',2)
    hold on
    subplot(2,5,10), plot(-IDloadsMean_r(:,5),'k', 'LineWidth',2)
    subplot(2,5,10), plot(-PassiveMomentsMean_r(:,5),'r--', 'LineWidth',2)
    title('Subtalar IE Moment Calibration Results')
    xlabel('Time Frames')
    ylabel('Moment (N-m)')
    axis([0 100 -50 100])
    legend({'Average Model Moment' 'Average ID Moment'});% 'Model Min/Max Moment' 'ID Min/Max Moment'})
    set(gca, 'FontSize', 10)
    saveas(fig, 'MomentMatchingAverageCurve.jpg')
    
end

function makeMuscleParameterPlots(Lmt, lmtilda, vmtilda, HipFEMA, HipAAMA, KneeMA, AnkleMA, SubtalarMA,Lmtorig,HipMAorig,HipAAMAorig,KneeMAorig,AnkleMAorig,SubtalarMAorig,optParams)

nMusc = optParams.nMusc;
OpensimMuscleLabels = optParams.MuscLabels;

Lmtmean = permute(mean(Lmt,2),[1 3 2]);
lmtildamean = permute(mean(lmtilda,2),[1 3 2]);
vmtildamean = permute(mean(vmtilda,2),[1 3 2]);
HipMAmean = permute(mean(HipFEMA,2),[1 3 2]);
HipAAMAmean = permute(mean(HipAAMA,2),[1 3 2]);
KneeMAmean = permute(mean(KneeMA,2),[1 3 2]);
AnkleMAmean = permute(mean(AnkleMA,2),[1 3 2]);
SubtalarMAmean = permute(mean(SubtalarMA,2),[1 3 2]);

Lmtorigmean = permute(mean(Lmtorig,3),[1 2 3]);
HipMAorigmean = permute(mean(HipMAorig,3),[1 2 3]);
HipAAMAorigmean = permute(mean(HipAAMAorig,3),[1 2 3]);
KneeMAorigmean = permute(mean(KneeMAorig,3),[1 2 3]);
AnkleMAorigmean = permute(mean(AnkleMAorig,3),[1 2 3]);
SubtalarMAorigmean = permute(mean(SubtalarMAorig,3),[1 2 3]);

Lmtmax = permute(max(Lmt,[],2),[1 2 3]);
lmtildamax = permute(max(lmtilda,[],2),[1 3 2]);
vmtildamax = permute(max(vmtilda,[],2),[1 3 2]);
HipMAmax = permute(max(HipFEMA,[],2),[1 3 2]);
HipAAMAmax = permute(max(HipAAMA,[],2),[1 3 2]);
KneeMAmax = permute(max(KneeMA,[],2),[1 3 2]);
AnkleMAmax = permute(max(AnkleMA,[],2),[1 3 2]);
SubtalarMAmax = permute(max(SubtalarMA,[],2),[1 3 2]);

Lmtmin = permute(min(Lmt,[],2),[1 3 2]);
lmtildamin = permute(min(lmtilda,[],2),[1 3 2]);
vmtildamin = permute(min(vmtilda,[],2),[1 3 2]);
HipMAmin = permute(min(HipFEMA,[],2),[1 3 2]);
HipAAMAmin = permute(min(HipAAMA,[],2),[1 3 2]);
KneeMAmin = permute(min(KneeMA,[],2),[1 3 2]);
AnkleMAmin = permute(min(AnkleMA,[],2),[1 3 2]);
SubtalarMAmin = permute(min(SubtalarMA,[],2),[1 3 2]);

fig = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:nMusc
    subplot(6,6,i), plot(0:optParams.SampleStep:100,Lmtmean(:,i))
    hold on
    subplot(6,6,i), plot(0:optParams.SampleStep:100,Lmtmax(:,i),'r-.')
    subplot(6,6,i), plot(0:optParams.SampleStep:100,Lmtmin(:,i),'r-.')
    subplot(6,6,i), plot(0:optParams.SampleStep:100,Lmtorigmean(:,i),'g')
    title([OpensimMuscleLabels{i} ' mean Lmt'])
    %         axis([0 100 .4 1.4]);
    set(gca, 'FontSize', 10)
end

fig = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:nMusc
    subplot(6,6,i), plot(0:optParams.SampleStep:100,lmtildamean(:,i))
    hold on
    subplot(6,6,i), plot(0:optParams.SampleStep:100,lmtildamax(:,i),'r-.')
    subplot(6,6,i), plot(0:optParams.SampleStep:100,lmtildamin(:,i),'r-.')
    title([OpensimMuscleLabels{i} ' mean lmtilda'])
    axis([0 100 .4 1.4]);
    set(gca, 'FontSize', 10)
end
%     saveas(fig, 'Musclmtilda.jpg')

fig = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:nMusc
    subplot(6,6,i), plot(0:optParams.SampleStep:100,vmtildamean(:,i))
    hold on
    subplot(6,6,i), plot(0:optParams.SampleStep:100,vmtildamax(:,i),'r-.')
    subplot(6,6,i), plot(0:optParams.SampleStep:100,vmtildamin(:,i),'r-.')
    title([OpensimMuscleLabels{i} ' mean vmtilda'])
    axis([0 100 -1 1]);
    set(gca, 'FontSize', 10)
end
%     saveas(fig, 'Muscvmtilda.jpg')

fig = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:nMusc
    subplot(6,6,i), plot(0:optParams.SampleStep:100,HipMAmean(:,i))
    hold on
    subplot(6,6,i), plot(0:optParams.SampleStep:100,HipMAorigmean(:,i),'g')
    subplot(6,6,i), plot(0:optParams.SampleStep:100,HipMAmax(:,i),'r-.')
    subplot(6,6,i), plot(0:optParams.SampleStep:100,HipMAmin(:,i),'r-.')
    title([OpensimMuscleLabels{i} ' mean Hip FE MA'])
    set(gca, 'FontSize', 10)
end
%     saveas(fig, 'HipFEMA.jpg')

fig = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:nMusc
    subplot(6,6,i), plot(0:optParams.SampleStep:100,HipAAMAmean(:,i))
    hold on
    subplot(6,6,i), plot(0:optParams.SampleStep:100,HipAAMAorigmean(:,i),'g')
    subplot(6,6,i), plot(0:optParams.SampleStep:100,HipAAMAmax(:,i),'r-.')
    subplot(6,6,i), plot(0:optParams.SampleStep:100,HipAAMAmin(:,i),'r-.')
    title([OpensimMuscleLabels{i} ' mean Hip AA MA'])
    set(gca, 'FontSize', 10)
end
%     saveas(fig, 'HipAAMA.jpg')

fig = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:nMusc
    subplot(6,6,i), plot(0:optParams.SampleStep:100,KneeMAmean(:,i))
    hold on
    subplot(6,6,i), plot(0:optParams.SampleStep:100,KneeMAorigmean(:,i),'g')
    subplot(6,6,i), plot(0:optParams.SampleStep:100,KneeMAmax(:,i),'r-.')
    subplot(6,6,i), plot(0:optParams.SampleStep:100,KneeMAmin(:,i),'r-.')
    title([OpensimMuscleLabels{i} ' mean Knee MA'])
    set(gca, 'FontSize', 10)
end
%     saveas(fig, 'KneeMA.jpg')

fig = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:nMusc
    subplot(6,6,i), plot(0:optParams.SampleStep:100,AnkleMAmean(:,i))
    hold on
    subplot(6,6,i), plot(0:optParams.SampleStep:100,AnkleMAorigmean(:,i),'g')
    subplot(6,6,i), plot(0:optParams.SampleStep:100,AnkleMAmax(:,i),'r-.')
    subplot(6,6,i), plot(0:optParams.SampleStep:100,AnkleMAmin(:,i),'r-.')
    title([OpensimMuscleLabels{i} ' mean Ankle MA'])
    set(gca, 'FontSize', 10)
end
%     saveas(fig, 'AnkleMA.jpg')

fig = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:nMusc
    subplot(6,6,i), plot(0:optParams.SampleStep:100,SubtalarMAmean(:,i))
    hold on
    subplot(6,6,i), plot(0:optParams.SampleStep:100,SubtalarMAorigmean(:,i),'g')
    subplot(6,6,i), plot(0:optParams.SampleStep:100,SubtalarMAmax(:,i),'r-.')
    subplot(6,6,i), plot(0:optParams.SampleStep:100,SubtalarMAmin(:,i),'r-.')
    title([OpensimMuscleLabels{i} ' mean Subtalar MA'])
    set(gca, 'FontSize', 10)
end

function makeExcitationPlots(Excitations,a,optParams)

nTrials_l = optParams.nTrials_l;
nTrialsPerSpeed = optParams.numTrialsPerSpeed;

MuscleLabels = optParams.MuscLabels;

a = reshape(a,optParams.nptsShort,optParams.nTrials,optParams.nMusc);

for k = 1:optParams.nSpeeds/2
    figure
    ExcitationMeans = permute(mean(Excitations(:,(k-1)*nTrialsPerSpeed+1:k*nTrialsPerSpeed,:),2),[1 3 2]);
    aMeans = permute(mean(a(:,(k-1)*nTrialsPerSpeed+1:k*nTrialsPerSpeed,:),2),[1 3 2]);
    for i = 1:35
        subplot(6,6,i), plot(0:optParams.SampleStep:100,ExcitationMeans(:,i),'LineWidth',2);
        hold on
        subplot(6,6,i), plot(0:optParams.SampleStep:100,aMeans(:,i),'r','LineWidth',2);
        axis([0 100 0 1])
        title(['Left ' MuscleLabels{i}]);
    end
    figure
    ExcitationMeans = permute(mean(Excitations(:,((k-1)*nTrialsPerSpeed+1:k*nTrialsPerSpeed)+nTrials_l,:),2),[1 3 2]);
    aMeans = permute(mean(a(:,((k-1)*nTrialsPerSpeed+1:k*nTrialsPerSpeed)+nTrials_l,:),2),[1 3 2]);
    for i = 1:35
        subplot(6,6,i), plot(0:optParams.SampleStep:100,ExcitationMeans(:,i),'LineWidth',2);
        hold on
        subplot(6,6,i), plot(0:optParams.SampleStep:100,aMeans(:,i),'r','LineWidth',2);
        axis([0 100 0 1])
        title(['Right ' MuscleLabels{i}]);
    end
end