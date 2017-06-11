function PrecalibrationOflmoandlts_v3

muscLabels = {'addbrev_l','addlong_l','addmagDist_l','addmagIsch_l','addmagMid_l','addmagProx_l','glmax1_l','glmax2_l','glmax3_l','glmed1_l','glmed2_l','glmed3_l','glmin1_l','glmin2_l','glmin3_l','iliacus_l','psoas_l','semimem_l','semiten_l','bflh_l','bfsh_l','recfem_l','vasmed_l','vaslat_l','vasint_l','gaslat_l','gasmed_l','tibant_l','tibpost_l','perbrev_l','perlong_l','pertert_l','soleus_l','edl_l','fdl_l';}

numOpt = 3;

close all

load etmData_refitted_both_outliersRemoved.mat
load MuscleFits_lhs.mat
load DigitizedPassiveMoments.mat

GaitData = [etmData_l.TMGait_0pt3 etmData_l.TMGait_0pt4 etmData_l.TMGait_0pt5 etmData_l.TMGait_0pt6 etmData_l.TMGait_0pt7 etmData_l.TMGait_0pt8...
    etmData_r.TMGait_0pt3 etmData_r.TMGait_0pt4 etmData_r.TMGait_0pt5 etmData_r.TMGait_0pt6 etmData_r.TMGait_0pt7 etmData_r.TMGait_0pt8];
DataLengths = [length(etmData_l.TMGait_0pt3) length(etmData_l.TMGait_0pt4) length(etmData_l.TMGait_0pt5) length(etmData_l.TMGait_0pt6) length(etmData_l.TMGait_0pt7) length(etmData_l.TMGait_0pt8)...
    length(etmData_r.TMGait_0pt3) length(etmData_r.TMGait_0pt4) length(etmData_r.TMGait_0pt5) length(etmData_r.TMGait_0pt6) length(etmData_r.TMGait_0pt7) length(etmData_r.TMGait_0pt8)];

Step = 5;
nTrials = length(GaitData);
nFramesGait = (101-1)/Step+1;

JAnglesGait = zeros(nFramesGait,6,nTrials);

for i = 1:nTrials
    
    JAnglesGait(:,:,i) = GaitData(i).JointAngles(21:Step:121,1:6)*pi/180;
    
end

% DataLengths = [0 DataLengths];
% color = {'b','r','g','k','c','y'};
% for i = 1:6
%     plot(mean(permute(JAnglesGait(:,1,sum(DataLengths(1:i))+1:sum(DataLengths(1:i))+DataLengths(i+1))*180/pi,[1 3 2]),2),color{i})
%     hold on
% end

% JAnglesGait = mean(JAnglesGait,3);

JAnglesGait = reshape(permute(JAnglesGait, [1 3 2]), nFramesGait*nTrials,6);

%Hip first, then knee, then ankle

%% Define muscle properties

nMusc = 35;
coefsFit = MuscleFits.Coefs;
coefsOrig = [];
for i = 1:nMusc
    coefsOrig = [coefsOrig; coefsFit{i}];
end

ncoefsFit = length(coefsOrig);

lmo =   [10.3 10.8 17.7 15.6 13.8 10.6 14.7 15.7 16.7 7.3 7.3 7.3 6.8 5.6 3.8 10.7 11.7 6.9 19.3 9.8 11 7.6 9.7 9.9 9.9 5.9 5.1 6.8 3.8 4.5 5.1 7.9 4.4 6.9 4.5]/100;
lts =   [3.6 13.0 9.0 22.1 4.8 4.3 5.0 7.3 7.0 5.7 6.6 4.6 1.6 2.6 5.1 9.4 9.7 37.8 24.5 32.2 10.4 34.6 11.2 13 10.6 38.2 40.1 24.1 28.2 14.9 33.3 10.0 28.2 36.7 37.8]/100;
alpha = [6.1 7.1 13.8 11.9 14.7 22.2 21.1 21.9 22.8 20.5 20.5 20.5 10.0 0 1.0 14.3 10.7 15.1 12.9 11.6 12.3 13.9 29.6 18.4 4.5 12.0 9.9 9.6 13.7 11.5 14.1 13.0 28.3 10.8 13.6]*pi/180;
Fmax = 2*[303.7 399.5 324.2 324.2 324.2 324.2 546.1 780.5 526.1 881.1 616.5 702.0 180.0 190.0 215.0 621.9 479.7 1162.7 301.9 705.2 315.8 848.8 1443.7 2255.4 1024.2 606.4 1308.0 673.7 905.6 305.9 653.3 90.0 3585.9 345.4 274.4];

lmoArnold = lmo;
ltsArnold = lts;

% load scaledlmolts_lmofmax350.mat

lmoDiff = (lmo-lmoArnold)*100
ltsDiff = (lts-ltsArnold)*100

% lmo(18) = lmo(18)*2;

sigma = 610e3;
LegVol = (47*1.7*80.5+1285)/100^3;
VolFraction = [1.47 2.26...
    1.97 1.97 1.97 1.97 ...
    3.52 5.02 3.39 ...
    1.82 1.27 1.45...
    0.45 0.48 0.54...
    2.48 3.80 3.46 2.60 2.92 1.40 3.79 6.06 11.66 3.84 2.11 3.62 1.91 1.49...
    0.53 1.14 0.16 6.21...
    0.97 0.43]/100;
Vmusc = LegVol*VolFraction;

optParams.lmo = lmo;
optParams.lts = lts;
optParams.coefs = coefsOrig;
optParams.nMusc = nMusc;
optParams.alpha = alpha;
optParams.Fmax = Fmax;
optParams.nTrials = nTrials;
optParams.nFrames = nFramesGait;
optParams.Vmusc = Vmusc;
optParams.sigma = sigma;

%% Assign joint angle ranges to data

PassiveM{1} = reshape(PassiveM{1},numel(PassiveM{1})/4,4);
PassiveM{1}(:,4) = [];
% PassiveM{1}(45:end,:) = [];
PassiveM{1} = PassiveM{1}(:);

PassiveM{2} = reshape(PassiveM{2},numel(PassiveM{2})/4,4);
PassiveM{2}(92:121,:) = [];
PassiveM{2} = PassiveM{2}(:);

optParams.MomentsPassive = [PassiveM{1}; PassiveM{2}; PassiveM{3}];

JAngles{1} = reshape(JAngles{1},numel(JAngles{1})/(4*6),4,6);
JAngles{1}(:,4,:) = [];
% JAngles{1}(45:end,:) = [];
JAngles{1} = reshape(JAngles{1},numel(JAngles{1})/(6),6);

JAngles{2} = reshape(JAngles{2},numel(JAngles{2})/(4*6),4,6);
JAngles{2}(92:121,:,:) = [];
JAngles{2} = reshape(JAngles{2},numel(JAngles{2})/(6),6);

nHip = length(JAngles{1});
nKnee = length(JAngles{2});
nAnkle = length(JAngles{3});

optParams.nHip = nHip;
optParams.nKnee = nKnee;
optParams.nAnkle = nAnkle;

JAngles = [JAngles{1}; JAngles{2}; JAngles{3}]*pi/180;

JAngles(:,1) = JAngles(:,1)-0*pi/180;
JAngles(:,2) = JAngles(:,2)-0*pi/180;
JAngles(:,4) = JAngles(:,4)-0*pi/180;

JAngles = [JAngles; JAnglesGait];

optParams.JAngles = JAngles;

%% Calculate passive forces for each joint angle set

JVels = zeros(size(JAngles));
nFrames = length(JAngles);
nFramesAll = nFrames;

onesCol = ones(nFrames,1);
zerosMat = zeros(nFrames,3);

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

MuscRef = zeros(1,nMusc);
MuscRef(1:17) = 1;
MuscRef([18:20 22]) = 2;
MuscRef([21 23:25]) = 3;
MuscRef([26:27]) = 4;
MuscRef([28:35]) = 5;

%% Combine all matrices

for i = 1:nMusc
    if MuscRef(i) == 1
        Mat{1} = [onesCol HipFEMat HipAAMat HipIEMat InteractionFEAA InteractionIEAA InteractionFEIE;...
            0*onesCol HipFEVelMat HipAAVelMat HipIEVelMat InteractionFEAAVel InteractionIEAAVel InteractionFEIEVel;...
            0*onesCol HipFEMAMat zerosMat zerosMat InteractionFEAAMA_FE zerosMat InteractionFEIEMA_FE;...
            0*onesCol zerosMat HipAAMAMat zerosMat InteractionFEAAMA_AA InteractionIEAAMA_AA zerosMat];
    elseif MuscRef(i) == 2
        Mat{2} = [onesCol HipFEMat HipAAMat HipIEMat KneeMat InteractionFEAA InteractionIEAA InteractionFEIE;...
            0*onesCol HipFEVelMat HipAAVelMat HipIEVelMat KneeVelMat InteractionFEAAVel InteractionIEAAVel InteractionFEIEVel;...
            0*onesCol HipFEMAMat zerosMat zerosMat zerosMat InteractionFEAAMA_FE zerosMat InteractionFEIEMA_FE;...
            0*onesCol zerosMat HipAAMAMat zerosMat zerosMat InteractionFEAAMA_AA InteractionIEAAMA_AA zerosMat;...
            0*onesCol zerosMat zerosMat zerosMat KneeMAMat zerosMat zerosMat zerosMat];
    elseif MuscRef(i) == 3
        Mat{3} = [onesCol KneeMat;...
            0*onesCol KneeVelMat;...
            0*onesCol KneeMAMat];
    elseif MuscRef(i) == 4
        Mat{4} = [onesCol KneeMat AnkleMat SubtalarMat InteractionAnkleSub;...
            0*onesCol KneeVelMat AnkleVelMat SubtalarVelMat InteractionAnkleSubVel;...
            0*onesCol KneeMAMat zerosMat zerosMat zerosMat;...
            0*onesCol zerosMat AnkleMAMat zerosMat InteractionAnkleMA;...
            0*onesCol zerosMat zerosMat SubtalarMAMat InteractionSubMA];
    elseif MuscRef(i) == 5
        Mat{5} = [onesCol AnkleMat SubtalarMat InteractionAnkleSub;...
            0*onesCol AnkleVelMat SubtalarVelMat InteractionAnkleSubVel;...
            0*onesCol AnkleMAMat zerosMat InteractionAnkleMA;...
            0*onesCol zerosMat SubtalarMAMat InteractionSubMA];
    end
end

optParams.Mat = Mat;
optParams.MuscRef = MuscRef;
optParams.nframesAll = nFramesAll;
optParams.onesCol = onesCol;

options = optimset('display', 'iter','maxIter',50);
% optionsPSO = psoset('NumParticles',40,'InitializationMethod','lhs','MaxFunEvals',1e6);

parameterChange = zeros(1,70);

for i = 1:numOpt
    
%     lowBound = min([lmo; lts])-0.01;
%     lowBound = -min(lowBound, 0.05);
%     upBound = 0.05*ones(1,35);
    lowBound = 0.5*ones(1,2*35);
    upBound = 2*ones(1,2*35);
    
    parameterChange = lsqnonlin(@Calibratelmolts, ones(1,2*35), lowBound, upBound, options, optParams);
%     parameterChange = pso_parallel(@Calibratelmolts, zeros(1,70), lowBound, upBound, optionsPSO, optParams);
    
    optParams.lmo = optParams.lmo.*parameterChange(1:35);
    optParams.lts = optParams.lts.*parameterChange(36:70);
    
    lmo = optParams.lmo;
    lts = optParams.lts;
    
%     Fmax = (optParams.Vmusc./lmo)*optParams.sigma;
    
    save scaledlmolts1.mat lmo lts
%     load scaledlmolts.mat
end

[moments, lmtilda, vmtilda, HipMA, HipAAMA, KneeMA, AnkleMA, SubtalarMA, muscMoments, passiveF, Lmt] = calcMoments(optParams, 1);

lmo = optParams.lmo;
lts = optParams.lts;

maxlm = max(lmtilda)
minlm = min(lmtilda)

JAnglesH = JAngles(1:nHip,1)*180/pi;
momentsH = moments(1:nHip,1);
PassiveMH = PassiveM{1};

JAnglesH = reshape(JAnglesH,nHip/3,3);
momentsH = reshape(momentsH,nHip/3,3);
PassiveMH = reshape(PassiveMH,nHip/3,3);

JAnglesK = JAngles(nHip+1:nHip+nKnee,3)*180/pi;
momentsK = moments(nHip+1:nHip+nKnee,3);
PassiveMK = PassiveM{2};

JAnglesK = reshape(JAnglesK,nKnee/4,4);
momentsK = reshape(momentsK,nKnee/4,4);
PassiveMK = reshape(PassiveMK,nKnee/4,4);

JAnglesA = JAngles(nHip+nKnee+1:nHip+nKnee+nAnkle,4)*180/pi;
momentsA = moments(nHip+nKnee+1:nHip+nKnee+nAnkle,4);
PassiveMA = PassiveM{3};

JAnglesA = reshape(JAnglesA,nAnkle/4,4);
momentsA = reshape(momentsA,nAnkle/4,4);
PassiveMA = reshape(PassiveMA,nAnkle/4,4);

figure

subplot(3,1,1), plot(JAnglesH,momentsH,'LineWidth',2);
hold on
subplot(3,1,1), plot(JAnglesH,PassiveMH,'--','LineWidth',2);
axis([-15 80 -60 100])
legend('K=15', 'K=60', 'K=90', 'K=110')
title('Hip Passive Moment Matching')

subplot(3,1,2), plot(JAnglesK, momentsK,'LineWidth',2);
hold on
subplot(3,1,2), plot(JAnglesK,PassiveMK,'--','LineWidth',2);
axis([0 120 -30 30])
legend('H=0 A=20', 'H=0 A=-15', 'H=15 A=0', 'H=-15 A=20')
title('Knee Passive Moment Matching')

subplot(3,1,3),plot(JAnglesA, momentsA,'LineWidth',2);
hold on
subplot(3,1,3), plot(JAnglesA,PassiveMA,'--','LineWidth',2);
axis([-30 32 -50 0])
legend('K=0', 'K=15', 'K=60', 'K=110')
title('Ankle Passive Moment Matching')

fig = figure;
bar([lmo-lmoArnold lts-ltsArnold])
set(gca, 'XTick', 1:2*nMusc, 'XTickLabel', [muscLabels muscLabels], 'FontSize', 10)
% rotateXLabels(gca,90)

fig = figure;
bar([(lmo-lmoArnold)./lmoArnold (lts-ltsArnold)./ltsArnold]*100)
set(gca, 'XTick', 1:2*nMusc, 'XTickLabel', [muscLabels muscLabels], 'FontSize', 10)
% rotateXLabels(gca,90)

lmtilda = lmtilda(nHip+nKnee+nAnkle+1:end,:);
lmtilda = reshape(lmtilda,nFramesGait,nTrials,35);
lmtildaMean = permute(mean(lmtilda,2),[1 3 2]);
% lmtildaMin = permute(min(lmtilda,[],2),[1 3 2]);
% lmtildaMax = permute(max(lmtilda,[],2),[1 3 2]);

figure
for i = 1:35
    subplot(6,6,i), plot(lmtildaMean(:,i))
    hold on
    plot(1.2*ones(1,101),'r--')
    plot(0.4*ones(1,101),'r--')
    axis([0 length(lmtilda(:,i))-1 0 1.5])
    title(muscLabels{i})
end

momentErrorH = sum(abs(momentsH(:)-PassiveMH(:)))/numel(momentsH)
momentErrorK = sum(abs(momentsK(:)-PassiveMK(:)))/numel(momentsK)
momentErrorA = sum(abs(momentsA(:)-PassiveMA(:)))/numel(momentsA)

momentsGait = moments(nHip+nKnee+nAnkle+1:end,:);

momentsGait = reshape(momentsGait,nFramesGait,numel(momentsGait)/(5*nFramesGait),5);
momentsGaitMean = permute(mean(momentsGait,2),[1 3 2]);
momentsGaitMax = permute(max(momentsGait,[],2),[1 3 2]);
momentsGaitMin = permute(min(momentsGait,[],2),[1 3 2]);

label = {'Hip FE' 'Hip AA' 'Knee FE' 'Ankle PFDF' 'Ankle IE'};
figure
for i = 1:5
    subplot(1,5,i), plot(0:nFramesGait-1, momentsGaitMean(:,i))
    hold on
    subplot(1,5,i), plot(0:nFramesGait-1, momentsGaitMax(:,i),'r--')
    subplot(1,5,i), plot(0:nFramesGait-1, momentsGaitMin(:,i),'r--')
    title(label{i})
end

%%keyboard

function cost = Calibratelmolts(parameterChange, optParams)

parameterChange = parameterChange(:)';

optParams.lmo = optParams.lmo.*parameterChange(1:35);
optParams.lts = optParams.lts.*parameterChange(36:70);

nHip = optParams.nHip;
nKnee = optParams.nKnee;
nAnkle = optParams.nAnkle;

[moments, lmtilda, ~, ~, ~, ~, ~, ~, ~, ~, Lmt] = calcMoments(optParams, 1);

moments = [moments(1:nHip,1); ...
    moments(nHip+1:nHip+nKnee,3);...
    moments(nHip+nKnee+1:nHip+nKnee+nAnkle,4)];

lmtildaGait = lmtilda(nHip+nKnee+nAnkle+1:end,:);
lmtildaGait = reshape(lmtildaGait,optParams.nFrames,optParams.nTrials,35);
lmtildaGaitMin = permute(min(min(lmtildaGait,[],1),[],2),[3 1 2]);
lmtildaGaitMax = permute(max(max(lmtildaGait,[],1),[],2),[3 1 2]);
% lmtildaGaitMean = permute(mean(lmtildaGait,2),[1 3 2]);
% lmtildaGaitMeanMean = permute(mean(mean(lmtildaGait,1),2),[3 1 2]);

% lmtildaGaitMin = min(lmtildaGait)';%
% lmtildaGaitMax = max(lmtildaGait)';%

% lmtildaPen = exp(-100*(.95*min(Lmt)-optParams.lts))'/10;
lmtildaminPen = exp(-10*(lmtildaGaitMin-.4));
lmtildamaxPen = exp(10*(lmtildaGaitMax-1.2));
% lmtildaminPen2 = exp(10*(lmtildaGaitMin-1));
% lmtildamaxPen2 = exp(-10*(lmtildaGaitMax-.7));

cost = [sqrt(10)*((moments(:)-optParams.MomentsPassive(:))/length(moments)^0.5); (parameterChange(:)-1)/0.25; lmtildaminPen; lmtildamaxPen; (lmtildaGaitMax(:)-1)/0.2; (parameterChange(1:35)'-parameterChange(36:70)')/0.25];

function [moments, lmtilda, vmtilda, HipMA, HipAAMA, KneeMA, AnkleMA, SubtalarMA, muscMoments, passiveF, Lmt, MuscForce] = calcMoments(optParams, DoMomentCalcs)

nMusc = optParams.nMusc;
nframesAll = optParams.nframesAll;
alpha = optParams.alpha;
% HipMA = optParams.HipMA;
% HipAAMA = optParams.HipAAMA;
% KneeMA = optParams.KneeMA;
% AnkleMA = optParams.AnkleMA;
% SubtalarMA = optParams.SubtalarMA;
% Lmt = optParams.Lmt;
% Vmt = optParams.Vmt;
onesCol = optParams.onesCol;
Mat = optParams.Mat;
MuscRef = optParams.MuscRef;
coefs = optParams.coefs;

% MACurvatureOffset = optParams.MACurvatureOffset;
lts = optParams.lts;
lmo = optParams.lmo;
Fmax = optParams.Fmax;
a = zeros(nframesAll,nMusc);

Lmt = zeros(nframesAll,nMusc);
Vmt = zeros(nframesAll,nMusc);
HipMA = zeros(nframesAll,nMusc);
HipAAMA = zeros(nframesAll,nMusc);
KneeMA = zeros(nframesAll,nMusc);
AnkleMA = zeros(nframesAll,nMusc);
SubtalarMA = zeros(nframesAll,nMusc);

k = 1;
for i = 1:nMusc
    if MuscRef(i) == 1
        Vec = Mat{1}*coefs(k:k+18);
        Lmt(:,i) = Vec(1:nframesAll);
        Vmt(:,i) = Vec(nframesAll+1:2*nframesAll);
        HipMA(:,i) = Vec(2*nframesAll+1:3*nframesAll);
        HipAAMA(:,i) = Vec(3*nframesAll+1:4*nframesAll);
        k = k+19;
    elseif MuscRef(i) == 2
        Vec = Mat{2}*coefs(k:k+21);
        Lmt(:,i) = Vec(1:nframesAll);
        Vmt(:,i) = Vec(nframesAll+1:2*nframesAll);
        HipMA(:,i) = Vec(2*nframesAll+1:3*nframesAll);
        HipAAMA(:,i) = Vec(3*nframesAll+1:4*nframesAll);
        KneeMA(:,i) = Vec(4*nframesAll+1:5*nframesAll);
        k = k+22;
    elseif MuscRef(i) == 3
        Vec = Mat{3}*coefs(k:k+3);
        Lmt(:,i) = Vec(1:nframesAll);
        Vmt(:,i) = Vec(nframesAll+1:2*nframesAll);
        KneeMA(:,i) = Vec(2*nframesAll+1:3*nframesAll);
        k = k+4;
    elseif MuscRef(i) == 4
        Vec = Mat{4}*coefs(k:k+12);
        Lmt(:,i) = Vec(1:nframesAll);
        Vmt(:,i) = Vec(nframesAll+1:2*nframesAll);
        KneeMA(:,i) = Vec(2*nframesAll+1:3*nframesAll);
        AnkleMA(:,i) = Vec(3*nframesAll+1:4*nframesAll);
        SubtalarMA(:,i) = Vec(4*nframesAll+1:5*nframesAll);
        k = k+13;
    elseif MuscRef(i) == 5
        Vec = Mat{5}*coefs(k:k+9);
        Lmt(:,i) = Vec(1:nframesAll);
        Vmt(:,i) = Vec(nframesAll+1:2*nframesAll);
        AnkleMA(:,i) = Vec(2*nframesAll+1:3*nframesAll);
        SubtalarMA(:,i) = Vec(3*nframesAll+1:4*nframesAll);
        k = k+10;
    end
end

lmtilda = (Lmt-onesCol*lts)./(onesCol*(lmo.*cos(alpha)));
vmtilda = Vmt./(10*onesCol*(lmo.*cos(alpha)));

Fmax = (optParams.Vmusc./lmo)*optParams.sigma;

if DoMomentCalcs
    passiveF = onesCol*(Fmax.*cos(alpha)).*passiveForce(lmtilda);
    MuscForce = passiveF;
    
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

function muscFLTilda = FLCurve(muscLength)

% muscFLTilda = (1+cos(2*pi*(muscLength-1)))/2;
%
% muscFLTilda(muscLength<.5) = 0;

tact1 = 9.404711698062630e-01;
tdeact1 = 1.073124043936821e+01;
c31 = 9.871922298750729e-01;
c41 = 1.967912928195565e+00;
tact2 = 3.553707328939712e-01;
tdeact2 = 4.446159579601148e+01;
c32 = 6.499995398513716e-01;
c42 = 2.061727424469695e+00;
tact3 = 3.319108437630654e-01;
tdeact3 = 2.990238558320256e+01;
c33 = 1.432634023963042e+00;
c43 = 2.524739871950831e+00;
c5 = 5.000000000000000e-02;

muscFLTilda = tact1.*exp(-tdeact1.*((muscLength-c31).^2).^(c41/2))+tact2.*exp(-tdeact2.*((muscLength-c32).^2).^(c42/2))+tact3.*exp(-tdeact3.*((muscLength-c33).^2).^(c43/2))+c5;

function muscFVTilda = FVCurve(muscVel)
% Coefficients determined manually
e1 = 0.708716646303580;
e2 = 5.65697837083118;
e3 = 0.191884073825674;
e4 = 0.896633651436123;
e5 = -0.165769379169328;
e6 = -0.000623692878996660;

muscFVTilda = -e1.*asinh(-e2.*muscVel-e3)+e4-e5.*(-e2.*muscVel-e3)+e6.*(-e2.*muscVel-e3).^3;

function passiveF = passiveForce(muscLength)

c1 = 3.3845;
c2 = -74.0906;
c3 = -2.4438;

passiveF = c1*exp(c2*(exp(muscLength)).^c3);
