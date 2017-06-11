function fitMuscleData

% close all
 
load MuscleFits_lhs.mat
Coefs = MuscleFits.Coefs;

load Patient4_etmData_right_reprocessed.mat

etmData_r = etmData;

load Patient4_etmData_left_reprocessed.mat

etmData_l = etmData;

TrialsStep = 1; % Reduces the number of trials used

Data = [etmData_r.TMGait_0pt2 etmData_l.TMGait_0pt2 ...
    etmData_r.TMGait_0pt3 etmData_l.TMGait_0pt3 ...
    etmData_r.TMGait_0pt4 etmData_l.TMGait_0pt4 ...
    etmData_r.TMGait_0pt5 etmData_l.TMGait_0pt5 ...
    etmData_r.TMGait_0pt6 etmData_l.TMGait_0pt6 ...
    etmData_r.TMGait_0pt7 etmData_l.TMGait_0pt7 ...
    etmData_r.TMGait_0pt8 etmData_l.TMGait_0pt8];

% Data = etmData.TMGait_0pt5;
nTrialsOrig = length(Data);
nTrials = ceil(nTrialsOrig/TrialsStep);
nMuscEMG = 16;
nMusc = 35;
plotResults = 1;

%% Extract original data

Time = zeros(141,nTrials);
JAngles = zeros(141,6,nTrials);
JVels = zeros(141,6,nTrials);
HipMA = zeros(141,nMusc,nTrials);
HipAAMA = zeros(141,nMusc,nTrials);
KneeMA = zeros(141,nMusc,nTrials);
AnkleMA = zeros(141,nMusc,nTrials);
SubtalarMA = zeros(141,nMusc,nTrials);
Lmt = zeros(141,nMusc,nTrials);
Vmt = zeros(141,nMusc,nTrials);
IDloads = zeros(141,5,nTrials);

% Store data in a more usable format
for i = 1:TrialsStep:nTrialsOrig
    
    j = ceil(i/TrialsStep);
    
    Time(:,j) = Data(i).Time;
    JAngles(:,:,j) = Data(i).JointAngles(:,:)*pi/180;
    JVels(:,:,j) = Data(i).JointVelocities(:,:)*pi/180;
    HipMA(:,:,j) = Data(i).HipFEMomentArms(:,:);
    HipAAMA(:,:,j) = Data(i).HipAAMomentArms(:,:);
    KneeMA(:,:,j) = Data(i).KneeMomentArms(:,:);
    AnkleMA(:,:,j) = Data(i).AnkleMomentArms(:,:);
    SubtalarMA(:,:,j) = Data(i).SubtalarMomentArms(:,:);
    Lmt(:,:,j) = Data(i).MuscleTendonLengths(:,:);
    Vmt(:,:,j) = Data(i).MuscleTendonVelocities(:,:);
    IDloads(:,:,j) = Data(i).InverseDynamicsLoads(:,:);
    
end

nTrials = ceil(nTrialsOrig/TrialsStep);

%% Set up matrices again

JAngles = reshape(permute(JAngles, [1 3 2]), 141*nTrials, 6);
JVels = reshape(permute(JVels, [1 3 2]), 141*nTrials, 6);
HipFEMA = reshape(permute(HipMA, [1 3 2]), 141*nTrials, nMusc);
HipAAMA = reshape(permute(HipAAMA, [1 3 2]), 141*nTrials, nMusc);
KneeMA = reshape(permute(KneeMA, [1 3 2]), 141*nTrials, nMusc);
AnkleMA = reshape(permute(AnkleMA, [1 3 2]), 141*nTrials, nMusc);
SubtalarMA = reshape(permute(SubtalarMA, [1 3 2]), 141*nTrials, nMusc);
Lmt = reshape(permute(Lmt, [1 3 2]), 141*nTrials, nMusc);
Vmt = reshape(permute(Vmt, [1 3 2]), 141*nTrials, nMusc);

JAngles(isnan(JAngles)) = 0;
JVels(isnan(JAngles)) = 0;
HipFEMA(isnan(HipFEMA)) = 0;
HipAAMA(isnan(HipAAMA)) = 0;
KneeMA(isnan(KneeMA)) = 0;
AnkleMA(isnan(AnkleMA)) = 0;
SubtalarMA(isnan(SubtalarMA)) = 0;
Lmt(isnan(Lmt)) = 0;
Vmt(isnan(Vmt)) = 0;

nFrames = 141*nTrials;

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

Lmtnew = zeros(size(Lmt));
Vmtnew = zeros(size(Vmt));

HipFEMAnew = zeros(size(HipFEMA));
HipAAMAnew = zeros(size(HipAAMA));
KneeMAnew = zeros(size(KneeMA));
AnkleMAnew = zeros(size(AnkleMA));
SubtalarMAnew = zeros(size(SubtalarMA));
Mat =[];

%% Calculate fits

for i = 1:nMusc
    if MuscRef(i) == 1
        Mat{i} = [onesCol HipFEMat HipAAMat HipIEMat InteractionFEAA InteractionIEAA InteractionFEIE;...
            0*onesCol HipFEVelMat HipAAVelMat HipIEVelMat InteractionFEAAVel InteractionIEAAVel InteractionFEIEVel;...
            0*onesCol HipFEMAMat zerosMat zerosMat InteractionFEAAMA_FE zerosMat InteractionFEIEMA_FE;...
            0*onesCol zerosMat HipAAMAMat zerosMat InteractionFEAAMA_AA InteractionIEAAMA_AA zerosMat];
        VecNew = Mat{i}*Coefs{i};
        Lmtnew(:,i) = VecNew(1:nFrames);
        Vmtnew(:,i) = VecNew(nFrames+1:2*nFrames);
        HipFEMAnew(:,i) = VecNew(2*nFrames+1:3*nFrames);
        HipAAMAnew(:,i) = VecNew(3*nFrames+1:4*nFrames);
    elseif MuscRef(i) == 2
        Mat{i} = [onesCol HipFEMat HipAAMat HipIEMat KneeMat InteractionFEAA InteractionIEAA InteractionFEIE;...
            0*onesCol HipFEVelMat HipAAVelMat HipIEVelMat KneeVelMat InteractionFEAAVel InteractionIEAAVel InteractionFEIEVel;...
            0*onesCol HipFEMAMat zerosMat zerosMat zerosMat InteractionFEAAMA_FE zerosMat InteractionFEIEMA_FE;...
            0*onesCol zerosMat HipAAMAMat zerosMat zerosMat InteractionFEAAMA_AA InteractionIEAAMA_AA zerosMat;...
            0*onesCol zerosMat zerosMat zerosMat KneeMAMat zerosMat zerosMat zerosMat];
        VecNew = Mat{i}*Coefs{i};
        Lmtnew(:,i) = VecNew(1:nFrames);
        Vmtnew(:,i) = VecNew(nFrames+1:2*nFrames);
        HipFEMAnew(:,i) = VecNew(2*nFrames+1:3*nFrames);
        HipAAMAnew(:,i) = VecNew(3*nFrames+1:4*nFrames);
        KneeMAnew(:,i) = VecNew(4*nFrames+1:5*nFrames);
    elseif MuscRef(i) == 3
        Mat{i} = [onesCol KneeMat;...
            0*onesCol KneeVelMat;...
            0*onesCol KneeMAMat];
        VecNew = Mat{i}*Coefs{i};
        Lmtnew(:,i) = VecNew(1:nFrames);
        Vmtnew(:,i) = VecNew(nFrames+1:2*nFrames);
        KneeMAnew(:,i) = VecNew(2*nFrames+1:3*nFrames);
    elseif MuscRef(i) == 4
        Mat{i} = [onesCol KneeMat AnkleMat SubtalarMat InteractionAnkleSub;...
            0*onesCol KneeVelMat AnkleVelMat SubtalarVelMat InteractionAnkleSubVel;...
            0*onesCol KneeMAMat zerosMat zerosMat zerosMat;...
            0*onesCol zerosMat AnkleMAMat zerosMat InteractionAnkleMA;...
            0*onesCol zerosMat zerosMat SubtalarMAMat InteractionSubMA];
        VecNew = Mat{i}*Coefs{i};
        Lmtnew(:,i) = VecNew(1:nFrames);
        Vmtnew(:,i) = VecNew(nFrames+1:2*nFrames);
        KneeMAnew(:,i) = VecNew(2*nFrames+1:3*nFrames);
        AnkleMAnew(:,i) = VecNew(3*nFrames+1:4*nFrames);
        SubtalarMAnew(:,i) = VecNew(4*nFrames+1:5*nFrames);
    elseif MuscRef(i) == 5
        Mat{i} = [onesCol AnkleMat SubtalarMat InteractionAnkleSub;...
            0*onesCol AnkleVelMat SubtalarVelMat InteractionAnkleSubVel;...
            0*onesCol AnkleMAMat zerosMat InteractionAnkleMA;...
            0*onesCol zerosMat SubtalarMAMat InteractionSubMA];
        VecNew = Mat{i}*Coefs{i};
        Lmtnew(:,i) = VecNew(1:nFrames);
        Vmtnew(:,i) = VecNew(nFrames+1:2*nFrames);
        AnkleMAnew(:,i) = VecNew(2*nFrames+1:3*nFrames);
        SubtalarMAnew(:,i) = VecNew(3*nFrames+1:4*nFrames);
    end
end

%% Plot fitting results

Lmt = reshape(Lmt,141,nTrials,nMusc);
Vmt = reshape(Vmt,141,nTrials,nMusc);
Lmtnew = reshape(Lmtnew,141,nTrials,nMusc);
Vmtnew = reshape(Vmtnew,141,nTrials,nMusc);

HipFEMA = reshape(HipFEMA,141,nTrials,nMusc);
HipAAMA = reshape(HipAAMA,141,nTrials,nMusc);
KneeMA = reshape(KneeMA,141,nTrials,nMusc);
AnkleMA = reshape(AnkleMA,141,nTrials,nMusc);
SubtalarMA = reshape(SubtalarMA,141,nTrials,nMusc);
HipFEMAnew = reshape(HipFEMAnew,141,nTrials,nMusc);
HipAAMAnew = reshape(HipAAMAnew,141,nTrials,nMusc);
KneeMAnew = reshape(KneeMAnew,141,nTrials,nMusc);
AnkleMAnew = reshape(AnkleMAnew,141,nTrials,nMusc);
SubtalarMAnew = reshape(SubtalarMAnew,141,nTrials,nMusc);
OpensimMuscleLabels = Data(1).OpensimMuscleLabels;

if plotResults
    
    Lmtmean = permute(mean(Lmt,2),[1 3 2]);
    Vmtmean = permute(mean(Vmt,2),[1 3 2]);
    HipFEMAmean = permute(mean(HipFEMA,2),[1 3 2]);
    HipAAMAmean = permute(mean(HipAAMA,2),[1 3 2]);
    KneeMAmean = permute(mean(KneeMA,2),[1 3 2]);
    AnkleMAmean = permute(mean(AnkleMA,2),[1 3 2]);
    SubtalarMAmean = permute(mean(SubtalarMA,2),[1 3 2]);
    
    Lmtmax = permute(max(Lmt,[],2),[1 3 2]);
    Vmtmax = permute(max(Vmt,[],2),[1 3 2]);
    HipFEMAmax = permute(max(HipFEMA,[],2),[1 3 2]);
    HipAAMAmax = permute(max(HipAAMA,[],2),[1 3 2]);
    KneeMAmax = permute(max(KneeMA,[],2),[1 3 2]);
    AnkleMAmax = permute(max(AnkleMA,[],2),[1 3 2]);
    SubtalarMAmax = permute(max(SubtalarMA,[],2),[1 3 2]);
    
    Lmtmin = permute(min(Lmt,[],2),[1 3 2]);
    Vmtmin = permute(min(Vmt,[],2),[1 3 2]);
    HipFEMAmin = permute(min(HipFEMA,[],2),[1 3 2]);
    HipAAMAmin = permute(min(HipAAMA,[],2),[1 3 2]);
    KneeMAmin = permute(min(KneeMA,[],2),[1 3 2]);
    AnkleMAmin = permute(min(AnkleMA,[],2),[1 3 2]);
    SubtalarMAmin = permute(min(SubtalarMA,[],2),[1 3 2]);
    
    Lmtnewmean = permute(mean(Lmtnew,2),[1 3 2]);
    Vmtnewmean = permute(mean(Vmtnew,2),[1 3 2]);
    HipFEMAnewmean = permute(mean(HipFEMAnew,2),[1 3 2]);
    HipAAMAnewmean = permute(mean(HipAAMAnew,2),[1 3 2]);
    KneeMAnewmean = permute(mean(KneeMAnew,2),[1 3 2]);
    AnkleMAnewmean = permute(mean(AnkleMAnew,2),[1 3 2]);
    SubtalarMAnewmean = permute(mean(SubtalarMAnew,2),[1 3 2]);
    
    Lmtnewmax = permute(max(Lmtnew,[],2),[1 3 2]);
    Vmtnewmax = permute(max(Vmtnew,[],2),[1 3 2]);
    HipFEMAnewmax = permute(max(HipFEMAnew,[],2),[1 3 2]);
    HipAAMAnewmax = permute(max(HipAAMAnew,[],2),[1 3 2]);
    KneeMAnewmax = permute(max(KneeMAnew,[],2),[1 3 2]);
    AnkleMAnewmax = permute(max(AnkleMAnew,[],2),[1 3 2]);
    SubtalarMAnewmax = permute(max(SubtalarMAnew,[],2),[1 3 2]);
    
    Lmtnewmin = permute(min(Lmtnew,[],2),[1 3 2]);
    Vmtnewmin = permute(min(Vmtnew,[],2),[1 3 2]);
    HipFEMAnewmin = permute(min(HipFEMAnew,[],2),[1 3 2]);
    HipAAMAnewmin = permute(min(HipAAMAnew,[],2),[1 3 2]);
    KneeMAnewmin = permute(min(KneeMAnew,[],2),[1 3 2]);
    AnkleMAnewmin = permute(min(AnkleMAnew,[],2),[1 3 2]);
    SubtalarMAnewmin = permute(min(SubtalarMAnew,[],2),[1 3 2]);
    
    for j = 1:25:nTrials
        
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        for i = 1:nMusc
            subplot(6,6,i), plot(0:140,Lmt(:,j,i))
            hold on
            subplot(6,6,i), plot(0:140,Lmtnew(:,j,i),'g')
            title([OpensimMuscleLabels{i} ' mean lmtilda'])
            %     axis([0 100 .4 1.4]);
            set(gca, 'FontSize', 10)
        end
        %pause
        % saveas(fig, 'Musclmtilda.jpg')
        
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        for i = 1:nMusc
            subplot(6,6,i), plot(0:140,Vmt(:,j,i))
            hold on
            subplot(6,6,i), plot(0:140,Vmtnew(:,j,i),'g')
            title([OpensimMuscleLabels{i} ' mean vmtilda'])
            %     axis([0 100 -1 1]);
            set(gca, 'FontSize', 10)
        end
        %pause
        % saveas(fig, 'Muscvmtilda.jpg')
        
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        for i = 1:nMusc
            subplot(6,6,i), plot(0:140,HipFEMA(:,j,i))
            hold on
            subplot(6,6,i), plot(0:140,HipFEMAnew(:,j,i),'g')
            title([OpensimMuscleLabels{i} ' mean Hip FE MA'])
            set(gca, 'FontSize', 10)
        end
        %pause
        % saveas(fig, 'HipFEMA.jpg')
        
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        for i = 1:nMusc
            subplot(6,6,i), plot(0:140,HipAAMA(:,j,i))
            hold on
            subplot(6,6,i), plot(0:140,HipAAMAnew(:,j,i),'g')
            title([OpensimMuscleLabels{i} ' mean Hip AA MA'])
            set(gca, 'FontSize', 10)
        end
        %pause
        % saveas(fig, 'HipAAMA.jpg')
        
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        for i = 1:nMusc
            subplot(6,6,i), plot(0:140,KneeMA(:,j,i))
            hold on
            subplot(6,6,i), plot(0:140,KneeMAnew(:,j,i),'g')
            title([OpensimMuscleLabels{i} ' mean Knee MA'])
            set(gca, 'FontSize', 10)
        end
        %pause
        % saveas(fig, 'KneeMA.jpg')
        
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        for i = 1:nMusc
            subplot(6,6,i), plot(0:140,AnkleMA(:,j,i))
            hold on
            subplot(6,6,i), plot(0:140,AnkleMAnew(:,j,i),'g')
            title([OpensimMuscleLabels{i} ' mean Ankle MA'])
            set(gca, 'FontSize', 10)
        end
        %pause
        % saveas(fig, 'AnkleMA.jpg')
        
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        for i = 1:nMusc
            subplot(6,6,i), plot(0:140,SubtalarMA(:,j,i))
            hold on
            subplot(6,6,i), plot(0:140,SubtalarMAnew(:,j,i),'g')
            title([OpensimMuscleLabels{i} ' mean Subtalar MA'])
            set(gca, 'FontSize', 10)
        end
        %pause
        % saveas(fig, 'SubtalarMA.jpg')
        
        close all
    end
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:nMusc
        subplot(6,6,i), plot(0:140,Lmtmean(:,i))
        hold on
        subplot(6,6,i), plot(0:140,Lmtnewmean(:,i),'g')
        subplot(6,6,i), plot(0:140,Lmtmax(:,i),'r-.')
        subplot(6,6,i), plot(0:140,Lmtmin(:,i),'r-.')
        title([OpensimMuscleLabels{i} ' mean lmtilda'])
        %     axis([0 100 .4 1.4]);
        set(gca, 'FontSize', 10)
    end
    % saveas(fig, 'Musclmtilda.jpg')
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:nMusc
        subplot(6,6,i), plot(0:140,Vmtmean(:,i))
        hold on
        subplot(6,6,i), plot(0:140,Vmtnewmean(:,i),'g')
        subplot(6,6,i), plot(0:140,Vmtmax(:,i),'r-.')
        subplot(6,6,i), plot(0:140,Vmtmin(:,i),'r-.')
        title([OpensimMuscleLabels{i} ' mean vmtilda'])
        %     axis([0 100 -1 1]);
        set(gca, 'FontSize', 10)
    end
    % saveas(fig, 'Muscvmtilda.jpg')
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:nMusc
        subplot(6,6,i), plot(0:140,HipFEMAmean(:,i))
        hold on
        subplot(6,6,i), plot(0:140,HipFEMAnewmean(:,i),'g')
        subplot(6,6,i), plot(0:140,HipFEMAmax(:,i),'r-.')
        subplot(6,6,i), plot(0:140,HipFEMAmin(:,i),'r-.')
        title([OpensimMuscleLabels{i} ' mean Hip FE MA'])
        set(gca, 'FontSize', 10)
    end
    % saveas(fig, 'HipFEMA.jpg')
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:nMusc
        subplot(6,6,i), plot(0:140,HipAAMAmean(:,i))
        hold on
        subplot(6,6,i), plot(0:140,HipAAMAnewmean(:,i),'g')
        subplot(6,6,i), plot(0:140,HipAAMAmax(:,i),'r-.')
        subplot(6,6,i), plot(0:140,HipAAMAmin(:,i),'r-.')
        title([OpensimMuscleLabels{i} ' mean Hip AA MA'])
        set(gca, 'FontSize', 10)
    end
    % saveas(fig, 'HipAAMA.jpg')
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:nMusc
        subplot(6,6,i), plot(0:140,KneeMAmean(:,i))
        hold on
        subplot(6,6,i), plot(0:140,KneeMAnewmean(:,i),'g')
        subplot(6,6,i), plot(0:140,KneeMAmax(:,i),'r-.')
        subplot(6,6,i), plot(0:140,KneeMAmin(:,i),'r-.')
        title([OpensimMuscleLabels{i} ' mean Knee MA'])
        set(gca, 'FontSize', 10)
    end
    % saveas(fig, 'KneeMA.jpg')
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:nMusc
        subplot(6,6,i), plot(0:140,AnkleMAmean(:,i))
        hold on
        subplot(6,6,i), plot(0:140,AnkleMAnewmean(:,i),'g')
        subplot(6,6,i), plot(0:140,AnkleMAmax(:,i),'r-.')
        subplot(6,6,i), plot(0:140,AnkleMAmin(:,i),'r-.')
        title([OpensimMuscleLabels{i} ' mean Ankle MA'])
        set(gca, 'FontSize', 10)
    end
    % saveas(fig, 'AnkleMA.jpg')
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:nMusc
        subplot(6,6,i), plot(0:140,SubtalarMAmean(:,i))
        hold on
        subplot(6,6,i), plot(0:140,SubtalarMAnewmean(:,i),'g')
        subplot(6,6,i), plot(0:140,SubtalarMAmax(:,i),'r-.')
        subplot(6,6,i), plot(0:140,SubtalarMAmin(:,i),'r-.')
        title([OpensimMuscleLabels{i} ' mean Subtalar MA'])
        set(gca, 'FontSize', 10)
    end
    % saveas(fig, 'SubtalarMA.jpg')
    
end

%% Save new moment arm curves

numTrialsGaitType = length(etmData_r.TMGait_0pt2);

for i = 1:numTrialsGaitType
    
    j = ceil(i/TrialsStep);
    
    etmData_r.TMGait_0pt2(i).HipFEMomentArms = permute(HipFEMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt2(i).HipAAMomentArms = permute(HipAAMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt2(i).KneeMomentArms = permute(KneeMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt2(i).AnkleMomentArms = permute(AnkleMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt2(i).SubtalarMomentArms = permute(SubtalarMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt2(i).MuscleTendonLengths = permute(Lmtnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt2(i).MuscleTendonVelocities = permute(Vmtnew(:,j,:),[1 3 2]);
    
end

HipFEMAnew(:,1:numTrialsGaitType,:) = [];
HipAAMAnew(:,1:numTrialsGaitType,:) = [];
KneeMAnew(:,1:numTrialsGaitType,:) = [];
AnkleMAnew(:,1:numTrialsGaitType,:) = [];
SubtalarMAnew(:,1:numTrialsGaitType,:) = [];
Lmtnew(:,1:numTrialsGaitType,:) = [];
Vmtnew(:,1:numTrialsGaitType,:) = [];

numTrialsGaitType = length(etmData_l.TMGait_0pt2);

for i = 1:numTrialsGaitType
    
    j = ceil(i/TrialsStep);
    
    etmData_l.TMGait_0pt2(i).HipFEMomentArms = permute(HipFEMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt2(i).HipAAMomentArms = permute(HipAAMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt2(i).KneeMomentArms = permute(KneeMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt2(i).AnkleMomentArms = permute(AnkleMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt2(i).SubtalarMomentArms = permute(SubtalarMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt2(i).MuscleTendonLengths = permute(Lmtnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt2(i).MuscleTendonVelocities = permute(Vmtnew(:,j,:),[1 3 2]);
    
end

HipFEMAnew(:,1:numTrialsGaitType,:) = [];
HipAAMAnew(:,1:numTrialsGaitType,:) = [];
KneeMAnew(:,1:numTrialsGaitType,:) = [];
AnkleMAnew(:,1:numTrialsGaitType,:) = [];
SubtalarMAnew(:,1:numTrialsGaitType,:) = [];
Lmtnew(:,1:numTrialsGaitType,:) = [];
Vmtnew(:,1:numTrialsGaitType,:) = [];

numTrialsGaitType = length(etmData_r.TMGait_0pt3);

for i = 1:numTrialsGaitType
    
    j = ceil(i/TrialsStep);
    
    etmData_r.TMGait_0pt3(i).HipFEMomentArms = permute(HipFEMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt3(i).HipAAMomentArms = permute(HipAAMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt3(i).KneeMomentArms = permute(KneeMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt3(i).AnkleMomentArms = permute(AnkleMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt3(i).SubtalarMomentArms = permute(SubtalarMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt3(i).MuscleTendonLengths = permute(Lmtnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt3(i).MuscleTendonVelocities = permute(Vmtnew(:,j,:),[1 3 2]);
    
end

HipFEMAnew(:,1:numTrialsGaitType,:) = [];
HipAAMAnew(:,1:numTrialsGaitType,:) = [];
KneeMAnew(:,1:numTrialsGaitType,:) = [];
AnkleMAnew(:,1:numTrialsGaitType,:) = [];
SubtalarMAnew(:,1:numTrialsGaitType,:) = [];
Lmtnew(:,1:numTrialsGaitType,:) = [];
Vmtnew(:,1:numTrialsGaitType,:) = [];

numTrialsGaitType = length(etmData_l.TMGait_0pt3);

for i = 1:numTrialsGaitType
    
    j = ceil(i/TrialsStep);
    
    etmData_l.TMGait_0pt3(i).HipFEMomentArms = permute(HipFEMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt3(i).HipAAMomentArms = permute(HipAAMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt3(i).KneeMomentArms = permute(KneeMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt3(i).AnkleMomentArms = permute(AnkleMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt3(i).SubtalarMomentArms = permute(SubtalarMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt3(i).MuscleTendonLengths = permute(Lmtnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt3(i).MuscleTendonVelocities = permute(Vmtnew(:,j,:),[1 3 2]);
    
end

HipFEMAnew(:,1:numTrialsGaitType,:) = [];
HipAAMAnew(:,1:numTrialsGaitType,:) = [];
KneeMAnew(:,1:numTrialsGaitType,:) = [];
AnkleMAnew(:,1:numTrialsGaitType,:) = [];
SubtalarMAnew(:,1:numTrialsGaitType,:) = [];
Lmtnew(:,1:numTrialsGaitType,:) = [];
Vmtnew(:,1:numTrialsGaitType,:) = [];

numTrialsGaitType = length(etmData_r.TMGait_0pt4);

for i = 1:numTrialsGaitType
    
    j = ceil(i/TrialsStep);
    
    etmData_r.TMGait_0pt4(i).HipFEMomentArms = permute(HipFEMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt4(i).HipAAMomentArms = permute(HipAAMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt4(i).KneeMomentArms = permute(KneeMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt4(i).AnkleMomentArms = permute(AnkleMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt4(i).SubtalarMomentArms = permute(SubtalarMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt4(i).MuscleTendonLengths = permute(Lmtnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt4(i).MuscleTendonVelocities = permute(Vmtnew(:,j,:),[1 3 2]);
    
end

HipFEMAnew(:,1:numTrialsGaitType,:) = [];
HipAAMAnew(:,1:numTrialsGaitType,:) = [];
KneeMAnew(:,1:numTrialsGaitType,:) = [];
AnkleMAnew(:,1:numTrialsGaitType,:) = [];
SubtalarMAnew(:,1:numTrialsGaitType,:) = [];
Lmtnew(:,1:numTrialsGaitType,:) = [];
Vmtnew(:,1:numTrialsGaitType,:) = [];

numTrialsGaitType = length(etmData_l.TMGait_0pt4);

for i = 1:numTrialsGaitType
    
    j = ceil(i/TrialsStep);
    
    etmData_l.TMGait_0pt4(i).HipFEMomentArms = permute(HipFEMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt4(i).HipAAMomentArms = permute(HipAAMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt4(i).KneeMomentArms = permute(KneeMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt4(i).AnkleMomentArms = permute(AnkleMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt4(i).SubtalarMomentArms = permute(SubtalarMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt4(i).MuscleTendonLengths = permute(Lmtnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt4(i).MuscleTendonVelocities = permute(Vmtnew(:,j,:),[1 3 2]);
    
end

HipFEMAnew(:,1:numTrialsGaitType,:) = [];
HipAAMAnew(:,1:numTrialsGaitType,:) = [];
KneeMAnew(:,1:numTrialsGaitType,:) = [];
AnkleMAnew(:,1:numTrialsGaitType,:) = [];
SubtalarMAnew(:,1:numTrialsGaitType,:) = [];
Lmtnew(:,1:numTrialsGaitType,:) = [];
Vmtnew(:,1:numTrialsGaitType,:) = [];

numTrialsGaitType = length(etmData_r.TMGait_0pt5);

for i = 1:numTrialsGaitType
    
    j = ceil(i/TrialsStep);
    
    etmData_r.TMGait_0pt5(i).HipFEMomentArms = permute(HipFEMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt5(i).HipAAMomentArms = permute(HipAAMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt5(i).KneeMomentArms = permute(KneeMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt5(i).AnkleMomentArms = permute(AnkleMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt5(i).SubtalarMomentArms = permute(SubtalarMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt5(i).MuscleTendonLengths = permute(Lmtnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt5(i).MuscleTendonVelocities = permute(Vmtnew(:,j,:),[1 3 2]);
    
end

HipFEMAnew(:,1:numTrialsGaitType,:) = [];
HipAAMAnew(:,1:numTrialsGaitType,:) = [];
KneeMAnew(:,1:numTrialsGaitType,:) = [];
AnkleMAnew(:,1:numTrialsGaitType,:) = [];
SubtalarMAnew(:,1:numTrialsGaitType,:) = [];
Lmtnew(:,1:numTrialsGaitType,:) = [];
Vmtnew(:,1:numTrialsGaitType,:) = [];

numTrialsGaitType = length(etmData_l.TMGait_0pt5);

for i = 1:numTrialsGaitType
    
    j = ceil(i/TrialsStep);
    
    etmData_l.TMGait_0pt5(i).HipFEMomentArms = permute(HipFEMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt5(i).HipAAMomentArms = permute(HipAAMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt5(i).KneeMomentArms = permute(KneeMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt5(i).AnkleMomentArms = permute(AnkleMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt5(i).SubtalarMomentArms = permute(SubtalarMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt5(i).MuscleTendonLengths = permute(Lmtnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt5(i).MuscleTendonVelocities = permute(Vmtnew(:,j,:),[1 3 2]);
    
end

HipFEMAnew(:,1:numTrialsGaitType,:) = [];
HipAAMAnew(:,1:numTrialsGaitType,:) = [];
KneeMAnew(:,1:numTrialsGaitType,:) = [];
AnkleMAnew(:,1:numTrialsGaitType,:) = [];
SubtalarMAnew(:,1:numTrialsGaitType,:) = [];
Lmtnew(:,1:numTrialsGaitType,:) = [];
Vmtnew(:,1:numTrialsGaitType,:) = [];

numTrialsGaitType = length(etmData_r.TMGait_0pt6);

for i = 1:numTrialsGaitType
    
    j = ceil(i/TrialsStep);
    
    etmData_r.TMGait_0pt6(i).HipFEMomentArms = permute(HipFEMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt6(i).HipAAMomentArms = permute(HipAAMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt6(i).KneeMomentArms = permute(KneeMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt6(i).AnkleMomentArms = permute(AnkleMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt6(i).SubtalarMomentArms = permute(SubtalarMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt6(i).MuscleTendonLengths = permute(Lmtnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt6(i).MuscleTendonVelocities = permute(Vmtnew(:,j,:),[1 3 2]);
    
end

HipFEMAnew(:,1:numTrialsGaitType,:) = [];
HipAAMAnew(:,1:numTrialsGaitType,:) = [];
KneeMAnew(:,1:numTrialsGaitType,:) = [];
AnkleMAnew(:,1:numTrialsGaitType,:) = [];
SubtalarMAnew(:,1:numTrialsGaitType,:) = [];
Lmtnew(:,1:numTrialsGaitType,:) = [];
Vmtnew(:,1:numTrialsGaitType,:) = [];

numTrialsGaitType = length(etmData_l.TMGait_0pt6);

for i = 1:numTrialsGaitType
    
    j = ceil(i/TrialsStep);
    
    etmData_l.TMGait_0pt6(i).HipFEMomentArms = permute(HipFEMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt6(i).HipAAMomentArms = permute(HipAAMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt6(i).KneeMomentArms = permute(KneeMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt6(i).AnkleMomentArms = permute(AnkleMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt6(i).SubtalarMomentArms = permute(SubtalarMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt6(i).MuscleTendonLengths = permute(Lmtnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt6(i).MuscleTendonVelocities = permute(Vmtnew(:,j,:),[1 3 2]);
    
end

HipFEMAnew(:,1:numTrialsGaitType,:) = [];
HipAAMAnew(:,1:numTrialsGaitType,:) = [];
KneeMAnew(:,1:numTrialsGaitType,:) = [];
AnkleMAnew(:,1:numTrialsGaitType,:) = [];
SubtalarMAnew(:,1:numTrialsGaitType,:) = [];
Lmtnew(:,1:numTrialsGaitType,:) = [];
Vmtnew(:,1:numTrialsGaitType,:) = [];

numTrialsGaitType = length(etmData_r.TMGait_0pt7);

for i = 1:numTrialsGaitType
    
    j = ceil(i/TrialsStep);
    
    etmData_r.TMGait_0pt7(i).HipFEMomentArms = permute(HipFEMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt7(i).HipAAMomentArms = permute(HipAAMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt7(i).KneeMomentArms = permute(KneeMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt7(i).AnkleMomentArms = permute(AnkleMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt7(i).SubtalarMomentArms = permute(SubtalarMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt7(i).MuscleTendonLengths = permute(Lmtnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt7(i).MuscleTendonVelocities = permute(Vmtnew(:,j,:),[1 3 2]);
    
end

HipFEMAnew(:,1:numTrialsGaitType,:) = [];
HipAAMAnew(:,1:numTrialsGaitType,:) = [];
KneeMAnew(:,1:numTrialsGaitType,:) = [];
AnkleMAnew(:,1:numTrialsGaitType,:) = [];
SubtalarMAnew(:,1:numTrialsGaitType,:) = [];
Lmtnew(:,1:numTrialsGaitType,:) = [];
Vmtnew(:,1:numTrialsGaitType,:) = [];

numTrialsGaitType = length(etmData_l.TMGait_0pt7);

for i = 1:numTrialsGaitType
    
    j = ceil(i/TrialsStep);
    
    etmData_l.TMGait_0pt7(i).HipFEMomentArms = permute(HipFEMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt7(i).HipAAMomentArms = permute(HipAAMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt7(i).KneeMomentArms = permute(KneeMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt7(i).AnkleMomentArms = permute(AnkleMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt7(i).SubtalarMomentArms = permute(SubtalarMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt7(i).MuscleTendonLengths = permute(Lmtnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt7(i).MuscleTendonVelocities = permute(Vmtnew(:,j,:),[1 3 2]);
    
end

HipFEMAnew(:,1:numTrialsGaitType,:) = [];
HipAAMAnew(:,1:numTrialsGaitType,:) = [];
KneeMAnew(:,1:numTrialsGaitType,:) = [];
AnkleMAnew(:,1:numTrialsGaitType,:) = [];
SubtalarMAnew(:,1:numTrialsGaitType,:) = [];
Lmtnew(:,1:numTrialsGaitType,:) = [];
Vmtnew(:,1:numTrialsGaitType,:) = [];

numTrialsGaitType = length(etmData_r.TMGait_0pt8);

for i = 1:numTrialsGaitType
    
    j = ceil(i/TrialsStep);
    
    try etmData_r.TMGait_0pt8(i).HipFEMomentArms = permute(HipFEMAnew(:,j,:),[1 3 2]);
    catch asdf = 0;
        keyboard
    end
    etmData_r.TMGait_0pt8(i).HipAAMomentArms = permute(HipAAMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt8(i).KneeMomentArms = permute(KneeMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt8(i).AnkleMomentArms = permute(AnkleMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt8(i).SubtalarMomentArms = permute(SubtalarMAnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt8(i).MuscleTendonLengths = permute(Lmtnew(:,j,:),[1 3 2]);
    etmData_r.TMGait_0pt8(i).MuscleTendonVelocities = permute(Vmtnew(:,j,:),[1 3 2]);
    
end

HipFEMAnew(:,1:numTrialsGaitType,:) = [];
HipAAMAnew(:,1:numTrialsGaitType,:) = [];
KneeMAnew(:,1:numTrialsGaitType,:) = [];
AnkleMAnew(:,1:numTrialsGaitType,:) = [];
SubtalarMAnew(:,1:numTrialsGaitType,:) = [];
Lmtnew(:,1:numTrialsGaitType,:) = [];
Vmtnew(:,1:numTrialsGaitType,:) = [];

numTrialsGaitType = length(etmData_l.TMGait_0pt8);

for i = 1:numTrialsGaitType
    
    j = ceil(i/TrialsStep);
    
    etmData_l.TMGait_0pt8(i).HipFEMomentArms = permute(HipFEMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt8(i).HipAAMomentArms = permute(HipAAMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt8(i).KneeMomentArms = permute(KneeMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt8(i).AnkleMomentArms = permute(AnkleMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt8(i).SubtalarMomentArms = permute(SubtalarMAnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt8(i).MuscleTendonLengths = permute(Lmtnew(:,j,:),[1 3 2]);
    etmData_l.TMGait_0pt8(i).MuscleTendonVelocities = permute(Vmtnew(:,j,:),[1 3 2]);
    
end

HipFEMAnew(:,1:numTrialsGaitType,:) = [];
HipAAMAnew(:,1:numTrialsGaitType,:) = [];
KneeMAnew(:,1:numTrialsGaitType,:) = [];
AnkleMAnew(:,1:numTrialsGaitType,:) = [];
SubtalarMAnew(:,1:numTrialsGaitType,:) = [];
Lmtnew(:,1:numTrialsGaitType,:) = [];
Vmtnew(:,1:numTrialsGaitType,:) = [];

save etmData_refitted_both_4.mat etmData_l etmData_r

