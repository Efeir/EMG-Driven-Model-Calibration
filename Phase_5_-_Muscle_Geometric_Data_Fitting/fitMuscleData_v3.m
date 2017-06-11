function fitMuscleData

close all

load Patient4_MuscleGeoFitData.mat

TrialsStep = 1; % Reduces the number of trials used

Data = MuscleGeoFitData;

nTrialsOrig = 1;
nTrials = ceil(nTrialsOrig/TrialsStep);
nMuscEMG = 16;
nMusc = 35;
nFrames = 1000;
plotResults = 0; % Switch to turn on plots

%% Load data

Time = zeros(nFrames,nTrials);
JAngles = zeros(nFrames,6,nTrials);
JVels = zeros(nFrames,6,nTrials);
HipMA = zeros(nFrames,nMusc,nTrials);
HipAAMA = zeros(nFrames,nMusc,nTrials);
KneeMA = zeros(nFrames,nMusc,nTrials);
AnkleMA = zeros(nFrames,nMusc,nTrials);
SubtalarMA = zeros(nFrames,nMusc,nTrials);
Lmt = zeros(nFrames,nMusc,nTrials);
Vmt = zeros(nFrames,nMusc,nTrials);
IDloads = zeros(nFrames,5,nTrials);

% Store data in a more usable format
for i = 1:TrialsStep:nTrialsOrig
    
    j = ceil(i/TrialsStep);
    
    Time(:,j) = Data(i).Time;
    JAngles(:,:,j) = Data(i).JointAngles(:,:);
    HipMA(:,:,j) = Data(i).HipFEMomentArms(:,:);
    HipAAMA(:,:,j) = Data(i).HipAAMomentArms(:,:);
    KneeMA(:,:,j) = Data(i).KneeMomentArms(:,:);
    AnkleMA(:,:,j) = Data(i).AnkleMomentArms(:,:);
    SubtalarMA(:,:,j) = Data(i).SubtalarMomentArms(:,:);
    Lmt(:,:,j) = Data(i).MuscleTendonLengths(:,:);
    
end

nTrials = 1

%% Set up matrices

JAngles = reshape(permute(JAngles, [1 3 2]), nFrames*nTrials, 6)*pi/180;
HipFEMA = reshape(permute(HipMA, [1 3 2]), nFrames*nTrials, nMusc);
HipAAMA = reshape(permute(HipAAMA, [1 3 2]), nFrames*nTrials, nMusc);
KneeMA = reshape(permute(KneeMA, [1 3 2]), nFrames*nTrials, nMusc);
AnkleMA = reshape(permute(AnkleMA, [1 3 2]), nFrames*nTrials, nMusc);
SubtalarMA = reshape(permute(SubtalarMA, [1 3 2]), nFrames*nTrials, nMusc);
Lmt = reshape(permute(Lmt, [1 3 2]), nFrames*nTrials, nMusc);

JAngles(isnan(JAngles)) = 0;
HipFEMA(isnan(HipFEMA)) = 0;
HipAAMA(isnan(HipAAMA)) = 0;
KneeMA(isnan(KneeMA)) = 0;
AnkleMA(isnan(AnkleMA)) = 0;
SubtalarMA(isnan(SubtalarMA)) = 0;
Lmt(isnan(Lmt)) = 0;

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

%% Perform regression to find coefficients

for i = 1:nMusc
    if MuscRef(i) == 1
        Mat{i} = [onesCol HipFEMat HipAAMat HipIEMat InteractionFEAA InteractionIEAA InteractionFEIE;...
            0*onesCol HipFEMAMat zerosMat zerosMat InteractionFEAAMA_FE zerosMat InteractionFEIEMA_FE;...
            0*onesCol zerosMat HipAAMAMat zerosMat InteractionFEAAMA_AA InteractionIEAAMA_AA zerosMat];
        Vec = [Lmt(:,i); HipFEMA(:,i);HipAAMA(:,i)];
        Coefs{i} = (Mat{i}\Vec);
    elseif MuscRef(i) == 2
        Mat{i} = [onesCol HipFEMat HipAAMat HipIEMat KneeMat InteractionFEAA InteractionIEAA InteractionFEIE;...
            0*onesCol HipFEMAMat zerosMat zerosMat zerosMat InteractionFEAAMA_FE zerosMat InteractionFEIEMA_FE;...
            0*onesCol zerosMat HipAAMAMat zerosMat zerosMat InteractionFEAAMA_AA InteractionIEAAMA_AA zerosMat;...
            0*onesCol zerosMat zerosMat zerosMat KneeMAMat zerosMat zerosMat zerosMat];
        Vec = [Lmt(:,i); HipFEMA(:,i);HipAAMA(:,i);KneeMA(:,i)];
        Coefs{i} = (Mat{i}\Vec);
    elseif MuscRef(i) == 3
        Mat{i} = [onesCol KneeMat;...
            0*onesCol KneeMAMat];
        Vec = [Lmt(:,i); KneeMA(:,i)];
        Coefs{i} = (Mat{i}\Vec);
    elseif MuscRef(i) == 4
        Mat{i} = [onesCol KneeMat AnkleMat SubtalarMat InteractionAnkleSub;...
            0*onesCol KneeMAMat zerosMat zerosMat zerosMat;...
            0*onesCol zerosMat AnkleMAMat zerosMat InteractionAnkleMA;...
            0*onesCol zerosMat zerosMat SubtalarMAMat InteractionSubMA];
        Vec = [Lmt(:,i); KneeMA(:,i);AnkleMA(:,i);SubtalarMA(:,i)];
        Coefs{i} = (Mat{i}\Vec);
    elseif MuscRef(i) == 5
        Mat{i} = [onesCol AnkleMat SubtalarMat InteractionAnkleSub;...
            0*onesCol AnkleMAMat zerosMat InteractionAnkleMA;...
            0*onesCol zerosMat SubtalarMAMat InteractionSubMA];
        Vec = [Lmt(:,i); AnkleMA(:,i);SubtalarMA(:,i)];
        Coefs{i} = (Mat{i}\Vec);
    end
end

%% Calculate fits

for i = 1:nMusc
    if MuscRef(i) == 1
        VecNew = Mat{i}*Coefs{i};
        Lmtnew(:,i) = VecNew(1:nFrames);
        HipFEMAnew(:,i) = VecNew(nFrames+1:2*nFrames);
        HipAAMAnew(:,i) = VecNew(2*nFrames+1:3*nFrames);
    elseif MuscRef(i) == 2
        VecNew = Mat{i}*Coefs{i};
        Lmtnew(:,i) = VecNew(1:nFrames);
        HipFEMAnew(:,i) = VecNew(nFrames+1:2*nFrames);
        HipAAMAnew(:,i) = VecNew(2*nFrames+1:3*nFrames);
        KneeMAnew(:,i) = VecNew(3*nFrames+1:4*nFrames);
    elseif MuscRef(i) == 3
        VecNew = Mat{i}*Coefs{i};
        Lmtnew(:,i) = VecNew(1:nFrames);
        KneeMAnew(:,i) = VecNew(nFrames+1:2*nFrames);
    elseif MuscRef(i) == 4
        VecNew = Mat{i}*Coefs{i};
        Lmtnew(:,i) = VecNew(1:nFrames);
        KneeMAnew(:,i) = VecNew(nFrames+1:2*nFrames);
        AnkleMAnew(:,i) = VecNew(2*nFrames+1:3*nFrames);
        SubtalarMAnew(:,i) = VecNew(3*nFrames+1:4*nFrames);
    elseif MuscRef(i) == 5
        VecNew = Mat{i}*Coefs{i};
        Lmtnew(:,i) = VecNew(1:nFrames);
        AnkleMAnew(:,i) = VecNew(nFrames+1:2*nFrames);
        SubtalarMAnew(:,i) = VecNew(2*nFrames+1:3*nFrames);
    end
end

%% Plot fitting results

Lmt = reshape(Lmt,nFrames,nMusc);
Lmtnew = reshape(Lmtnew,nFrames,nMusc);
Vmtnew = reshape(Vmtnew,nFrames,nMusc);

HipFEMA = reshape(HipFEMA,nFrames,nMusc);
HipAAMA = reshape(HipAAMA,nFrames,nMusc);
KneeMA = reshape(KneeMA,nFrames,nMusc);
AnkleMA = reshape(AnkleMA,nFrames,nMusc);
SubtalarMA = reshape(SubtalarMA,nFrames,nMusc);
HipFEMAnew = reshape(HipFEMAnew,nFrames,nMusc);
HipAAMAnew = reshape(HipAAMAnew,nFrames,nMusc);
KneeMAnew = reshape(KneeMAnew,nFrames,nMusc);
AnkleMAnew = reshape(AnkleMAnew,nFrames,nMusc);
SubtalarMAnew = reshape(SubtalarMAnew,nFrames,nMusc);
OpensimMuscleLabels = Data(1).OpensimMuscleLabels;

LmtError = abs(Lmt-Lmtnew);
HipFEError = abs(HipMA-HipFEMAnew);
HipAAError = abs(HipAAMA-HipAAMAnew);
KneeError = abs(KneeMA-KneeMAnew);
AnkleError = abs(AnkleMA-AnkleMAnew);
SubtalarError = abs(SubtalarMA-SubtalarMAnew);

LmtErrorMed = (median(LmtError))
HipFEErrorMed = (median(HipFEError))
HipAAErrorMed = (median(HipAAError))
KneeErrorMed = (median(KneeError))
AnkleErrorMed = (median(AnkleError))
SubtalarErrorMed = (median(SubtalarError))

LmtError = LmtError(:);
HipFEError = HipFEError(:);
HipAAError = HipAAError(:);
KneeError = KneeError(:);
AnkleError = AnkleError(:);
SubtalarError = SubtalarError(:);

HipFEError(HipFEError==0) = [];
HipAAError(HipAAError==0) = [];
KneeError(KneeError==0) = [];
AnkleError(AnkleError==0) = [];
SubtalarError(SubtalarError==0) = [];

% Make plots of errors at all sample points
figure
subplot(1,5,1), plot(sort(HipFEError))
title('Hip FE Moment Arm Errors')
ylabel('Errors (m)')
% axis([0 1000 0 0.001])
subplot(1,5,2), plot(sort(HipAAError))
title('Hip AA Moment Arm Errors')
% axis([0 1000 0 0.001])
subplot(1,5,3), plot(sort(KneeError))
title('Knee FE Moment Arm Errors')
% axis([0 1000 0 0.001])
subplot(1,5,4), plot(sort(AnkleError))
title('Ankle Moment Arm Errors')
% axis([0 1000 0 0.001])
subplot(1,5,5), plot(sort(SubtalarError))
title('Subtalar Moment Arm Errors')
% axis([0 1000 0 0.001])

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
            subplot(6,6,i), plot(0:999,Lmt(:,j,i))
            hold on
            subplot(6,6,i), plot(0:999,Lmtnew(:,j,i),'g')
            title([OpensimMuscleLabels{i} ' mean lmtilda'])
            %     axis([0 100 .4 1.4]);
            set(gca, 'FontSize', 10)
        end
        pause
        % saveas(fig, 'Musclmtilda.jpg')
        
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        for i = 1:nMusc
            subplot(6,6,i), plot(0:999,HipFEMA(:,j,i))
            hold on
            subplot(6,6,i), plot(0:999,HipFEMAnew(:,j,i),'g')
            title([OpensimMuscleLabels{i} ' mean Hip FE MA'])
            set(gca, 'FontSize', 10)
        end
        pause
        % saveas(fig, 'HipFEMA.jpg')
        
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        for i = 1:nMusc
            subplot(6,6,i), plot(0:999,HipAAMA(:,j,i))
            hold on
            subplot(6,6,i), plot(0:999,HipAAMAnew(:,j,i),'g')
            title([OpensimMuscleLabels{i} ' mean Hip AA MA'])
            set(gca, 'FontSize', 10)
        end
        pause
        % saveas(fig, 'HipAAMA.jpg')
        
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        for i = 1:nMusc
            subplot(6,6,i), plot(0:999,KneeMA(:,j,i))
            hold on
            subplot(6,6,i), plot(0:999,KneeMAnew(:,j,i),'g')
            title([OpensimMuscleLabels{i} ' mean Knee MA'])
            set(gca, 'FontSize', 10)
        end
        pause
        % saveas(fig, 'KneeMA.jpg')
        
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        for i = 1:nMusc
            subplot(6,6,i), plot(0:999,AnkleMA(:,j,i))
            hold on
            subplot(6,6,i), plot(0:999,AnkleMAnew(:,j,i),'g')
            title([OpensimMuscleLabels{i} ' mean Ankle MA'])
            set(gca, 'FontSize', 10)
        end
        pause
        % saveas(fig, 'AnkleMA.jpg')
        
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        for i = 1:nMusc
            subplot(6,6,i), plot(0:999,SubtalarMA(:,j,i))
            hold on
            subplot(6,6,i), plot(0:999,SubtalarMAnew(:,j,i),'g')
            title([OpensimMuscleLabels{i} ' mean Subtalar MA'])
            set(gca, 'FontSize', 10)
        end
        pause
        % saveas(fig, 'SubtalarMA.jpg')
        
        close all
    end
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:nMusc
        subplot(6,6,i), plot(0:999,Lmtmean(:,i))
        hold on
        subplot(6,6,i), plot(0:999,Lmtnewmean(:,i),'g')
        subplot(6,6,i), plot(0:999,Lmtmax(:,i),'r-.')
        subplot(6,6,i), plot(0:999,Lmtmin(:,i),'r-.')
        title([OpensimMuscleLabels{i} ' mean lmtilda'])
        %     axis([0 100 .4 1.4]);
        set(gca, 'FontSize', 10)
    end
    % saveas(fig, 'Musclmtilda.jpg')
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:nMusc
        subplot(6,6,i), plot(0:999,HipFEMAmean(:,i))
        hold on
        subplot(6,6,i), plot(0:999,HipFEMAnewmean(:,i),'g')
        subplot(6,6,i), plot(0:999,HipFEMAmax(:,i),'r-.')
        subplot(6,6,i), plot(0:999,HipFEMAmin(:,i),'r-.')
        title([OpensimMuscleLabels{i} ' mean Hip FE MA'])
        set(gca, 'FontSize', 10)
    end
    % saveas(fig, 'HipFEMA.jpg')
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:nMusc
        subplot(6,6,i), plot(0:999,HipAAMAmean(:,i))
        hold on
        subplot(6,6,i), plot(0:999,HipAAMAnewmean(:,i),'g')
        subplot(6,6,i), plot(0:999,HipAAMAmax(:,i),'r-.')
        subplot(6,6,i), plot(0:999,HipAAMAmin(:,i),'r-.')
        title([OpensimMuscleLabels{i} ' mean Hip AA MA'])
        set(gca, 'FontSize', 10)
    end
    % saveas(fig, 'HipAAMA.jpg')
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:nMusc
        subplot(6,6,i), plot(0:999,KneeMAmean(:,i))
        hold on
        subplot(6,6,i), plot(0:999,KneeMAnewmean(:,i),'g')
        subplot(6,6,i), plot(0:999,KneeMAmax(:,i),'r-.')
        subplot(6,6,i), plot(0:999,KneeMAmin(:,i),'r-.')
        title([OpensimMuscleLabels{i} ' mean Knee MA'])
        set(gca, 'FontSize', 10)
    end
    % saveas(fig, 'KneeMA.jpg')
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:nMusc
        subplot(6,6,i), plot(0:999,AnkleMAmean(:,i))
        hold on
        subplot(6,6,i), plot(0:999,AnkleMAnewmean(:,i),'g')
        subplot(6,6,i), plot(0:999,AnkleMAmax(:,i),'r-.')
        subplot(6,6,i), plot(0:999,AnkleMAmin(:,i),'r-.')
        title([OpensimMuscleLabels{i} ' mean Ankle MA'])
        set(gca, 'FontSize', 10)
    end
    % saveas(fig, 'AnkleMA.jpg')
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:nMusc
        subplot(6,6,i), plot(0:999,SubtalarMAmean(:,i))
        hold on
        subplot(6,6,i), plot(0:999,SubtalarMAnewmean(:,i),'g')
        subplot(6,6,i), plot(0:999,SubtalarMAmax(:,i),'r-.')
        subplot(6,6,i), plot(0:999,SubtalarMAmin(:,i),'r-.')
        title([OpensimMuscleLabels{i} ' mean Subtalar MA'])
        set(gca, 'FontSize', 10)
    end
    % saveas(fig, 'SubtalarMA.jpg')
    
end

%% Save new moment arm curves

% MuscleFits.Mats = Mat;
MuscleFits.Coefs = Coefs;
MuscleFits.MuscRef = MuscRef;

save MuscleFits_lhs.mat MuscleFits

