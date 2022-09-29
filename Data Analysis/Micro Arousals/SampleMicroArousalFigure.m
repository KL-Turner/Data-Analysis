clear; clc;
microArousalFileID = 'MicroArousalData.mat';
load(microArousalFileID)
animalIDs = {'T115','T117','T118'};
EMG = []; normEMG = []; catEMG = []; normCatEMG = [];
diameter = []; normDiameter = []; catDiameter = []; normCatDiameter = [];
delta = []; normDelta = []; catDelta = []; normCatDelta = [];
for aa = 1:length(animalIDs)
    vIDs = fieldnames(Data.(animalIDs{1,aa}));
    for bb = 1:length(vIDs)
        for cc = 1:size(Data.(animalIDs{1,aa}).(vIDs{bb,1}).vesselData,1)
            Data.(animalIDs{1,aa}).(vIDs{bb,1}).normVesselData_long(cc,:) = sgolayfilt(Data.(animalIDs{1,aa}).(vIDs{bb,1}).vesselData_long(cc,:) - mean(Data.(animalIDs{1,aa}).(vIDs{bb,1}).vesselData_long(cc,1:50)),3,17);
            Data.(animalIDs{1,aa}).(vIDs{bb,1}).normEmgData_long(cc,:) = Data.(animalIDs{1,aa}).(vIDs{bb,1}).emgData_long(cc,:) - mean(Data.(animalIDs{1,aa}).(vIDs{bb,1}).emgData_long(cc,1:300));
            Data.(animalIDs{1,aa}).(vIDs{bb,1}).normDeltaData_long(cc,:) = Data.(animalIDs{1,aa}).(vIDs{bb,1}).deltaData_long(cc,:) - mean(Data.(animalIDs{1,aa}).(vIDs{bb,1}).deltaData_long(cc,1:300));
        end
    end
end
for aa = 1:length(animalIDs)
    vIDs = fieldnames(Data.(animalIDs{1,aa}));
    for bb = 1:length(vIDs)
        EMG = cat(1,EMG,mean(Data.(animalIDs{1,aa}).(vIDs{bb,1}).emgData_long,1));
        diameter = cat(1,diameter,mean(Data.(animalIDs{1,aa}).(vIDs{bb,1}).vesselData_long,1));
        delta = cat(1,delta,mean(Data.(animalIDs{1,aa}).(vIDs{bb,1}).deltaData_long,1));
        catEMG = cat(1,catEMG,Data.(animalIDs{1,aa}).(vIDs{bb,1}).emgData_long);
        catDiameter = cat(1,catDiameter,Data.(animalIDs{1,aa}).(vIDs{bb,1}).vesselData_long);
        catDelta = cat(1,catDelta,Data.(animalIDs{1,aa}).(vIDs{bb,1}).deltaData_long);
        normEMG = cat(1,normEMG,mean(Data.(animalIDs{1,aa}).(vIDs{bb,1}).normEmgData_long,1));
        try
            normDiameter = cat(1,normDiameter,mean(Data.(animalIDs{1,aa}).(vIDs{bb,1}).normVesselData_long,1));
        catch
            keyboard
        end
        normDelta = cat(1,normDelta,mean(Data.(animalIDs{1,aa}).(vIDs{bb,1}).normDeltaData_long,1));
        normCatEMG = cat(1,normCatEMG,Data.(animalIDs{1,aa}).(vIDs{bb,1}).normEmgData_long);
        normCatDiameter = cat(1,normCatDiameter,Data.(animalIDs{1,aa}).(vIDs{bb,1}).normVesselData_long);
        normCatDelta = cat(1,normCatDelta,Data.(animalIDs{1,aa}).(vIDs{bb,1}).normDeltaData_long);
    end
end
%%
meanEMG = mean(EMG,1);
stdEMG = std(EMG,0,1);
meanDiameter = mean(diameter,1);
stdDiameter = std(diameter,0,1);
meanDelta = mean(delta,1);
stdDelta = std(delta,0,1);
meanNormEMG = mean(normEMG,1);
stdNormEMG = std(normEMG,0,1);
meanNormDiameter = mean(normDiameter,1);
stdNormDiameter = std(normDiameter,0,1);
meanNormDelta = mean(normDelta,1);
stdNormDelta = std(normDelta,0,1);
meanCatEMG = mean(catEMG,1);
stdCatEMG = std(catEMG,0,1);
meanCatDiameter = mean(catDiameter,1);
stdCatDiameter = std(catDiameter,0,1);
meanCatDelta = mean(catDelta,1);
stdCatDelta = std(catDelta,0,1);
meanNormCatEMG = mean(normCatEMG,1);
stdNormCatEMG = std(normCatEMG,0,1);
meanNormCatDiameter = mean(normCatDiameter,1);
stdNormCatDiameter = std(normCatDiameter,0,1);
meanNormCatDelta = mean(normCatDelta,1);
stdNormCatDelta = std(normCatDelta,0,1);
%%
figure;
sgtitle('vessel average, n = 8 arterioles, n = 3 mice')
subplot(3,1,1);
plot((1:length(meanDiameter))/5,meanDiameter,'r-','LineWidth',2)
hold on
plot((1:length(meanDiameter))/5,meanDiameter + stdDiameter,'r-','LineWidth',0.5)
plot((1:length(meanDiameter))/5,meanDiameter - stdDiameter,'r-','LineWidth',0.5)
xlabel('Time (sec)')
ylabel('\DeltaD/D (%)')
xticks([0,5,10,15,20,25,30,35,40])
xticklabels({'-20','-15','-10','-5','0','5','10','15','20'})
set(gca,'box','off')
axis tight
subplot(3,1,2);
plot((1:length(meanEMG))/30,meanEMG,'k-','LineWidth',2)
hold on
plot((1:length(meanEMG))/30,meanEMG + stdEMG,'k-','LineWidth',0.5)
plot((1:length(meanEMG))/30,meanEMG - stdEMG,'k-','LineWidth',0.5)
xlabel('Time (sec)')
ylabel('Power (a.u.)')
xticks([0,5,10,15,20,25,30,35,40])
xticklabels({'-20','-15','-10','-5','0','5','10','15','20'})
set(gca,'box','off')
axis tight
subplot(3,1,3);
plot((1:length(meanDelta))/30,meanDelta,'b-','LineWidth',2)
hold on
plot((1:length(meanDelta))/30,meanDelta + stdDelta,'b-','LineWidth',0.5)
plot((1:length(meanDelta))/30,meanDelta - stdDelta,'b-','LineWidth',0.5)
xlabel('Time (sec)')
ylabel('\DeltaP/P (%)')
xticks([0,5,10,15,20,25,30,35,40])
xticklabels({'-20','-15','-10','-5','0','5','10','15','20'})
set(gca,'box','off')
axis tight
%%
figure;
sgtitle('vessel average, n = 8 arterioles, n = 3 mice')
for aa = 1:size(catDelta,1)
    hold on
    plot(catDelta(aa,:))
end
xlabel('Time (sec)')
ylabel('\DeltaP/P (%)')
xticks([0,5,10,15,20,25,30,35,40])
xticklabels({'-20','-15','-10','-5','0','5','10','15','20'})
set(gca,'box','off')
axis tight
%%
figure;
sgtitle('-20:0 mean subtracted, vessel average, n = 10 arterioles, n = 4 mice')
subplot(3,1,1);
plot((1:length(meanNormDiameter))/5,meanNormDiameter,'r-','LineWidth',2)
hold on
plot((1:length(meanNormDiameter))/5,meanNormDiameter + stdNormDiameter,'r-','LineWidth',0.5)
plot((1:length(meanNormDiameter))/5,meanNormDiameter - stdNormDiameter,'r-','LineWidth',0.5)
xlabel('Time (sec)')
ylabel('\DeltaD/D (%)')
xticks([0,5,10,15,20,25,30,35,40])
xticklabels({'-20','-15','-10','-5','0','5','10','15','20'})
set(gca,'box','off')
axis tight
subplot(3,1,2);
plot((1:length(meanNormEMG))/30,meanNormEMG,'k-','LineWidth',2)
hold on
plot((1:length(meanNormEMG))/30,meanNormEMG + stdNormEMG,'k-','LineWidth',0.5)
plot((1:length(meanNormEMG))/30,meanNormEMG - stdNormEMG,'k-','LineWidth',0.5)
xlabel('Time (sec)')
ylabel('Power (a.u.)')
xticks([0,5,10,15,20,25,30,35,40])
xticklabels({'-20','-15','-10','-5','0','5','10','15','20'})
set(gca,'box','off')
axis tight
subplot(3,1,3);
plot((1:length(meanNormDelta))/30,meanNormDelta,'b-','LineWidth',2)
hold on
plot((1:length(meanNormDelta))/30,meanNormDelta + stdNormDelta,'b-','LineWidth',0.5)
plot((1:length(meanNormDelta))/30,meanNormDelta - stdNormDelta,'b-','LineWidth',0.5)
xlabel('Time (sec)')
ylabel('\DeltaP/P (%)')
xticks([0,5,10,15,20,25,30,35,40])
xticklabels({'-20','-15','-10','-5','0','5','10','15','20'})
set(gca,'box','off')
axis tight
%%
figure;
sgtitle('event average, n = 44 events')
subplot(3,1,1);
plot((1:length(meanCatDiameter))/5,meanCatDiameter,'r-','LineWidth',2)
hold on
plot((1:length(meanCatDiameter))/5,meanCatDiameter + stdCatDiameter,'r-','LineWidth',0.5)
plot((1:length(meanCatDiameter))/5,meanCatDiameter - stdCatDiameter,'r-','LineWidth',0.5)
xlabel('Time (sec)')
ylabel('\DeltaD/D (%)')
xticks([0,5,10,15,20,25,30,35,40])
xticklabels({'-20','-15','-10','-5','0','5','10','15','20'})
set(gca,'box','off')
axis tight
subplot(3,1,2);
plot((1:length(meanCatEMG))/30,meanCatEMG,'k-','LineWidth',2)
hold on
plot((1:length(meanCatEMG))/30,meanCatEMG + stdCatEMG,'k-','LineWidth',0.5)
plot((1:length(meanCatEMG))/30,meanCatEMG - stdCatEMG,'k-','LineWidth',0.5)
xlabel('Time (sec)')
ylabel('Power (a.u.)')
xticks([0,5,10,15,20,25,30,35,40])
xticklabels({'-20','-15','-10','-5','0','5','10','15','20'})
set(gca,'box','off')
axis tight
subplot(3,1,3);
plot((1:length(meanCatDelta))/30,meanCatDelta,'b-','LineWidth',2)
hold on
plot((1:length(meanCatDelta))/30,meanCatDelta + stdCatDelta,'b-','LineWidth',0.5)
plot((1:length(meanCatDelta))/30,meanCatDelta - stdCatDelta,'b-','LineWidth',0.5)
xlabel('Time (sec)')
ylabel('\DeltaP/P (%)')
xticks([0,5,10,15,20,25,30,35,40])
xticklabels({'-20','-15','-10','-5','0','5','10','15','20'})
set(gca,'box','off')
axis tight
%%
figure;
sgtitle('-20:0 mean subtracted, event average, n = 44 events')
subplot(3,1,1);
plot((1:length(meanNormCatDiameter))/5,meanNormCatDiameter,'r-','LineWidth',2)
hold on
plot((1:length(meanNormCatDiameter))/5,meanNormCatDiameter + stdNormCatDiameter,'r-','LineWidth',0.5)
plot((1:length(meanNormCatDiameter))/5,meanNormCatDiameter - stdNormCatDiameter,'r-','LineWidth',0.5)
xlabel('Time (sec)')
ylabel('\DeltaD/D (%)')
xticks([0,5,10,15,20,25,30,35,40])
xticklabels({'-20','-15','-10','-5','0','5','10','15','20'})
set(gca,'box','off')
axis tight
subplot(3,1,2);
plot((1:length(meanNormCatEMG))/30,meanNormCatEMG,'k-','LineWidth',2)
hold on
plot((1:length(meanNormCatEMG))/30,meanNormCatEMG + stdNormCatEMG,'k-','LineWidth',0.5)
plot((1:length(meanNormCatEMG))/30,meanNormCatEMG - stdNormCatEMG,'k-','LineWidth',0.5)
xlabel('Time (sec)')
ylabel('Power (a.u.)')
xticks([0,5,10,15,20,25,30,35,40])
xticklabels({'-20','-15','-10','-5','0','5','10','15','20'})
set(gca,'box','off')
axis tight
subplot(3,1,3);
plot((1:length(meanDelta))/30,meanNormCatDelta,'b-','LineWidth',2)
hold on
plot((1:length(meanNormCatDelta))/30,meanNormCatDelta + stdNormCatDelta,'b-','LineWidth',0.5)
plot((1:length(meanNormCatDelta))/30,meanNormCatDelta - stdNormCatDelta,'b-','LineWidth',0.5)
xlabel('Time (sec)')
ylabel('\DeltaP/P (%)')
xticks([0,5,10,15,20,25,30,35,40])
xticklabels({'-20','-15','-10','-5','0','5','10','15','20'})
set(gca,'box','off')
axis tight
%%
figure;
subplot(3,1,1);
plot((1:length(meanNormCatDiameter))/5,meanNormCatDiameter,'r-','LineWidth',2)
hold on
plot((1:length(meanNormCatDiameter))/5,meanNormCatDiameter + stdNormCatDiameter,'r-','LineWidth',0.5)
plot((1:length(meanNormCatDiameter))/5,meanNormCatDiameter - stdNormCatDiameter,'r-','LineWidth',0.5)
xlabel('Time (sec)')
ylabel('\DeltaD/D (%)')
xticks([0,5,10,15,20,25,30,35,40])
xticklabels({'-20','-15','-10','-5','0','5','10','15','20'})
set(gca,'box','off')
axis tight
subplot(3,1,2);
plot((1:length(meanNormCatEMG))/30,meanNormCatEMG,'k-','LineWidth',2)
hold on
plot((1:length(meanNormCatEMG))/30,meanNormCatEMG + stdNormCatEMG,'k-','LineWidth',0.5)
plot((1:length(meanNormCatEMG))/30,meanNormCatEMG - stdNormCatEMG,'k-','LineWidth',0.5)
xlabel('Time (sec)')
ylabel('Power (a.u.)')
xticks([0,5,10,15,20,25,30,35,40])
xticklabels({'-20','-15','-10','-5','0','5','10','15','20'})
set(gca,'box','off')
axis tight
subplot(3,1,3);
plot((1:length(meanNormCatDelta))/30,meanNormCatDelta,'b-','LineWidth',2)
hold on
plot((1:length(meanNormCatDelta))/30,meanNormCatDelta + stdNormCatDelta,'b-','LineWidth',0.5)
plot((1:length(meanNormCatDelta))/30,meanNormCatDelta - stdNormCatDelta,'b-','LineWidth',0.5)
xlabel('Time (sec)')
ylabel('\DeltaP/P (%)')
xticks([0,5,10,15,20,25,30,35,40])
xticklabels({'-20','-15','-10','-5','0','5','10','15','20'})
set(gca,'box','off')
axis tight

