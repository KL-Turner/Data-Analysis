function [] = BlinkResponses_Pupil(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

resultsStruct = 'Results_BlinkResponses';
load(resultsStruct);
animalIDs = fieldnames(Results_BlinkResponses);
timeVector = (0:20*30)/30 - 10;
data.Awake.zDiameter = []; data.Awake.whisk = []; data.Awake.HbT = []; data.Awake.cort = []; data.Awake.hip = []; data.Awake.EMG = [];
data.Asleep.zDiameter = []; data.Asleep.whisk = []; data.Asleep.HbT = []; data.Asleep.cort = []; data.Asleep.hip = []; data.Asleep.EMG = [];
data.Awake.zDiameter_T = []; data.Awake.whisk_T = []; data.Awake.HbT_T = []; data.Awake.cort_T = []; data.Awake.hip_T = []; data.Awake.EMG_T = [];
data.Asleep.zDiameter_T = []; data.Asleep.whisk_T = []; data.Asleep.HbT_T = []; data.Asleep.cort_T = []; data.Asleep.hip_T = []; data.Asleep.EMG_T = [];
data.Awake.zDiameter_F = []; data.Awake.whisk_F = []; data.Awake.HbT_F = []; data.Awake.cort_F = []; data.Awake.hip_F = []; data.Awake.EMG_F = [];
data.Asleep.zDiameter_F = []; data.Asleep.whisk_F = []; data.Asleep.HbT_F = []; data.Asleep.cort_F = []; data.Asleep.hip_F = []; data.Asleep.EMG_F = [];
blinkStates = {'Awake','Asleep'};
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    for bb = 1:length(blinkStates)
        blinkState = blinkStates{1,bb};
        %%
        if isempty(Results_BlinkResponses.(animalID).(blinkState).zDiameter) == false
            data.(blinkState).zDiameter  = cat(1,data.(blinkState).zDiameter,Results_BlinkResponses.(animalID).(blinkState).zDiameter);
        end
        data.(blinkState).whisk  = cat(1,data.(blinkState).whisk,Results_BlinkResponses.(animalID).(blinkState).whisk);
        data.(blinkState).HbT  = cat(1,data.(blinkState).HbT,Results_BlinkResponses.(animalID).(blinkState).LH_HbT,Results_BlinkResponses.(animalID).(blinkState).RH_HbT);
        data.(blinkState).cort = cat(3,data.(blinkState).cort,Results_BlinkResponses.(animalID).(blinkState).LH_cort,Results_BlinkResponses.(animalID).(blinkState).RH_cort);
        data.(blinkState).hip = cat(3,data.(blinkState).hip,Results_BlinkResponses.(animalID).(blinkState).hip);
        data.(blinkState).EMG = cat(1,data.(blinkState).EMG,Results_BlinkResponses.(animalID).(blinkState).EMG);
        %%
        if isempty(Results_BlinkResponses.(animalID).(blinkState).zDiameter_T) == false
            data.(blinkState).zDiameter_T  = cat(1,data.(blinkState).zDiameter_T,Results_BlinkResponses.(animalID).(blinkState).zDiameter_T);
        end
        if isempty(Results_BlinkResponses.(animalID).(blinkState).whisk_T) == false
            data.(blinkState).whisk_T  = cat(1,data.(blinkState).whisk_T,Results_BlinkResponses.(animalID).(blinkState).whisk_T);
            data.(blinkState).HbT_T  = cat(1,data.(blinkState).HbT_T,Results_BlinkResponses.(animalID).(blinkState).LH_HbT_T,Results_BlinkResponses.(animalID).(blinkState).RH_HbT_T);
            data.(blinkState).cort_T = cat(3,data.(blinkState).cort_T,Results_BlinkResponses.(animalID).(blinkState).LH_cort_T,Results_BlinkResponses.(animalID).(blinkState).RH_cort_T);
            data.(blinkState).hip_T = cat(3,data.(blinkState).hip_T,Results_BlinkResponses.(animalID).(blinkState).hip_T);
            data.(blinkState).EMG_T = cat(1,data.(blinkState).EMG_T,Results_BlinkResponses.(animalID).(blinkState).EMG_T);
        end
        %%
        if isempty(Results_BlinkResponses.(animalID).(blinkState).zDiameter_F) == false
            data.(blinkState).zDiameter_F  = cat(1,data.(blinkState).zDiameter_F,Results_BlinkResponses.(animalID).(blinkState).zDiameter_F);
        end
        if isempty(Results_BlinkResponses.(animalID).(blinkState).whisk_F) == false
            data.(blinkState).whisk_F  = cat(1,data.(blinkState).whisk_F,Results_BlinkResponses.(animalID).(blinkState).whisk_F);
            data.(blinkState).HbT_F  = cat(1,data.(blinkState).HbT_F,Results_BlinkResponses.(animalID).(blinkState).LH_HbT_F,Results_BlinkResponses.(animalID).(blinkState).RH_HbT_F);
            data.(blinkState).cort_F = cat(3,data.(blinkState).cort_F,Results_BlinkResponses.(animalID).(blinkState).LH_cort_F,Results_BlinkResponses.(animalID).(blinkState).RH_cort_F);
            data.(blinkState).hip_F = cat(3,data.(blinkState).hip_F,Results_BlinkResponses.(animalID).(blinkState).hip_F);
            data.(blinkState).EMG_F = cat(1,data.(blinkState).EMG_F,Results_BlinkResponses.(animalID).(blinkState).EMG_F);
        end
        T = Results_BlinkResponses.(animalID).(blinkState).T;
        F = Results_BlinkResponses.(animalID).(blinkState).F;
    end
end
%
for bb = 1:length(blinkStates)
    blinkState = blinkStates{1,bb};
    data.(blinkState).meanDiameter = mean(data.(blinkState).zDiameter,1);
    data.(blinkState).stdDiameter = std(data.(blinkState).zDiameter,0,1)./sqrt(size(data.(blinkState).zDiameter,1));
    data.(blinkState).meanHbT = mean(data.(blinkState).HbT,1);
    data.(blinkState).stdHbT = std(data.(blinkState).HbT,0,1)./sqrt(size(data.(blinkState).HbT,1));
    data.(blinkState).meanCort = mean(data.(blinkState).cort,3);
    data.(blinkState).meanHip = mean(data.(blinkState).hip,3);
    data.(blinkState).meanEMG = mean(data.(blinkState).EMG,1);
    data.(blinkState).stdEMG = std(data.(blinkState).EMG,0,1)./sqrt(size(data.(blinkState).EMG,1));
    data.(blinkState).meanWhisk = mean(data.(blinkState).whisk*100,1);
    data.(blinkState).stdWhisk = std(data.(blinkState).whisk*100,0,1)./sqrt(size(data.(blinkState).whisk,1));
    
    data.(blinkState).meanDiameter_T = mean(data.(blinkState).zDiameter_T,1);
    data.(blinkState).stdDiameter_T = std(data.(blinkState).zDiameter_T,0,1)./sqrt(size(data.(blinkState).zDiameter_T,1));
    data.(blinkState).meanHbT_T = mean(data.(blinkState).HbT_T,1);
    data.(blinkState).stdHbT_T = std(data.(blinkState).HbT_T,0,1)./sqrt(size(data.(blinkState).HbT_T,1));
    data.(blinkState).meanCort_T = mean(data.(blinkState).cort_T,3);
    data.(blinkState).meanHip_T = mean(data.(blinkState).hip_T,3);
    data.(blinkState).meanEMG_T = mean(data.(blinkState).EMG_T,1);
    data.(blinkState).stdEMG_T = std(data.(blinkState).EMG_T,0,1)./sqrt(size(data.(blinkState).EMG_T,1));
    data.(blinkState).meanWhisk_T = mean(data.(blinkState).whisk_T*100,1);
    data.(blinkState).stdWhisk_T = std(data.(blinkState).whisk_T*100,0,1)./sqrt(size(data.(blinkState).whisk_T,1));
    
    data.(blinkState).meanDiameter_F = mean(data.(blinkState).zDiameter_F,1);
    data.(blinkState).stdDiameter_F = std(data.(blinkState).zDiameter_F,0,1)./sqrt(size(data.(blinkState).zDiameter_F,1));
    data.(blinkState).meanHbT_F = mean(data.(blinkState).HbT_F,1);
    data.(blinkState).stdHbT_F = std(data.(blinkState).HbT_F,0,1)./sqrt(size(data.(blinkState).HbT_F,1));
    data.(blinkState).meanCort_F = mean(data.(blinkState).cort_F,3);
    data.(blinkState).meanHip_F = mean(data.(blinkState).hip_F,3);
    data.(blinkState).meanEMG_F = mean(data.(blinkState).EMG_F,1);
    data.(blinkState).stdEMG_F = std(data.(blinkState).EMG_F,0,1)./sqrt(size(data.(blinkState).EMG_F,1));
    data.(blinkState).meanWhisk_F = mean(data.(blinkState).whisk_F*100,1);
    data.(blinkState).stdWhisk_F = std(data.(blinkState).whisk_F*100,0,1)./sqrt(size(data.(blinkState).whisk_F,1));
end
%
%% WHISK BEHAVIOR
summaryFigure = figure;
subplot(2,3,1);
p1 = plot(timeVector,data.Awake.meanWhisk,'LineWidth',2);
hold on
p2 = plot(timeVector,data.Awake.meanWhisk_T,'LineWidth',2);
p3 = plot(timeVector,data.Awake.meanWhisk_F,'LineWidth',2);
title('Awake Whisk Probability')
ylabel('Percentage (%)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
legend([p1,p2,p3],'All blinks','Low whisk','High Whisk')
axis square
%%
subplot(2,3,2);
plot(timeVector,data.Awake.meanDiameter,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.Awake.meanDiameter + data.Awake.stdDiameter,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.Awake.meanDiameter - data.Awake.stdDiameter,'color',colors('smoky black'),'LineWidth',0.5)
title('Awake Pupil Diameter')
ylabel('Z Units')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
axis square

figure;
for aa = 1:size(data.Awake.HbT_T,1)
    hold on
    plot(data.Awake.HbT_T(aa,:))
end
xline(300)
xline(330)

figure;
for aa = 1:size(data.Awake.HbT_T,1)
    hold on
    plot(data.Awake.HbT_T(aa,:) - mean(data.Awake.HbT(aa,1))) 
end
xline(300)
xline(330)

 
%%
ax8 = subplot(6,6,8);
plot(timeVector,data.Awake.meanDiameter_T,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.Awake.meanDiameter_T + data.Awake.stdDiameter_T,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.Awake.meanDiameter_T - data.Awake.stdDiameter_T,'color',colors('smoky black'),'LineWidth',0.5)
title('Awake Pupil Diameter')
ylabel('Z Units')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
axis square
ax9 = subplot(6,6,9);
plot(timeVector,data.Awake.meanDiameter_F,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.Awake.meanDiameter_F + data.Awake.stdDiameter_F,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.Awake.meanDiameter_F - data.Awake.stdDiameter_F,'color',colors('smoky black'),'LineWidth',0.5)
title('Awake Pupil Diameter')
ylabel('Z Units')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
axis square
ax10 = subplot(6,6,10);
plot(timeVector,data.Asleep.meanDiameter,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.Asleep.meanDiameter + data.Asleep.stdDiameter,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.Asleep.meanDiameter - data.Asleep.stdDiameter,'color',colors('smoky black'),'LineWidth',0.5)
title('Asleep Pupil Diameter')
ylabel('Z Units')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
axis square
ax11 = subplot(6,6,11);
plot(timeVector,data.Asleep.meanDiameter_T,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.Asleep.meanDiameter_T + data.Asleep.stdDiameter_T,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.Asleep.meanDiameter_T - data.Asleep.stdDiameter_T,'color',colors('smoky black'),'LineWidth',0.5)
title('Asleep Pupil Diameter')
ylabel('Z Units')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
axis square
ax12 = subplot(6,6,12);
plot(timeVector,data.Asleep.meanDiameter_F,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.Asleep.meanDiameter_F + data.Asleep.stdDiameter_F,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.Asleep.meanDiameter_F - data.Asleep.stdDiameter_F,'color',colors('smoky black'),'LineWidth',0.5)
title('Asleep Pupil Diameter')
ylabel('Z Units')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
axis square
%% HbT
summaryFigure = figure;
ax13 = subplot(6,6,13);
plot(timeVector,data.Awake.meanHbT,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.Awake.meanHbT + data.Awake.stdHbT,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.Awake.meanHbT - data.Awake.stdHbT,'color',colors('smoky black'),'LineWidth',0.5)
title('Awake HbT')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
axis square
ax14 = subplot(6,6,14);
plot(timeVector,data.Awake.meanHbT_T,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.Awake.meanHbT_T + data.Awake.stdHbT_T,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.Awake.meanHbT_T - data.Awake.stdHbT_T,'color',colors('smoky black'),'LineWidth',0.5)
title('Awake HbT')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
axis square
ax15 = subplot(6,6,15);
plot(timeVector,data.Awake.meanHbT_F,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.Awake.meanHbT_F + data.Awake.stdHbT_F,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.Awake.meanHbT_F - data.Awake.stdHbT_F,'color',colors('smoky black'),'LineWidth',0.5)
title('Awake HbT')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
axis square
ax16 = subplot(6,6,16);
plot(timeVector,data.Asleep.meanHbT,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.Asleep.meanHbT + data.Asleep.stdHbT,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.Asleep.meanHbT - data.Asleep.stdHbT,'color',colors('smoky black'),'LineWidth',0.5)
title('Asleep HbT')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
axis square
ax17 = subplot(6,6,17);
plot(timeVector,data.Asleep.meanHbT_T,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.Asleep.meanHbT_T + data.Asleep.stdHbT_T,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.Asleep.meanHbT_T - data.Asleep.stdHbT_T,'color',colors('smoky black'),'LineWidth',0.5)
title('Asleep HbT')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
axis square
ax18 = subplot(6,6,18);
plot(timeVector,data.Asleep.meanHbT_F,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.Asleep.meanHbT_F + data.Asleep.stdHbT_F,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.Asleep.meanHbT_F - data.Asleep.stdHbT_F,'color',colors('smoky black'),'LineWidth',0.5)
title('Asleep HbT')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
axis square
%% EMG
summaryFigure = figure;
ax19 = subplot(6,6,19);
plot(timeVector,data.Awake.meanEMG,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.Awake.meanEMG + data.Awake.stdEMG,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.Awake.meanEMG - data.Awake.stdEMG,'color',colors('smoky black'),'LineWidth',0.5)
title('Awake EMG')
ylabel('Power (a.u.)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
axis square
ax20 = subplot(6,6,20);
plot(timeVector,data.Awake.meanEMG_T,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.Awake.meanEMG_T + data.Awake.stdEMG_T,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.Awake.meanEMG_T - data.Awake.stdEMG_T,'color',colors('smoky black'),'LineWidth',0.5)
title('Awake EMG')
ylabel('Power (a.u.)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
axis square
ax21 = subplot(6,6,21);
plot(timeVector,data.Awake.meanEMG_F,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.Awake.meanEMG_F + data.Awake.stdEMG_F,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.Awake.meanEMG_F - data.Awake.stdEMG_F,'color',colors('smoky black'),'LineWidth',0.5)
title('Awake EMG')
ylabel('Power (a.u.)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
axis square
ax22 = subplot(6,6,22);
plot(timeVector,data.Asleep.meanEMG,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.Asleep.meanEMG + data.Asleep.stdEMG,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.Asleep.meanEMG - data.Asleep.stdEMG,'color',colors('smoky black'),'LineWidth',0.5)
title('Asleep EMG')
ylabel('Power (a.u.)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
axis square
ax23 = subplot(6,6,23);
plot(timeVector,data.Asleep.meanEMG_T,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.Asleep.meanEMG_T + data.Asleep.stdEMG_T,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.Asleep.meanEMG_T - data.Asleep.stdEMG_T,'color',colors('smoky black'),'LineWidth',0.5)
title('Asleep EMG')
ylabel('Power (a.u.)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
axis square
ax24 = subplot(6,6,24);
plot(timeVector,data.Asleep.meanEMG_F,'color',colors('smoky black'),'LineWidth',2);
hold on
plot(timeVector,data.Asleep.meanEMG_F + data.Asleep.stdEMG_F,'color',colors('smoky black'),'LineWidth',0.5)
plot(timeVector,data.Asleep.meanEMG_F - data.Asleep.stdEMG_F,'color',colors('smoky black'),'LineWidth',0.5)
title('Asleep EMG')
ylabel('Power (a.u.)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
axis square
%% CORT
summaryFigure = figure;
ax25 = subplot(6,6,25);
imagesc(T,F,data.Awake.meanCort)
title('Awake Cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
caxis([-0.2,0.2])
axis square
axis xy
ax26 = subplot(6,6,26);
imagesc(T,F,data.Awake.meanCort_T)
title('Awake Cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
caxis([-0.2,0.2])
axis square
axis xy
ax27 = subplot(6,6,27);
imagesc(T,F,data.Awake.meanCort_F)
title('Awake Cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
caxis([-0.2,0.2])
axis square
axis xy
ax28 = subplot(6,6,28);
imagesc(T,F,data.Asleep.meanCort)
title('Asleep Cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
caxis([-0.2,0.2])
axis square
axis xy
ax29 = subplot(6,6,29);
imagesc(T,F,data.Asleep.meanCort_T)
title('Asleep Cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
caxis([-0.2,0.2])
axis square
axis xy
ax30 = subplot(6,6,30);
imagesc(T,F,data.Asleep.meanCort_F)
title('Asleep Cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
caxis([-0.2,0.2])
axis square
axis xy
%% HIP
summaryFigure = figure;
ax31 = subplot(6,6,31);
imagesc(T,F,data.Awake.meanHip)
title('Awake Hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
caxis([-0.2,0.2])
axis square
axis xy
ax32 = subplot(6,6,32);
imagesc(T,F,data.Awake.meanHip_T)
title('Awake Hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
caxis([-0.2,0.2])
axis square
axis xy
ax33 = subplot(6,6,33);
imagesc(T,F,data.Awake.meanHip_F)
title('Awake Hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
caxis([-0.2,0.2])
axis square
axis xy
ax34 = subplot(6,6,34);
imagesc(T,F,data.Asleep.meanHip)
title('Asleep Hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
caxis([-0.2,0.2])
axis square
axis xy
ax35 = subplot(6,6,35);
imagesc(T,F,data.Asleep.meanHip_T)
title('Asleep Hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
caxis([-0.2,0.2])
axis square
axis xy
ax36 = subplot(6,6,36);
imagesc(T,F,data.Asleep.meanHip_F)
title('Asleep Hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
caxis([-0.2,0.2])
axis square
axis xy
%% link figure axis
% linkaxes([ax1,ax2,ax3])
% linkaxes([ax4,ax5,ax6])
% linkaxes([ax7,ax8,ax9])
% linkaxes([ax10,ax11,ax12])
% linkaxes([ax13,ax14,ax15])
% linkaxes([ax16,ax17,ax18])
% linkaxes([ax19,ax20,ax21])
% linkaxes([ax22,ax23,ax24])
% linkaxes([ax25,ax26,ax27])
% linkaxes([ax28,ax29,ax30])
% linkaxes([ax31,ax32,ax33])
% linkaxes([ax34,ax35,ax36])
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Blink Evoked Data' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'BlinkEvoked']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'BlinkEvoked'])
end

end
