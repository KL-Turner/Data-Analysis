function [] = Fig2_TBD(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

% behavior colors
colorRest = [(0/256),(166/256),(81/256)];
colorWhisk = [(31/256),(120/256),(179/256)];
colorStim = [(255/256),(28/256),(206/256)];
colorNREM = [(191/256),(0/256),(255/256)];
colorREM = [(254/256),(139/256),(0/256)];
colorAlert = [(255/256),(191/256),(0/256)];
colorAsleep = [(0/256),(128/256),(255/256)];
colorAll = [(183/256),(115/256),(51/256)];
%% set-up and process data
resultsStruct = 'Results_BehavData.mat';
load(resultsStruct);
animalIDs = fieldnames(Results_BehavData);
behavFields = {'Rest','Whisk','Stim','NREM','REM'};
% mean HbT comparison between behaviors
% pre-allocate the date for each day
% pre-allocate
for cc = 1:length(behavFields)
    behavField = behavFields{1,cc};
    data.(behavField).indMeanDiameter = [];
    data.(behavField).indDiameter = [];
    data.(behavField).indMeanzDiameter = [];
    data.(behavField).indzDiameter = [];
end
% concatenate
for cc = 1:length(animalIDs)
    animalID = animalIDs{cc,1};
    for dd = 1:length(behavFields)
        behavField = behavFields{1,dd};
        data.(behavField).indMeanDiameter = cat(1,data.(behavField).indMeanDiameter,mean(Results_BehavData.(animalID).(behavField).mmDiameter.eventMeans,'omitnan'));
        data.(behavField).indMeanzDiameter = cat(1,data.(behavField).indMeanzDiameter,mean(Results_BehavData.(animalID).(behavField).zDiameter.eventMeans,'omitnan'));
        indDiameter = []; indzDiameter = [];
        for ee = 1:length(Results_BehavData.(animalID).(behavField).mmArea.indData)
            indDiameter = cat(2,indDiameter,Results_BehavData.(animalID).(behavField).mmDiameter.indData{ee,1});
            indzDiameter = cat(2,indzDiameter,Results_BehavData.(animalID).(behavField).zDiameter.indData{ee,1});
        end
        data.(behavField).indDiameter = cat(2,data.(behavField).indDiameter,indDiameter);
        data.(behavField).indzDiameter = cat(2,data.(behavField).indzDiameter,indzDiameter);
    end
end
% mean/std
for ee = 1:length(behavFields)
    behavField = behavFields{1,ee};
    % diameter
    data.(behavField).meanDiameter = mean(data.(behavField).indMeanDiameter,1,'omitnan');
    data.(behavField).stdDiameter = std(data.(behavField).indMeanDiameter,0,1,'omitnan');
    realIndex = ~isnan(data.(behavField).indDiameter);
    data.(behavField).indDiameter = data.(behavField).indDiameter(realIndex);
    % z diameter
    data.(behavField).meanzDiameter = mean(data.(behavField).indMeanzDiameter,1,'omitnan');
    data.(behavField).stdzDiameter = std(data.(behavField).indMeanzDiameter,0,1,'omitnan');
    realIndex = ~isnan(data.(behavField).indzDiameter);
    data.(behavField).indzDiameter = data.(behavField).indzDiameter(realIndex);
end
%%
resultsStruct = 'Results_Evoked';
load(resultsStruct);
whiskDataTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
solenoidNames = {'LPadSol','RPadSol','AudSol'};
dataTypes = {'mmArea','mmDiameter','zArea','zDiameter'};
timeVector = (0:12*30)/30 - 2;
animalIDs = fieldnames(Results_Evoked);
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    data.controlSolenoid.(dataType) = [];
    data.stimSolenoid.(dataType) = [];
    data.briefWhisk.(dataType) = [];
    data.interWhisk.(dataType) = [];
    data.extendWhisk.(dataType) = [];
end
for cc = 1:length(animalIDs)
    animalID = animalIDs{cc,1};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        % whisking
        data.briefWhisk.(dataType) = cat(1,data.briefWhisk.(dataType),Results_Evoked.(animalID).Whisk.(dataType).ShortWhisks.mean);
        data.interWhisk.(dataType) = cat(1,data.interWhisk.(dataType),Results_Evoked.(animalID).Whisk.(dataType).IntermediateWhisks.mean);
        data.extendWhisk.(dataType) = cat(1,data.extendWhisk.(dataType),Results_Evoked.(animalID).Whisk.(dataType).LongWhisks.mean);
        % solenoids
        data.stimSolenoid.(dataType) = cat(1,data.stimSolenoid.(dataType),Results_Evoked.(animalID).Stim.(dataType).LPadSol.mean,Results_Evoked.(animalID).Stim.(dataType).RPadSol.mean);
        data.controlSolenoid.(dataType) = cat(1,data.controlSolenoid.(dataType),Results_Evoked.(animalID).Stim.(dataType).AudSol.mean);
    end
end
%
for dd = 1:length(dataTypes)
    dataType = dataTypes{1,dd};
    procData.briefWhisk.(dataType).mean = mean(data.briefWhisk.(dataType),1);
    procData.briefWhisk.(dataType).std = std(data.briefWhisk.(dataType),0,1)./sqrt(size(data.briefWhisk.(dataType),1));
    procData.interWhisk.(dataType).mean = mean(data.interWhisk.(dataType),1);
    procData.interWhisk.(dataType).std = std(data.interWhisk.(dataType),0,1)./sqrt(size(data.interWhisk.(dataType),1));
    procData.extendWhisk.(dataType).mean = mean(data.extendWhisk.(dataType),1);
    procData.extendWhisk.(dataType).std = std(data.extendWhisk.(dataType),0,1)./sqrt(size(data.extendWhisk.(dataType),1));
    procData.stimSolenoid.(dataType).mean = mean(data.stimSolenoid.(dataType),1);
    procData.stimSolenoid.(dataType).std = std(data.stimSolenoid.(dataType),0,1)./sqrt(size(data.stimSolenoid.(dataType),1));
    procData.controlSolenoid.(dataType).mean = mean(data.controlSolenoid.(dataType),1);
    procData.controlSolenoid.(dataType).std = std(data.controlSolenoid.(dataType),0,1)./sqrt(size(data.controlSolenoid.(dataType),1));
end
%% set-up and process data
resultsStruct = 'Results_CrossCorrelation';
load(resultsStruct);
animalIDs = fieldnames(Results_CrossCorrelation);
behavFields = {'Rest','NREM','REM','Alert','Asleep','All'};
dataTypes = {'mmArea','mmDiameter','zArea','zDiameter'};
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            % pre-allocate necessary variable fields
            data.(behavField).(dataType).dummCheck = 1;
            if isfield(data.(behavField).(dataType),'LH_xcVals') == false
                data.(behavField).(dataType).LH_xcVals = [];
                data.(behavField).(dataType).RH_xcVals = [];
                data.(behavField).(dataType).lags = [];
            end
            if isfield(Results_CrossCorrelation.(animalID),behavField) == true
                data.(behavField).(dataType).LH_xcVals = cat(1,data.(behavField).(dataType).LH_xcVals,Results_CrossCorrelation.(animalID).(behavField).LH_HbT.(dataType).xcVals);
                data.(behavField).(dataType).RH_xcVals = cat(1,data.(behavField).(dataType).RH_xcVals,Results_CrossCorrelation.(animalID).(behavField).RH_HbT.(dataType).xcVals);
                data.(behavField).(dataType).lags = cat(1,data.(behavField).(dataType).lags,Results_CrossCorrelation.(animalID).(behavField).LH_HbT.(dataType).lags,Results_CrossCorrelation.(animalID).(behavField).RH_HbT.(dataType).lags);
                data.(behavField).animalID{aa,1} = animalID;
                data.(behavField).behavior{aa,1} = behavField;
                data.(behavField).LH{aa,1} = 'LH';
                data.(behavField).RH{aa,1} = 'RH';
            end
        end
    end
end
% take the averages of each field through the proper dimension
for dd = 1:length(behavFields)
    behavField = behavFields{1,dd};
    for ff = 1:length(dataTypes)
        dataType = dataTypes{1,ff};
        data.(behavField).(dataType).xcVals = cat(1,data.(behavField).(dataType).LH_xcVals,data.(behavField).(dataType).RH_xcVals);
        data.(behavField).(dataType).meanXcVals = mean(data.(behavField).(dataType).xcVals,1);
        data.(behavField).(dataType).stdXcVals = std(data.(behavField).(dataType).xcVals,0,1);
        data.(behavField).(dataType).meanLags = mean(data.(behavField).(dataType).lags,1);
    end
end
% % find max/time to peak
% for gg = 1:length(behavFields)
%     behavField = behavFields{1,gg};
%     for hh = 1:size(data.(behavField).(dataType).LH_xcVals,1)
%         % hbt
%         gammaArray = data.(behavField).(dataType).LH_xcVals(hh,:);
%         [gammaMax,gammaIndex] = max(gammaArray);
%         data.(behavField).gammaPeak(hh,1) = gammaMax;
%         data.(behavField).gammaTTP(hh,1) = data.(behavField).meanLags(gammaIndex)/30;
%     end
% end
% % mean/std
% for ii = 1:length(behavFields)
%     behavField = behavFields{1,ii};
%     data.(behavField).meanMuaPeak = mean(data.(behavField).muaPeak,1);
%     data.(behavField).stdMuaPeak = std(data.(behavField).muaPeak,0,1);
%     data.(behavField).meanMuaTTP = mean(data.(behavField).muaTTP,1);
%     data.(behavField).stdMuaTTP = std(data.(behavField).muaTTP,0,1);
% end
%% variables for loops
resultsStruct = 'Results_Coherence';
load(resultsStruct);
animalIDs = fieldnames(Results_Coherence);
behavFields = {'Rest','NREM','REM','Awake','Asleep','All'};
dataTypes = {'mmArea','mmDiameter','zArea','zDiameter'};
%% pre-allocate data structure
for aa = 1:length(behavFields)
    behavField = behavFields{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        data.(behavField).(dataType).HbTC = [];
        data.(behavField).(dataType).HbTf = [];
        
        data.(behavField).(dataType).gammaC = [];
        data.(behavField).(dataType).gammaf = [];
    end
end
% power spectra during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(Results_Coherence.(animalID).(behavField).(dataType).LH_HbT.C) == false
                data.(behavField).(dataType).HbTC = cat(2,data.(behavField).(dataType).HbTC,Results_Coherence.(animalID).(behavField).(dataType).LH_HbT.C,Results_Coherence.(animalID).(behavField).(dataType).RH_HbT.C);
                data.(behavField).(dataType).HbTf = cat(1,data.(behavField).(dataType).HbTf,Results_Coherence.(animalID).(behavField).(dataType).LH_HbT.f,Results_Coherence.(animalID).(behavField).(dataType).RH_HbT.f);
                
                data.(behavField).(dataType).gammaC = cat(2,data.(behavField).(dataType).gammaC,Results_Coherence.(animalID).(behavField).(dataType).LH_gammaBandPower.C,Results_Coherence.(animalID).(behavField).(dataType).RH_gammaBandPower.C);
                data.(behavField).(dataType).gammaf = cat(1,data.(behavField).(dataType).gammaf,Results_Coherence.(animalID).(behavField).(dataType).LH_gammaBandPower.f,Results_Coherence.(animalID).(behavField).(dataType).RH_gammaBandPower.f);
            end
        end
    end
end
% take mean/StD of S/f
for aa = 1:length(behavFields)
    behavField = behavFields{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        data.(behavField).(dataType).meanHbTC = mean(data.(behavField).(dataType).HbTC,2);
        data.(behavField).(dataType).stdHbTC = std(data.(behavField).(dataType).HbTC,0,2);
        data.(behavField).(dataType).meanHbTf = mean(data.(behavField).(dataType).HbTf,1);
        
        data.(behavField).(dataType).meangammaC = mean(data.(behavField).(dataType).gammaC,2);
        data.(behavField).(dataType).stdgammaC = std(data.(behavField).(dataType).gammaC,0,2);
        data.(behavField).(dataType).meangammaf = mean(data.(behavField).(dataType).gammaf,1);
    end
end
%% figures
summaryFigure = figure;
sgtitle('Figure 2')
%%
subplot(3,3,1);
%
p1 = plot(timeVector,procData.interWhisk.zDiameter.mean,'color',colors('vegas gold'),'LineWidth',2);
hold on
plot(timeVector,procData.interWhisk.zDiameter.mean + procData.interWhisk.zDiameter.std,'color',colors('vegas gold'),'LineWidth',0.5)
plot(timeVector,procData.interWhisk.zDiameter.mean - procData.interWhisk.zDiameter.std,'color',colors('vegas gold'),'LineWidth',0.5)
%
p2 = plot(timeVector,procData.stimSolenoid.zDiameter.mean,'color',colors('dark candy apple red'),'LineWidth',2);
plot(timeVector,procData.stimSolenoid.zDiameter.mean + procData.stimSolenoid.zDiameter.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
plot(timeVector,procData.stimSolenoid.zDiameter.mean - procData.stimSolenoid.zDiameter.std,'color',colors('dark candy apple red'),'LineWidth',0.5)
%
p3 = plot(timeVector,procData.controlSolenoid.zDiameter.mean,'color',colors('deep carrot orange'),'LineWidth',2);
hold on
plot(timeVector,procData.controlSolenoid.zDiameter.mean + procData.controlSolenoid.zDiameter.std,'color',colors('deep carrot orange'),'LineWidth',0.5)
plot(timeVector,procData.controlSolenoid.zDiameter.mean - procData.controlSolenoid.zDiameter.std,'color',colors('deep carrot orange'),'LineWidth',0.5)
ylabel('\DeltaZ Units')
xlabel('Time (s)')
legend([p1,p2,p3],'Whisk','Stim','Aud')
set(gca,'box','off')
xlim([-2,10])
axis square
%% mm pupil diameter scatter
ax2 = subplot(3,3,2);
scatter(ones(1,length(data.Rest.indMeanDiameter))*1,data.Rest.indMeanDiameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Rest.meanDiameter,data.Rest.stdDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.Whisk.indMeanDiameter))*2,data.Whisk.indMeanDiameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Whisk.meanDiameter,data.Whisk.stdDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.Stim.indMeanDiameter))*3,data.Stim.indMeanDiameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.Stim.meanDiameter,data.Stim.stdDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.NREM.indMeanDiameter))*4,data.NREM.indMeanDiameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.NREM.meanDiameter,data.NREM.stdDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.REM.indMeanDiameter))*5,data.REM.indMeanDiameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.REM.meanDiameter,data.REM.stdDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
title('mm Diameter')
ylabel('Diameter (mm)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields) + 1])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% mm pupil diameter scatter
ax4 = subplot(3,3,3);
scatter(ones(1,length(data.Rest.indMeanzDiameter))*1,data.Rest.indMeanzDiameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Rest.meanzDiameter,data.Rest.stdzDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.Whisk.indMeanzDiameter))*2,data.Whisk.indMeanzDiameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Whisk.meanzDiameter,data.Whisk.stdzDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.Stim.indMeanzDiameter))*3,data.Stim.indMeanzDiameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.Stim.meanzDiameter,data.Stim.stdzDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.NREM.indMeanzDiameter))*4,data.NREM.indMeanzDiameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.NREM.meanzDiameter,data.NREM.stdzDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.REM.indMeanzDiameter))*5,data.REM.indMeanzDiameter,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.REM.meanzDiameter,data.REM.stdzDiameter,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
title('Z-scored diameter')
ylabel('Z units')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields) + 1])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%%
subplot(3,4,5);
L1 = semilogx(data.Rest.zDiameter.meanHbTf,data.Rest.zDiameter.meanHbTC,'color',colorRest,'LineWidth',2);
hold on
rectangle('Position',[0.005,0.1,0.1 - 0.005,1],'FaceColor','w','EdgeColor','w')
L2 = semilogx(data.NREM.zDiameter.meanHbTf,data.NREM.zDiameter.meanHbTC,'color',colorNREM,'LineWidth',2);
rectangle('Position',[0.005,0.1,1/30 - 0.005,1],'FaceColor','w','EdgeColor','w')
L3 = semilogx(data.REM.zDiameter.meanHbTf,data.REM.zDiameter.meanHbTC,'color',colorREM,'LineWidth',2);
rectangle('Position',[0.005,0.1,1/60 - 0.005,1],'FaceColor','w','EdgeColor','w')
L4 = semilogx(data.Awake.zDiameter.meanHbTf,data.Awake.zDiameter.meanHbTC,'color',colorAlert,'LineWidth',2);
L5 = semilogx(data.Asleep.zDiameter.meanHbTf,data.Asleep.zDiameter.meanHbTC,'color',colorAsleep,'LineWidth',2);
L6 = semilogx(data.All.zDiameter.meanHbTf,data.All.zDiameter.meanHbTC,'color',colorAll,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
title('HbT and zDiameter')
ylabel('Coherence')
xlabel('Freq (Hz)')
legend([L1,L2,L3,L4,L5,L6],'Rest','NREM','REM','Alert','Asleep','All','Location','NorthEast')
% axis square
xlim([0.003,1])
ylim([0,1])
set(gca,'box','off')
%% HbT:Pupil Stats
% 0.003:0.01
% 0.3:0.5
%%
subplot(3,4,7);
L1 = semilogx(data.Rest.zDiameter.meangammaf,data.Rest.zDiameter.meangammaC,'color',colorRest,'LineWidth',2);
hold on
rectangle('Position',[0.005,0.1,0.1 - 0.005,1],'FaceColor','w','EdgeColor','w')
L2 = semilogx(data.NREM.zDiameter.meangammaf,data.NREM.zDiameter.meangammaC,'color',colorNREM,'LineWidth',2);
rectangle('Position',[0.005,0.1,1/30 - 0.005,1],'FaceColor','w','EdgeColor','w')
L3 = semilogx(data.REM.zDiameter.meangammaf,data.REM.zDiameter.meangammaC,'color',colorREM,'LineWidth',2);
rectangle('Position',[0.005,0.1,1/60 - 0.005,1],'FaceColor','w','EdgeColor','w')
L4 = semilogx(data.Awake.zDiameter.meangammaf,data.Awake.zDiameter.meangammaC,'color',colorAlert,'LineWidth',2);
L5 = semilogx(data.Asleep.zDiameter.meangammaf,data.Asleep.zDiameter.meangammaC,'color',colorAsleep,'LineWidth',2);
L6 = semilogx(data.All.zDiameter.meangammaf,data.All.zDiameter.meangammaC,'color',colorAll,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
title('gamma-band and zDiameter')
ylabel('Coherence')
xlabel('Freq (Hz)')
% axis square
xlim([0.003,1])
ylim([0,1])
set(gca,'box','off')
%% Gamma:Pupil Stats
% 0.003:0.01
% 0.3:0.5
%%
subplot(3,6,13);
freq = 30;
lagSec = 30;
plot(data.Rest.zDiameter.meanLags,data.Rest.zDiameter.meanXcVals,'color',colorRest);
hold on
plot(data.NREM.zDiameter.meanLags,data.NREM.zDiameter.meanXcVals,'color',colorNREM);
plot(data.REM.zDiameter.meanLags,data.REM.zDiameter.meanXcVals,'color',colorREM);
plot(data.Alert.zDiameter.meanLags,data.Alert.zDiameter.meanXcVals,'color',colorAlert);
plot(data.Asleep.zDiameter.meanLags,data.Asleep.zDiameter.meanXcVals,'color',colorAsleep);
plot(data.All.zDiameter.meanLags,data.All.zDiameter.meanXcVals,'color',colorAll);
title({'Blank-SAP treated RH REM','MUA-[HbT] XCorr'})
xticks([-lagSec*freq,-lagSec*freq/2,0,lagSec*freq/2,lagSec*freq])
xticklabels({'-30','-15','0','15','30'})
xlim([-lagSec*freq,lagSec*freq])
xlabel('Lags (s)')
ylabel('Correlation')
title('zDiameter')
axis square
set(gca,'box','off')
%% HbT XCorr stats
% peak (awake,NREM,REM)
% time-to-peak (awake, NREM, REM)
%% Gamma XCorr
%% Gamma XCorr
% peak stats (Alert, Asleep, All)
% time-to-peak (alert, asleep, all)
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig2_TBD']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Fig2_TBD'])
end

end
