function [] = CrossCorrelation_Pupil(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose:
%________________________________________________________________________________________________________________________

%% set-up and process data
resultsStruct = 'Results_CrossCorrelation';
load(resultsStruct);
animalIDs = fieldnames(Results_CrossCorrelation);
behavFields = {'Rest','NREM','REM','Alert','Asleep','All'};
dataTypes = {'mmArea','mmDiameter','zArea','zDiameter'};
colorRest = [(0/256),(166/256),(81/256)];
colorNREM = [(191/256),(0/256),(255/256)];
colorREM = [(254/256),(139/256),(0/256)];
colorAlert = [(255/256),(191/256),(0/256)];
colorAsleep = [(0/256),(128/256),(255/256)];
colorAll = [(183/256),(115/256),(51/256)];
%% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            % pre-allocate necessary variable fields
            data.(behavField).(dataType).dummCheck = 1;
            if isfield(data.(behavField).(dataType),'xcVals') == false
                data.(behavField).(dataType).xcVals = [];
                data.(behavField).(dataType).lags = [];
            end
            if isfield(Results_CrossCorrelation.(animalID),behavField) == true
                data.(behavField).(dataType).xcVals = cat(1,data.(behavField).(dataType).xcVals,Results_CrossCorrelation.(animalID).(behavField).LH_HbT.(dataType).xcVals,Results_CrossCorrelation.(animalID).(behavField).RH_HbT.(dataType).xcVals);
                data.(behavField).(dataType).lags = cat(1,data.(behavField).(dataType).lags,Results_CrossCorrelation.(animalID).(behavField).LH_HbT.(dataType).lags,Results_CrossCorrelation.(animalID).(behavField).RH_HbT.(dataType).lags);
            end
        end
    end
end
%% take the averages of each field through the proper dimension
for dd = 1:length(behavFields)
    behavField = behavFields{1,dd};
    for ff = 1:length(dataTypes)
        dataType = dataTypes{1,ff};
        data.(behavField).(dataType).meanXcVals = mean(data.(behavField).(dataType).xcVals,1);
        data.(behavField).(dataType).stdXcVals = std(data.(behavField).(dataType).xcVals,0,1);
        data.(behavField).(dataType).meanLags = mean(data.(behavField).(dataType).lags,1);
    end
end
%% Awake Rest XCorr
freq = 30;
lagSec = 30;
summaryFigure = figure;
sgtitle('Awake Rest LFP-HbT Cross-correlation')
%% RH MUA-HbT XCorr [Blank-SAP] REM
ax1 = subplot(2,2,1);
p1 = plot(data.Rest.mmArea.meanLags,data.Rest.mmArea.meanXcVals,'color',colorRest);
hold on
p2 = plot(data.NREM.mmArea.meanLags,data.NREM.mmArea.meanXcVals,'color',colorNREM);
p3 = plot(data.REM.mmArea.meanLags,data.REM.mmArea.meanXcVals,'color',colorREM);
p4 = plot(data.Alert.mmArea.meanLags,data.Alert.mmArea.meanXcVals,'color',colorAlert);
p5 = plot(data.Asleep.mmArea.meanLags,data.Asleep.mmArea.meanXcVals,'color',colorAsleep);
p6 = plot(data.All.mmArea.meanLags,data.All.mmArea.meanXcVals,'color',colorAll);
title({'Blank-SAP treated RH REM','MUA-[HbT] XCorr'})
xticks([-lagSec*freq,-lagSec*freq/2,0,lagSec*freq/2,lagSec*freq])
xticklabels({'-30','-15','0','15','30'})
xlim([-lagSec*freq,lagSec*freq])
xlabel('Lags (s)')
ylabel('Correlation')
title('mmArea')
axis square
set(gca,'box','off')
legend([p1,p2,p3,p4,p5,p6],'Rest','NREM','REM','Alert','Asleep','All')
ax2 = subplot(2,2,2);
plot(data.Rest.mmDiameter.meanLags,data.Rest.mmDiameter.meanXcVals,'color',colorRest);
hold on
plot(data.NREM.mmDiameter.meanLags,data.NREM.mmDiameter.meanXcVals,'color',colorNREM);
plot(data.REM.mmDiameter.meanLags,data.REM.mmDiameter.meanXcVals,'color',colorREM);
plot(data.Alert.mmDiameter.meanLags,data.Alert.mmDiameter.meanXcVals,'color',colorAlert);
plot(data.Asleep.mmDiameter.meanLags,data.Asleep.mmDiameter.meanXcVals,'color',colorAsleep);
plot(data.All.mmDiameter.meanLags,data.All.mmDiameter.meanXcVals,'color',colorAll);
title({'Blank-SAP treated RH REM','MUA-[HbT] XCorr'})
xticks([-lagSec*freq,-lagSec*freq/2,0,lagSec*freq/2,lagSec*freq])
xticklabels({'-30','-15','0','15','30'})
xlim([-lagSec*freq,lagSec*freq])
xlabel('Lags (s)')
ylabel('Correlation')
title('mmDiameter')
axis square
set(gca,'box','off')
ax3 = subplot(2,2,3);
plot(data.Rest.zArea.meanLags,data.Rest.zArea.meanXcVals,'color',colorRest);
hold on
plot(data.NREM.zArea.meanLags,data.NREM.zArea.meanXcVals,'color',colorNREM);
plot(data.REM.zArea.meanLags,data.REM.zArea.meanXcVals,'color',colorREM);
plot(data.Alert.zArea.meanLags,data.Alert.zArea.meanXcVals,'color',colorAlert);
plot(data.Asleep.zArea.meanLags,data.Asleep.zArea.meanXcVals,'color',colorAsleep);
plot(data.All.zArea.meanLags,data.All.zArea.meanXcVals,'color',colorAll);
title({'Blank-SAP treated RH REM','MUA-[HbT] XCorr'})
xticks([-lagSec*freq,-lagSec*freq/2,0,lagSec*freq/2,lagSec*freq])
xticklabels({'-30','-15','0','15','30'})
xlim([-lagSec*freq,lagSec*freq])
xlabel('Lags (s)')
ylabel('Correlation')
title('zArea')
axis square
set(gca,'box','off')
ax4 = subplot(2,2,4);
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
% figure characteristics
linkaxes([ax1,ax2,ax3,ax4],'xy')
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Pupil-HbT Cross Correlation' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Pupil_CrossCorrelation']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Pupil_CrossCorrelation'])
end

end
