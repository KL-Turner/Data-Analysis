clear; clc;
currentFolder = pwd;
fileparts = strsplit(currentFolder,filesep);
rootFolder = fullfile(fileparts{1:end});
animalIDs = {'T115','T116','T117','T118','T125','T126'};
modelType = 'Manual';
baselineType = 'manualSelection';
leadTime = 10;
lagTime = 10;
vesselSamplingRate = 5;
emgSamplingRate = 30;
[z,p,k] = butter(4,1/(vesselSamplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    Data.(animalID) = [];
    dataLocation = [rootFolder '/' animalID '/2P Data/'];
    cd(dataLocation)
    % find and load SleepData.mat strut
    sleepDataFileStruct = dir('*_SleepData.mat');
    sleepDataFile = {sleepDataFileStruct.name}';
    sleepDataFileID = char(sleepDataFile);
    load(sleepDataFileID)
    % load the baseline structure
    baselinesFileStruct = dir('*_RestingBaselines.mat');
    baselinesFile = {baselinesFileStruct.name}';
    baselinesFileID = char(baselinesFile);
    load(baselinesFileID)
    % NREM file IDs
    sleepFileIDs = unique(SleepData.Manual.NREM.FileIDs);
    for bb = 1:length(sleepFileIDs)
        sleepFileID = sleepFileIDs(bb,1);
        sleepFileID = sleepFileID{1,1};
        directory = pwd;      % Full path of the directory to be searched in
        filesAndFolders = dir(directory);     % Returns all the files and folders in the directory
        filesInDir = filesAndFolders(~([filesAndFolders.isdir]));  % Returns only the files in the directory
        fileNames = {filesInDir.name};
        stringToBeFound = sleepFileID;
        numOfFiles = length(filesInDir);
        fileFound = false;
        while fileFound == false
            for cc = 1:numOfFiles
                fileName = fileNames{1,cc};
                found = strfind(fileName,stringToBeFound);
                if isempty(found) == false
                    fileFound = true;
                    mergedDataFileID = fileName;
                    [~,~,fileDate,~,~,vesselID] = GetFileInfo2_2P(fileName);
                    strDay = ConvertDate_2P_eLife2020(fileDate);
                    load(mergedDataFileID,'-mat')
                    % vessel diameter
                    vesselDiameter = MergedData.data.vesselDiameter.data;
                    normVesselDiameter = (vesselDiameter - RestingBaselines.(baselineType).vesselDiameter.data.(vesselID).(strDay))./(RestingBaselines.(baselineType).vesselDiameter.data.(vesselID).(strDay));
                    filtVesselDiameter = filtfilt(sos,g,normVesselDiameter)*100;
                    % EMG
                    EMG = MergedData.data.EMG.data;
                    if isfield(Data.(animalID),vesselID) == false
                        Data.(animalID).(vesselID).timePoint = [];
                        Data.(animalID).(vesselID).vesselData = [];
                        Data.(animalID).(vesselID).emgData = [];
                        Data.(animalID).(vesselID).fileID ={};
                    end
                    break
                end
            end
        end
        [figHandle] = GenerateSingleFigures_2P_eLife2020(mergedDataFileID,baselineType,'n',RestingBaselines);
        useFile = input('Use file for micro-arousal transitions? (y/n): ','s');
        if strcmp(useFile,'y') == true
            keepSearching = 'y';
            while strcmp(keepSearching,'y') == true
                timePoint = input('Input time of micro-arousal (s): ');
                vesselVector = filtVesselDiameter(floor((timePoint - leadTime)*vesselSamplingRate):floor((timePoint + lagTime)*vesselSamplingRate));
                emgVector = EMG(floor((timePoint - leadTime)*emgSamplingRate):floor((timePoint + lagTime)*emgSamplingRate));
                Data.(animalID).(vesselID).timePoint = cat(1,Data.(animalID).(vesselID).timePoint,timePoint);
                Data.(animalID).(vesselID).vesselData = cat(1,Data.(animalID).(vesselID).vesselData,vesselVector);
                Data.(animalID).(vesselID).emgData = cat(1,Data.(animalID).(vesselID).emgData,emgVector);
                Data.(animalID).(vesselID).fileID = cat(1,Data.(animalID).(vesselID).fileID,mergedDataFileID);
                keepSearching = input('Keep searching? (y/n)','s');
            end
            close(figHandle)
        else
            close(figHandle)
        end
    end
end
%% Generate Figure
animalIDs = {'T115','T116','T117','T118'};
EMG_all = []; normEMG_all = []; catEMG_all = []; normCatEMG_all = [];
Diameter_all = []; normDiameter_all = []; catDiameter_all = []; normCatDiameter_all = [];
for aa = 1:length(animalIDs)
    vIDs = fieldnames(Data.(animalIDs{1,aa}));
    for bb = 1:length(vIDs)
        for cc = 1:size(Data.(animalIDs{1,aa}).(vIDs{bb,1}).vesselData,1)
            Data.(animalIDs{1,aa}).(vIDs{bb,1}).normVesselData(cc,:) = Data.(animalIDs{1,aa}).(vIDs{bb,1}).vesselData(cc,:) - mean(Data.(animalIDs{1,aa}).(vIDs{bb,1}).vesselData(cc,1:51));
            Data.(animalIDs{1,aa}).(vIDs{bb,1}).normEmgData(cc,:) = Data.(animalIDs{1,aa}).(vIDs{bb,1}).emgData(cc,:) - mean(Data.(animalIDs{1,aa}).(vIDs{bb,1}).emgData(cc,1:301));
        end
    end
end
for aa = 1:length(animalIDs)
    vIDs = fieldnames(Data.(animalIDs{1,aa}));
    for bb = 1:length(vIDs)
        EMG_all = cat(1,EMG_all,mean(Data.(animalIDs{1,aa}).(vIDs{bb,1}).emgData,1));
        Diameter_all = cat(1,Diameter_all,mean(Data.(animalIDs{1,aa}).(vIDs{bb,1}).vesselData,1));
        normEMG_all = cat(1,normEMG_all,mean(Data.(animalIDs{1,aa}).(vIDs{bb,1}).normEmgData,1));
        normDiameter_all = cat(1,normDiameter_all,mean(Data.(animalIDs{1,aa}).(vIDs{bb,1}).normVesselData,1));
        catEMG_all = cat(1,catEMG_all,Data.(animalIDs{1,aa}).(vIDs{bb,1}).emgData);
        catDiameter_all = cat(1,catDiameter_all,Data.(animalIDs{1,aa}).(vIDs{bb,1}).vesselData);
        normCatEMG_all = cat(1,normCatEMG_all,Data.(animalIDs{1,aa}).(vIDs{bb,1}).normEmgData);
        normCatDiameter_all = cat(1,normCatDiameter_all,Data.(animalIDs{1,aa}).(vIDs{bb,1}).normVesselData);
    end
end
meanEMG = mean(EMG_all,1);
stdEMG = std(EMG_all,0,1);
meanDiameter = mean(Diameter_all,1);
stdDiameter = std(Diameter_all,0,1);
meanNormEMG = mean(normEMG_all,1);
stdNormEMG = std(normEMG_all,0,1);
meanNormDiameter = mean(normDiameter_all,1);
stdNormDiameter = std(normDiameter_all,0,1);
meanCatEMG = mean(catEMG_all,1);
stdCatEMG = std(catEMG_all,0,1);
meanCatDiameter = mean(catDiameter_all,1);
stdCatDiameter = std(catDiameter_all,0,1);
meanNormCatEMG = mean(normCatEMG_all,1);
stdNormCatEMG = std(normCatEMG_all,0,1);
meanNormCatDiameter = mean(normCatDiameter_all,1);
stdNormCatDiameter = std(normCatDiameter_all,0,1);
figure;
ax1 = subplot(2,2,1);
yyaxis left
plot((1:length(meanDiameter))/5,meanDiameter,'r-','LineWidth',2)
hold on
plot((1:length(meanDiameter))/5,meanDiameter + stdDiameter,'r-','LineWidth',0.5)
plot((1:length(meanDiameter))/5,meanDiameter - stdDiameter,'r-','LineWidth',0.5)
title('vessel average, n = 10 arterioles, n = 4 mice')
xlabel('Time (sec)')
ylabel('\DeltaD/D (%)')
xticks([0,5,10,15,20])
xticklabels({'-10','-5','0','5','10'})
yyaxis right
plot((1:length(meanEMG))/30,meanEMG,'k-','LineWidth',2)
hold on
plot((1:length(meanEMG))/30,meanEMG + stdEMG,'k-','LineWidth',0.5)
plot((1:length(meanEMG))/30,meanEMG - stdEMG,'k-','LineWidth',0.5)
ylabel({'EMG power','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
axis tight
ax1.YAxis(1).Color = 'r';
ax1.YAxis(2).Color = 'k';
ax2 = subplot(2,2,2);
yyaxis left
plot((1:length(meanNormDiameter))/5,meanNormDiameter,'r-','LineWidth',2)
hold on
plot((1:length(meanNormDiameter))/5,meanNormDiameter + stdNormDiameter,'r-','LineWidth',0.5)
plot((1:length(meanNormDiameter))/5,meanNormDiameter - stdNormDiameter,'r-','LineWidth',0.5)
title('-10:0 mean subtracted, vessel average, n = 10 arterioles, n = 4 mice')
xlabel('Time (sec)')
ylabel('\DeltaD/D (%)')
xticks([0,5,10,15,20])
xticklabels({'-10','-5','0','5','10'})
yyaxis right
plot((1:length(meanNormEMG))/30,meanNormEMG,'k-','LineWidth',2)
hold on
plot((1:length(meanNormEMG))/30,meanNormEMG + stdNormEMG,'k-','LineWidth',0.5)
plot((1:length(meanNormEMG))/30,meanNormEMG - stdNormEMG,'k-','LineWidth',0.5)
ylabel({'EMG power','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
axis tight
ax2.YAxis(1).Color = 'r';
ax2.YAxis(2).Color = 'k';
ax3 = subplot(2,2,3);
yyaxis left
plot((1:length(meanCatDiameter))/5,meanCatDiameter,'r-','LineWidth',2)
hold on
plot((1:length(meanCatDiameter))/5,meanCatDiameter + stdCatDiameter,'r-','LineWidth',0.5)
plot((1:length(meanCatDiameter))/5,meanCatDiameter - stdCatDiameter,'r-','LineWidth',0.5)
title('event average, n = 79 events')
xlabel('Time (sec)')
ylabel('\DeltaD/D (%)')
xticks([0,5,10,15,20])
xticklabels({'-10','-5','0','5','10'})
yyaxis right
plot((1:length(meanCatEMG))/30,meanCatEMG,'k-','LineWidth',2)
hold on
plot((1:length(meanCatEMG))/30,meanCatEMG + stdCatEMG,'k-','LineWidth',0.5)
plot((1:length(meanCatEMG))/30,meanCatEMG - stdCatEMG,'k-','LineWidth',0.5)
ylabel({'EMG power','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
axis tight
ax3.YAxis(1).Color = 'r';
ax3.YAxis(2).Color = 'k';
ax4 = subplot(2,2,4);
yyaxis left
plot((1:length(meanNormCatDiameter))/5,meanNormCatDiameter,'r-','LineWidth',2)
hold on
plot((1:length(meanNormCatDiameter))/5,meanNormCatDiameter + stdNormCatDiameter,'r-','LineWidth',0.5)
plot((1:length(meanNormCatDiameter))/5,meanNormCatDiameter - stdNormCatDiameter,'r-','LineWidth',0.5)
title('-10:0 mean subtracted, event average, n = 79 events')
xlabel('Time (sec)')
ylabel('\DeltaD/D (%)')
xticks([0,5,10,15,20])
xticklabels({'-10','-5','0','5','10'})
yyaxis right
plot((1:length(meanNormCatEMG))/30,meanNormCatEMG,'k-','LineWidth',2)
hold on
plot((1:length(meanNormCatEMG))/30,meanNormCatEMG + stdNormCatEMG,'k-','LineWidth',0.5)
plot((1:length(meanNormCatEMG))/30,meanNormCatEMG - stdNormCatEMG,'k-','LineWidth',0.5)
ylabel({'EMG power','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
axis tight
ax4.YAxis(1).Color = 'r';
ax4.YAxis(2).Color = 'k';
