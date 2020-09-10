function [] = ElectromyographyWhiskingTransition()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Determine the transition of each animal from awake rest to whisking
%________________________________________________________________________________________________________________________

clear; clc;
%% function parameters
IOS_animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
params.minTime.Rest = 10;
currentFolder = pwd;
fileparts = strsplit(currentFolder,filesep);
rootFolder = fullfile(fileparts{1:end});
%% only run analysis for valid animal IDs
for aa = 1:length(IOS_animalIDs)
    disp(num2str(aa))
    clearvars -except rootFolder params IOS_animalIDs aa dataMeans_A dataMeans_B
    animalID = IOS_animalIDs{1,aa};
    dataLocation = [rootFolder '\' animalID '\Bilateral Imaging\'];
    cd(dataLocation)
    % find and load RestData.mat struct
    restDataFileStruct = dir('*_RestData.mat');
    restDataFile = {restDataFileStruct.name}';
    restDataFileID = char(restDataFile);
    load(restDataFileID)
    % find and load manual baseline event information
    manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
    manualBaselineFile = {manualBaselineFileStruct.name}';
    manualBaselineFileID = char(manualBaselineFile);
    load(manualBaselineFileID)
    % find and load RestingBaselines.mat strut
    baselineDataFileStruct = dir('*_RestingBaselines.mat');
    baselineDataFile = {baselineDataFileStruct.name}';
    baselineDataFileID = char(baselineDataFile);
    load(baselineDataFileID)
    % sampling rate
    samplingRate = RestData.CBV.adjLH.CBVCamSamplingRate;
    % criteria for the resting
    RestCriteria.Fieldname = {'durations'};
    RestCriteria.Comparison = {'gt'};
    RestCriteria.Value = {params.minTime.Rest};
    RestPuffCriteria.Fieldname = {'puffDistances'};
    RestPuffCriteria.Comparison = {'gt'};
    RestPuffCriteria.Value = {5};
    %% analyze [EMG] after periods of rest
    % pull data from RestData.mat structure
    [restLogical] = FilterEvents_IOS_Manuscript2020(RestData.CBV_HbT.adjLH,RestCriteria);
    [puffLogical] = FilterEvents_IOS_Manuscript2020(RestData.CBV_HbT.adjLH,RestPuffCriteria);
    combRestLogical = logical(restLogical.*puffLogical);
    restFileIDs = RestData.CBV_HbT.adjLH.fileIDs(combRestLogical,:);
    restEventTimes = RestData.CBV_HbT.adjLH.eventTimes(combRestLogical,:);
    restDurations = RestData.CBV_HbT.adjLH.durations(combRestLogical,:);
    LH_RestingData = RestData.CBV_HbT.adjLH.data(combRestLogical,:);
    % keep only the data that occurs within the manually-approved awake regions
    [~,finalRestFileIDs,finalRestDurations,finalRestEventTimes] = RemoveInvalidData_IOS_Manuscript2020(LH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    % pull data from each procData fileID
    cc = 1;
    for bb = 1:length(finalRestFileIDs)
        fileID_A = [animalID '_' finalRestFileIDs{bb,1} '_ProcData.mat'];
        fileID_B = [animalID '_' finalRestFileIDs{bb,1} '_RawData.mat'];
        strDay = ConvertDate_IOS_Manuscript2020(finalRestFileIDs{bb,1}(1:6));
        if finalRestEventTimes(bb,1) + finalRestDurations(bb,1) <= 895
            load(fileID_A,'-mat')
            load(fileID_B,'-mat')
            fpass = [300,3000];
            analogSamplingRate = 20000;
            analogExpectedLength = 900*analogSamplingRate;
            trimmedEMG = RawData.data.EMG(1:min(analogExpectedLength,length(RawData.data.EMG)));
            [z,p,k] = butter(3,fpass/(analogSamplingRate/2));
            [sos,g] = zp2sos(z,p,k);
            filtEMG = filtfilt(sos,g,trimmedEMG - mean(trimmedEMG));
            % kernelWidth = 0.5;
            % smoothingKernel = gausswin(kernelWidth*ProcData.notes.analogSamplingRate)/sum(gausswin(kernelWidth*ProcData.notes.analogSamplingRate));
            % EMGPwr = log10(conv(filtEMG.^2,smoothingKernel,'same'));
            EMGPwr = log10(filtEMG.^2);
            resampEMG = resample(EMGPwr,ProcData.notes.dsFs,ProcData.notes.analogSamplingRate);
            emg_A = resampEMG - RestingBaselines.manualSelection.EMG.emg.(strDay);
            emg_B = ProcData.data.EMG.emg - RestingBaselines.manualSelection.EMG.emg.(strDay);
            startIndex = round(((finalRestEventTimes(bb,1) + finalRestDurations(bb,1)) - 10)*samplingRate);
            endIndex = round(((finalRestEventTimes(bb,1) + finalRestDurations(bb,1)) + 5)*samplingRate);
            data_A(cc,:) = emg_A(startIndex:endIndex);
            data_B(cc,:) = emg_B(startIndex:endIndex); %#ok<*AGROW>
            cc = cc + 1;
        end
    end
    dataMeans_A(aa,:) = mean(data_A,1);
    dataMeans_B(aa,:) = mean(data_B,1);
    cd(rootFolder)
end
timeVec = -10:(1/samplingRate):5;
meanEMG_A = mean(dataMeans_A,1);
stdEMG_A = std(dataMeans_A,0,1);
meanEMG_B = mean(dataMeans_B,1);
stdEMG_B = std(dataMeans_B,0,1);
%% figure
figure
p1 = plot(timeVec,meanEMG_A - mean(meanEMG_A(1:10*ProcData.notes.dsFs)),'color',colors_Manuscript2020('rich black'),'LineWidth',2);
hold on
% plot(timeVec,meanEMG_A + stdEMG_A,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
% plot(timeVec,meanEMG_A - stdEMG_A,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
p2 = plot(timeVec,meanEMG_B - mean(meanEMG_B(1:10*ProcData.notes.dsFs)),'color',colors_Manuscript2020('sapphire'),'LineWidth',2);
% plot(timeVec,meanEMG_B + stdEMG_B,'color',colors_Manuscript2020('cyan'),'LineWidth',0.5)
% plot(timeVec,meanEMG_B - stdEMG_B,'color',colors_Manuscript2020('cyan'),'LineWidth',0.5)
xline(0,'color','r','LineWidth',2)
title('rest-whisk EMG transition, n = 14')
ylabel('EMG power (a.u.)')
xlabel('Peri-whisk time (s)')
legend([p1,p2],'Unfiltered','Filtered')
xlim([-10,5])
set(gca,'box','off')
%% figure
figure
plot(timeVec,meanEMG_A - mean(meanEMG_A(1:10*ProcData.notes.dsFs)),'color',colors_Manuscript2020('rich black'),'LineWidth',2);
hold on
for dd = 1:size(dataMeans_A,1)
    plot(timeVec,dataMeans_A(dd,:) - mean(dataMeans_A(dd,1:10*ProcData.notes.dsFs)),'LineWidth',0.5)
end
xline(0,'color','r','LineWidth',2)
title('rest-whisk EMG transition, Unfiltered n = 14')
ylabel('EMG power (a.u.)')
xlabel('Peri-whisk time (s)')
xlim([-10,5])
set(gca,'box','off')

end
