function [AnalysisResults] = AnalyzeVesselTransitionalAverages_eLife2020(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the transitions between different arousal-states (2PLSM)
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs = {'T115','T116','T117','T118','T125','T126'};
modelType = 'Manual';
%% only run analysis for valid animal IDs
if any(strcmp(animalIDs,animalID))
    dataLocation = [rootFolder '/' animalID '/2P Data/'];
    cd(dataLocation)
    % character list of all MergedData files
    mergedDirectory = dir('*_MergedData.mat');
    mergedDataFiles = {mergedDirectory.name}';
    mergedDataFileIDs = char(mergedDataFiles);
    % find and load EventData.mat struct
    eventDataFileStruct = dir('*_EventData.mat');
    eventDataFile = {eventDataFileStruct.name}';
    eventDataFileID = char(eventDataFile);
    load(eventDataFileID)
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
    % find and load SleepData.mat strut
    sleepDataFileStruct = dir('*_SleepData.mat');
    sleepDataFile = {sleepDataFileStruct.name}';
    sleepDataFileID = char(sleepDataFile);
    load(sleepDataFileID)
    % lowpass filter
    samplingRate = RestData.vesselDiameter.data.samplingRate;
    [z,p,k] = butter(4,1/(samplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    %% NREM to REM transition
    if isfield(SleepData.(modelType),'REM') == true
        nremTransition = [];
        % pull data from SleepData.mat structure
        remVesselIDs = SleepData.(modelType).REM.VesselIDs;
        remBinTimes = SleepData.(modelType).REM.BinTimes;
        remFileIDs = SleepData.(modelType).REM.FileIDs;
        timeVector = -30:(1/samplingRate):70;
        for aa = 1:length(remBinTimes)
            remVID = remVesselIDs{aa,1};
            binTimes = remBinTimes{aa,1};
            remFileID = remFileIDs{aa,1};
            if binTimes(1) >= 30 && binTimes(end) <= 870 && strcmp(remVID(1),'V') == false
                startTime = (binTimes(1) - 30)*samplingRate;
                endTime = (binTimes(1) + 70)*samplingRate;
                for bb = 1:size(mergedDataFileIDs,1)
                    mergedDataFileID = mergedDataFileIDs(bb,:);
                    [~,~,fileDate,fileID,~,~] = GetFileInfo2_2P_eLife2020(mergedDataFileID);
                    strDay = ConvertDate_2P_eLife2020(fileDate);
                    if strcmp(remFileID,fileID) == true
                        load(mergedDataFileID,'-mat')
                        vesselData = MergedData.data.vesselDiameter.data(startTime:endTime);
                        normVesselData = (vesselData - RestingBaselines.manualSelection.vesselDiameter.data.(remVID).(strDay))./RestingBaselines.manualSelection.vesselDiameter.data.(remVID).(strDay);
                        filtVesselData = filtfilt(sos,g,normVesselData)*100;
                    end
                end
                if isfield(nremTransition,remVID) == false
                    nremTransition.(remVID) = filtVesselData;
                else
                    nremTransition.(remVID) = vertcat(nremTransition.(remVID),filtVesselData);
                end
            end
        end
        remTransitionVesselIDs = fieldnames(nremTransition);
        for cc = 1:length(remTransitionVesselIDs)
            vID = remTransitionVesselIDs{cc,1};
            % save results
            AnalysisResults.(animalID).Transitions.NREMtoREM.(vID).mean = mean(nremTransition.(vID),1);
            AnalysisResults.(animalID).Transitions.NREMtoREM.(vID).StD = std(nremTransition.(vID),0,1);
            AnalysisResults.(animalID).Transitions.NREMtoREM.(vID).timeVector = timeVector;
        end
        %% REM to awake transition
        remTransition = [];
        % pull data from SleepData.mat structure
        remVesselIDs = SleepData.(modelType).REM.VesselIDs;
        remBinTimes = SleepData.(modelType).REM.BinTimes;
        remFileIDs = SleepData.(modelType).REM.FileIDs;
        timeVector = -30:(1/samplingRate):30;
        for aa = 1:length(remBinTimes)
            remVID = remVesselIDs{aa,1};
            binTimes = remBinTimes{aa,1};
            remFileID = remFileIDs{aa,1};
            if binTimes(end) <= 870 && strcmp(remVID(1),'V') == false
                startTime = binTimes(end - 6)*samplingRate;
                endTime = (binTimes(end) + 30)*samplingRate;
                for bb = 1:size(mergedDataFileIDs,1)
                    mergedDataFileID = mergedDataFileIDs(bb,:);
                    [~,~,fileDate,fileID,~,~] = GetFileInfo2_2P_eLife2020(mergedDataFileID);
                    strDay = ConvertDate_2P_eLife2020(fileDate);
                    if strcmp(remFileID,fileID) == true
                        load(mergedDataFileID,'-mat')
                        vesselData = MergedData.data.vesselDiameter.data(startTime:endTime);
                        normVesselData = (vesselData - RestingBaselines.manualSelection.vesselDiameter.data.(remVID).(strDay))./RestingBaselines.manualSelection.vesselDiameter.data.(remVID).(strDay);
                        filtVesselData = filtfilt(sos,g,normVesselData)*100;
                    end
                end
                if isfield(remTransition,remVID) == false
                    remTransition.(remVID) = filtVesselData;
                else
                    remTransition.(remVID) = vertcat(remTransition.(remVID),filtVesselData);
                end
            end
        end
        remTransitionVesselIDs = fieldnames(remTransition);
        for cc = 1:length(remTransitionVesselIDs)
            vID = remTransitionVesselIDs{cc,1};
            % save results
            AnalysisResults.(animalID).Transitions.REMtoAwake.(vID).mean = mean(remTransition.(vID),1);
            AnalysisResults.(animalID).Transitions.REMtoAwake.(vID).StD = std(remTransition.(vID),0,1);
            AnalysisResults.(animalID).Transitions.REMtoAwake.(vID).timeVector = timeVector;
        end
    else
        % save results
        AnalysisResults.(animalID).Transitions.NREMtoREM = [];
        AnalysisResults.(animalID).Transitions.REMtoAwake = [];
    end
end
% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults')

end
