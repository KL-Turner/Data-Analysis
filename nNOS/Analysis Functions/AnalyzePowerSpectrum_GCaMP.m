function [Results_PowerSpec_GCaMP] = AnalyzePowerSpectrum_GCaMP(animalID,group,set,rootFolder,delim,Results_PowerSpec_GCaMP)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
modelType = 'Forest';
params.minTime.Rest = 10;
params.minTime.NREM = 30;
params.minTime.REM = 60;
% only run analysis for valid animal IDs
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Bilateral Imaging'];
cd(dataLocation)
% character list of all ProcData file IDs
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% find and load RestData.mat struct
restDataFileStruct = dir('*_RestData.mat');
restDataFile = {restDataFileStruct.name}';
restDataFileID = char(restDataFile);
load(restDataFileID,'-mat')
% find and load manual baseline event information
manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
manualBaselineFile = {manualBaselineFileStruct.name}';
manualBaselineFileID = char(manualBaselineFile);
load(manualBaselineFileID,'-mat')
% find and load RestingBaselines.mat struct
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
% find and load SleepData.mat struct
sleepDataFileStruct = dir('*_SleepData.mat');
sleepDataFile = {sleepDataFileStruct.name}';
sleepDataFileID = char(sleepDataFile);
load(sleepDataFileID,'-mat')
% find and load Forest_ScoringResults.mat struct
forestScoringResultsFileID = [animalID '_Forest_ScoringResults.mat'];
load(forestScoringResultsFileID,'-mat')
% criteria for resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
RestPuffCriteria.Fieldname = {'puffDistances'};
RestPuffCriteria.Comparison = {'gt'};
RestPuffCriteria.Value = {5};
% loop variables
hemispheres = {'LH','RH','fLH','fRH'};
dataTypes = {'HbT','HbO','HbR','GCaMP7s'};
% go through each valid data type for arousal based coherence analysis
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        %% Rest
        clear restingData procRestData restData
        samplingRate = RestData.(dataType).(hemisphere).CBVCamSamplingRate;
        [restLogical] = FilterEvents_IOS(RestData.(dataType).(hemisphere),RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.(dataType).(hemisphere),RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.(dataType).(hemisphere).fileIDs(combRestLogical,:);
        restEventTimes = RestData.(dataType).(hemisphere).eventTimes(combRestLogical,:);
        restDurations = RestData.(dataType).(hemisphere).durations(combRestLogical,:);
        restingData = RestData.(dataType).(hemisphere).data(combRestLogical,:);
        % keep only the data that occurs within the manually-approved alert regions
        [restData,~,~,~] = RemoveInvalidData_IOS(restingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        % detrend and truncate data to minimum length to match events
        finalRestData = zeros(length(restData{1,1}),length(restData));
        for cc = 1:length(restData)
            if length(restData{cc,1}) < params.minTime.Rest*samplingRate
                restChunkSampleDiff = params.minTime.Rest*samplingRate - length(restData{cc,1});
                restPad = (ones(1,restChunkSampleDiff))*restData{cc,1}(end);
                procRestData{cc,1} = horzcat(restData{cc,1},restPad);
                procRestData{cc,1} = detrend(procRestData{cc,1},'constant');
                finalRestData(:,cc) = procRestData{cc,1};
            else
                procRestData{cc,1} = detrend(restData{cc,1}(1:(params.minTime.Rest*samplingRate)),'constant');
                finalRestData(:,cc) = procRestData{cc,1};
            end
        end
        % parameters for coherencyc - information available in function
        params.tapers = [1,1]; % Tapers [n, 2n - 1]
        params.pad = 1;
        params.Fs = samplingRate;
        params.fpass = [0,1]; % Pass band [0, nyquist]
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the power spectra of the desired signals
        [restS,restf,restsErr] = mtspectrumc(finalRestData,params);
        % save results
        Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).Rest.S = restS;
        Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).Rest.f = restf;
        Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).Rest.sErr = restsErr;
        %% Alert
        clear alertData procAlertData finalAlertData scoringLabels
        zz = 1;
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,~,fileID] = GetFileInfo_IOS(procDataFileID);
            for dd = 1:length(ScoringResults.fileIDs)
                if strcmp(fileID,ScoringResults.fileIDs{dd,1}) == true
                    scoringLabels = ScoringResults.labels{dd,1};
                end
            end
            % check labels to match arousal state
            if sum(strcmp(scoringLabels,'Not Sleep')) > 144 % 36 bins (180 total) or 3 minutes of asleep
                load(procDataFileID,'-mat')
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    alertData{zz,1} = ProcData.data.(dataType).(hemisphere);
                    zz = zz + 1;
                end
            end
        end
        % detrend data
        if isempty(alertData) == false
            finalAlertData = zeros(length(alertData{1,1}),length(alertData));
            for cc = 1:length(alertData)
                procAlertData{cc,1} = detrend(alertData{cc,1},'constant');
                finalAlertData(:,cc) = procAlertData{cc,1};
            end
            % calculate the power spectra of the desired signals
            [alertS,alertf,alertsErr] = mtspectrumc(finalAlertData,params);
            % save results
            Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).Alert.S = alertS;
            Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).Alert.f = alertf;
            Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).Alert.sErr = alertsErr;
        else
            % save results
            Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).Alert.S = [];
            Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).Alert.f = [];
            Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).Alert.sErr = [];
        end
        %% Asleep
        clear asleepData procAsleepData finalAsleepData scoringLabels
        zz = 1;
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,~,fileID] = GetFileInfo_IOS(procDataFileID);
            for dd = 1:length(ScoringResults.fileIDs)
                if strcmp(fileID,ScoringResults.fileIDs{dd,1}) == true
                    scoringLabels = ScoringResults.labels{dd,1};
                end
            end
            % check labels to match arousal state
            if sum(strcmp(scoringLabels,'Not Sleep')) < 36   % 36 bins (180 total) or 3 minutes of alert
                load(procDataFileID,'-mat')
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    asleepData{zz,1} = ProcData.data.(dataType).(hemisphere);
                    zz = zz + 1;
                end
            end
        end
        % detrend data
        if isempty(asleepData) == false
            finalAsleepData = zeros(length(asleepData{1,1}),length(asleepData));
            for cc = 1:length(asleepData)
                procAsleepData{cc,1} = detrend(asleepData{cc,1},'constant');
                finalAsleepData(:,cc) = procAsleepData{cc,1};
            end
            % calculate the power spectra of the desired signals
            params.tapers = [10,19]; % Tapers [n, 2n - 1]
            [asleepS,asleepf,asleepsErr] = mtspectrumc(finalAsleepData,params);
            % save results
            Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).Asleep.S = asleepS;
            Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).Asleep.f = asleepf;
            Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).Asleep.sErr = asleepsErr;
        else
            % save results
            Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).Asleep.S = [];
            Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).Asleep.f = [];
            Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).Asleep.sErr = [];
        end
        %% All
        clear allData procAllData finalAllData
        zz = 1;
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            load(procDataFileID,'-mat')
            puffs = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(puffs) == true
                allData{zz,1} = ProcData.data.(dataType).(hemisphere);
                zz = zz + 1;
            end
        end
        % detrend data
        finalAllData = zeros(length(allData{1,1}),length(allData));
        for cc = 1:length(allData)
            procAllData{cc,1} = detrend(allData{cc,1},'constant');
            finalAllData(:,cc) = procAllData{cc,1};
        end
        % calculate the power spectra of the desired signals
        params.tapers = [10,19]; % Tapers [n, 2n - 1]
        [allS,allf,allsErr] = mtspectrumc(finalAllData,params);
        % save results
        Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).All.S = allS;
        Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).All.f = allf;
        Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).All.sErr = allsErr;
        %% NREM
        clear nremData procNREMData finalNREMData
        [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.(dataType).(hemisphere),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        % detrend and truncate data to minimum length to match events
        finalNREM = zeros(length(nremData{1,1}),length(nremData));
        for ee = 1:length(nremData)
            procNREMData{ee,1} = detrend(nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
            finalNREM(:,ee) = procNREMData{ee,1};
        end
        % calculate the power spectra of the desired signals
        params.tapers = [3,5]; % Tapers [n, 2n - 1]
        [nremS,nremf,nremsErr] = mtspectrumc(finalNREM,params);
        % save results
        Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).NREM.S = nremS;
        Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).NREM.f = nremf;
        Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).NREM.sErr = nremsErr;
        %% REM
        clear remData procREMData finalREMData
        [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.(dataType).(hemisphere),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        % detrend and truncate data to minimum length to match events
        finalREM = zeros(length(remData{1,1}),length(remData));
        for ee = 1:length(remData)
            procREMData{ee,1} = detrend(remData{ee,1}(1:(params.minTime.REM*samplingRate)),'constant');
            finalREM(:,ee) = procREMData{ee,1};
        end
        % calculate the power spectra of the desired signals
        params.tapers = [5,9]; % Tapers [n, 2n - 1]
        [remS,remf,remsErr] = mtspectrumc(finalREM,params);
        % save results
        Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).REM.S = remS;
        Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).REM.f = remf;
        Results_PowerSpec_GCaMP.(group).(animalID).(hemisphere).(dataType).REM.sErr = remsErr;
    end
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_PowerSpec_GCaMP.mat','Results_PowerSpec_GCaMP')
cd([rootFolder delim 'Data'])