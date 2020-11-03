function [AnalysisResults] = AnalyzeMeanVesselDiameter_eLife2020(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the arteriole diameter D/D during different arousal states (2PLSM)
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs = {'T115','T116','T117','T118','T125','T126'};
modelType = 'Manual';
params.Offset = 2;
params.minTime.Rest = 10;
params.minTime.Whisk = params.Offset + 5;
params.minTime.NREM = 30;
params.minTime.REM = 60;
%% only run analysis for valid animal IDs
if any(strcmp(animalIDs,animalID))
    dataLocation = [rootFolder '\' animalID '\2P Data\'];
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
    % find and load EventData.mat struct
    eventDataFileStruct = dir('*_EventData.mat');
    eventDataFile = {eventDataFileStruct.name}';
    eventDataFileID = char(eventDataFile);
    load(eventDataFileID)
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
    % criteria for whisking
    WhiskCriteria.Fieldname = {'duration','duration'};
    WhiskCriteria.Comparison = {'gt','lt'};
    WhiskCriteria.Value = {2,5};
    % criteria for resting
    RestCriteria.Fieldname = {'durations'};
    RestCriteria.Comparison = {'gt'};
    RestCriteria.Value = {params.minTime.Rest};
    %% analyze arteriole D/D during periods of rest
    % pull data from RestData.mat structure
    [restLogical] = FilterEvents_2P_eLife2020(RestData.vesselDiameter.data,RestCriteria);
    combRestLogical = logical(restLogical);
    restVesselData = RestData.vesselDiameter.data.data(combRestLogical,:);
    restFileIDs = RestData.vesselDiameter.data.fileIDs(combRestLogical,:);
    restVesselIDs = RestData.vesselDiameter.data.vesselIDs(combRestLogical,:);
    restDurations = RestData.vesselDiameter.data.durations(combRestLogical,:);
    restEventTimes = RestData.vesselDiameter.data.eventTimes(combRestLogical,:);
    % keep only the data that occurs within the manually-approved awake regions
    [finalRestVesselData,finalRestFileIDs,finalRestVesselIDs,~,~] = RemoveInvalidData_2P_eLife2020(restVesselData,restFileIDs,restVesselIDs,restDurations,restEventTimes,ManualDecisions);
    % go through the data and normalize + filter each rest epoch based on individual vessels
    uniqueRestVesselIDs = unique(finalRestVesselIDs);
    for aa = 1:length(uniqueRestVesselIDs)
        cc = 1;
        for bb = 1:length(finalRestVesselIDs)
            if strcmp(uniqueRestVesselIDs{aa,1},finalRestVesselIDs{bb,1})
                strDay = ConvertDate_2P_eLife2020(finalRestFileIDs{bb,1}(1:6));
                tempRestData.(uniqueRestVesselIDs{aa,1}){cc,1} = filtfilt(sos,g,((finalRestVesselData{bb,1} - RestingBaselines.manualSelection.vesselDiameter.data.(uniqueRestVesselIDs{aa,1}).(strDay))/RestingBaselines.manualSelection.vesselDiameter.data.(uniqueRestVesselIDs{aa,1}).(strDay)));
                tempRestFileIDs.(uniqueRestVesselIDs{aa,1}){cc,1} = finalRestFileIDs{bb,1};
                cc = cc + 1;
            end
        end
    end
    % take the average of each vessel's individual resting event
    for dd = 1:length(uniqueRestVesselIDs)
        for ee = 1:length(tempRestData.(uniqueRestVesselIDs{dd,1}))
            tempRestDataMeans.(uniqueRestVesselIDs{dd,1})(ee,1) = mean(tempRestData.(uniqueRestVesselIDs{dd,1}){ee,1})*100;
            tempRestDataMaxs.(uniqueRestVesselIDs{dd,1})(ee,1) = max(tempRestData.(uniqueRestVesselIDs{dd,1}){ee,1})*100;
            tempRestDataInd.(uniqueRestVesselIDs{dd,1}){ee,1} = tempRestData.(uniqueRestVesselIDs{dd,1}){ee,1}*100;
        end
    end
    % take the average of each vessel's total resting events
    for ff = 1:length(uniqueRestVesselIDs)
        % save results
        AnalysisResults.(animalID).MeanVesselDiameter.Rest.(uniqueRestVesselIDs{ff,1}).indEvents = tempRestDataInd.(uniqueRestVesselIDs{ff,1});
        AnalysisResults.(animalID).MeanVesselDiameter.Rest.(uniqueRestVesselIDs{ff,1}).mean = tempRestDataMeans.(uniqueRestVesselIDs{ff,1});
        AnalysisResults.(animalID).MeanVesselDiameter.Rest.(uniqueRestVesselIDs{ff,1}).max = tempRestDataMaxs.(uniqueRestVesselIDs{ff,1});
        AnalysisResults.(animalID).MeanVesselDiameter.Rest.(uniqueRestVesselIDs{ff,1}).fileIDs = tempRestFileIDs.(uniqueRestVesselIDs{ff,1});
    end
    %% analyze arteriole D/D during periods of moderate whisking (2-5 seconds)
    % pull data from EventData.mat structure
    [whiskLogical] = FilterEvents_2P_eLife2020(EventData.vesselDiameter.data.whisk,WhiskCriteria);
    combWhiskLogical = logical(whiskLogical);
    whiskVesselData = EventData.vesselDiameter.data.whisk.data(combWhiskLogical,:);
    whiskFileIDs = EventData.vesselDiameter.data.whisk.fileIDs(combWhiskLogical,:);
    whiskVesselIDs = EventData.vesselDiameter.data.whisk.vesselIDs(combWhiskLogical,:);
    whiskDurations = EventData.vesselDiameter.data.whisk.duration(combWhiskLogical,:);
    whiskEventTimes = EventData.vesselDiameter.data.whisk.eventTime(combWhiskLogical,:);
    % keep only the data that occurs within the manually-approved awake regions
    [finalWhiskVesselData,finalWhiskFileIDs,finalWhiskVesselIDs,~,~] = RemoveInvalidData_2P_eLife2020(whiskVesselData,whiskFileIDs,whiskVesselIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    % go through the data and normalize + filter each whisk event based on individual vessels
    uniqueWhiskVesselIDs = unique(finalWhiskVesselIDs);
    for aa = 1:length(uniqueWhiskVesselIDs)
        cc = 1;
        for bb = 1:length(finalWhiskVesselIDs)
            if strcmp(uniqueWhiskVesselIDs{aa,1},finalWhiskVesselIDs{bb,1})
                strDay = ConvertDate_2P_eLife2020(finalWhiskFileIDs{bb,1}(1:6));
                tempWhiskData.(uniqueWhiskVesselIDs{aa,1}){cc,1} = filtfilt(sos,g,((finalWhiskVesselData(bb,:) - RestingBaselines.manualSelection.vesselDiameter.data.(uniqueWhiskVesselIDs{aa,1}).(strDay))/RestingBaselines.manualSelection.vesselDiameter.data.(uniqueWhiskVesselIDs{aa,1}).(strDay)));
                tempWhiskFileIDs.(uniqueWhiskVesselIDs{aa,1}){cc,1} = finalWhiskFileIDs{bb,1};
                cc = cc + 1;
            end
        end
    end
    % mean-subtract 2 seconds prior to whisk
    for dd = 1:length(uniqueWhiskVesselIDs)
        for ee = 1:length(tempWhiskData.(uniqueWhiskVesselIDs{dd,1}))
            tempWhiskDataB.(uniqueWhiskVesselIDs{dd,1}){ee,1} = tempWhiskData.(uniqueWhiskVesselIDs{dd,1}){ee,1} - mean(tempWhiskData.(uniqueWhiskVesselIDs{dd,1}){ee,1}(1:params.Offset*samplingRate));
        end
    end
    % take the average of each vessel's individual whisking event from onset through 5 seconds
    for dd = 1:length(uniqueWhiskVesselIDs)
        for ee = 1:length(tempWhiskData.(uniqueWhiskVesselIDs{dd,1}))
            tempWhiskDataMeans.(uniqueWhiskVesselIDs{dd,1})(ee,1) = mean(tempWhiskDataB.(uniqueWhiskVesselIDs{dd,1}){ee,1}(params.Offset*samplingRate:params.minTime.Whisk*samplingRate))*100;
            tempWhiskDataMaxs.(uniqueWhiskVesselIDs{dd,1})(ee,1) = max(tempWhiskDataB.(uniqueWhiskVesselIDs{dd,1}){ee,1}(params.Offset*samplingRate:params.minTime.Whisk*samplingRate))*100;
            tempWhiskDataInd.(uniqueWhiskVesselIDs{dd,1}){ee,1} = tempWhiskDataB.(uniqueWhiskVesselIDs{dd,1}){ee,1}(params.Offset*samplingRate:params.minTime.Whisk*samplingRate)*100;
        end
    end
    % take the average of each vessel's total whisking events
    for ff = 1:length(uniqueWhiskVesselIDs)
        % save results
        AnalysisResults.(animalID).MeanVesselDiameter.Whisk.(uniqueWhiskVesselIDs{ff,1}).indEvents = tempWhiskDataInd.(uniqueWhiskVesselIDs{ff,1});
        AnalysisResults.(animalID).MeanVesselDiameter.Whisk.(uniqueWhiskVesselIDs{ff,1}).mean = tempWhiskDataMeans.(uniqueWhiskVesselIDs{ff,1});
        AnalysisResults.(animalID).MeanVesselDiameter.Whisk.(uniqueWhiskVesselIDs{ff,1}).max = tempWhiskDataMaxs.(uniqueWhiskVesselIDs{ff,1});
        AnalysisResults.(animalID).MeanVesselDiameter.Whisk.(uniqueWhiskVesselIDs{ff,1}).fileIDs = tempWhiskFileIDs.(uniqueWhiskVesselIDs{ff,1});
    end
    AnalysisResults.(animalID).MeanVesselDiameter.Whisk.allFileIDs = finalWhiskFileIDs;
    %% analyze arteriole D/D during periods of NREM sleep
    if isfield(SleepData.(modelType),'NREM') == true
        % pull data from SleepData.mat structure
        nremVesselData = SleepData.(modelType).NREM.data.vesselDiameter.data;
        nremFileIDs = SleepData.(modelType).NREM.FileIDs;
        nremVesselIDs = SleepData.(modelType).NREM.VesselIDs;
        % filter each event based on individual vessels
        uniqueNREMVesselIDs = unique(nremVesselIDs);
        for aa = 1:length(uniqueNREMVesselIDs)
            cc = 1;
            for bb = 1:length(nremVesselIDs)
                if strcmp(uniqueNREMVesselIDs{aa,1},nremVesselIDs{bb,1})
                    tempNREMData.(uniqueNREMVesselIDs{aa,1}){cc,1} = filtfilt(sos,g,nremVesselData{bb,1});
                    tempNREMFileIDs.(uniqueNREMVesselIDs{aa,1}){cc,1} = nremFileIDs{bb,1};
                    cc = cc + 1;
                end
            end
        end
        % take the average of each vessel's individual NREM events
        for dd = 1:length(uniqueNREMVesselIDs)
            for ee = 1:length(tempNREMData.(uniqueNREMVesselIDs{dd,1}))
                tempNREMDataMeans.(uniqueNREMVesselIDs{dd,1})(ee,1) = mean(tempNREMData.(uniqueNREMVesselIDs{dd,1}){ee,1})*100;
                tempNREMDataMaxs.(uniqueNREMVesselIDs{dd,1})(ee,1) = max(tempNREMData.(uniqueNREMVesselIDs{dd,1}){ee,1})*100;
                tempNREMDataInd.(uniqueNREMVesselIDs{dd,1}){ee,1} = tempNREMData.(uniqueNREMVesselIDs{dd,1}){ee,1}*100;
            end
        end
        % take the average of each vessel's total NREM events
        for ff = 1:length(uniqueNREMVesselIDs)
            % save results
            AnalysisResults.(animalID).MeanVesselDiameter.NREM.(uniqueNREMVesselIDs{ff,1}).indEvents = tempNREMDataInd.(uniqueNREMVesselIDs{ff,1});
            AnalysisResults.(animalID).MeanVesselDiameter.NREM.(uniqueNREMVesselIDs{ff,1}).mean = tempNREMDataMeans.(uniqueNREMVesselIDs{ff,1});
            AnalysisResults.(animalID).MeanVesselDiameter.NREM.(uniqueNREMVesselIDs{ff,1}).max = tempNREMDataMaxs.(uniqueNREMVesselIDs{ff,1});
            AnalysisResults.(animalID).MeanVesselDiameter.NREM.(uniqueNREMVesselIDs{ff,1}).fileIDs = tempNREMFileIDs.(uniqueNREMVesselIDs{ff,1});
        end
    end
    %% analyze arteriole D/D during periods of REM sleep
    if isfield(SleepData.(modelType),'REM') == true
        % pull data from SleepData.mat structure
        remVesselData = SleepData.(modelType).REM.data.vesselDiameter.data;
        remFileIDs = SleepData.(modelType).REM.FileIDs;
        remVesselIDs = SleepData.(modelType).REM.VesselIDs;
        % filter each event based on individual vessels
        uniqueREMVesselIDs = unique(remVesselIDs);
        for aa = 1:length(uniqueREMVesselIDs)
            cc = 1;
            for bb = 1:length(remVesselIDs)
                if strcmp(uniqueREMVesselIDs{aa,1},remVesselIDs{bb,1})
                    tempREMData.(uniqueREMVesselIDs{aa,1}){cc,1} = filtfilt(sos,g,remVesselData{bb,1});
                    tempREMFileIDs.(uniqueREMVesselIDs{aa,1}){cc,1} = remFileIDs{bb,1};
                    cc = cc + 1;
                end
            end
        end
        % take the average of each vessel's individual REM events
        for dd = 1:length(uniqueREMVesselIDs)
            for ee = 1:length(tempREMData.(uniqueREMVesselIDs{dd,1}))
                tempREMDataMeans.(uniqueREMVesselIDs{dd,1})(ee,1) = mean(tempREMData.(uniqueREMVesselIDs{dd,1}){ee,1})*100;
                tempREMDataMaxs.(uniqueREMVesselIDs{dd,1})(ee,1) = max(tempREMData.(uniqueREMVesselIDs{dd,1}){ee,1})*100;
                tempREMDataInd.(uniqueREMVesselIDs{dd,1}){ee,1} = tempREMData.(uniqueREMVesselIDs{dd,1}){ee,1}*100;
            end
        end
        % take the average of each vessel's total resting events
        for ff = 1:length(uniqueREMVesselIDs)
            % save results
            AnalysisResults.(animalID).MeanVesselDiameter.REM.(uniqueREMVesselIDs{ff,1}).indEvents = tempREMDataInd.(uniqueREMVesselIDs{ff,1});
            AnalysisResults.(animalID).MeanVesselDiameter.REM.(uniqueREMVesselIDs{ff,1}).mean = tempREMDataMeans.(uniqueREMVesselIDs{ff,1});
            AnalysisResults.(animalID).MeanVesselDiameter.REM.(uniqueREMVesselIDs{ff,1}).max = tempREMDataMaxs.(uniqueREMVesselIDs{ff,1});
            AnalysisResults.(animalID).MeanVesselDiameter.REM.(uniqueREMVesselIDs{ff,1}).fileIDs = tempREMFileIDs.(uniqueREMVesselIDs{ff,1});
        end
    end
    % save data
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end


end
