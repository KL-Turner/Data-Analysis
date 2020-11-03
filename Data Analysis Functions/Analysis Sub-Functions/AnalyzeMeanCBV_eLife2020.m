function [AnalysisResults] = AnalyzeMeanCBV_eLife2020(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the hemodynamic signal [HbT] during different arousal states (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
IOS_animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
Iso_animalIDs = {'T108','T109','T110','T111','T119','T120','T121','T122','T123'};
modelType = 'Forest';
params.minTime.Rest = 10;
params.Offset = 2;
params.minTime.Whisk = params.Offset + 5;
params.minTime.Stim = params.Offset + 2;
params.minTime.NREM = 30;
params.minTime.REM = 60;
%% only run analysis for valid animal IDs
if any(strcmp(IOS_animalIDs,animalID))
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
    samplingRate = RestData.CBV.adjLH.CBVCamSamplingRate;
    [z,p,k] = butter(4,1/(samplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    % criteria for the whisking
    WhiskCriteria.Fieldname = {'duration','duration','puffDistance'};
    WhiskCriteria.Comparison = {'gt','lt','gt'};
    WhiskCriteria.Value = {2,5,5};
    WhiskPuffCriteria.Fieldname = {'puffDistance'};
    WhiskPuffCriteria.Comparison = {'gt'};
    WhiskPuffCriteria.Value = {5};
    % criteria for the resting
    RestCriteria.Fieldname = {'durations'};
    RestCriteria.Comparison = {'gt'};
    RestCriteria.Value = {params.minTime.Rest};
    RestPuffCriteria.Fieldname = {'puffDistances'};
    RestPuffCriteria.Comparison = {'gt'};
    RestPuffCriteria.Value = {5};
    % criteria for the stimulation
    StimCriteriaA.Value = {'RPadSol'};
    StimCriteriaA.Fieldname = {'solenoidName'};
    StimCriteriaA.Comparison = {'equal'};
    StimCriteriaB.Value = {'LPadSol'};
    StimCriteriaB.Fieldname = {'solenoidName'};
    StimCriteriaB.Comparison = {'equal'};
    %% analyze [HbT] during periods of rest
    % pull data from RestData.mat structure
    [restLogical] = FilterEvents_IOS_eLife2020(RestData.CBV_HbT.adjLH,RestCriteria);
    [puffLogical] = FilterEvents_IOS_eLife2020(RestData.CBV_HbT.adjLH,RestPuffCriteria);
    combRestLogical = logical(restLogical.*puffLogical);
    restFileIDs = RestData.CBV_HbT.adjLH.fileIDs(combRestLogical,:);
    restEventTimes = RestData.CBV_HbT.adjLH.eventTimes(combRestLogical,:);
    restDurations = RestData.CBV_HbT.adjLH.durations(combRestLogical,:);
    LH_RestingData = RestData.CBV_HbT.adjLH.data(combRestLogical,:);
    RH_RestingData = RestData.CBV_HbT.adjRH.data(combRestLogical,:);
    % keep only the data that occurs within the manually-approved awake regions
    [LH_finalRestData,finalRestFileIDs,~,~] = RemoveInvalidData_IOS_eLife2020(LH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    [RH_finalRestData,~,~,~] = RemoveInvalidData_IOS_eLife2020(RH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    % filter [HbT]
    for gg = 1:length(LH_finalRestData)
        LH_ProcRestData{gg,1} = filtfilt(sos,g,LH_finalRestData{gg,1}); %#ok<*AGROW>
        RH_ProcRestData{gg,1} = filtfilt(sos,g,RH_finalRestData{gg,1});
    end
    % take mean [HbT] during resting epochs
    for nn = 1:length(LH_ProcRestData)
        LH_restCBVMean(nn,1) = mean(LH_ProcRestData{nn,1}(1:end));
        RH_restCBVMean(nn,1) = mean(RH_ProcRestData{nn,1}(1:end));
    end
    % save results
    AnalysisResults.(animalID).MeanCBV.Rest.CBV_HbT.MeanAdjLH = LH_restCBVMean;
    AnalysisResults.(animalID).MeanCBV.Rest.CBV_HbT.MeanAdjRH = RH_restCBVMean;
    AnalysisResults.(animalID).MeanCBV.Rest.CBV_HbT.IndAdjLH = LH_ProcRestData;
    AnalysisResults.(animalID).MeanCBV.Rest.CBV_HbT.IndAdjRH = RH_ProcRestData;
    AnalysisResults.(animalID).MeanCBV.Rest.CBV_HbT.FileIDs = finalRestFileIDs;
    %% analyze [HbT] during periods of moderate whisking (2-5 seconds)
    % pull data from EventData.mat structure
    [whiskLogical] = FilterEvents_IOS_eLife2020(EventData.CBV_HbT.adjLH.whisk,WhiskCriteria);
    [puffLogical] = FilterEvents_IOS_eLife2020(EventData.CBV_HbT.adjLH.whisk,WhiskPuffCriteria);
    combWhiskLogical = logical(whiskLogical.*puffLogical);
    whiskFileIDs = EventData.CBV_HbT.adjLH.whisk.fileIDs(combWhiskLogical,:);
    whiskEventTimes = EventData.CBV_HbT.adjLH.whisk.eventTime(combWhiskLogical,:);
    whiskDurations = EventData.CBV_HbT.adjLH.whisk.duration(combWhiskLogical,:);
    LH_whiskData = EventData.CBV_HbT.adjLH.whisk.data(combWhiskLogical,:);
    RH_whiskData = EventData.CBV_HbT.adjRH.whisk.data(combWhiskLogical,:);
    % keep only the data that occurs within the manually-approved awake regions
    [LH_finalWhiskData,finalWhiskFileIDs,~,~] = RemoveInvalidData_IOS_eLife2020(LH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    [RH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS_eLife2020(RH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    % filter [HbT] and mean-subtract 2 seconds prior to whisk
    for gg = 1:size(LH_finalWhiskData,1)
        LH_ProcWhiskData_temp = filtfilt(sos,g,LH_finalWhiskData(gg,:));
        LH_ProcWhiskData(gg,:) = LH_ProcWhiskData_temp - mean(LH_ProcWhiskData_temp(1:params.Offset*samplingRate));
        RH_ProcWhiskData_temp = filtfilt(sos,g,RH_finalWhiskData(gg,:));
        RH_ProcWhiskData(gg,:) = RH_ProcWhiskData_temp - mean(RH_ProcWhiskData_temp(1:params.Offset*samplingRate));
    end
    % take mean [HbT] during whisking epochs from onset through 5 seconds
    for nn = 1:size(LH_ProcWhiskData,1)
        LH_whiskCBVMean{nn,1} = mean(LH_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
        RH_whiskCBVMean{nn,1} = mean(RH_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
        LH_whiskCBV{nn,1} = LH_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
        RH_whiskCBV{nn,1} = RH_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
    end
    % save results
    AnalysisResults.(animalID).MeanCBV.Whisk.CBV_HbT.MeanAdjLH = cell2mat(LH_whiskCBVMean);
    AnalysisResults.(animalID).MeanCBV.Whisk.CBV_HbT.MeanAdjRH = cell2mat(RH_whiskCBVMean);
    AnalysisResults.(animalID).MeanCBV.Whisk.CBV_HbT.IndAdjLH = LH_whiskCBV;
    AnalysisResults.(animalID).MeanCBV.Whisk.CBV_HbT.IndAdjRH = RH_whiskCBV;
    AnalysisResults.(animalID).MeanCBV.Whisk.CBV_HbT.FileIDs = finalWhiskFileIDs;
    %% analyze [HbT] during periods of stimulation
    % pull data from EventData.mat structure
    LH_stimFilter = FilterEvents_IOS_eLife2020(EventData.CBV_HbT.adjLH.stim,StimCriteriaA);
    RH_stimFilter = FilterEvents_IOS_eLife2020(EventData.CBV_HbT.adjRH.stim,StimCriteriaB);
    [LH_stimHbTData] = EventData.CBV_HbT.adjLH.stim.data(LH_stimFilter,:);
    [RH_stimHbTData] = EventData.CBV_HbT.adjRH.stim.data(RH_stimFilter,:);
    [LH_stimFileIDs] = EventData.CBV_HbT.adjLH.stim.fileIDs(LH_stimFilter,:);
    [RH_stimFileIDs] = EventData.CBV_HbT.adjRH.stim.fileIDs(RH_stimFilter,:);
    [LH_stimEventTimes] = EventData.CBV_HbT.adjLH.stim.eventTime(LH_stimFilter,:);
    [RH_stimEventTimes] = EventData.CBV_HbT.adjRH.stim.eventTime(RH_stimFilter,:);
    LH_stimDurations = zeros(length(LH_stimEventTimes),1);
    RH_stimDurations = zeros(length(RH_stimEventTimes),1);
    % keep only the data that occurs within the manually-approved awake regions
    [LH_finalStimData,LH_finalStimFileIDs,~,~] = RemoveInvalidData_IOS_eLife2020(LH_stimHbTData,LH_stimFileIDs,LH_stimDurations,LH_stimEventTimes,ManualDecisions);
    [RH_finalStimData,RH_finalStimFileIDs,~,~] = RemoveInvalidData_IOS_eLife2020(RH_stimHbTData,RH_stimFileIDs,RH_stimDurations,RH_stimEventTimes,ManualDecisions);
    % filter [HbT] and mean-subtract 2 seconds prior to stimulus (left hem)
    for gg = 1:size(LH_finalStimData,1)
        LH_ProcStimData_temp = filtfilt(sos,g,LH_finalStimData(gg,:));
        LH_ProcStimData(gg,:) = LH_ProcStimData_temp - mean(LH_ProcStimData_temp(1:params.Offset*samplingRate));
    end
    % filter [HbT] and mean-subtract 2 seconds prior to stimulus (right hem)
    for gg = 1:size(RH_finalStimData,1)
        RH_ProcStimData_temp = filtfilt(sos,g,RH_finalStimData(gg,:));
        RH_ProcStimData(gg,:) = RH_ProcStimData_temp - mean(RH_ProcStimData_temp(1:params.Offset*samplingRate));
    end
    % take mean [HbT] 1-2 seconds after stimulation (left hem)
    for nn = 1:size(LH_ProcStimData,1)
        LH_stimCBVMean{nn,1} = mean(LH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate),2);
        LH_stimCBV{nn,1} = LH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate);
    end
    % take mean [HbT] 1-2 seconds after stimulation (right hem)
    for nn = 1:size(RH_ProcStimData,1)
        RH_stimCBVMean{nn,1} = mean(RH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate),2);
        RH_stimCBV{nn,1} = RH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate);
    end
    % save results
    AnalysisResults.(animalID).MeanCBV.Stim.CBV_HbT.MeanAdjLH = cell2mat(LH_stimCBVMean);
    AnalysisResults.(animalID).MeanCBV.Stim.CBV_HbT.MeanAdjRH = cell2mat(RH_stimCBVMean);
    AnalysisResults.(animalID).MeanCBV.Stim.CBV_HbT.IndAdjLH = LH_stimCBV;
    AnalysisResults.(animalID).MeanCBV.Stim.CBV_HbT.IndAdjRH = RH_stimCBV;
    AnalysisResults.(animalID).MeanCBV.Stim.CBV_HbT.LH_FileIDs = LH_finalStimFileIDs;
    AnalysisResults.(animalID).MeanCBV.Stim.CBV_HbT.RH_FileIDs = RH_finalStimFileIDs;
    %% analyze [HbT] during periods of NREM sleep
    % pull data from SleepData.mat structure
    [LH_nremData,nremFileIDs,~] = RemoveStimSleepData_IOS_eLife2020(animalID,SleepData.(modelType).NREM.data.CBV_HbT.LH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    [RH_nremData,~,~] = RemoveStimSleepData_IOS_eLife2020(animalID,SleepData.(modelType).NREM.data.CBV_HbT.RH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    % filter and take mean [HbT] during NREM epochs
    for nn = 1:length(LH_nremData)
        LH_nremCBVMean(nn,1) = mean(filtfilt(sos,g,LH_nremData{nn,1}(1:end)));
        RH_nremCBVMean(nn,1) = mean(filtfilt(sos,g,RH_nremData{nn,1}(1:end)));
    end
    % save results
    AnalysisResults.(animalID).MeanCBV.NREM.CBV_HbT.MeanAdjLH = LH_nremCBVMean;
    AnalysisResults.(animalID).MeanCBV.NREM.CBV_HbT.MeanAdjRH = RH_nremCBVMean;
    AnalysisResults.(animalID).MeanCBV.NREM.CBV_HbT.IndAdjLH = LH_nremData;
    AnalysisResults.(animalID).MeanCBV.NREM.CBV_HbT.IndAdjRH = RH_nremData;
    AnalysisResults.(animalID).MeanCBV.NREM.CBV_HbT.FileIDs = nremFileIDs;
    %% analyze [HbT] during periods of REM sleep
    % pull data from SleepData.mat structure
    [LH_remData,remFileIDs,~] = RemoveStimSleepData_IOS_eLife2020(animalID,SleepData.(modelType).REM.data.CBV_HbT.LH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    [RH_remData,~,~] = RemoveStimSleepData_IOS_eLife2020(animalID,SleepData.(modelType).REM.data.CBV_HbT.RH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    % filter and take mean [HbT] during REM epochs
    for nn = 1:length(LH_remData)
        LH_remCBVMean(nn,1) = mean(filtfilt(sos,g,LH_remData{nn,1}(1:end)));
        RH_remCBVMean(nn,1) = mean(filtfilt(sos,g,RH_remData{nn,1}(1:end)));
    end
    % save results
    AnalysisResults.(animalID).MeanCBV.REM.CBV_HbT.MeanAdjLH = LH_remCBVMean;
    AnalysisResults.(animalID).MeanCBV.REM.CBV_HbT.MeanAdjRH = RH_remCBVMean;
    AnalysisResults.(animalID).MeanCBV.REM.CBV_HbT.IndAdjLH = LH_remData;
    AnalysisResults.(animalID).MeanCBV.REM.CBV_HbT.IndAdjRH = RH_remData;
    AnalysisResults.(animalID).MeanCBV.REM.CBV_HbT.FileIDs = remFileIDs;
    %% analyze [HbT] during periods of isolfurane
    if any(strcmp(Iso_animalIDs,animalID))
        dataLocation = [rootFolder '\' animalID '\Isoflurane Trials\'];
        cd(dataLocation)
        % pull ProcData.mat file associated with isoflurane administration
        procDataFileStruct = dir('*_ProcData.mat');
        procDataFile = {procDataFileStruct.name}';
        procDataFileID = char(procDataFile);
        load(procDataFileID)
        % extract left and right [HbT] changes during the last 100 seconds of data
        isoLH_HbT = ProcData.data.CBV_HbT.adjLH((end - samplingRate*100):end);
        filtIsoLH_HbT = filtfilt(sos,g,isoLH_HbT);
        isoRH_HbT = ProcData.data.CBV_HbT.adjRH((end - samplingRate*100):end);
        filtIsoRH_HbT = filtfilt(sos,g,isoRH_HbT);
        % save results
        AnalysisResults.(animalID).MeanCBV.Iso.CBV_HbT.adjLH = mean(filtIsoLH_HbT);
        AnalysisResults.(animalID).MeanCBV.Iso.CBV_HbT.adjRH = mean(filtIsoRH_HbT);
        AnalysisResults.(animalID).MeanCBV.Iso.CBV_HbT.FileIDs = procDataFileID;
    end
    % save data
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end

end
