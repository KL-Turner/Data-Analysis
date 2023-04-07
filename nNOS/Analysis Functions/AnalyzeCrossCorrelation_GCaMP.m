function [Results_CrossCorr_GCaMP] = AnalyzeCrossCorrelation_GCaMP(animalID,group,set,rootFolder,delim,Results_CrossCorr_GCaMP)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataTypes = {'LH','RH'};
modelType = 'Forest';
params.minTime.Rest = 10;
params.minTime.NREM = 30;
params.minTime.REM = 60;
% only run analysis for valid animal IDs
dataLocation = [rootFolder delim 'Data' delim set delim group delim animalID delim 'Bilateral Imaging'];
cd(dataLocation)
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
% find and load RestingBaselines.mat strut
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
% find and load SleepData.mat strut
sleepDataFileStruct = dir('*_SleepData.mat');
sleepDataFile = {sleepDataFileStruct.name}';
sleepDataFileID = char(sleepDataFile);
load(sleepDataFileID,'-mat')
% lowpass filter
samplingRate = RestData.CBV_HbT.LH.CBVCamSamplingRate;
[z,p,k] = butter(4,1/(samplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
% criteria for resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
RestPuffCriteria.Fieldname = {'puffDistances'};
RestPuffCriteria.Comparison = {'gt'};
RestPuffCriteria.Value = {5};
% go through each valid data type for arousal-based cross-correlation analysis
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    % pull a few necessary numbers from the RestData.mat struct such as trial duration and sampling rate
    trialDuration_sec = RestData.CBV_HbT.LH.trialDuration_sec;
    % cross-correlation analysis for resting data
    % pull data from RestData.mat structure
    [restLogical] = FilterEvents_IOS(RestData.CBV_HbT.(dataType),RestCriteria);
    [puffLogical] = FilterEvents_IOS(RestData.CBV_HbT.(dataType),RestPuffCriteria);
    combRestLogical = logical(restLogical.*puffLogical);
    restFileIDs = RestData.CBV_HbT.(dataType).fileIDs(combRestLogical,:);
    restDurations = RestData.CBV_HbT.(dataType).durations(combRestLogical,:);
    restEventTimes = RestData.CBV_HbT.(dataType).eventTimes(combRestLogical,:);
    restingHbTData = RestData.CBV_HbT.(dataType).data(combRestLogical,:);
    restingGCaMPData = RestData.GCaMP7s.(['cor' dataType]).data(combRestLogical,:); % normdata is nan
    % keep only the data that occurs within the manually-approved awake regions
    [restFinalRestHbTData,restFinalFileIDs,restFinalDurations,restFinalEventTimes] = RemoveInvalidData_IOS(restingHbTData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    [restFinalRestGCaMPData,~,~,~] = RemoveInvalidData_IOS(restingGCaMPData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    cc = 1;
    for bb = 1:length(restFinalFileIDs)
        % check whether the event occurs in the appropriate time frame
        restStartTime = ceil(restFinalEventTimes(bb,1)*10)/10; % *10/10 used to round to first decimal place in a floor/ceil fashion.
        restDuration = floor(restFinalDurations(bb,1)*10)/10;
        if restStartTime >= 0.5 && (restStartTime + restDuration) <= (trialDuration_sec - 0.5)
            % remove the number of samples due to rounding up to start and rounding down to end. This is done to keep the HbT/GCaMP vectores aligned positionally with the upcoming
            % spectral analysis which is at 10 Hz
            leadSamples = round((restStartTime - restFinalEventTimes(bb,1))*samplingRate);
            lagSamples = round((restFinalDurations(bb,1) - restDuration)*samplingRate);
            % load in CBV_HbT from rest period
            restHbT = restFinalRestHbTData{bb,1};
            restGCaMP = restFinalRestGCaMPData{bb,1};
            % remove leading/lag samples due to rounding to nearest 0.1 up/0.1 down
            restSnipHbT = restHbT(1 + leadSamples:end - lagSamples);
            restSnipGCaMP = restGCaMP(1 + leadSamples:end - lagSamples);
            restFiltHbT = filtfilt(sos,g,detrend(restSnipHbT,'constant'));
            try
            restFiltGCaMP = filtfilt(sos,g,detrend(restSnipGCaMP,'constant'));
            catch
                keyboard
            end
            % only take the first 10 seconds of the epoch. occassionally a sample gets lost from rounding during the
            % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
            if length(restFiltHbT) < params.minTime.Rest*samplingRate
                restChunkSampleDiff = params.minTime.Rest*samplingRate - length(restFiltHbT);
                restPadHbT = (ones(1,restChunkSampleDiff))*restFiltHbT(end);
                restPadGCaMP = (ones(1,restChunkSampleDiff))*restFiltGCaMP(end);
                restShortHbT = horzcat(restFiltHbT,restPadHbT);
                restShortGCaMP = horzcat(restFiltGCaMP,restPadGCaMP);
            else
                restShortHbT = restFiltHbT(1:params.minTime.Rest*samplingRate);
                restShortGCaMP = restFiltGCaMP(1:params.minTime.Rest*samplingRate);
            end
            restProcData.HbT{cc,1} = restShortHbT;
            restProcData.GCaMP{cc,1} = restShortGCaMP;
            cc = cc + 1;
        end
        % set parameters for cross-correlation analysis
        restLagTime = 5; % seconds
        restFrequency = 10;
        restMaxLag = restLagTime*restFrequency;
        % run cross-correlation analysis - average through time
        for dd = 1:length(restProcData.HbT)
            restHbTarray = restProcData.HbT{dd,1};
            restGCaMParray = restProcData.GCaMP{dd,1};
            [restHbTvGCaMPxcVals(dd,:),restGCaMP_lags] = xcorr(restHbTarray,restGCaMParray,restMaxLag,'coeff');
        end
        restMeanHbTvGCaMPxcVals = mean(restHbTvGCaMPxcVals,1);
        restStdHbTvGCaMPxcVals = std(restHbTvGCaMPxcVals,0,1);
    end
    % save results
    Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).Rest.GCaMP_lags = restGCaMP_lags;
    Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).Rest.HbTvGCaMPxcVals = restMeanHbTvGCaMPxcVals;
    Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).Rest.HbTvGCaMPxcVals_std = restStdHbTvGCaMPxcVals;
    % cross-correlation analysis for NREM
    NREM_sleepTime = params.minTime.NREM; % seconds
    [NREM_finalHbT,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.CBV_HbT.(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    [NREM_finalGCaMP,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GCaMP7s.(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    % adjust[HbT] and GCaMP events to match the edits made to the length of each spectrogram
    mm = 1;
    for ll = 1:length(NREM_finalHbT)
        if isempty(NREM_finalHbT) == false
            NREM_HbTVals = NREM_finalHbT{ll,1}(1:NREM_sleepTime*samplingRate);
            NREM_GCaMPVals = NREM_finalGCaMP{ll,1}(1:NREM_sleepTime*samplingRate);
            NREM_finalHbTVals{mm,1} = filtfilt(sos,g,detrend(NREM_HbTVals,'constant'));
            NREM_finalGCaMPVals{mm,1} = filtfilt(sos,g,detrend(NREM_GCaMPVals,'constant'));
            mm = mm + 1;
        end
    end
    % run cross-correlation analysis - average through time
    NREM_lagTime = 5;   % Seconds
    NREM_frequency = 10;
    NREM_maxLag = NREM_lagTime*NREM_frequency;
    for nn = 1:length(NREM_finalHbTVals)
        NREM_HbT_array = NREM_finalHbTVals{nn,1};
        NREM_GCaMP_array = NREM_finalGCaMPVals{nn,1};
        [NREM_HbTvGCaMPxcVals(nn,:),NREM_GCaMP_lags] = xcorr(NREM_HbT_array,NREM_GCaMP_array,NREM_maxLag,'coeff');
    end
    NREM_meanHbTvGCaMPxcVals = mean(NREM_HbTvGCaMPxcVals,1);
    NREM_stdHbTvGCaMPxcVals = std(NREM_HbTvGCaMPxcVals,0,1);
    % save results
    Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).Rest.GCaMP_lags = NREM_GCaMP_lags;
    Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).Rest.HbTvGCaMPxcVals = NREM_meanHbTvGCaMPxcVals;
    Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).Rest.HbTvGCaMPxcVals_std = NREM_stdHbTvGCaMPxcVals;
    % cross-correlation analysis for REM
    REM_sleepTime = params.minTime.REM; % seconds
    [REM_finalHbT,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.CBV_HbT.(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    [REM_finalGCaMP,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GCaMP7s.(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    % adjust[HbT] and GCaMP events to match the edits made to the length of each spectrogram
    mm = 1;
    for ll = 1:length(REM_finalHbT)
        if isempty(REM_finalHbT) == false
            REM_HbTVals = REM_finalHbT{ll,1}(1:REM_sleepTime*samplingRate);
            REM_GCaMPVals = REM_finalGCaMP{ll,1}(1:REM_sleepTime*samplingRate);
            REM_finalHbTVals{mm,1} = filtfilt(sos,g,detrend(REM_HbTVals,'constant'));
            REM_finalGCaMPVals{mm,1} = filtfilt(sos,g,detrend(REM_GCaMPVals,'constant'));
            mm = mm + 1;
        end
    end
    % run cross-correlation analysis - average through time
    REM_lagTime = 5; % Seconds
    REM_frequency = 10;
    REM_maxLag = REM_lagTime*REM_frequency;
    for nn = 1:length(REM_finalHbTVals)
        REM_HbT_array = REM_finalHbTVals{nn,1};
        REM_GCaMP_array = REM_finalGCaMPVals{nn,1};
        [REM_HbTvGCaMPxcVals(nn,:),REM_GCaMP_lags] = xcorr(REM_HbT_array,REM_GCaMP_array,REM_maxLag,'coeff');
    end
    REM_meanHbTvGCaMPxcVals = mean(REM_HbTvGCaMPxcVals,1);
    REM_stdHbTvGCaMPxcVals = std(REM_HbTvGCaMPxcVals,0,1);
    % save results
    Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).Rest.GCaMP_lags = REM_GCaMP_lags;
    Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).Rest.HbTvGCaMPxcVals = REM_meanHbTvGCaMPxcVals;
    Results_CrossCorr_GCaMP.(group).(animalID).(hemisphere).Rest.HbTvGCaMPxcVals_std = REM_stdHbTvGCaMPxcVals;
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_CrossCorrGCaMP.mat','Results_CrossCorr_GCaMP')
cd([rootFolder delim 'Data'])

end
