function [Results_CrossCorrelation] = AnalyzeCrossCorrelation_Pupil(animalID,rootFolder,delim,Results_CrossCorrelation)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the cross-correlation between neural activity and hemodynamics [HbT] (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
dataTypes = {'mmArea','mmDiameter','zArea','zDiameter'};
hemispheres = {'LH_HbT','RH_HbT'};
modelType = 'Forest';
params.minTime.Rest = 10;
params.minTime.NREM = 30;
params.minTime.REM = 60;
%% only run analysis for valid animal IDs
dataLocation = [rootFolder delim 'Data' delim animalID delim 'Bilateral Imaging'];
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
% character list of all ProcData file IDs
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
%
scoringResultsFileStruct = dir('*Forest_ScoringResults.mat');
scoringResultsFile = {scoringResultsFileStruct.name}';
scoringResultsFileID = char(scoringResultsFile);
load(scoringResultsFileID,'-mat')
% go through each valid data type for arousal-based cross-correlation analysis
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        % pull a few necessary numbers from the RestData.mat struct such as trial duration and sampling rate
        oneSecSpecFs = 30;   % Hz
        %% cross-correlation analysis for resting data
        % pull data from RestData.mat structure
        [restLogical] = FilterEvents_IOS(RestData.Pupil.(dataType),RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.Pupil.(dataType),RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.Pupil.(dataType).fileIDs(combRestLogical,:);
        restDurations = RestData.Pupil.(dataType).durations(combRestLogical,:);
        restEventTimes = RestData.Pupil.(dataType).eventTimes(combRestLogical,:);
        HbT_restData = RestData.Pupil.(hemisphere).data(combRestLogical,:);
        Pupil_restData = RestData.Pupil.(dataType).data(combRestLogical,:);
        % keep only the data that occurs within the manually-approved awake regions
        [HbT_finalRestData,~,~,~] = RemoveInvalidData_IOS(HbT_restData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [Pupil_finalRestData,~,~,~] = RemoveInvalidData_IOS(Pupil_restData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        cc = 1;
        restProcData = [];
        for dd = 1:length(HbT_finalRestData)
            if sum(isnan(Pupil_finalRestData{dd,1})) == 0
                if length(HbT_finalRestData{bb,1}) < params.minTime.Rest*samplingRate
                    restChunkSampleDiff = params.minTime.Rest*samplingRate - length(HbT_finalRestData{bb,1});
                    HbT_restPad = (ones(1,restChunkSampleDiff))*HbT_finalRestData{bb,1}(end);
                    Pupil_restPad = (ones(1,restChunkSampleDiff))*Pupil_finalRestData{bb,1}(end);
                    HbT_ProcRestData = horzcat(HbT_finalRestData{bb,1},HbT_restPad); %#ok<*AGROW>
                    Pupil_ProcRestData = horzcat(Pupil_finalRestData{bb,1},Pupil_restPad);
                    HbT_ProcRestData = filtfilt(sos,g,detrend(HbT_ProcRestData{bb,1},'constant'));
                    Pupil_ProcRestData = filtfilt(sos,g,detrend(Pupil_ProcRestData{bb,1},'constant'));
                else
                    HbT_ProcRestData = filtfilt(sos,g,detrend(HbT_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant'));
                    Pupil_ProcRestData = filtfilt(sos,g,detrend(Pupil_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant'));
                end
                restProcData.HbT{cc,1} = HbT_ProcRestData;
                restProcData.Pupil{cc,1} = Pupil_ProcRestData;
                cc = cc + 1;
            end
        end
        if isempty(restProcData) == false
            % set parameters for cross-correlation analysis
            restLagTime = 5;   % seconds
            restFrequency = oneSecSpecFs;   % Hz
            restMaxLag = restLagTime*restFrequency;
            % run cross-correlation analysis - average through time
            for dd = 1:length(restProcData.HbT)
                restHbTarray = restProcData.HbT{dd,1};
                restPupilarray = restProcData.Pupil{dd,1};
                [restXcVals(dd,:),restPupil_lags] = xcorr(restHbTarray,restPupilarray,restMaxLag,'coeff'); %#ok<*AGROW>
            end
            restMeanXcVals = mean(restXcVals,1);
            % save results
            Results_CrossCorrelation.(animalID).Rest.(hemisphere).(dataType).lags = restPupil_lags;
            Results_CrossCorrelation.(animalID).Rest.(hemisphere).(dataType).xcVals = restMeanXcVals;
        else
            % save results
            Results_CrossCorrelation.(animalID).Rest.(hemisphere).(dataType).lags = [];
            Results_CrossCorrelation.(animalID).Rest.(hemisphere).(dataType).xcVals = [];
        end
        %% analyze neural-hemo coherence during periods of alert
        zz = 1;
        clear HbT_awakeData Pupil_awakeData
        HbT_awakeData = []; Pupil_awakeData = [];
        for dd = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(dd,:);
            [~,~,awakeDataFileID] = GetFileInfo_IOS(procDataFileID);
            scoringLabels = [];
            for cc = 1:length(ScoringResults.fileIDs)
                if strcmp(awakeDataFileID,ScoringResults.fileIDs{cc,1}) == true
                    scoringLabels = ScoringResults.labels{cc,1};
                end
            end
            % check labels to match arousal state
            if sum(strcmp(scoringLabels,'Not Sleep')) > 144   % 36 bins (180 total) or 3 minutes of asleep
                load(procDataFileID,'-mat')
                % don't include trials with stimulation
                if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
                    try
                        puffs = ProcData.data.stimulations.LPadSol;
                    catch
                        puffs = ProcData.data.solenoids.LPadSol;
                    end
                    if isempty(puffs) == true
                        if sum(isnan(ProcData.data.Pupil.(dataType))) == 0
                            if strcmp(hemisphere,'LH_HbT') == true
                                HbT_awakeData{zz,1} = ProcData.data.CBV_HbT.adjLH;
                            elseif strcmp(hemisphere,'RH_HbT') == true
                                HbT_awakeData{zz,1} = ProcData.data.CBV_HbT.adjRH;
                            end
                            Pupil_awakeData{zz,1} = ProcData.data.Pupil.(dataType);
                            zz = zz + 1;
                        end
                    end
                end
            end
        end
        % filter and detrend data
        if isempty(HbT_awakeData) == false
            awakeProcData = [];
            for dd = 1:length(HbT_awakeData)
                awakeProcData.HbT{dd,1} = filtfilt(sos,g,detrend(HbT_awakeData{dd,1},'constant'));
                awakeProcData.Pupil{dd,1} = filtfilt(sos,g,detrend(Pupil_awakeData{dd,1},'constant'));
            end
            % set parameters for cross-correlation analysis
            awakeLagTime = 30;   % seconds
            awakeFrequency = oneSecSpecFs;   % Hz
            awakeMaxLag = awakeLagTime*awakeFrequency;
            % run cross-correlation analysis - average through time
            for dd = 1:length(awakeProcData.HbT)
                awakeHbTarray = awakeProcData.HbT{dd,1};
                awakePupilarray = awakeProcData.Pupil{dd,1};
                [awakeXcVals(dd,:),awakePupil_lags] = xcorr(awakeHbTarray,awakePupilarray,awakeMaxLag,'coeff'); %#ok<*AGROW>
            end
            awakeMeanXcVals = mean(awakeXcVals,1);
            % save results
            Results_CrossCorrelation.(animalID).Alert.(hemisphere).(dataType).lags = awakePupil_lags;
            Results_CrossCorrelation.(animalID).Alert.(hemisphere).(dataType).xcVals = awakeMeanXcVals;
        else
            % save results
            Results_CrossCorrelation.(animalID).Alert.(hemisphere).(dataType).lags = [];
            Results_CrossCorrelation.(animalID).Alert.(hemisphere).(dataType).xcVals = [];
        end
        %% analyze neural-hemo coherence during periods of asleep
        zz = 1;
        clear HbT_asleepData Pupil_asleepData
        HbT_asleepData = []; Pupil_asleepData = [];
        for dd = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(dd,:);
            [~,~,asleepDataFileID] = GetFileInfo_IOS(procDataFileID);
            scoringLabels = [];
            for cc = 1:length(ScoringResults.fileIDs)
                if strcmp(asleepDataFileID,ScoringResults.fileIDs{cc,1}) == true
                    scoringLabels = ScoringResults.labels{cc,1};
                end
            end
            % check labels to match arousal state
            if sum(strcmp(scoringLabels,'Not Sleep')) < 36   % 36 bins (180 total) or 3 minutes of asleep
                load(procDataFileID,'-mat')
                % don't include trials with stimulation
                if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
                    try
                        puffs = ProcData.data.stimulations.LPadSol;
                    catch
                        puffs = ProcData.data.solenoids.LPadSol;
                    end
                    if isempty(puffs) == true
                        if sum(isnan(ProcData.data.Pupil.(dataType))) == 0
                            if strcmp(hemisphere,'LH_HbT') == true
                                HbT_asleepData{zz,1} = ProcData.data.CBV_HbT.adjLH;
                            elseif strcmp(hemisphere,'RH_HbT') == true
                                HbT_asleepData{zz,1} = ProcData.data.CBV_HbT.adjRH;
                            end
                            Pupil_asleepData{zz,1} = ProcData.data.Pupil.(dataType);
                            zz = zz + 1;
                        end
                    end
                end
            end
        end
        % filter and detrend data
        if isempty(HbT_asleepData) == false
            asleepProcData = [];
            for dd = 1:length(HbT_asleepData)
                asleepProcData.HbT{dd,1} = filtfilt(sos,g,detrend(HbT_asleepData{dd,1},'constant'));
                asleepProcData.Pupil{dd,1} = filtfilt(sos,g,detrend(Pupil_asleepData{dd,1},'constant'));
            end
            % set parameters for cross-correlation analysis
            asleepLagTime = 30;   % seconds
            asleepFrequency = oneSecSpecFs;   % Hz
            asleepMaxLag = asleepLagTime*asleepFrequency;
            % run cross-correlation analysis - average through time
            for dd = 1:length(asleepProcData.HbT)
                asleepHbTarray = asleepProcData.HbT{dd,1};
                asleepPupilarray = asleepProcData.Pupil{dd,1};
                [asleepXcVals(dd,:),asleepPupil_lags] = xcorr(asleepHbTarray,asleepPupilarray,asleepMaxLag,'coeff'); %#ok<*AGROW>
            end
            asleepMeanXcVals = mean(asleepXcVals,1);
            % save results
            Results_CrossCorrelation.(animalID).Asleep.(hemisphere).(dataType).lags = asleepPupil_lags;
            Results_CrossCorrelation.(animalID).Asleep.(hemisphere).(dataType).xcVals = asleepMeanXcVals;
        else
            % save results
            Results_CrossCorrelation.(animalID).Asleep.(hemisphere).(dataType).lags = [];
            Results_CrossCorrelation.(animalID).Asleep.(hemisphere).(dataType).xcVals = [];
        end
        %% analyze neural-hemo coherence during periods of all data
        zz = 1;
        clear HbT_allData Pupil_allData
        HbT_allData = []; Pupil_allData = [];
        for dd = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(dd,:);
            [~,~,~] = GetFileInfo_IOS(procDataFileID);
            load(procDataFileID,'-mat')
            if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
                try
                    puffs = ProcData.data.stimulations.LPadSol;
                catch
                    puffs = ProcData.data.solenoids.LPadSol;
                end
                % don't include trials with stimulation
                if isempty(puffs) == true
                    if sum(isnan(ProcData.data.Pupil.(dataType))) == 0
                        if strcmp(hemisphere,'LH_HbT') == true
                            HbT_allData{zz,1} = ProcData.data.CBV_HbT.adjLH;
                        elseif strcmp(hemisphere,'RH_HbT') == true
                            HbT_allData{zz,1} = ProcData.data.CBV_HbT.adjRH;
                        end
                        Pupil_allData{zz,1} = ProcData.data.Pupil.(dataType);
                        zz = zz + 1;
                    end
                end
            end
        end
        % filter and detrend data
        if isempty(HbT_allData) == false
            allProcData = [];
            for dd = 1:length(HbT_allData)
                allProcData.HbT{dd,1} = filtfilt(sos,g,detrend(HbT_allData{dd,1},'constant'));
                allProcData.Pupil{dd,1} = filtfilt(sos,g,detrend(Pupil_allData{dd,1},'constant'));
            end
            % set parameters for cross-correlation analysis
            allLagTime = 30;   % seconds
            allFrequency = oneSecSpecFs;   % Hz
            allMaxLag = allLagTime*allFrequency;
            % run cross-correlation analysis - average through time
            for dd = 1:length(allProcData.HbT)
                allHbTarray = allProcData.HbT{dd,1};
                allPupilarray = allProcData.Pupil{dd,1};
                [allXcVals(dd,:),allPupil_lags] = xcorr(allHbTarray,allPupilarray,allMaxLag,'coeff'); %#ok<*AGROW>
            end
            allMeanXcVals = mean(allXcVals,1);
            % save results
            Results_CrossCorrelation.(animalID).All.(hemisphere).(dataType).lags = allPupil_lags;
            Results_CrossCorrelation.(animalID).All.(hemisphere).(dataType).xcVals = allMeanXcVals;
        else
            % save results
            Results_CrossCorrelation.(animalID).All.(hemisphere).(dataType).lags = [];
            Results_CrossCorrelation.(animalID).All.(hemisphere).(dataType).xcVals = [];
        end
        %% cross-correlation analysis for NREM
        if isempty(SleepData.(modelType).NREM.data.Pupil) == false
            NREM_sleepTime = params.minTime.NREM;   % seconds
            [NREM_finalHbT,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.Pupil.(hemisphere).data,SleepData.(modelType).NREM.data.Pupil.fileIDs,SleepData.(modelType).NREM.data.Pupil.binTimes);
            [NREM_finalPupil,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.Pupil.(dataType).data,SleepData.(modelType).NREM.data.Pupil.fileIDs,SleepData.(modelType).NREM.data.Pupil.binTimes);
            NREM_finalHbTVals = []; NREM_finalPupilVals = [];
            if isempty(NREM_finalHbT) == false
                % adjust[HbT] and Pupil events to match the edits made to the length of each spectrogram
                mm = 1;
                for ll = 1:length(NREM_finalHbT)
                    if sum(isnan(NREM_finalPupil{ll,1})) == 0
                        NREM_HbTVals = NREM_finalHbT{ll,1}(1:NREM_sleepTime*samplingRate);
                        NREM_PupilVals = NREM_finalPupil{ll,1}(1:NREM_sleepTime*samplingRate);
                        NREM_finalHbTVals{mm,1} = filtfilt(sos,g,detrend(NREM_HbTVals,'constant'));
                        NREM_finalPupilVals{mm,1} = filtfilt(sos,g,detrend(NREM_PupilVals,'constant'));
                        mm = mm + 1;
                    end
                end
                if isempty(NREM_finalHbTVals) == false
                    % run cross-correlation analysis - average through time
                    NREM_lagTime = 10;   % Seconds
                    NREM_frequency = oneSecSpecFs;   % Hz
                    NREM_maxLag = NREM_lagTime*NREM_frequency;
                    for nn = 1:length(NREM_finalHbTVals)
                        NREM_HbT_array = NREM_finalHbTVals{nn,1};
                        NREM_Pupil_array = NREM_finalPupilVals{nn,1};
                        [NREM_xcVals(nn,:),NREM_Pupil_lags] = xcorr(NREM_HbT_array,NREM_Pupil_array,NREM_maxLag,'coeff');
                    end
                    NREM_meanXcVals = mean(NREM_xcVals,1);
                    % save results
                    Results_CrossCorrelation.(animalID).NREM.(hemisphere).(dataType).lags = NREM_Pupil_lags;
                    Results_CrossCorrelation.(animalID).NREM.(hemisphere).(dataType).xcVals = NREM_meanXcVals;
                else
                    % save results
                    Results_CrossCorrelation.(animalID).NREM.(hemisphere).(dataType).lags = [];
                    Results_CrossCorrelation.(animalID).NREM.(hemisphere).(dataType).xcVals = [];
                end
            end
        else
            % save results
            Results_CrossCorrelation.(animalID).NREM.(hemisphere).(dataType).lags = [];
            Results_CrossCorrelation.(animalID).NREM.(hemisphere).(dataType).xcVals = [];
        end
        %% cross-correlation analysis for REM
        if isempty(SleepData.(modelType).REM.data.Pupil) == false
            REM_sleepTime = params.minTime.REM;   % seconds
            [REM_finalHbT,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.Pupil.(hemisphere).data,SleepData.(modelType).REM.data.Pupil.fileIDs,SleepData.(modelType).REM.data.Pupil.binTimes);
            [REM_finalPupil,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.Pupil.(dataType).data,SleepData.(modelType).REM.data.Pupil.fileIDs,SleepData.(modelType).REM.data.Pupil.binTimes);
            REM_finalHbTVals = []; REM_finalPupilVals = [];
            if isempty(REM_finalHbT) == false
                % adjust[HbT] and Pupil events to match the edits made to the length of each spectrogram
                mm = 1;
                for ll = 1:length(REM_finalHbT)
                    if sum(isnan(REM_finalPupil{ll,1})) == 0
                        REM_HbTVals = REM_finalHbT{ll,1}(1:REM_sleepTime*samplingRate);
                        REM_PupilVals = REM_finalPupil{ll,1}(1:REM_sleepTime*samplingRate);
                        REM_finalHbTVals{mm,1} = filtfilt(sos,g,detrend(REM_HbTVals,'constant'));
                        REM_finalPupilVals{mm,1} = filtfilt(sos,g,detrend(REM_PupilVals,'constant'));
                        mm = mm + 1;
                    end
                end
                if isempty(REM_finalHbTVals) == false
                    % run cross-correlation analysis - average through time
                    REM_lagTime = 10;   % Seconds
                    REM_frequency = oneSecSpecFs;   % Hz
                    REM_maxLag = REM_lagTime*REM_frequency;
                    for nn = 1:length(REM_finalHbTVals)
                        REM_HbT_array = REM_finalHbTVals{nn,1};
                        REM_Pupil_array = REM_finalPupilVals{nn,1};
                        [REM_xcVals(nn,:),REM_Pupil_lags] = xcorr(REM_HbT_array,REM_Pupil_array,REM_maxLag,'coeff');
                    end
                    REM_meanXcVals = mean(REM_xcVals,1);
                    % save results
                    Results_CrossCorrelation.(animalID).REM.(hemisphere).(dataType).lags = REM_Pupil_lags;
                    Results_CrossCorrelation.(animalID).REM.(hemisphere).(dataType).xcVals = REM_meanXcVals;
                else
                    % save results
                    Results_CrossCorrelation.(animalID).REM.(hemisphere).(dataType).lags = [];
                    Results_CrossCorrelation.(animalID).REM.(hemisphere).(dataType).xcVals = [];
                end
            end
        else
            % save results
            Results_CrossCorrelation.(animalID).REM.(hemisphere).(dataType).lags = [];
            Results_CrossCorrelation.(animalID).REM.(hemisphere).(dataType).xcVals = [];
        end
    end
end
% save data
cd([rootFolder delim])
save('Results_CrossCorrelation.mat','Results_CrossCorrelation')

end
