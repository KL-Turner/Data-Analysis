function [Results_NeuralHemoCoherGCaMP] = AnalyzeNeuralHemoCoherence_GCaMP(animalID,group,rootFolder,delim,Results_NeuralHemoCoherGCaMP)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Analyze the spectral coherence between neural-hemodynamic [HbT] signals (IOS)
%________________________________________________________________________________________________________________________

hemDataTypes = {'LH','RH'};
modelType = 'Forest';
params.minTime.Rest = 10;
params.minTime.NREM = 30;
params.minTime.REM = 60;
% only run analysis for valid animal IDs
dataLocation = [rootFolder delim 'Data' delim group delim animalID delim 'Bilateral Imaging'];
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
% go through each valid data type for arousal based coherence analysis
for zzz = 1:length(hemDataTypes)
    hemDataType = hemDataTypes{1,zzz};
    %% analyze neural-hemo coherence during periods of rest
    % pull data from RestData.mat structure
    samplingRate = RestData.CBV_HbT.LH.CBVCamSamplingRate;
    [restLogical] = FilterEvents_IOS(RestData.CBV_HbT.(hemDataType),RestCriteria);
    [puffLogical] = FilterEvents_IOS(RestData.CBV_HbT.(hemDataType),RestPuffCriteria);
    combRestLogical = logical(restLogical.*puffLogical);
    restFileIDs = RestData.CBV_HbT.(hemDataType).fileIDs(combRestLogical,:);
    restEventTimes = RestData.CBV_HbT.(hemDataType).eventTimes(combRestLogical,:);
    restDurations = RestData.CBV_HbT.(hemDataType).durations(combRestLogical,:);
    dataType = 'GCaMP7s';
    HbT_RestingData = RestData.CBV_HbT.(hemDataType).data(combRestLogical,:);
    Neural_RestingData = RestData.GCaMP7s.(['cor' hemDataType]).NormData(combRestLogical,:);
    fHbT_RestingData = RestData.CBV_HbT.(['frontal' hemDataType]).data(combRestLogical,:);
    fNeural_RestingData = RestData.GCaMP7s.(['corFrontal' hemDataType]).NormData(combRestLogical,:);
    % keep only the data that occurs within the manually-approved alert regions
    [HbT_finalRestData,~,~,~] = RemoveInvalidData_IOS(HbT_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    [Neural_finalRestData,~,~,~] = RemoveInvalidData_IOS(Neural_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    [fHbT_finalRestData,~,~,~] = RemoveInvalidData_IOS(fHbT_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    [fNeural_finalRestData,~,~,~] = RemoveInvalidData_IOS(fNeural_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    clear HbT_ProcRestData Neural_ProcRestData
    clear fHbT_ProcRestData fNeural_ProcRestData
    % filter, detrend, and truncate data to minimum length to match events
    for bb = 1:length(HbT_finalRestData)
        if length(HbT_finalRestData{bb,1}) < params.minTime.Rest*samplingRate
            restChunkSampleDiff = params.minTime.Rest*samplingRate - length(HbT_finalRestData{bb,1});
            HbT_restPad = (ones(1,restChunkSampleDiff))*HbT_finalRestData{bb,1}(end);
            Neural_restPad = (ones(1,restChunkSampleDiff))*Neural_finalRestData{bb,1}(end);
            HbT_ProcRestData{bb,1} = horzcat(HbT_finalRestData{bb,1},HbT_restPad);
            Neural_ProcRestData{bb,1} = horzcat(Neural_finalRestData{bb,1},Neural_restPad);
            HbT_ProcRestData{bb,1} = detrend(HbT_ProcRestData{bb,1},'constant');
            Neural_ProcRestData{bb,1} = detrend(Neural_ProcRestData{bb,1},'constant');
            fHbT_restPad = (ones(1,restChunkSampleDiff))*fHbT_finalRestData{bb,1}(end);
            fNeural_restPad = (ones(1,restChunkSampleDiff))*fNeural_finalRestData{bb,1}(end);
            fHbT_ProcRestData{bb,1} = horzcat(fHbT_finalRestData{bb,1},fHbT_restPad);
            fNeural_ProcRestData{bb,1} = horzcat(fNeural_finalRestData{bb,1},fNeural_restPad);
            fHbT_ProcRestData{bb,1} = detrend(fHbT_ProcRestData{bb,1},'constant');
            fNeural_ProcRestData{bb,1} = detrend(fNeural_ProcRestData{bb,1},'constant');
        else
            HbT_ProcRestData{bb,1} = detrend(HbT_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
            Neural_ProcRestData{bb,1} = detrend(Neural_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
            fHbT_ProcRestData{bb,1} = detrend(fHbT_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
            fNeural_ProcRestData{bb,1} = detrend(fNeural_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
        end
    end
    % pre-allocate coherence matrix
    HbT_restData = zeros(length(HbT_ProcRestData{1,1}),length(HbT_ProcRestData));
    Neural_restData = zeros(length(Neural_ProcRestData{1,1}),length(Neural_ProcRestData));
    fHbT_restData = zeros(length(fHbT_ProcRestData{1,1}),length(fHbT_ProcRestData));
    fNeural_restData = zeros(length(fNeural_ProcRestData{1,1}),length(fNeural_ProcRestData));
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
    for cc = 1:length(HbT_ProcRestData)
        HbT_restData(:,cc) = HbT_ProcRestData{cc,1};
        Neural_restData(:,cc) = Neural_ProcRestData{cc,1};
        fHbT_restData(:,cc) = fHbT_ProcRestData{cc,1};
        fNeural_restData(:,cc) = fNeural_ProcRestData{cc,1};
    end
    % parameters for coherencyc - information available in function
    params.tapers = [1,1]; % Tapers [n, 2n - 1]
    params.pad = 1;
    params.Fs = samplingRate;
    params.fpass = [0,1]; % Pass band [0, nyquist]
    params.trialave = 1;
    params.err = [2,0.05];
    % calculate the coherence between desired signals
    [C_RestData,~,~,~,~,f_RestData,confC_RestData,~,cErr_RestData] = coherencyc(HbT_restData,Neural_restData,params);
    % save results
    Results_NeuralHemoCoherGCaMP.(animalID).Rest.(dataType).(hemDataType).C = C_RestData;
    Results_NeuralHemoCoherGCaMP.(animalID).Rest.(dataType).(hemDataType).f = f_RestData;
    Results_NeuralHemoCoherGCaMP.(animalID).Rest.(dataType).(hemDataType).confC = confC_RestData;
    Results_NeuralHemoCoherGCaMP.(animalID).Rest.(dataType).(hemDataType).cErr = cErr_RestData;
    % calculate the coherence between desired signals
    [fC_RestData,~,~,~,~,ff_RestData,fconfC_RestData,~,fcErr_RestData] = coherencyc(fHbT_restData,fNeural_restData,params);
    % save results
    Results_NeuralHemoCoherGCaMP.(animalID).Rest.(dataType).(hemDataType).fC = fC_RestData;
    Results_NeuralHemoCoherGCaMP.(animalID).Rest.(dataType).(hemDataType).ff = ff_RestData;
    Results_NeuralHemoCoherGCaMP.(animalID).Rest.(dataType).(hemDataType).fconfC = fconfC_RestData;
    Results_NeuralHemoCoherGCaMP.(animalID).Rest.(dataType).(hemDataType).fcErr = fcErr_RestData;
    %% analyze neural-hemo coherence during periods of alert
    zz = 1;
    clear HbT_AlertData Neural_AlertData HbT_ProcAlertData Neural_ProcAlertData
    clear fHbT_AlertData fNeural_AlertData fHbT_ProcAlertData fNeural_ProcAlertData
    HbT_AlertData = [];
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        scoringLabels = [];
        for cc = 1:length(ScoringResults.fileIDs)
            if strcmp(procDataFileID,ScoringResults.fileIDs{cc,1}) == true
                scoringLabels = ScoringResults.labels{cc,1};
            end
        end
        % check labels to match arousal state
        if sum(strcmp(scoringLabels,'Not Sleep')) > 144 % 36 bins (180 total) or 3 minutes of sleep
            load(procDataFileID,'-mat')
            puffs = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(puffs) == true
                HbT_AlertData{zz,1} = ProcData.data.CBV_HbT.(hemDataType);
                Neural_AlertData{zz,1} = ProcData.data.(dataType).(['cor' hemDataType]);
                fHbT_AlertData{zz,1} = ProcData.data.CBV_HbT.(['frontal' hemDataType]);
                fNeural_AlertData{zz,1} = ProcData.data.(dataType).(['corFrontal' hemDataType]);
                zz = zz + 1;
            end
        end
    end
    % filter and detrend data
    if isempty(HbT_AlertData) == false
        for bb = 1:length(HbT_AlertData)
            HbT_ProcAlertData{bb,1} = detrend(HbT_AlertData{bb,1},'constant');
            Neural_ProcAlertData{bb,1} = detrend(Neural_AlertData{bb,1},'constant');
            fHbT_ProcAlertData{bb,1} = detrend(fHbT_AlertData{bb,1},'constant');
            fNeural_ProcAlertData{bb,1} = detrend(fNeural_AlertData{bb,1},'constant');
        end
        % pre-allocate coherence matrix
        HbT_alertData = zeros(length(HbT_ProcAlertData{1,1}),length(HbT_ProcAlertData));
        Neural_alertData = zeros(length(Neural_ProcAlertData{1,1}),length(Neural_ProcAlertData));
        fHbT_alertData = zeros(length(fHbT_ProcAlertData{1,1}),length(fHbT_ProcAlertData));
        fNeural_alertData = zeros(length(fNeural_ProcAlertData{1,1}),length(fNeural_ProcAlertData));
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
        for cc = 1:length(HbT_ProcAlertData)
            HbT_alertData(:,cc) = HbT_ProcAlertData{cc,1};
            Neural_alertData(:,cc) = Neural_ProcAlertData{cc,1};
            fHbT_alertData(:,cc) = fHbT_ProcAlertData{cc,1};
            fNeural_alertData(:,cc) = fNeural_ProcAlertData{cc,1};
        end
        % calculate the coherence between desired signals
        params.tapers = [10,19]; % Tapers [n, 2n - 1]
        [C_AlertData,~,~,~,~,f_AlertData,confC_AlertData,~,cErr_AlertData] = coherencyc(HbT_alertData,Neural_alertData,params);
        % save results
        Results_NeuralHemoCoherGCaMP.(animalID).Alert.(dataType).(hemDataType).C = C_AlertData;
        Results_NeuralHemoCoherGCaMP.(animalID).Alert.(dataType).(hemDataType).f = f_AlertData;
        Results_NeuralHemoCoherGCaMP.(animalID).Alert.(dataType).(hemDataType).confC = confC_AlertData;
        Results_NeuralHemoCoherGCaMP.(animalID).Alert.(dataType).(hemDataType).cErr = cErr_AlertData;
        % calculate the coherence between desired signals
        [fC_AlertData,~,~,~,~,ff_AlertData,fconfC_AlertData,~,fcErr_AlertData] = coherencyc(fHbT_alertData,fNeural_alertData,params);
        % save results
        Results_NeuralHemoCoherGCaMP.(animalID).Alert.(dataType).(hemDataType).fC = fC_AlertData;
        Results_NeuralHemoCoherGCaMP.(animalID).Alert.(dataType).(hemDataType).ff = ff_AlertData;
        Results_NeuralHemoCoherGCaMP.(animalID).Alert.(dataType).(hemDataType).fconfC = fconfC_AlertData;
        Results_NeuralHemoCoherGCaMP.(animalID).Alert.(dataType).(hemDataType).fcErr = fcErr_AlertData;
    else
        % save results
        Results_NeuralHemoCoherGCaMP.(animalID).Alert.(dataType).(hemDataType).C = [];
        Results_NeuralHemoCoherGCaMP.(animalID).Alert.(dataType).(hemDataType).f = [];
        Results_NeuralHemoCoherGCaMP.(animalID).Alert.(dataType).(hemDataType).confC = [];
        Results_NeuralHemoCoherGCaMP.(animalID).Alert.(dataType).(hemDataType).cErr = [];
        % save results
        Results_NeuralHemoCoherGCaMP.(animalID).Alert.(dataType).(hemDataType).fC = [];
        Results_NeuralHemoCoherGCaMP.(animalID).Alert.(dataType).(hemDataType).ff = [];
        Results_NeuralHemoCoherGCaMP.(animalID).Alert.(dataType).(hemDataType).fconfC = [];
        Results_NeuralHemoCoherGCaMP.(animalID).Alert.(dataType).(hemDataType).fcErr = [];
    end
    %% analyze neural-hemo coherence during periods of asleep
    zz = 1;
    clear HbT_AsleepData Neural_AsleepData HbT_ProcAsleepData Neural_ProcAsleepData
    clear fHbT_AsleepData fNeural_AsleepData fHbT_ProcAsleepData fNeural_ProcAsleepData
    HbT_AsleepData = [];
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        scoringLabels = [];
        for cc = 1:length(ScoringResults.fileIDs)
            if strcmp(procDataFileID,ScoringResults.fileIDs{cc,1}) == true
                scoringLabels = ScoringResults.labels{cc,1};
            end
        end
        % check labels to match arousal state
        if sum(strcmp(scoringLabels,'Not Sleep')) < 36 % 36 bins (180 total) or 3 minutes of alert
            load(procDataFileID,'-mat')
            puffs = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(puffs) == true
                HbT_AsleepData{zz,1} = ProcData.data.CBV_HbT.(hemDataType);
                Neural_AsleepData{zz,1} = ProcData.data.(dataType).(['cor' hemDataType]);
                fHbT_AsleepData{zz,1} = ProcData.data.CBV_HbT.(['frontal' hemDataType]);
                fNeural_AsleepData{zz,1} = ProcData.data.(dataType).(['corFrontal' hemDataType]);
                zz = zz + 1;
            end
        end
    end
    % filter and detrend data
    if isempty(HbT_AsleepData) == false
        for bb = 1:length(HbT_AsleepData)
            HbT_ProcAsleepData{bb,1} = detrend(HbT_AsleepData{bb,1},'constant');
            Neural_ProcAsleepData{bb,1} = detrend(Neural_AsleepData{bb,1},'constant');
            fHbT_ProcAsleepData{bb,1} = detrend(fHbT_AsleepData{bb,1},'constant');
            fNeural_ProcAsleepData{bb,1} = detrend(fNeural_AsleepData{bb,1},'constant');
        end
        % pre-allocate coherence matrix
        HbT_asleepData = zeros(length(HbT_ProcAsleepData{1,1}),length(HbT_ProcAsleepData));
        Neural_asleepData = zeros(length(Neural_ProcAsleepData{1,1}),length(Neural_ProcAsleepData));
        fHbT_asleepData = zeros(length(fHbT_ProcAsleepData{1,1}),length(fHbT_ProcAsleepData));
        fNeural_asleepData = zeros(length(fNeural_ProcAsleepData{1,1}),length(fNeural_ProcAsleepData));
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
        for cc = 1:length(HbT_ProcAsleepData)
            HbT_asleepData(:,cc) = HbT_ProcAsleepData{cc,1};
            Neural_asleepData(:,cc) = Neural_ProcAsleepData{cc,1};
            fHbT_asleepData(:,cc) = fHbT_ProcAsleepData{cc,1};
            fNeural_asleepData(:,cc) = fNeural_ProcAsleepData{cc,1};
        end
        % calculate the coherence between desired signals
        params.tapers = [10,19]; % Tapers [n, 2n - 1]
        [C_AsleepData,~,~,~,~,f_AsleepData,confC_AsleepData,~,cErr_AsleepData] = coherencyc(HbT_asleepData,Neural_asleepData,params);
        % save results
        Results_NeuralHemoCoherGCaMP.(animalID).Asleep.(dataType).(hemDataType).C = C_AsleepData;
        Results_NeuralHemoCoherGCaMP.(animalID).Asleep.(dataType).(hemDataType).f = f_AsleepData;
        Results_NeuralHemoCoherGCaMP.(animalID).Asleep.(dataType).(hemDataType).confC = confC_AsleepData;
        Results_NeuralHemoCoherGCaMP.(animalID).Asleep.(dataType).(hemDataType).cErr = cErr_AsleepData;
        % calculate the coherence between desired signals
        [fC_AsleepData,~,~,~,~,ff_AsleepData,fconfC_AsleepData,~,fcErr_AsleepData] = coherencyc(fHbT_asleepData,fNeural_asleepData,params);
        % save results
        Results_NeuralHemoCoherGCaMP.(animalID).Asleep.(dataType).(hemDataType).fC = fC_AsleepData;
        Results_NeuralHemoCoherGCaMP.(animalID).Asleep.(dataType).(hemDataType).ff = ff_AsleepData;
        Results_NeuralHemoCoherGCaMP.(animalID).Asleep.(dataType).(hemDataType).fconfC = fconfC_AsleepData;
        Results_NeuralHemoCoherGCaMP.(animalID).Asleep.(dataType).(hemDataType).fcErr = fcErr_AsleepData;
    else
        % save results
        Results_NeuralHemoCoherGCaMP.(animalID).Asleep.(dataType).(hemDataType).C = [];
        Results_NeuralHemoCoherGCaMP.(animalID).Asleep.(dataType).(hemDataType).f = [];
        Results_NeuralHemoCoherGCaMP.(animalID).Asleep.(dataType).(hemDataType).confC = [];
        Results_NeuralHemoCoherGCaMP.(animalID).Asleep.(dataType).(hemDataType).cErr = [];
        % save results
        Results_NeuralHemoCoherGCaMP.(animalID).Asleep.(dataType).(hemDataType).fC = [];
        Results_NeuralHemoCoherGCaMP.(animalID).Asleep.(dataType).(hemDataType).ff = [];
        Results_NeuralHemoCoherGCaMP.(animalID).Asleep.(dataType).(hemDataType).fconfC = [];
        Results_NeuralHemoCoherGCaMP.(animalID).Asleep.(dataType).(hemDataType).fcErr = [];
    end
    %% analyze neural-hemo coherence during periods of all data
    zz = 1;
    clear HbT_AllData Neural_AllData HbT_ProcAllData Neural_ProcAllData
    clear fHbT_AllData fNeural_AllData fHbT_ProcAllData fNeural_ProcAllData
    HbT_AllData = [];
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        load(procDataFileID,'-mat')
        puffs = ProcData.data.stimulations.LPadSol;
        % don't include trials with stimulation
        if isempty(puffs) == true
            HbT_AllData{zz,1} = ProcData.data.CBV_HbT.(hemDataType);
            Neural_AllData{zz,1} = ProcData.data.(dataType).(['cor' hemDataType]);
            fHbT_AllData{zz,1} = ProcData.data.CBV_HbT.(['frontal' hemDataType]);
            fNeural_AllData{zz,1} = ProcData.data.(dataType).(['corFrontal' hemDataType]);
            zz = zz + 1;
        end
    end
    % filter and detrend data
    if isempty(HbT_AllData) == false
        for bb = 1:length(HbT_AllData)
            HbT_ProcAllUnstimData{bb,1} = detrend(HbT_AllData{bb,1},'constant');
            Neural_ProcAllUnstimData{bb,1} = detrend(Neural_AllData{bb,1},'constant');
            fHbT_ProcAllUnstimData{bb,1} = detrend(fHbT_AllData{bb,1},'constant');
            fNeural_ProcAllUnstimData{bb,1} = detrend(fNeural_AllData{bb,1},'constant');
        end
        % pre-allocate coherence matrix
        HbT_allData = zeros(length(HbT_ProcAllUnstimData{1,1}),length(HbT_ProcAllUnstimData));
        Neural_allData = zeros(length(Neural_ProcAllUnstimData{1,1}),length(Neural_ProcAllUnstimData));
        fHbT_allData = zeros(length(fHbT_ProcAllUnstimData{1,1}),length(fHbT_ProcAllUnstimData));
        fNeural_allData = zeros(length(fNeural_ProcAllUnstimData{1,1}),length(fNeural_ProcAllUnstimData));
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
        for cc = 1:length(HbT_ProcAllUnstimData)
            HbT_allData(:,cc) = HbT_ProcAllUnstimData{cc,1};
            Neural_allData(:,cc) = Neural_ProcAllUnstimData{cc,1};
            fHbT_allData(:,cc) = fHbT_ProcAllUnstimData{cc,1};
            fNeural_allData(:,cc) = fNeural_ProcAllUnstimData{cc,1};
        end
        % calculate the coherence between desired signals
        params.tapers = [10,19]; % Tapers [n, 2n - 1]
        [C_AllUnstimData,~,~,~,~,f_AllUnstimData,confC_AllUnstimData,~,cErr_AllUnstimData] = coherencyc(HbT_allData,Neural_allData,params);
        % save results
        Results_NeuralHemoCoherGCaMP.(animalID).All.(dataType).(hemDataType).C = C_AllUnstimData;
        Results_NeuralHemoCoherGCaMP.(animalID).All.(dataType).(hemDataType).f = f_AllUnstimData;
        Results_NeuralHemoCoherGCaMP.(animalID).All.(dataType).(hemDataType).confC = confC_AllUnstimData;
        Results_NeuralHemoCoherGCaMP.(animalID).All.(dataType).(hemDataType).cErr = cErr_AllUnstimData;
        % calculate the coherence between desired signals
        [fC_AllUnstimData,~,~,~,~,ff_AllUnstimData,fconfC_AllUnstimData,~,fcErr_AllUnstimData] = coherencyc(fHbT_allData,fNeural_allData,params);
        % save results
        Results_NeuralHemoCoherGCaMP.(animalID).All.(dataType).(hemDataType).fC = fC_AllUnstimData;
        Results_NeuralHemoCoherGCaMP.(animalID).All.(dataType).(hemDataType).ff = ff_AllUnstimData;
        Results_NeuralHemoCoherGCaMP.(animalID).All.(dataType).(hemDataType).fconfC = fconfC_AllUnstimData;
        Results_NeuralHemoCoherGCaMP.(animalID).All.(dataType).(hemDataType).fcErr = fcErr_AllUnstimData;
    end
    %% analyze neural-hemo coherence during periods of NREM
    % pull data from SleepData.mat structure
    [HbT_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.CBV_HbT.(hemDataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    [Neural_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GCaMP7s.(hemDataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    [fHbT_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.CBV_HbT.(['frontal' hemDataType]),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    [fNeural_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GCaMP7s.(['frontal' hemDataType]),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    % filter, detrend, and truncate data to minimum length to match events
    for ee = 1:length(HbT_nremData)
        HbT_nremData{ee,1} = detrend(HbT_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        Neural_nremData{ee,1} = detrend(Neural_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        fHbT_nremData{ee,1} = detrend(fHbT_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        fNeural_nremData{ee,1} = detrend(fNeural_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
    end
    % pre-allocate coherence matrix
    HbT_nrem = zeros(length(HbT_nremData{1,1}),length(HbT_nremData));
    Neural_nrem = zeros(length(Neural_nremData{1,1}),length(Neural_nremData));
    fHbT_nrem = zeros(length(fHbT_nremData{1,1}),length(fHbT_nremData));
    fNeural_nrem = zeros(length(fNeural_nremData{1,1}),length(fNeural_nremData));
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
    for ff = 1:length(HbT_nremData)
        HbT_nrem(:,ff) = HbT_nremData{ff,1};
        Neural_nrem(:,ff) = Neural_nremData{ff,1};
        fHbT_nrem(:,ff) = fHbT_nremData{ff,1};
        fNeural_nrem(:,ff) = fNeural_nremData{ff,1};
    end
    % calculate the coherence between desired signals
    params.tapers = [3,5]; % Tapers [n, 2n - 1]
    [C_nrem,~,~,~,~,f_nrem,confC_nrem,~,cErr_nrem] = coherencyc(HbT_nrem,Neural_nrem,params);
    % save results
    Results_NeuralHemoCoherGCaMP.(animalID).NREM.(dataType).(hemDataType).C = C_nrem;
    Results_NeuralHemoCoherGCaMP.(animalID).NREM.(dataType).(hemDataType).f = f_nrem;
    Results_NeuralHemoCoherGCaMP.(animalID).NREM.(dataType).(hemDataType).confC = confC_nrem;
    Results_NeuralHemoCoherGCaMP.(animalID).NREM.(dataType).(hemDataType).cErr = cErr_nrem;
    % calculate the coherence between desired signals
    [fC_nrem,~,~,~,~,ff_nrem,fconfC_nrem,~,fcErr_nrem] = coherencyc(fHbT_nrem,fNeural_nrem,params);
    % save results
    Results_NeuralHemoCoherGCaMP.(animalID).NREM.(dataType).(hemDataType).fC = fC_nrem;
    Results_NeuralHemoCoherGCaMP.(animalID).NREM.(dataType).(hemDataType).ff = ff_nrem;
    Results_NeuralHemoCoherGCaMP.(animalID).NREM.(dataType).(hemDataType).fconfC = fconfC_nrem;
    Results_NeuralHemoCoherGCaMP.(animalID).NREM.(dataType).(hemDataType).fcErr = fcErr_nrem;
    %% analyze neural-hemo coherence during periods of REM
    % pull data from SleepData.mat structure
    [HbT_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.CBV_HbT.(hemDataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    [Neural_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GCaMP7s.(hemDataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    [fHbT_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.CBV_HbT.(['frontal' hemDataType]),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    [fNeural_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GCaMP7s.(['frontal' hemDataType]),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    % filter, detrend, and truncate data to minimum length to match events
    for gg = 1:length(HbT_remData)
        HbT_remData{gg,1} = detrend(HbT_remData{gg,1}(1:(params.minTime.REM*samplingRate)),'constant');
        Neural_remData{gg,1} = detrend(Neural_remData{gg,1}(1:(params.minTime.REM*samplingRate)),'constant');
        fHbT_remData{gg,1} = detrend(fHbT_remData{gg,1}(1:(params.minTime.REM*samplingRate)),'constant');
        fNeural_remData{gg,1} = detrend(fNeural_remData{gg,1}(1:(params.minTime.REM*samplingRate)),'constant');
    end
    % pre-allocate coherence matrix
    HbT_rem = zeros(length(HbT_remData{1,1}),length(HbT_remData));
    Neural_rem = zeros(length(Neural_remData{1,1}),length(Neural_remData));
    fHbT_rem = zeros(length(fHbT_remData{1,1}),length(fHbT_remData));
    fNeural_rem = zeros(length(fNeural_remData{1,1}),length(fNeural_remData));
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
    for hh = 1:length(HbT_remData)
        HbT_rem(:,hh) = HbT_remData{hh,1};
        Neural_rem(:,hh) = Neural_remData{hh,1};
        fHbT_rem(:,hh) = fHbT_remData{hh,1};
        fNeural_rem(:,hh) = fNeural_remData{hh,1};
    end
    % calculate the coherence between desired signals
    params.tapers = [5,9]; % Tapers [n, 2n - 1]
    [C_rem,~,~,~,~,f_rem,confC_rem,~,cErr_rem] = coherencyc(HbT_rem,Neural_rem,params);
    % save results
    Results_NeuralHemoCoherGCaMP.(animalID).REM.(dataType).(hemDataType).C = C_rem;
    Results_NeuralHemoCoherGCaMP.(animalID).REM.(dataType).(hemDataType).f = f_rem;
    Results_NeuralHemoCoherGCaMP.(animalID).REM.(dataType).(hemDataType).confC = confC_rem;
    Results_NeuralHemoCoherGCaMP.(animalID).REM.(dataType).(hemDataType).cErr = cErr_rem;
    % calculate the coherence between desired signals
    [fC_rem,~,~,~,~,ff_rem,fconfC_rem,~,fcErr_rem] = coherencyc(fHbT_rem,fNeural_rem,params);
    % save results
    Results_NeuralHemoCoherGCaMP.(animalID).REM.(dataType).(hemDataType).fC = fC_rem;
    Results_NeuralHemoCoherGCaMP.(animalID).REM.(dataType).(hemDataType).ff = ff_rem;
    Results_NeuralHemoCoherGCaMP.(animalID).REM.(dataType).(hemDataType).fconfC = fconfC_rem;
    Results_NeuralHemoCoherGCaMP.(animalID).REM.(dataType).(hemDataType).fcErr = fcErr_rem;
end
% save data
cd([rootFolder delim 'Results_Turner'])
save('Results_NeuralHemoCoherGCaMP.mat','Results_NeuralHemoCoherGCaMP')
cd([rootFolder delim 'Data'])

end
