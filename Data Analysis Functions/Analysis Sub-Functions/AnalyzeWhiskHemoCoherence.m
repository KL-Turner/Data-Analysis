function [AnalysisResults] = AnalyzeWhiskHemoCoherence(animalID,group,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the spectral coherence between neural-hemodynamic [HbT] signals (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
hemDataTypes = {'adjLH','adjRH'};
%% only run analysis for valid animal IDs
dataLocation = [rootFolder '\' group '\' animalID '\Bilateral Imaging\'];
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
% lowpass filter
samplingRate = RestData.CBV_HbT.adjLH.CBVCamSamplingRate;
[z,p,k] = butter(4,1/(samplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
% go through each valid data type for arousal-based coherence analysis
for zzz = 1:length(hemDataTypes)
    hemDataType = hemDataTypes{1,zzz};
    %% analyze neural-hemo coherence during periods of alert
    zz = 1;
    clear HbT_AwakeData Gamma_AwakeData HbT_ProcAwakeData Gamma_ProcAwakeData
    HbT_AwakeData = [];
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,~,allDataFileID] = GetFileInfo_IOS(procDataFileID);
        scoringLabels = [];
        for cc = 1:length(ScoringResults.fileIDs)
            if strcmp(allDataFileID,ScoringResults.fileIDs{cc,1}) == true
                scoringLabels = ScoringResults.labels{cc,1};
            end
        end
        % check labels to match arousal state
        if sum(strcmp(scoringLabels,'Not Sleep')) > 144   % 36 bins (180 total) or 3 minutes of sleep
            load(procDataFileID,'-mat')
            puffs = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(puffs) == true
                HbT_AwakeData{zz,1} = ProcData.data.CBV_HbT.(hemDataType)(2:end-1);
                whisk_AwakeData{zz,1} = ProcData.data.binWhiskerAngle(1:26998);
                zz = zz + 1;
            end
        end
    end
    % filter and detrend data
    if isempty(HbT_AwakeData) == false
        for bb = 1:length(HbT_AwakeData)
            HbT_ProcAwakeData{bb,1} = filtfilt(sos,g,detrend(HbT_AwakeData{bb,1},'constant'));
            whisk_ProcAwakeData{bb,1} = detrend(whisk_AwakeData{bb,1},'constant');
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        HbT_awakeData = zeros(length(HbT_ProcAwakeData{1,1}),length(HbT_ProcAwakeData));
        whisk_awakeData = zeros(length(whisk_ProcAwakeData{1,1}),length(whisk_ProcAwakeData));
        for cc = 1:length(HbT_ProcAwakeData)
            HbT_awakeData(:,cc) = HbT_ProcAwakeData{cc,1};
            whisk_awakeData(:,cc) = whisk_ProcAwakeData{cc,1};
        end
        % parameters for coherencyc - information available in function
        params.tapers = [1,1];   % Tapers [n, 2n - 1]
        params.pad = 1;
        params.Fs = samplingRate;
        params.fpass = [0,15];   % Pass band [0, nyquist]
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the coherence between desired signals
        [C_AwakeData,~,~,~,~,f_AwakeData,confC_AwakeData,~,cErr_AwakeData] = coherencyc(HbT_awakeData,whisk_awakeData,params);
        % save results
        AnalysisResults.(animalID).WhiskHemoCoherence.Awake.(hemDataType).C = C_AwakeData;
        AnalysisResults.(animalID).WhiskHemoCoherence.Awake.(hemDataType).f = f_AwakeData;
        AnalysisResults.(animalID).WhiskHemoCoherence.Awake.(hemDataType).confC = confC_AwakeData;
        AnalysisResults.(animalID).WhiskHemoCoherence.Awake.(hemDataType).cErr = cErr_AwakeData;
    else
        % save results
        AnalysisResults.(animalID).NeuralHemoCoherence.Awake.(hemDataType).C = [];
        AnalysisResults.(animalID).NeuralHemoCoherence.Awake.(hemDataType).f = [];
        AnalysisResults.(animalID).NeuralHemoCoherence.Awake.(hemDataType).confC = [];
        AnalysisResults.(animalID).NeuralHemoCoherence.Awake.(hemDataType).cErr = [];
    end
    %% analyze neural-hemo coherence during periods of asleep
    zz = 1;
    clear HbT_SleepData Gamma_SleepData HbT_ProcSleepData Gamma_ProcSleepData
    HbT_SleepData = [];
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,~,allDataFileID] = GetFileInfo_IOS(procDataFileID);
        scoringLabels = [];
        for cc = 1:length(ScoringResults.fileIDs)
            if strcmp(allDataFileID,ScoringResults.fileIDs{cc,1}) == true
                scoringLabels = ScoringResults.labels{cc,1};
            end
        end
        % check labels to match arousal state
        if sum(strcmp(scoringLabels,'Not Sleep')) < 36   % 36 bins (180 total) or 3 minutes of awake
            load(procDataFileID,'-mat')
            puffs = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(puffs) == true
                HbT_SleepData{zz,1} = ProcData.data.CBV_HbT.(hemDataType)(2:end-1);
                whisk_SleepData{zz,1} = ProcData.data.binWhiskerAngle(1:26998);
                zz = zz + 1;
            end
        end
    end
    % filter and detrend data
    if isempty(HbT_SleepData) == false
        for bb = 1:length(HbT_SleepData)
            HbT_ProcSleepData{bb,1} = filtfilt(sos,g,detrend(HbT_SleepData{bb,1},'constant'));
            whisk_ProcSleepData{bb,1} = detrend(whisk_SleepData{bb,1},'constant');
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        HbT_sleepData = zeros(length(HbT_ProcSleepData{1,1}),length(HbT_ProcSleepData));
        whisk_sleepData = zeros(length(whisk_ProcSleepData{1,1}),length(whisk_ProcSleepData));
        for cc = 1:length(HbT_ProcSleepData)
            HbT_sleepData(:,cc) = HbT_ProcSleepData{cc,1};
            whisk_sleepData(:,cc) = whisk_ProcSleepData{cc,1};
        end
        % parameters for coherencyc - information available in function
        params.tapers = [1,1];   % Tapers [n, 2n - 1]
        params.pad = 1;
        params.Fs = samplingRate;
        params.fpass = [0,15];   % Pass band [0, nyquist]
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the coherence between desired signals
        [C_SleepData,~,~,~,~,f_SleepData,confC_SleepData,~,cErr_SleepData] = coherencyc(HbT_sleepData,whisk_sleepData,params);
        % save results
        AnalysisResults.(animalID).WhiskHemoCoherence.Sleep.(hemDataType).C = C_SleepData;
        AnalysisResults.(animalID).WhiskHemoCoherence.Sleep.(hemDataType).f = f_SleepData;
        AnalysisResults.(animalID).WhiskHemoCoherence.Sleep.(hemDataType).confC = confC_SleepData;
        AnalysisResults.(animalID).WhiskHemoCoherence.Sleep.(hemDataType).cErr = cErr_SleepData;
    else
        % save results
        AnalysisResults.(animalID).WhiskHemoCoherence.Sleep.(hemDataType).C = [];
        AnalysisResults.(animalID).WhiskHemoCoherence.Sleep.(hemDataType).f = [];
        AnalysisResults.(animalID).WhiskHemoCoherence.Sleep.(hemDataType).confC = [];
        AnalysisResults.(animalID).WhiskHemoCoherence.Sleep.(hemDataType).cErr = [];
    end
    %% analyze neural-hemo coherence during periods of all data
    zz = 1;
    clear HbT_AllUnstimData Gamma_AllUnstimData HbT_ProcAllUnstimData Gamma_ProcAllUnstimData
    HbT_AllUnstimData = [];
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        load(procDataFileID,'-mat')
        puffs = ProcData.data.stimulations.LPadSol;
        % don't include trials with stimulation
        if isempty(puffs) == true
            HbT_AllUnstimData{zz,1} = ProcData.data.CBV_HbT.(hemDataType)(2:end-1);
            whisk_AllUnstimData{zz,1} = ProcData.data.binWhiskerAngle(1:26998);
            zz = zz + 1;
        end
    end
    % filter and detrend data
    if isempty(HbT_AllUnstimData) == false
        for bb = 1:length(HbT_AllUnstimData)
            HbT_ProcAllUnstimData{bb,1} = filtfilt(sos,g,detrend(HbT_AllUnstimData{bb,1},'constant'));
            whisk_ProcAllUnstimData{bb,1} = detrend(whisk_AllUnstimData{bb,1},'constant');
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        HbT_allUnstimData = zeros(length(HbT_ProcAllUnstimData{1,1}),length(HbT_ProcAllUnstimData));
        whisk_allUnstimData = zeros(length(whisk_ProcAllUnstimData{1,1}),length(whisk_ProcAllUnstimData));
        for cc = 1:length(HbT_ProcAllUnstimData)
            HbT_allUnstimData(:,cc) = HbT_ProcAllUnstimData{cc,1};
            whisk_allUnstimData(:,cc) = whisk_ProcAllUnstimData{cc,1};
        end
        % parameters for coherencyc - information available in function
        params.tapers = [1,1];   % Tapers [n, 2n - 1]
        params.pad = 1;
        params.Fs = samplingRate;
        params.fpass = [0,15];   % Pass band [0, nyquist]
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the coherence between desired signals
        [C_AllUnstimData,~,~,~,~,f_AllUnstimData,confC_AllUnstimData,~,cErr_AllUnstimData] = coherencyc(HbT_allUnstimData,whisk_allUnstimData,params);
        % save results
        AnalysisResults.(animalID).WhiskHemoCoherence.All.(hemDataType).C = C_AllUnstimData;
        AnalysisResults.(animalID).WhiskHemoCoherence.All.(hemDataType).f = f_AllUnstimData;
        AnalysisResults.(animalID).WhiskHemoCoherence.All.(hemDataType).confC = confC_AllUnstimData;
        AnalysisResults.(animalID).WhiskHemoCoherence.All.(hemDataType).cErr = cErr_AllUnstimData;
    end
end
% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults','-v7.3')

end
