function [AnalysisResults] = AnalyzePowerSpectrum2(animalID,group,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the spectral power of hemodynamic [HbT] and neural signals (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
paramsA.minTime.Rest = 10;
paramsA.minTime.NREM = 30;
paramsA.minTime.REM = 60;
%% only run analysis for valid animal IDs
dataLocation = [rootFolder '\' group '\' animalID '\Bilateral Imaging\'];
cd(dataLocation)
% character list of all RawData file IDs
rawDataFileStruct = dir('*_RawData.mat');
rawDataFiles = {rawDataFileStruct.name}';
rawDataFileIDs = char(rawDataFiles);
% character list of all ProcData file IDs
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% find and load Forest_ScoringResults.mat struct
forestScoringResultsFileID = [animalID '_Forest_ScoringResults.mat'];
load(forestScoringResultsFileID,'-mat')
% lowpass filter
samplingRate = 20000;
dsFs = 1000;
% parameters for mtspectrumc - information available in function
paramsA.tapers = [3,5];   % Tapers [n, 2n - 1]
paramsA.pad = 1;
paramsA.Fs = dsFs;
paramsA.fpass = [1,100];   % Pass band [0, nyquist]
paramsA.trialave = 1;
paramsA.err = [2,0.05];
paramsB.tapers = [3,5];   % Tapers [n, 2n - 1]
paramsB.pad = 1;
paramsB.Fs = dsFs;
paramsB.fpass = [300,3000];   % Pass band [0, nyquist]
paramsB.trialave = 1;
paramsB.err = [2,0.05];
%% analyze power spectra during periods of alert
xx = 1; yy = 1; zz = 1;
LH_AlertData = []; RH_AlertData = []; LH_AsleepData = []; RH_AsleepData = []; LH_AllData = []; RH_AllData = [];
for bb = 1:size(rawDataFileIDs,1)
    rawDataFileID = rawDataFileIDs(bb,:);
    procDataFileID = procDataFileIDs(bb,:);
    [~,~,allDataFileID] = GetFileInfo_IOS(procDataFileID);
    scoringLabels = [];
    for cc = 1:length(ScoringResults.fileIDs)
        if strcmp(allDataFileID,ScoringResults.fileIDs{cc,1}) == true
            scoringLabels = ScoringResults.labels{cc,1};
        end
    end
    load(procDataFileID,'-mat')
    puffs = ProcData.data.stimulations.LPadSol;
    % don't include trials with stimulation
    if isempty(puffs) == true
        load(rawDataFileID,'-mat')
        LH_AllData{xx,1} = resample(RawData.data.cortical_LH,dsFs,samplingRate);
        RH_AllData{xx,1} = resample(RawData.data.cortical_RH,dsFs,samplingRate);
        xx = xx + 1;
        % check labels to match arousal state
        if sum(strcmp(scoringLabels,'Not Sleep')) > 144   % 36 bins (180 total) or 3 minutes of sleep
            LH_AlertData{yy,1} = resample(RawData.data.cortical_LH,dsFs,samplingRate);
            RH_AlertData{yy,1} = resample(RawData.data.cortical_RH,dsFs,samplingRate);
            yy = yy + 1;
        elseif sum(strcmp(scoringLabels,'Not Sleep')) < 36   % 36 bins (180 total) or 3 minutes of awake
            LH_AsleepData{zz,1} = resample(RawData.data.cortical_LH,dsFs,samplingRate);
            RH_AsleepData{zz,1} = resample(RawData.data.cortical_RH,dsFs,samplingRate);
            zz = zz + 1;
        end
    end
end
%% alert
if isempty(LH_AlertData) == false
    % detrend data
    for bb = 1:length(LH_AlertData)
        LH_ProcAlertData{bb,1} = detrend(LH_AlertData{bb,1},'constant');
        RH_ProcAlertData{bb,1} = detrend(RH_AlertData{bb,1},'constant');
    end
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
    LH_alertData = zeros(length(LH_ProcAlertData{1,1}),length(LH_ProcAlertData));
    RH_alertData = zeros(length(RH_ProcAlertData{1,1}),length(RH_ProcAlertData));
    for cc = 1:length(LH_ProcAlertData)
        LH_alertData(:,cc) = LH_ProcAlertData{cc,1};
        RH_alertData(:,cc) = RH_ProcAlertData{cc,1};
    end
    % calculate the power spectra of the desired signals
    [LH_alert_S_gamma,LH_alert_f_gamma,LH_alert_sErr_gamma] = mtspectrumc(LH_alertData,paramsA);
    [RH_alert_S_gamma,RH_alert_f_gamma,RH_alert_sErr_gamma] = mtspectrumc(RH_alertData,paramsA);
    % calculate the power spectra of the desired signals
    [LH_alert_S_mua,LH_alert_f_mua,LH_alert_sErr_mua] = mtspectrumc(LH_alertData,paramsB);
    [RH_alert_S_mua,RH_alert_f_mua,RH_alert_sErr_mua] = mtspectrumc(RH_alertData,paramsB);
    % save results
    AnalysisResults.(animalID).PowerSpectra2.Alert.gammaBandPower.LH.S = LH_alert_S_gamma;
    AnalysisResults.(animalID).PowerSpectra2.Alert.gammaBandPower.LH.f = LH_alert_f_gamma;
    AnalysisResults.(animalID).PowerSpectra2.Alert.gammaBandPower.LH.sErr = LH_alert_sErr_gamma;
    AnalysisResults.(animalID).PowerSpectra2.Alert.gammaBandPower.RH.S = RH_alert_S_gamma;
    AnalysisResults.(animalID).PowerSpectra2.Alert.gammaBandPower.RH.f = RH_alert_f_gamma;
    AnalysisResults.(animalID).PowerSpectra2.Alert.gammaBandPower.RH.sErr = RH_alert_sErr_gamma;
    AnalysisResults.(animalID).PowerSpectra2.Alert.multiUnitActivity.LH.S = LH_alert_S_mua;
    AnalysisResults.(animalID).PowerSpectra2.Alert.multiUnitActivity.LH.f = LH_alert_f_mua;
    AnalysisResults.(animalID).PowerSpectra2.Alert.multiUnitActivity.LH.sErr = LH_alert_sErr_mua;
    AnalysisResults.(animalID).PowerSpectra2.Alert.multiUnitActivity.RH.S = RH_alert_S_mua;
    AnalysisResults.(animalID).PowerSpectra2.Alert.multiUnitActivity.RH.f = RH_alert_f_mua;
    AnalysisResults.(animalID).PowerSpectra2.Alert.multiUnitActivity.RH.sErr = RH_alert_sErr_mua;
else
    % save results
    AnalysisResults.(animalID).PowerSpectra2.Alert.gammaBandPower.LH.S = [];
    AnalysisResults.(animalID).PowerSpectra2.Alert.gammaBandPower.LH.f = [];
    AnalysisResults.(animalID).PowerSpectra2.Alert.gammaBandPower.LH.sErr = [];
    AnalysisResults.(animalID).PowerSpectra2.Alert.gammaBandPower.RH.S = [];
    AnalysisResults.(animalID).PowerSpectra2.Alert.gammaBandPower.RH.f = [];
    AnalysisResults.(animalID).PowerSpectra2.Alert.gammaBandPower.RH.sErr = [];
    AnalysisResults.(animalID).PowerSpectra2.Alert.multiUnitActivity.LH.S = [];
    AnalysisResults.(animalID).PowerSpectra2.Alert.multiUnitActivity.LH.f = [];
    AnalysisResults.(animalID).PowerSpectra2.Alert.multiUnitActivity.LH.sErr = [];
    AnalysisResults.(animalID).PowerSpectra2.Alert.multiUnitActivity.RH.S = [];
    AnalysisResults.(animalID).PowerSpectra2.Alert.multiUnitActivity.RH.f = [];
    AnalysisResults.(animalID).PowerSpectra2.Alert.multiUnitActivity.RH.sErr = [];
end
%% Asleep
if isempty(LH_AsleepData) == false
    % detrend data
    for bb = 1:length(LH_AsleepData)
        LH_ProcAsleepData{bb,1} = detrend(LH_AsleepData{bb,1},'constant');
        RH_ProcAsleepData{bb,1} = detrend(RH_AsleepData{bb,1},'constant');
    end
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
    LH_asleepData = zeros(length(LH_ProcAsleepData{1,1}),length(LH_ProcAsleepData));
    RH_asleepData = zeros(length(RH_ProcAsleepData{1,1}),length(RH_ProcAsleepData));
    for cc = 1:length(LH_ProcAsleepData)
        LH_asleepData(:,cc) = LH_ProcAsleepData{cc,1};
        RH_asleepData(:,cc) = RH_ProcAsleepData{cc,1};
    end
    % calculate the power spectra of the desired signals
    [LH_asleep_S_gamma,LH_asleep_f_gamma,LH_asleep_sErr_gamma] = mtspectrumc(LH_asleepData,paramsA);
    [RH_asleep_S_gamma,RH_asleep_f_gamma,RH_asleep_sErr_gamma] = mtspectrumc(RH_asleepData,paramsA);
    % calculate the power spectra of the desired signals
    [LH_asleep_S_mua,LH_asleep_f_mua,LH_asleep_sErr_mua] = mtspectrumc(LH_asleepData,paramsB);
    [RH_asleep_S_mua,RH_asleep_f_mua,RH_asleep_sErr_mua] = mtspectrumc(RH_asleepData,paramsB);
    % save results
    AnalysisResults.(animalID).PowerSpectra2.Asleep.gammaBandPower.LH.S = LH_asleep_S_gamma;
    AnalysisResults.(animalID).PowerSpectra2.Asleep.gammaBandPower.LH.f = LH_asleep_f_gamma;
    AnalysisResults.(animalID).PowerSpectra2.Asleep.gammaBandPower.LH.sErr = LH_asleep_sErr_gamma;
    AnalysisResults.(animalID).PowerSpectra2.Asleep.gammaBandPower.RH.S = RH_asleep_S_gamma;
    AnalysisResults.(animalID).PowerSpectra2.Asleep.gammaBandPower.RH.f = RH_asleep_f_gamma;
    AnalysisResults.(animalID).PowerSpectra2.Asleep.gammaBandPower.RH.sErr = RH_asleep_sErr_gamma;
    AnalysisResults.(animalID).PowerSpectra2.Asleep.multiUnitActivity.LH.S = LH_asleep_S_mua;
    AnalysisResults.(animalID).PowerSpectra2.Asleep.multiUnitActivity.LH.f = LH_asleep_f_mua;
    AnalysisResults.(animalID).PowerSpectra2.Asleep.multiUnitActivity.LH.sErr = LH_asleep_sErr_mua;
    AnalysisResults.(animalID).PowerSpectra2.Asleep.multiUnitActivity.RH.S = RH_asleep_S_mua;
    AnalysisResults.(animalID).PowerSpectra2.Asleep.multiUnitActivity.RH.f = RH_asleep_f_mua;
    AnalysisResults.(animalID).PowerSpectra2.Asleep.multiUnitActivity.RH.sErr = RH_asleep_sErr_mua;
else
    % save results
    AnalysisResults.(animalID).PowerSpectra2.Asleep.gammaBandPower.LH.S = [];
    AnalysisResults.(animalID).PowerSpectra2.Asleep.gammaBandPower.LH.f = [];
    AnalysisResults.(animalID).PowerSpectra2.Asleep.gammaBandPower.LH.sErr = [];
    AnalysisResults.(animalID).PowerSpectra2.Asleep.gammaBandPower.RH.S = [];
    AnalysisResults.(animalID).PowerSpectra2.Asleep.gammaBandPower.RH.f = [];
    AnalysisResults.(animalID).PowerSpectra2.Asleep.gammaBandPower.RH.sErr = [];
    AnalysisResults.(animalID).PowerSpectra2.Asleep.multiUnitActivity.LH.S = [];
    AnalysisResults.(animalID).PowerSpectra2.Asleep.multiUnitActivity.LH.f = [];
    AnalysisResults.(animalID).PowerSpectra2.Asleep.multiUnitActivity.LH.sErr = [];
    AnalysisResults.(animalID).PowerSpectra2.Asleep.multiUnitActivity.RH.S = [];
    AnalysisResults.(animalID).PowerSpectra2.Asleep.multiUnitActivity.RH.f = [];
    AnalysisResults.(animalID).PowerSpectra2.Asleep.multiUnitActivity.RH.sErr = [];
end
%% All
if isempty(LH_AllData) == false
    % detrend data
    for bb = 1:length(LH_AllData)
        LH_ProcAllData{bb,1} = detrend(LH_AllData{bb,1},'constant');
        RH_ProcAllData{bb,1} = detrend(RH_AllData{bb,1},'constant');
    end
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
    LH_allData = zeros(length(LH_ProcAllData{1,1}),length(LH_ProcAllData));
    RH_allData = zeros(length(RH_ProcAllData{1,1}),length(RH_ProcAllData));
    for cc = 1:length(LH_ProcAllData)
        LH_allData(:,cc) = LH_ProcAllData{cc,1};
        RH_allData(:,cc) = RH_ProcAllData{cc,1};
    end
    % calculate the power spectra of the desired signals
    [LH_all_S_gamma,LH_all_f_gamma,LH_all_sErr_gamma] = mtspectrumc(LH_allData,paramsA);
    [RH_all_S_gamma,RH_all_f_gamma,RH_all_sErr_gamma] = mtspectrumc(RH_allData,paramsA);
    % calculate the power spectra of the desired signals
    [LH_all_S_mua,LH_all_f_mua,LH_all_sErr_mua] = mtspectrumc(LH_allData,paramsB);
    [RH_all_S_mua,RH_all_f_mua,RH_all_sErr_mua] = mtspectrumc(RH_allData,paramsB);
    % save results
    AnalysisResults.(animalID).PowerSpectra2.All.gammaBandPower.LH.S = LH_all_S_gamma;
    AnalysisResults.(animalID).PowerSpectra2.All.gammaBandPower.LH.f = LH_all_f_gamma;
    AnalysisResults.(animalID).PowerSpectra2.All.gammaBandPower.LH.sErr = LH_all_sErr_gamma;
    AnalysisResults.(animalID).PowerSpectra2.All.gammaBandPower.RH.S = RH_all_S_gamma;
    AnalysisResults.(animalID).PowerSpectra2.All.gammaBandPower.RH.f = RH_all_f_gamma;
    AnalysisResults.(animalID).PowerSpectra2.All.gammaBandPower.RH.sErr = RH_all_sErr_gamma;
    AnalysisResults.(animalID).PowerSpectra2.All.multiUnitActivity.LH.S = LH_all_S_mua;
    AnalysisResults.(animalID).PowerSpectra2.All.multiUnitActivity.LH.f = LH_all_f_mua;
    AnalysisResults.(animalID).PowerSpectra2.All.multiUnitActivity.LH.sErr = LH_all_sErr_mua;
    AnalysisResults.(animalID).PowerSpectra2.All.multiUnitActivity.RH.S = RH_all_S_mua;
    AnalysisResults.(animalID).PowerSpectra2.All.multiUnitActivity.RH.f = RH_all_f_mua;
    AnalysisResults.(animalID).PowerSpectra2.All.multiUnitActivity.RH.sErr = RH_all_sErr_mua;
else
    % save results
    AnalysisResults.(animalID).PowerSpectra2.All.gammaBandPower.LH.S = [];
    AnalysisResults.(animalID).PowerSpectra2.All.gammaBandPower.LH.f = [];
    AnalysisResults.(animalID).PowerSpectra2.All.gammaBandPower.LH.sErr = [];
    AnalysisResults.(animalID).PowerSpectra2.All.gammaBandPower.RH.S = [];
    AnalysisResults.(animalID).PowerSpectra2.All.gammaBandPower.RH.f = [];
    AnalysisResults.(animalID).PowerSpectra2.All.gammaBandPower.RH.sErr = [];
    AnalysisResults.(animalID).PowerSpectra2.All.multiUnitActivity.LH.S = [];
    AnalysisResults.(animalID).PowerSpectra2.All.multiUnitActivity.LH.f = [];
    AnalysisResults.(animalID).PowerSpectra2.All.multiUnitActivity.LH.sErr = [];
    AnalysisResults.(animalID).PowerSpectra2.All.multiUnitActivity.RH.S = [];
    AnalysisResults.(animalID).PowerSpectra2.All.multiUnitActivity.RH.f = [];
    AnalysisResults.(animalID).PowerSpectra2.All.multiUnitActivity.RH.sErr = [];
end
%% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults','-v7.3')

end