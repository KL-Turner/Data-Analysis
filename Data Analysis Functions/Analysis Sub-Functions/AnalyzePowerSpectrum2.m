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

%% analyze power spectra during periods of alert
behavFields = {'Alert','Asleep','All'};
dataTypes = {'LH','RH','Hip'};
xx = 1; yy = 1; zz = 1;
analogFs = 20000;
dsFs = 1000;
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
        Data.All.LH{xx,1} = resample(RawData.data.cortical_LH,dsFs,analogFs);
        Data.All.RH{xx,1} = resample(RawData.data.cortical_RH,dsFs,analogFs);
        Data.All.Hip{yy,1} = resample(RawData.data.hippocampus,dsFs,analogFs);
        xx = xx + 1;
        % check labels to match arousal state
        if sum(strcmp(scoringLabels,'Not Sleep')) > 144   % 36 bins (180 total) or 3 minutes of sleep
            Data.Alert.LH{yy,1} = resample(RawData.data.cortical_LH,dsFs,analogFs);
            Data.Alert.RH{yy,1} = resample(RawData.data.cortical_RH,dsFs,analogFs);
            Data.Alert.Hip{yy,1} = resample(RawData.data.hippocampus,dsFs,analogFs);
            yy = yy + 1;
        elseif sum(strcmp(scoringLabels,'Not Sleep')) < 36   % 36 bins (180 total) or 3 minutes of awake
            Data.Asleep.LH{zz,1} = resample(RawData.data.cortical_LH,dsFs,analogFs);
            Data.Asleep.RH{zz,1} = resample(RawData.data.cortical_RH,dsFs,analogFs);
            Data.Asleep.Hip{yy,1} = resample(RawData.data.hippocampus,dsFs,analogFs);
            zz = zz + 1;
        end
    end
end
%% Calculate LFP power spectrum
for aa = 1:length(behavFields)
    behavField = behavFields{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        if isempty(Data.(behavField).(dataType)) == false
            % detrend data
            procData = {};
            for cc = 1:length(Data.(behavField).(dataType))
                procData{cc,1} = detrend(Data.(behavField).(dataType){cc,1},'constant');
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            data = zeros(length(procData{1,1}),length(procData));
            for dd = 1:length(procData)
                data(:,dd) = procData{dd,1};
            end
            % calculate the power spectra of the desired signals
            [S,f,sErr] = mtspectrumc(data,params);
            % save results
            AnalysisResults.(animalID).LFP.(behavField).(dataType).S = S;
            AnalysisResults.(animalID).LFP.(behavField).(dataType).f = f;
            AnalysisResults.(animalID).LFP.(behavField).(dataType).sErr = sErr;
        else
            % save results
            AnalysisResults.(animalID).LFP.(behavField).(dataType).S = [];
            AnalysisResults.(animalID).LFP.(behavField).(dataType).f = [];
            AnalysisResults.(animalID).LFP.(behavField).(dataType).sErr = [];
        end
    end
end
end
