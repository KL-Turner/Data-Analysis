function [Results_PowerSpecLFP] = AnalyzePowerSpectrum_LFP(animalID,group,rootFolder,delim,Results_PowerSpecLFP)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the spectral power of hemodynamic [HbT] and neural signals (IOS)
%________________________________________________________________________________________________________________________

%% only run analysis for valid animal IDs
dataLocation = [rootFolder delim group delim animalID delim 'Bilateral Imaging'];
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
%% analyze power spectra during periods of alert/asleep/all
behavFields = {'Alert','Asleep','All'};
dataTypes = {'LH','RH','Hip'};
xx = 1; yy = 1; zz = 1;
analogFs = 20000;
dsFs = 1000;
params.tapers = [5,9];   % Tapers [n, 2n - 1]
params.pad = 1;
params.Fs = dsFs;
params.fpass = [1,100];   % Pass band [0, nyquist]
params.trialave = 1;
params.err = [2,0.05];
Data.Alert = []; Data.Asleep = []; Data.All = [];
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
        motionArtifact = ProcData.notes.motionArtifact;
        if motionArtifact == false
            Data.All.LH{xx,1} = resample(RawData.data.cortical_LH,dsFs,analogFs);
            Data.All.RH{xx,1} = resample(RawData.data.cortical_RH,dsFs,analogFs);
            Data.All.Hip{xx,1} = resample(RawData.data.hippocampus,dsFs,analogFs);
            xx = xx + 1;
        end
        % check labels to match arousal state
        if sum(strcmp(scoringLabels,'Not Sleep')) > 144   % 36 bins (180 total) or 3 minutes of sleep
            if motionArtifact == false
                Data.Alert.LH{yy,1} = resample(RawData.data.cortical_LH,dsFs,analogFs);
                Data.Alert.RH{yy,1} = resample(RawData.data.cortical_RH,dsFs,analogFs);
                Data.Alert.Hip{yy,1} = resample(RawData.data.hippocampus,dsFs,analogFs);
                yy = yy + 1;
            end
        elseif sum(strcmp(scoringLabels,'Not Sleep')) < 36   % 36 bins (180 total) or 3 minutes of awake
            if motionArtifact == false
                Data.Asleep.LH{zz,1} = resample(RawData.data.cortical_LH,dsFs,analogFs);
                Data.Asleep.RH{zz,1} = resample(RawData.data.cortical_RH,dsFs,analogFs);
                Data.Asleep.Hip{zz,1} = resample(RawData.data.hippocampus,dsFs,analogFs);
                zz = zz + 1;
            end
        end
    end
end
%% Calculate LFP power spectrum
for aa = 1:length(behavFields)
    behavField = behavFields{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        if isempty(Data.(behavField)) == false
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
            %             dataGPU = gpuArray(data);
            [S,f,sErr] = mtspectrumc(data,params);
            % save results
            Results_PowerSpecLFP.(animalID).(behavField).(dataType).S = S;
            Results_PowerSpecLFP.(animalID).(behavField).(dataType).f = f;
            Results_PowerSpecLFP.(animalID).(behavField).(dataType).sErr = sErr;
        else
            % save results
            Results_PowerSpecLFP.(animalID).(behavField).(dataType).S = [];
            Results_PowerSpecLFP.(animalID).(behavField).(dataType).f = [];
            Results_PowerSpecLFP.(animalID).(behavField).(dataType).sErr = [];
        end
    end
end
% save data
cd(rootFolder)
save('Results_PowerSpecLFP.mat','Results_PowerSpecLFP')

end
