function [AnalysisResults] = AnalyzeCoherence_IOS(dataTypes,params,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%
%   Purpose: Analyze the spectral coherence between bilateral hemodynamic and neural signals.
%________________________________________________________________________________________________________________________

% find and load RestData.mat struct
restDataFileStruct = dir('*_RestData.mat');
restDataFile = {restDataFileStruct.name}';
restDataFileID = char(restDataFile);
load(restDataFileID)

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

% identify animal's ID and pull important infortmat
fileBreaks = strfind(restDataFileID, '_');
animalID = restDataFileID(1:fileBreaks(1)-1);
manualFileIDs = unique(RestingBaselines.manualSelection.baselineFileInfo.fileIDs);
samplingRate = RestData.CBV_HbT.adjLH.CBVCamSamplingRate;

RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};

PuffCriteria.Fieldname = {'puffDistances'};
PuffCriteria.Comparison = {'gt'};
PuffCriteria.Value = {5};

for a = 1:length(dataTypes)
    dataType = dataTypes{1,a};
    %% Analyze coherence during periods of rest
    disp(['AnalyzeCoherence: ' dataType ' during awake rest']); disp(' ')
    % use the RestCriteria we specified earlier to find unstim resting events that are greater than the criteria
    if strcmp(dataType,'CBV_HbT') == true
        [restLogical] = FilterEvents_IOS(RestData.(dataType).adjLH,RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.(dataType).adjLH,PuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        unstimRestFiles = RestData.(dataType).adjLH.fileIDs(combRestLogical,:);
        LH_unstimRestingData = RestData.(dataType).adjLH.data(combRestLogical,:);
        RH_unstimRestingData = RestData.(dataType).adjRH.data(combRestLogical,:);
    else
        [restLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),PuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        unstimRestFiles = RestData.cortical_LH.(dataType).fileIDs(combRestLogical,:);
        LH_unstimRestingData = RestData.cortical_LH.(dataType).NormData(combRestLogical,:);
        RH_unstimRestingData = RestData.cortical_RH.(dataType).NormData(combRestLogical,:);
    end
    
    % identify the unique days and the unique number of files from the list of unstim resting events
    restUniqueDays = GetUniqueDays_IOS(unstimRestFiles);
    restUniqueFiles = unique(unstimRestFiles);
    restNumberOfFiles = length(unique(unstimRestFiles));
    
    % decimate the file list to only include those files that occur within the desired number of target minutes
    clear restFiltLogical
    for c = 1:length(restUniqueDays)
        restDay = restUniqueDays(c);
        d = 1;
        for e = 1:restNumberOfFiles
            restFile = restUniqueFiles(e);
            restFileID = restFile{1}(1:6);
            if strcmp(restDay,restFileID) && sum(strcmp(restFile,manualFileIDs)) == 1
                restFiltLogical{c,1}(e,1) = 1; %#ok<*AGROW>
                d = d + 1;
            else
                restFiltLogical{c,1}(e,1) = 0;
            end
        end
    end
    restFinalLogical = any(sum(cell2mat(restFiltLogical'),2),2);
    
    % extract unstim the resting events that correspond to the acceptable file list and the acceptable resting criteria
    clear restFileFilter
    filtRestFiles = restUniqueFiles(restFinalLogical,:);
    for f = 1:length(unstimRestFiles)
        restLogic = strcmp(unstimRestFiles{f},filtRestFiles);
        restLogicSum = sum(restLogic);
        if restLogicSum == 1
            restFileFilter(f,1) = 1;
        else
            restFileFilter(f,1) = 0;
        end
    end
    restFinalFileFilter = logical(restFileFilter);
    LH_finalRestData = LH_unstimRestingData(restFinalFileFilter,:);
    RH_finalRestData = RH_unstimRestingData(restFinalFileFilter,:);
    
    % only take the first 10 seconds of the epoch. occassionunstimy a sample gets lost from rounding during the
    % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
    % lowpass filter and detrend each segment
    [B, A] = butter(3,1/(samplingRate/2),'low');
    clear LH_ProcRestData
    clear RH_ProcRestData
    for g = 1:length(LH_finalRestData)
        if length(LH_finalRestData{g,1}) < params.minTime.Rest*samplingRate
            restChunkSampleDiff = params.minTime.Rest*samplingRate - length(LH_finalRestData{g,1});
            LH_restPad = (ones(1,restChunkSampleDiff))*LH_finalRestData{g,1}(end);
            RH_restPad = (ones(1,restChunkSampleDiff))*RH_finalRestData{g,1}(end);
            LH_ProcRestData{g,1} = horzcat(LH_finalRestData{g,1},LH_restPad);
            RH_ProcRestData{g,1} = horzcat(RH_finalRestData{g,1},RH_restPad);
            LH_ProcRestData{g,1} = detrend(filtfilt(B,A,LH_ProcRestData{g,1}),'constant');
            RH_ProcRestData{g,1} = detrend(filtfilt(B,A,RH_ProcRestData{g,1}),'constant');
        else
            LH_ProcRestData{g,1} = detrend(filtfilt(B,A,LH_finalRestData{g,1}(1:(params.minTime.Rest*samplingRate))),'constant');
            RH_ProcRestData{g,1} = detrend(filtfilt(B,A,RH_finalRestData{g,1}(1:(params.minTime.Rest*samplingRate))),'constant');
        end
    end
    
    % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
    LH_restData = zeros(length(LH_ProcRestData{1,1}),length(LH_ProcRestData));
    RH_restData = zeros(length(RH_ProcRestData{1,1}),length(RH_ProcRestData));
    for h = 1:length(LH_ProcRestData)
        LH_restData(:,h) = LH_ProcRestData{h,1};
        RH_restData(:,h) = RH_ProcRestData{h,1};
    end
    
    % parameters for coherencyc_IOS - information available in function
    params.tapers = [3 5];   % Tapers [n, 2n - 1]
    params.pad = 1;
    params.Fs = samplingRate;   % Sampling Rate
    params.fpass = [0 1];   % Pass band [0, nyquist]
    params.trialave = 1;
    params.err = [2 0.05];
    
    % calculate the coherence between desired signals
    [C_RestData,~, ~, ~, ~,f_RestData,confC_RestData,~,cErr_RestData] = coherencyc_IOS(LH_restData,RH_restData,params);
    
    % nboot = 1000;
    % rest_CI = bootci(nboot, @coherencyc2, LH_restData', RH_restData', params);
    
    % summary figure
    restCoherence = figure;
    plot(f_RestData,C_RestData,'k')
    hold on;
    plot(f_RestData,cErr_RestData,'color',colors_IOS('battleship grey'))
    xlabel('Freq (Hz)');
    ylabel('Coherence');
    title([animalID  ' ' dataType ' coherence for resting data']);
    set(gca,'Ticklength',[0 0]);
    legend('Coherence','Jackknife Lower','Jackknife Upper','Location','Southeast');
    set(legend,'FontSize',6);
    ylim([0 1])
    xlim([0 1])
    axis square
    
    % save results
    AnalysisResults.Coherence.Rest.(dataType).C = C_RestData;
    AnalysisResults.Coherence.Rest.(dataType).f = f_RestData;
    AnalysisResults.Coherence.Rest.(dataType).confC = confC_RestData;
    AnalysisResults.Coherence.Rest.(dataType).cErr = cErr_RestData;
    
    % save figure
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Combined Imaging/Figures/Analysis Coherence/'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(restCoherence, [dirpath animalID '_Rest_' dataType '_Coherence']);
    close(restCoherence)
    
    %% Analyze coherence during periods of NREM sleep
    % pull data from SleepData.mat structure
    disp(['AnalyzeCoherence: ' dataType ' during NREM']); disp(' ')
    if strcmp(dataType,'CBV_HbT') == true
        LH_nremData = SleepData.NREM.data.(dataType).LH;
        RH_nremData = SleepData.NREM.data.(dataType).RH;
    else
        LH_nremData = SleepData.NREM.data.cortical_LH.(dataType);
        RH_nremData = SleepData.NREM.data.cortical_RH.(dataType);
    end
    
    % detrend - data is already lowpass filtered
    for j = 1:length(LH_nremData)
        LH_nremData{j,1} = detrend(LH_nremData{j,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        RH_nremData{j,1} = detrend(RH_nremData{j,1}(1:(params.minTime.NREM*samplingRate)),'constant');
    end
    
    % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
    LH_nrem = zeros(length(LH_nremData{1,1}),length(LH_nremData));
    RH_nrem = zeros(length(RH_nremData{1,1}),length(RH_nremData));
    for k = 1:length(LH_nremData)
        LH_nrem(:,k) = LH_nremData{k,1};
        RH_nrem(:,k) = RH_nremData{k,1};
    end
    
    % parameters for coherencyc_IOS - information available in function
    params.tapers = [3 5];   % Tapers [n, 2n - 1]
    params.pad = 1;
    params.Fs = samplingRate;   % Sampling Rate
    params.fpass = [0 1];   % Pass band [0, nyquist]
    params.trialave = 1;
    params.err = [2 0.05];
    
    % calculate the coherence between desired signals
    [C_nrem,~,~,~,~,f_nrem,confC_nrem,~,cErr_nrem] = coherencyc_IOS(LH_nrem,RH_nrem,params);
    
    % nboot = 1000;
    % nrem_CI = bootci(nboot, @coherencyc2, LH_nrem', RH_nrem', params);
    
    % summary figure
    nremCoherence = figure;
    plot(f_nrem,C_nrem,'k')
    hold on;
    plot(f_nrem,cErr_nrem,'color',colors_IOS('battleship grey'))
    xlabel('Freq (Hz)');
    ylabel('Coherence');
    title([animalID  ' ' dataType ' coherence for NREM data']);
    set(gca,'Ticklength',[0 0]);
    legend('Coherence','Jackknife Lower','Jackknife Upper','Location','Southeast');
    set(legend,'FontSize',6);
    ylim([0 1])
    xlim([0 1])
    axis square
    
    % save results
    AnalysisResults.Coherence.NREM.(dataType).C = C_nrem;
    AnalysisResults.Coherence.NREM.(dataType).f = f_nrem;
    AnalysisResults.Coherence.NREM.(dataType).confC = confC_nrem;
    AnalysisResults.Coherence.NREM.(dataType).cErr = cErr_nrem;
    
    % save figure
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Combined Imaging/Figures/Analysis Coherence/'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(nremCoherence, [dirpath animalID '_NREM_' dataType '_Coherence']);
    close(nremCoherence)
    
    %% Analyze coherence during periods of REM sleep
    % pull data from SleepData.mat structure
    disp(['AnalyzeCoherence: ' dataType ' during REM']); disp(' ')
    if strcmp(dataType,'CBV_HbT') == true
        LH_remData = SleepData.REM.data.(dataType).LH;
        RH_remData = SleepData.REM.data.(dataType).RH;
    else
        LH_remData = SleepData.REM.data.cortical_LH.(dataType);
        RH_remData = SleepData.REM.data.cortical_RH.(dataType);
    end
    
    % detrend - data is already lowpass filtered
    for m = 1:length(LH_remData)
        LH_remData{m,1} = detrend(LH_remData{m,1}(1:(params.minTime.REM*samplingRate)),'constant');
        RH_remData{m,1} = detrend(RH_remData{m,1}(1:(params.minTime.REM*samplingRate)),'constant');
    end
    
    % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
    LH_rem = zeros(length(LH_remData{1,1}),length(LH_remData));
    RH_rem = zeros(length(RH_remData{1,1}),length(RH_remData));
    for n = 1:length(LH_remData)
        LH_rem(:,n) = LH_remData{n,1};
        RH_rem(:,n) = RH_remData{n,1};
    end
    
    % parameters for coherencyc_IOS - information available in function
    params.tapers = [3 5];   % Tapers [n, 2n - 1]
    params.pad = 1;
    params.Fs = samplingRate;   % Sampling Rate
    params.fpass = [0 1];   % Pass band [0, nyquist]
    params.trialave = 1;
    params.err = [2 0.05];
    
    % calculate the coherence between desired signals
    [C_rem,~, ~, ~, ~,f_rem, confC_rem, ~,cErr_rem] = coherencyc_IOS(LH_rem,RH_rem,params);
    
    % nboot = 1000;
    % rem_CI = bootci(nboot, @coherencyc2, LH_rem', RH_rem', params);
    
    % summary figure
    remCoherence = figure;
    plot(f_rem,C_rem,'k')
    hold on;
    plot(f_rem,cErr_rem,'color',colors_IOS('battleship grey'))
    xlabel('Freq (Hz)');
    ylabel('Coherence');
    title([animalID  ' ' dataType ' coherence for REM data']);
    set(gca,'Ticklength',[0 0]);
    legend('Coherence','Jackknife Lower','Jackknife Upper','Location','Southeast');
    set(legend,'FontSize',6);
    ylim([0 1])
    xlim([0 1])
    axis square
    
    % save results
    AnalysisResults.Coherence.REM.(dataType).C = C_rem;
    AnalysisResults.Coherence.REM.(dataType).f = f_rem;
    AnalysisResults.Coherence.REM.(dataType).confC = confC_rem;
    AnalysisResults.Coherence.REM.(dataType).cErr = cErr_rem;
    
    % save figure
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Combined Imaging/Figures/Analysis Coherence/'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(remCoherence, [dirpath animalID '_REM_' dataType '_Coherence']);
    close(remCoherence)
end

%% save results strucure
save([animalID '_AnalysisResults.mat'],'AnalysisResults');

end