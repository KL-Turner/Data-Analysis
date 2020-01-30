function [AnalysisResults] = AnalyzePowerSpectrum_IOS(dataTypes,params,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%
%   Purpose: Analyze the spectral power between of hemodynamic and neural signals.
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

%% Analyze coherence during periods of rest
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};

PuffCriteria.Fieldname = {'puffDistances'};
PuffCriteria.Comparison = {'gt'};
PuffCriteria.Value = {5};

for a = 1:length(dataTypes)
    dataType = dataTypes{1,a};
    %% Analyze power spectra during periods of rest
    disp(['AnalyzePowerSpectra: ' dataType ' during Rest']); disp(' ')
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
        LH_unstimRestingData =RestData.cortical_LH.(dataType).NormData(combRestLogical,:);
        RH_unstimRestingData = RestData.cortical_RH.(dataType).NormData(combRestLogical,:);
        Hip_unstimRestingData = RestData.hippocampus.(dataType).NormData(combRestLogical,:);
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
    if strcmp(dataType,'CBV_HbT') == false
        Hip_finalRestData = Hip_unstimRestingData(restFinalFileFilter,:);
    end
    
    % only take the first 10 seconds of the epoch. occassionunstimy a sample gets lost from rounding during the
    % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
    % lowpass filter and detrend each segment
    [B, A] = butter(3,1/(samplingRate/2),'low');
    clear LH_ProcRestData
    clear RH_ProcRestData
    clear Hip_ProcRestData
    for g = 1:length(LH_finalRestData)
        if length(LH_finalRestData{g,1}) < params.minTime.Rest*samplingRate
            restChunkSampleDiff = params.minTime.Rest*samplingRate - length(LH_finalRestData{g,1});
            LH_restPad = (ones(1,restChunkSampleDiff))*LH_finalRestData{g,1}(end);
            RH_restPad = (ones(1,restChunkSampleDiff))*RH_finalRestData{g,1}(end);
            LH_ProcRestData{g,1} = horzcat(LH_finalRestData{g,1},LH_restPad);
            RH_ProcRestData{g,1} = horzcat(RH_finalRestData{g,1},RH_restPad);
            LH_ProcRestData{g,1} = detrend(filtfilt(B,A,LH_ProcRestData{g,1}),'constant');
            RH_ProcRestData{g,1} = detrend(filtfilt(B,A,RH_ProcRestData{g,1}),'constant');
            if strcmp(dataType,'CBV_HbT') == false
                Hip_restPad = (ones(1,restChunkSampleDiff))*Hip_finalRestData{g,1}(end);
                Hip_ProcRestData{g,1} = horzcat(Hip_finalRestData{g,1},Hip_restPad);
                Hip_ProcRestData{g,1} = detrend(filtfilt(B,A,Hip_ProcRestData{g,1}),'constant');
            end
        else
            LH_ProcRestData{g,1} = detrend(filtfilt(B,A,LH_finalRestData{g,1}(1:(params.minTime.Rest*samplingRate))),'constant');
            RH_ProcRestData{g,1} = detrend(filtfilt(B,A,RH_finalRestData{g,1}(1:(params.minTime.Rest*samplingRate))),'constant');
            if strcmp(dataType,'CBV_HbT') == false
                Hip_ProcRestData{g,1} = detrend(filtfilt(B,A,Hip_finalRestData{g,1}(1:(params.minTime.Rest*samplingRate))),'constant');
            end
        end
    end
    
    % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
    LH_restData = zeros(length(LH_ProcRestData{1,1}),length(LH_ProcRestData));
    RH_restData = zeros(length(RH_ProcRestData{1,1}),length(RH_ProcRestData));
    if strcmp(dataType,'CBV_HbT') == false
        Hip_restData = zeros(length(Hip_ProcRestData{1,1}),length(Hip_ProcRestData));
    end
    for n = 1:length(LH_ProcRestData)
        LH_restData(:,n) = LH_ProcRestData{n,1};
        RH_restData(:,n) = RH_ProcRestData{n,1};
        if strcmp(dataType,'CBV_HbT') == false
            Hip_restData(:,n) = Hip_ProcRestData{n,1};
        end
    end
    
    % parameters for coherencyc_IOS - information available in function
    params.tapers = [3 5];   % Tapers [n, 2n - 1]
    params.pad = 1;
    params.Fs = samplingRate;   % Sampling Rate
    params.fpass = [0 1];   % Pass band [0, nyquist]
    params.trialave = 1;
    params.err = [2 0.05];
    
    % calculate the power spectra of the desired signals
    [LH_rest_S,LH_rest_f,LH_rest_sErr] = mtspectrumc_IOS(LH_restData,params);
    [RH_rest_S,RH_rest_f,RH_rest_sErr] = mtspectrumc_IOS(RH_restData,params);
    if strcmp(dataType,'CBV_HbT') == false
        [Hip_rest_S,Hip_rest_f,Hip_rest_sErr] = mtspectrumc_IOS(Hip_restData,params);
    end
    
    % nboot = 1000;
    % LH_restCI = bootci(nboot,@mtspectrumc2,LH_restData',params);
    % RH_restCI = bootci(nboot,@mtspectrumc2,RH_restData',params);
    % if strcmp(dataType,'CBV_HbT') == false
    %     Hip_restCI = bootci(nboot, @mtspectrumc2,Hip_restData',params);
    % end
    
    % summary figures
    LH_RestPower = figure;
    loglog(LH_rest_f,LH_rest_S,'k')
    hold on;
    loglog(LH_rest_f,LH_rest_sErr,'color',colors_IOS('battleship grey'))
    xlabel('Freq (Hz)');
    ylabel('Power');
    title([animalID  ' adjLH ' dataType ' Power during awake rest']);
    set(gca,'Ticklength',[0 0]);
    legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
    set(legend,'FontSize',6);
    xlim([0 1])
    axis square
    
    RH_RestPower = figure;
    loglog(RH_rest_f,RH_rest_S,'k')
    hold on;
    loglog(RH_rest_f,RH_rest_sErr,'color',colors_IOS('battleship grey'))
    xlabel('Freq (Hz)');
    ylabel('Power');
    title([animalID  ' adjRH ' dataType ' Power during awake rest']);
    set(gca,'Ticklength',[0 0]);
    legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
    set(legend,'FontSize',6);
    xlim([0 1])
    axis square
    
    if strcmp(dataType,'CBV_HbT') == false
        Hip_RestPower = figure;
        loglog(Hip_rest_f,Hip_rest_S,'k')
        hold on;
        loglog(Hip_rest_f,Hip_rest_sErr,'color',colors_IOS('battleship grey'))
        xlabel('Freq (Hz)');
        ylabel('Power');
        title([animalID  ' Hippocampal ' dataType ' Power during awake rest']);
        set(gca,'Ticklength',[0 0]);
        legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
        set(legend,'FontSize',6);
        xlim([0 1])
        axis square
    end
    
    % save results
    AnalysisResults.PowerSpectra.Rest.(dataType).adjLH.S = LH_rest_S;
    AnalysisResults.PowerSpectra.Rest.(dataType).adjLH.f = LH_rest_f;
    AnalysisResults.PowerSpectra.Rest.(dataType).adjLH.sErr = LH_rest_sErr;
    AnalysisResults.PowerSpectra.Rest.(dataType).adjRH.S = RH_rest_S;
    AnalysisResults.PowerSpectra.Rest.(dataType).adjRH.f = RH_rest_f;
    AnalysisResults.PowerSpectra.Rest.(dataType).adjRH.sErr = RH_rest_sErr;
    if strcmp(dataType,'CBV_HbT') == false
        AnalysisResults.PowerSpectra.Rest.(dataType).Hip.S = Hip_rest_S;
        AnalysisResults.PowerSpectra.Rest.(dataType).Hip.f = Hip_rest_f;
        AnalysisResults.PowerSpectra.Rest.(dataType).Hip.sErr = Hip_rest_sErr;
    end
    
    % save figures
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Combined Imaging/Figures/Analysis Power Spectra/'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(LH_RestPower,[dirpath animalID '_Rest_LH_' dataType '_PowerSpectra']);
    close(LH_RestPower)
    savefig(RH_RestPower,[dirpath animalID '_Rest_RH_' dataType '_PowerSpectra']);
    close(RH_RestPower)
    if strcmp(dataType,'CBV_HbT') == false
        savefig(Hip_RestPower,[dirpath animalID '_Rest_Hippocampal_' dataType '_PowerSpectra']);
        close(Hip_RestPower)
    end
    
    %% Analyze power spectra during periods of NREM sleep
    % pull data from SleepData.mat structure
    disp(['AnalyzePowerSpectra: ' dataType ' during NREM']); disp(' ')
    if strcmp(dataType,'CBV_HbT') == true
        LH_nremData = SleepData.NREM.data.(dataType).LH;
        RH_nremData = SleepData.NREM.data.(dataType).RH;
    else
        LH_nremData = SleepData.NREM.data.cortical_LH.(dataType);
        RH_nremData = SleepData.NREM.data.cortical_RH.(dataType);
        Hip_nremData = SleepData.NREM.data.hippocampus.(dataType);
    end
    
    % detrend - data is already lowpass filtered
    for j = 1:length(LH_nremData)
        LH_nremData{j,1} = detrend(LH_nremData{j,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        RH_nremData{j,1} = detrend(RH_nremData{j,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        if strcmp(dataType,'CBV_HbT') == false
            Hip_nremData{j,1} = detrend(Hip_nremData{j,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        end
    end
    
    % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
    LH_nrem = zeros(length(LH_nremData{1,1}),length(LH_nremData));
    RH_nrem = zeros(length(RH_nremData{1,1}),length(RH_nremData));
    if strcmp(dataType,'CBV_HbT') == false
        Hip_nrem = zeros(length(Hip_nremData{1,1}),length(Hip_nremData));
    end
    for k = 1:length(LH_nremData)
        LH_nrem(:,k) = LH_nremData{k,1};
        RH_nrem(:,k) = RH_nremData{k,1};
        if strcmp(dataType,'CBV_HbT') == false
            Hip_nrem(:,k) = Hip_nremData{k,1};
        end
    end
    
    % parameters for coherencyc_IOS - information available in function
    params.tapers = [3 5];   % Tapers [n, 2n - 1]
    params.pad = 1;
    params.Fs = samplingRate;   % Sampling Rate
    params.fpass = [0 1];   % Pass band [0, nyquist]
    params.trialave = 1;
    params.err = [2 0.05];
    
    % calculate the power spectra of the desired signals
    [LH_nrem_S,LH_nrem_f,LH_nrem_sErr] = mtspectrumc_IOS(LH_nrem,params);
    [RH_nrem_S,RH_nrem_f,RH_nrem_sErr] = mtspectrumc_IOS(RH_nrem,params);
    if strcmp(dataType,'CBV_HbT') == false
        [Hip_nrem_S,Hip_nrem_f,Hip_nrem_sErr] = mtspectrumc_IOS(Hip_nrem,params);
    end
    
    % nboot = 1000;
    % LH_nremCI = bootci(nboot,@mtspectrumc2,LH_nremData',params);
    % RH_nremCI = bootci(nboot,@mtspectrumc2,RH_nremData',params);
    % if strcmp(dataType,'CBV_HbT') == false
    %     Hip_nremCI = bootci(nboot,@mtspectrumc2,Hip_nremData',params);
    % end
    
    % summary figures
    LH_nremPower = figure;
    loglog(LH_nrem_f,LH_nrem_S,'k')
    hold on;
    loglog(LH_nrem_f,LH_nrem_sErr,'color',colors_IOS('battleship grey'))
    xlabel('Freq (Hz)');
    ylabel('Power');
    title([animalID  ' adjLH ' dataType ' Power during NREM']);
    set(gca,'Ticklength',[0 0]);
    legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
    set(legend,'FontSize',6);
    xlim([0 1])
    axis square
    
    RH_nremPower = figure;
    loglog(RH_nrem_f,RH_nrem_S,'k')
    hold on;
    loglog(RH_nrem_f,RH_nrem_sErr,'color',colors_IOS('battleship grey'))
    xlabel('Freq (Hz)');
    ylabel('Power');
    title([animalID  ' adjRH ' dataType ' Power during NREM']);
    set(gca,'Ticklength',[0 0]);
    legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
    set(legend,'FontSize',6);
    xlim([0 1])
    axis square
    
    if strcmp(dataType,'CBV_HbT') == false
        Hip_nremPower = figure;
        loglog(Hip_nrem_f,Hip_nrem_S,'k')
        hold on;
        loglog(Hip_nrem_f,Hip_nrem_sErr,'color',colors_IOS('battleship grey'))
        xlabel('Freq (Hz)');
        ylabel('Power');
        title([animalID  ' Hippocampal ' dataType ' Power during NREM']);
        set(gca,'Ticklength',[0 0]);
        legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
        set(legend,'FontSize',6);
        xlim([0 1])
        axis square
    end
    
    % save results
    AnalysisResults.PowerSpectra.NREM.(dataType).adjLH.S = LH_nrem_S;
    AnalysisResults.PowerSpectra.NREM.(dataType).adjLH.f = LH_nrem_f;
    AnalysisResults.PowerSpectra.NREM.(dataType).adjLH.sErr = LH_nrem_sErr;
    AnalysisResults.PowerSpectra.NREM.(dataType).adjRH.S = RH_nrem_S;
    AnalysisResults.PowerSpectra.NREM.(dataType).adjRH.f = RH_nrem_f;
    AnalysisResults.PowerSpectra.NREM.(dataType).adjRH.sErr = RH_nrem_sErr;
    if strcmp(dataType,'CBV_HbT') == false
        AnalysisResults.PowerSpectra.NREM.(dataType).Hip.S = Hip_nrem_S;
        AnalysisResults.PowerSpectra.NREM.(dataType).Hip.f = Hip_nrem_f;
        AnalysisResults.PowerSpectra.NREM.(dataType).Hip.sErr = Hip_nrem_sErr;
    end
    
    % save figures
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Combined Imaging/Figures/Analysis Power Spectra/'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(LH_nremPower,[dirpath animalID '_NREM_LH_' dataType '_PowerSpectra']);
    close(LH_nremPower)
    savefig(RH_nremPower,[dirpath animalID '_NREM_RH_' dataType '_PowerSpectra']);
    close(RH_nremPower)
    if strcmp(dataType,'CBV_HbT') == false
        savefig(Hip_nremPower,[dirpath animalID '_NREM_Hippocampal_' dataType '_PowerSpectra']);
        close(Hip_nremPower)
    end
    
    %% Analyze power spectra during periods of REM sleep
    % pull data from SleepData.mat structure
    disp(['AnalyzePowerSpectra: ' dataType ' during REM']); disp(' ')
    if strcmp(dataType,'CBV_HbT') == true
        LH_remData = SleepData.REM.data.(dataType).LH;
        RH_remData = SleepData.REM.data.(dataType).RH;
    else
        LH_remData = SleepData.REM.data.cortical_LH.(dataType);
        RH_remData = SleepData.REM.data.cortical_RH.(dataType);
        Hip_remData = SleepData.REM.data.hippocampus.(dataType);
    end
    
    % detrend - data is already lowpass filtered
    for j = 1:length(LH_remData)
        LH_remData{j,1} = detrend(LH_remData{j,1}(1:(params.minTime.REM*samplingRate)),'constant');
        RH_remData{j,1} = detrend(RH_remData{j,1}(1:(params.minTime.REM*samplingRate)),'constant');
        if strcmp(dataType,'CBV_HbT') == false
            Hip_remData{j,1} = detrend(Hip_remData{j,1}(1:(params.minTime.REM*samplingRate)),'constant');
        end
    end
    
    % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
    LH_rem = zeros(length(LH_remData{1,1}),length(LH_remData));
    RH_rem = zeros(length(RH_remData{1,1}),length(RH_remData));
    if strcmp(dataType,'CBV_HbT') == false
        Hip_rem = zeros(length(Hip_remData{1,1}),length(Hip_remData));
    end
    for k = 1:length(LH_remData)
        LH_rem(:,k) = LH_remData{k,1};
        RH_rem(:,k) = RH_remData{k,1};
        if strcmp(dataType,'CBV_HbT') == false
            Hip_rem(:,k) = Hip_remData{k,1};
        end
    end
    
    % parameters for coherencyc_IOS - information available in function
    params.tapers = [3 5];   % Tapers [n, 2n - 1]
    params.pad = 1;
    params.Fs = samplingRate;   % Sampling Rate
    params.fpass = [0 1];   % Pass band [0, nyquist]
    params.trialave = 1;
    params.err = [2 0.05];
    
    % calculate the power spectra of the desired signals
    [LH_rem_S,LH_rem_f,LH_rem_sErr] = mtspectrumc_IOS(LH_rem,params);
    [RH_rem_S,RH_rem_f,RH_rem_sErr] = mtspectrumc_IOS(RH_rem,params);
    if strcmp(dataType,'CBV_HbT') == false
        [Hip_rem_S,Hip_rem_f,Hip_rem_sErr] = mtspectrumc_IOS(Hip_rem,params);
    end
    
    % nboot = 1000;
    % LH_remCI = bootci(nboot,@mtspectrumc2,LH_remData',params);
    % RH_remCI = bootci(nboot,@mtspectrumc2,RH_remData',params);
    % if strcmp(dataType,'CBV_HbT') == false
    %     Hip_remCI = bootci(nboot,@mtspectrumc2,Hip_remData',params);
    % end
    
    % summary figures
    LH_remPower = figure;
    loglog(LH_rem_f,LH_rem_S,'k')
    hold on;
    loglog(LH_rem_f,LH_rem_sErr,'color',colors_IOS('battleship grey'))
    xlabel('Freq (Hz)');
    ylabel('Power');
    title([animalID  ' adjLH ' dataType ' Power during REM']);
    set(gca,'Ticklength',[0 0]);
    legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
    set(legend,'FontSize',6);
    xlim([0 1])
    axis square
    
    RH_remPower = figure;
    loglog(RH_rem_f,RH_rem_S,'k')
    hold on;
    loglog(RH_rem_f,RH_rem_sErr,'color',colors_IOS('battleship grey'))
    xlabel('Freq (Hz)');
    ylabel('Power');
    title([animalID  ' adjRH ' dataType ' Power during REM']);
    set(gca,'Ticklength',[0 0]);
    legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
    set(legend,'FontSize',6);
    xlim([0 1])
    axis square
    
    if strcmp(dataType,'CBV_HbT') == false
        Hip_remPower = figure;
        loglog(Hip_rem_f,Hip_rem_S,'k')
        hold on;
        loglog(Hip_rem_f,Hip_rem_sErr,'color',colors_IOS('battleship grey'))
        xlabel('Freq (Hz)');
        ylabel('Power');
        title([animalID  ' Hippocampal ' dataType ' Power during REM']);
        set(gca,'Ticklength',[0 0]);
        legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
        set(legend,'FontSize',6);
        xlim([0 1])
        axis square
    end
    
    % save results
    AnalysisResults.PowerSpectra.REM.(dataType).adjLH.S = LH_rem_S;
    AnalysisResults.PowerSpectra.REM.(dataType).adjLH.f = LH_rem_f;
    AnalysisResults.PowerSpectra.REM.(dataType).adjLH.sErr = LH_rem_sErr;
    AnalysisResults.PowerSpectra.REM.(dataType).adjRH.S = RH_rem_S;
    AnalysisResults.PowerSpectra.REM.(dataType).adjRH.f = RH_rem_f;
    AnalysisResults.PowerSpectra.REM.(dataType).adjRH.sErr = RH_rem_sErr;
    if strcmp(dataType,'CBV_HbT') == false
        AnalysisResults.PowerSpectra.REM.(dataType).Hip.S = Hip_rem_S;
        AnalysisResults.PowerSpectra.REM.(dataType).Hip.f = Hip_rem_f;
        AnalysisResults.PowerSpectra.REM.(dataType).Hip.sErr = Hip_rem_sErr;
    end
    
    % save figures
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Combined Imaging/Figures/Analysis Power Spectra/'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(LH_remPower,[dirpath animalID '_REM_LH_' dataType '_PowerSpectra']);
    close(LH_remPower)
    savefig(RH_remPower,[dirpath animalID '_REM_RH_' dataType '_PowerSpectra']);
    close(RH_remPower)
    if strcmp(dataType,'CBV_HbT') == false
        savefig(Hip_remPower,[dirpath animalID '_REM_Hippocampal_' dataType '_PowerSpectra']);
        close(Hip_remPower)
    end
end

%% save results struct
save([animalID '_AnalysisResults.mat'], 'AnalysisResults');

end