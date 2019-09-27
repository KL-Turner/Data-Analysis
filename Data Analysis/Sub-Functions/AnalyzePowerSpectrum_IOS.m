function [AnalysisResults] = AnalyzePowerSpectrum_IOS(dataTypes,baselineType,params,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%________________________________________________________________________________________________________________________

% list of all Procdata.mat files
procDataFileStruct = dir('*_Procdata.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);

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

% parameters for coherencyc_IOS - information available in function
samplingRate = RestData.CBV.LH.CBVCamSamplingRate;
params.tapers = [3 5];   % Tapers [n, 2n - 1]
params.pad = 1;
params.Fs = samplingRate;   % Sampling Rate
params.fpass = [0 1];   % Pass band [0, nyquist]
params.trialave = 1;
params.err = [2 0.05];

% identify animal's ID and pull important infortmat
fileBreaks = strfind(restDataFileID, '_');
animalID = restDataFileID(1:fileBreaks(1)-1);
trialDuration_min = RestData.CBV.LH.trialDuration_sec/60;   % min
manualFileIDs = unique(RestingBaselines.manualSelection.baselineFileInfo.fileIDs);
fileTarget = params.targetMinutes/trialDuration_min;
filterSets = {'manualSelection','setDuration','entireDuration'};

%% Analyze coherence during periods of rest
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};

PuffCriteria.Fieldname = {'puffDistances'};
PuffCriteria.Comparison = {'gt'};
PuffCriteria.Value = {5};

for a = 1:length(dataTypes)
    dataType = dataTypes{1,a};
    for b = 1:length(filterSets)
        filterSet = filterSets{1,b};        
        %% Analyze power spectra during periods of rest
        % use the RestCriteria we specified earlier to find all resting events that are greater than the criteria
        if strcmp(dataType, 'CBV') == true || strcmp(dataType,'CBV_HbT') == true
            [restLogical] = FilterEvents_IOS(RestData.(dataType).LH,RestCriteria);
            [puffLogical] = FilterEvents_IOS(RestData.(dataType).LH,PuffCriteria);
            combRestLogical = logical(restLogical.*puffLogical);
            allRestFiles = RestData.(dataType).LH.fileIDs(combRestLogical,:);
            if strcmp(dataType,'CBV') == true
                LH_allRestingData = RestData.(dataType).LH.NormData(combRestLogical,:);
                RH_allRestingData = RestData.(dataType).RH.NormData(combRestLogical,:);
            else
                LH_allRestingData = RestData.(dataType).LH.data(combRestLogical,:);
                RH_allRestingData = RestData.(dataType).RH.data(combRestLogical,:);
            end
        else
            [restLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestCriteria);
            [puffLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),PuffCriteria);
            combRestLogical = logical(restLogical.*puffLogical);
            allRestFiles = RestData.cortical_LH.(dataType).fileIDs(combRestLogical,:);
            LH_allRestingData =RestData.cortical_LH.(dataType).NormData(combRestLogical,:);
            RH_allRestingData = RestData.cortical_RH.(dataType).NormData(combRestLogical,:);
            Hip_allRestingData = RestData.hippocampus.(dataType).NormData(combRestLogical,:);
        end
        
        % identify the unique days and the unique number of files from the list of all resting events
        restUniqueDays = GetUniqueDays_IOS(allRestFiles);
        restUniqueFiles = unique(allRestFiles);
        restNumberOfFiles = length(unique(allRestFiles));
        
        % decimate the file list to only include those files that occur within the desired number of target minutes
        clear restFiltLogical
        for c = 1:length(restUniqueDays)
            restDay = restUniqueDays(c);
            d = 1;
            for e = 1:restNumberOfFiles
                restFile = restUniqueFiles(e);
                restFileID = restFile{1}(1:6);
                if strcmp(filterSet,'manualSelection') == true
                    if strcmp(restDay,restFileID) && sum(strcmp(restFile,manualFileIDs)) == 1
                        restFiltLogical{c,1}(e,1) = 1; %#ok<*AGROW>
                        d = d + 1;
                    else
                        restFiltLogical{c,1}(e,1) = 0;
                    end
                elseif strcmp(filterSet,'setDuration') == true
                    if strcmp(restDay,restFileID) && d <= fileTarget
                        restFiltLogical{c,1}(e,1) = 1;
                        d = d + 1;
                    else
                        restFiltLogical{c,1}(e,1) = 0;
                    end
                elseif strcmp(filterSet,'entireDuration') == true
                    if strcmp(restDay,restFileID)
                        restFiltLogical{c,1}(e,1) = 1;
                        d = d + 1;
                    else
                        restFiltLogical{c,1}(e,1) = 0;
                    end
                end
            end
        end
        restFinalLogical = any(sum(cell2mat(restFiltLogical'),2),2);
        
        % extract all the resting events that correspond to the acceptable file list and the acceptable resting criteria
        clear restFileFilter
        filtRestFiles = restUniqueFiles(restFinalLogical,:);
        for f = 1:length(allRestFiles)
            restLogic = strcmp(allRestFiles{f},filtRestFiles);
            restLogicSum = sum(restLogic);
            if restLogicSum == 1
                restFileFilter(f,1) = 1;
            else
                restFileFilter(f,1) = 0;
            end
        end
        restFinalFileFilter = logical(restFinalFileFilter);
        LH_finalRestData = LH_allRestingData(restFinalFileFilter,:);
        RH_finalRestData = RH_allRestingData(restFinalFileFilter,:);
        if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
            Hip_finalRestData = Hip_allRestingData(finalFileFilter,:);
        end
        
        % only take the first 10 seconds of the epoch. occassionally a sample gets lost from rounding during the
        % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
        % lowpass filter and detrend each segment
        [B, A] = butter(4,1/(samplingRate/2),'low');
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
                if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
                    Hip_restPad = (ones(1,restChunkSampleDiff))*Hip_finalRestData{g,1}(end);
                    Hip_ProcRestData{g,1} = horzcat(Hip_finalRestData{g,1},Hip_restPad);
                    Hip_ProcRestData{g,1} = detrend(filtfilt(B,A,Hip_ProcRestData{g,1}),'constant');
                end
            else
                LH_ProcRestData{g,1} = detrend(filtfilt(B,A,LH_finalRestData{g,1}(1:(params.minTime.Rest*samplingRate))),'constant');
                RH_ProcRestData{g,1} = detrend(filtfilt(B,A,RH_finalRestData{g,1}(1:(params.minTime.Rest*samplingRate))),'constant');
                if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
                    Hip_ProcRestData{g,1} = detrend(filtfilt(B,A,Hip_finalRestData{g,1}(1:(params.minTime.Rest*samplingRate))),'constant');
                end
            end
        end
        
        % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontally)
        LH_restData = zeros(length(LH_finalRestData{1,1}),length(LH_finalRestData));
        RH_restData = zeros(length(RH_finalRestData{1,1}),length(RH_finalRestData));
        if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
            Hip_restData = zeros(length(Hip_finalRestData{1,1}),length(Hip_finalRestData));
        end
        for n = 1:length(LH_finalRestData)
            LH_restData(:,n) = LH_finalRestData{n,1};
            RH_restData(:,n) = RH_finalRestData{n,1};
            if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
                Hip_restData(:, n) = Hip_finalRestData{n,1};
            end
        end
        
        % calculate the power spectra of the desired signals
        disp(['Analyzing the power spectrum of the LH RestData (' filterSet ') ' dataType ' signal power...']); disp(' ')
        [LH_rest_S,LH_rest_f,LH_rest_sErr] = mtspectrumc_IOS(LH_restData,params);
        disp(['Analyzing the power spectrum of the RH RestData ' dataType ' signal power...']); disp(' ')
        [RH_rest_S,RH_rest_f,RH_rest_sErr] = mtspectrumc_IOS(RH_restData,params);
        if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
            disp(['Analyzing the power spectrum of the Hippocampal RestData (' filterSet ') ' dataType ' signal power...']); disp(' ')
            [Hip_rest_S,Hip_rest_f,Hip_rest_sErr] = mtspectrumc_IOS(Hip_restData,params);
        end
        
        % nboot = 1000;
        % LH_restCI = bootci(nboot,@mtspectrumc2,LH_restData',params);
        % RH_restCI = bootci(nboot,@mtspectrumc2,RH_restData',params);
        % if strcmp(dataType, 'CBV') == false && strcmp(dataType,'CBV_HbT') == false
        %     Hip_restCI = bootci(nboot, @mtspectrumc2,Hip_restData',params);
        % end
        
        % summary figures
        LH_RestPower = figure;
        loglog(LH_rest_f,LH_rest_S,'k')
        hold on;
        loglog(LH_rest_f,LH_rest_sErr,'color',colors_IOS('battleship grey'))
        xlabel('Freq (Hz)');
        ylabel('Power');
        title([animalID  ' LH ' dataType ' Power during awake rest - ' filterSet]);
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
        title([animalID  ' RH ' dataType ' Power during awake rest - ' filterSet]);
        set(gca,'Ticklength',[0 0]);
        legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
        set(legend,'FontSize',6);
        xlim([0 1])
        axis square
        
        if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
            Hip_RestPower = figure;
            loglog(Hip_rest_f,Hip_rest_S,'k')
            hold on;
            loglog(Hip_rest_f,Hip_rest_sErr,'color',colors_IOS('battleship grey'))
            xlabel('Freq (Hz)');
            ylabel('Power');
            title([animalID  ' Hippocampal ' dataType ' Power during awake rest - ' filterSet]);
            set(gca,'Ticklength',[0 0]);
            legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
            set(legend,'FontSize',6);
            xlim([0 1])
            axis square
        end
        
        % save results
        AnalysisResults.PowerSpectra.Rest.(dataType).LH.S = LH_rest_S;
        AnalysisResults.PowerSpectra.Rest.(dataType).LH.f = LH_rest_f;
        AnalysisResults.PowerSpectra.Rest.(dataType).LH.sErr = LH_rest_sErr;
        AnalysisResults.PowerSpectra.Rest.(dataType).RH.S = RH_rest_S;
        AnalysisResults.PowerSpectra.Rest.(dataType).RH.f = RH_rest_f;
        AnalysisResults.PowerSpectra.Rest.(dataType).RH.sErr = RH_rest_sErr;
        if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
            AnalysisResults.PowerSpectra.Rest.(dataType).Hip.S = Hip_rest_S;
            AnalysisResults.PowerSpectra.Rest.(dataType).Hip.f = Hip_rest_f;
            AnalysisResults.PowerSpectra.Rest.(dataType).Hip.sErr = Hip_rest_sErr;
        end
        
        % save figures
        [pathstr, ~, ~] = fileparts(cd);
        dirpath = [pathstr '/Figures/Analysis Power Spectra/'];
        if ~exist(dirpath, 'dir')
            mkdir(dirpath);
        end
        savefig(LH_RestPower, [dirpath animalID '_Rest_LH_' filterSet '_' dataType '_PowerSpectra']);
        savefig(RH_RestPower, [dirpath animalID '_Rest_RH_' filterSet '_' dataType '_PowerSpectra']);
        if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
            savefig(Hip_RestPower, [dirpath animalID '_Rest_Hippocampal_' filterSet '_' dataType '_PowerSpectra']);
        end
    end
    
    %% Analyze power spectra during periods of NREM sleep
    % pull data from SleepData.mat structure
    if strcmp(dataType,'CBV') == true || strcmp(dataType,'CBV_HbT') == true
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
        if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
            Hip_nremData{j,1} = detrend(Hip_nremData{j,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        end
    end
    
    % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontally)
    LH_nrem = zeros(length(LH_nremData{1,1}),length(LH_nremData));
    RH_nrem = zeros(length(RH_nremData{1,1}),length(RH_nremData));
    if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
        Hip_nrem = zeros(length(Hip_nremData{1,1}),length(Hip_nremData));
    end
    for k = 1:length(LH_nremData)
        LH_nrem(:,k) = LH_nremData{k,1};
        RH_nrem(:,k) = RH_nremData{k,1};
        if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
            Hip_nrem(:,k) = Hip_nremData{k,1};
        end
    end
    
    % calculate the power spectra of the desired signals
    disp(['Analyzing the power spectrum of the LH NREM ' dataType ' signal power...']); disp(' ')
    [LH_nrem_S,LH_nrem_f,LH_nrem_sErr] = mtspectrumc_IOS(LH_nremData,params);
    disp(['Analyzing the power spectrum of the RH NREM ' dataType ' signal power...']); disp(' ')
    [RH_nrem_S,RH_nrem_f,RH_nrem_sErr] = mtspectrumc_IOS(RH_nremData,params);
    if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
        disp(['Analyzing the power spectrum of the Hippocampal NREM ' dataType ' signal power...']); disp(' ')
        [Hip_nrem_S,Hip_nrem_f,Hip_nrem_sErr] = mtspectrumc_IOS(Hip_nremData,params);
    end
    
    % nboot = 1000;
    % LH_nremCI = bootci(nboot,@mtspectrumc2,LH_nremData',params);
    % RH_nremCI = bootci(nboot,@mtspectrumc2,RH_nremData',params);
    % if strcmp(dataType, 'CBV') == false && strcmp(dataType,'CBV_HbT') == false
    %     Hip_nremCI = bootci(nboot,@mtspectrumc2,Hip_nremData',params);
    % end
    
    % summary figures
    LH_nremPower = figure;
    loglog(LH_nrem_f,LH_nrem_S,'k')
    hold on;
    loglog(LH_nrem_f,LH_nrem_sErr,'color',colors_IOS('battleship grey'))
    xlabel('Freq (Hz)');
    ylabel('Power');
    title([animalID  ' LH ' dataType ' Power during NREM']);
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
    title([animalID  ' RH ' dataType ' Power during NREM']);
    set(gca,'Ticklength',[0 0]);
    legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
    set(legend,'FontSize',6);
    xlim([0 1])
    axis square
    
    if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
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
    AnalysisResults.PowerSpectra.NREM.(dataType).LH.S = LH_nrem_S;
    AnalysisResults.PowerSpectra.NREM.(dataType).LH.f = LH_nrem_f;
    AnalysisResults.PowerSpectra.NREM.(dataType).LH.sErr = LH_nrem_sErr;
    AnalysisResults.PowerSpectra.NREM.(dataType).RH.S = RH_nrem_S;
    AnalysisResults.PowerSpectra.NREM.(dataType).RH.f = RH_nrem_f;
    AnalysisResults.PowerSpectra.NREM.(dataType).RH.sErr = RH_nrem_sErr;
    if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
        AnalysisResults.PowerSpectra.NREM.(dataType).Hip.S = Hip_nrem_S;
        AnalysisResults.PowerSpectra.NREM.(dataType).Hip.f = Hip_nrem_f;
        AnalysisResults.PowerSpectra.NREM.(dataType).Hip.sErr = Hip_nrem_sErr;
    end
    
    % save figures
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Analysis Power Spectra/'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(LH_nremPower, [dirpath animalID '_NREM_LH_' filterSet '_' dataType '_PowerSpectra']);
    savefig(RH_nremPower, [dirpath animalID '_NREM_RH_' filterSet '_' dataType '_PowerSpectra']);
    if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
        savefig(Hip_nremPower, [dirpath animalID '_NREM_Hippocampal_' filterSet '_' dataType '_PowerSpectra']);
    end
    
    %% Analyze power spectra during periods of REM sleep
    % pull data from SleepData.mat structure
    if strcmp(dataType,'CBV') == true || strcmp(dataType,'CBV_HbT') == true
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
        if strcmp(dataType, 'CBV') == false && strcmp(dataType,'CBV_HbT') == false
            Hip_remData{j,1} = detrend(Hip_remData{j,1}(1:(params.minTime.REM*samplingRate)),'constant');
        end
    end
    
    % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontally)
    LH_rem = zeros(length(LH_remData{1,1}),length(LH_remData));
    RH_rem = zeros(length(RH_remData{1,1}),length(RH_remData));
    if strcmp(dataType, 'CBV') == false && strcmp(dataType,'CBV_HbT') == false
        Hip_rem = zeros(length(Hip_remData{1,1}),length(Hip_remData));
    end
    for k = 1:length(LH_remData)
        LH_rem(:,k) = LH_remData{k,1};
        RH_rem(:,k) = RH_remData{k,1};
        if strcmp(dataType, 'CBV') == false && strcmp(dataType,'CBV_HbT') == false
            Hip_rem(:,k) = Hip_remData{k,1};
        end
    end
    
    % calculate the power spectra of the desired signals
    disp(['Analyzing the power spectrum of the LH REM ' dataType ' signal power...']); disp(' ')
    [LH_rem_S,LH_rem_f,LH_rem_sErr] = mtspectrumc_IOS(LH_remData,params);
    disp(['Analyzing the power spectrum of the RH REM ' dataType ' signal power...']); disp(' ')
    [RH_rem_S,RH_rem_f,RH_rem_sErr] = mtspectrumc_IOS(RH_remData,params);
    if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
        disp(['Analyzing the power spectrum of the Hippocampal REM ' dataType ' signal power...']); disp(' ')
        [Hip_rem_S,Hip_rem_f,Hip_rem_sErr] = mtspectrumc_IOS(Hip_remData,params);
    end
    
    % nboot = 1000;
    % LH_remCI = bootci(nboot,@mtspectrumc2,LH_remData',params);
    % RH_remCI = bootci(nboot,@mtspectrumc2,RH_remData',params);
    % if strcmp(dataType, 'CBV') == false && strcmp(dataType,'CBV_HbT') == false
    %     Hip_remCI = bootci(nboot,@mtspectrumc2,Hip_remData',params);
    % end
    
    % summary figures
    LH_remPower = figure;
    loglog(LH_rem_f,LH_rem_S,'k')
    hold on;
    loglog(LH_rem_f,LH_rem_sErr,'color',colors_IOS('battleship grey'))
    xlabel('Freq (Hz)');
    ylabel('Power');
    title([animalID  ' LH ' dataType ' Power during REM']);
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
    title([animalID  ' RH ' dataType ' Power during REM']);
    set(gca,'Ticklength',[0 0]);
    legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
    set(legend,'FontSize',6);
    xlim([0 1])
    axis square
    
    if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
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
    AnalysisResults.PowerSpectra.REM.(dataType).LH.S = LH_rem_S;
    AnalysisResults.PowerSpectra.REM.(dataType).LH.f = LH_rem_f;
    AnalysisResults.PowerSpectra.REM.(dataType).LH.sErr = LH_rem_sErr;
    AnalysisResults.PowerSpectra.REM.(dataType).RH.S = RH_rem_S;
    AnalysisResults.PowerSpectra.REM.(dataType).RH.f = RH_rem_f;
    AnalysisResults.PowerSpectra.REM.(dataType).RH.sErr = RH_rem_sErr;
    if strcmp(dataType, 'CBV') == false && strcmp(dataType,'CBV_HbT') == false
        AnalysisResults.PowerSpectra.REM.(dataType).Hip.S = Hip_rem_S;
        AnalysisResults.PowerSpectra.REM.(dataType).Hip.f = Hip_rem_f;
        AnalysisResults.PowerSpectra.REM.(dataType).Hip.sErr = Hip_rem_sErr;
    end
    
    % save figures
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Analysis Power Spectra/'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(LH_remPower, [dirpath animalID '_REM_LH_' filterSet '_' dataType '_PowerSpectra']);
    savefig(RH_remPower, [dirpath animalID '_REM_RH_' filterSet '_' dataType '_PowerSpectra']);
    if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
        savefig(Hip_remPower, [dirpath animalID '_REM_Hippocampal_' filterSet '_' dataType '_PowerSpectra']);
    end
    
    %% Analyze power spectra during all data
    for o = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(o,:);
        load(procDataFileID);
        [~,fileDate,~] = GetFileInfo_IOS(procDataFileID);
        AD_strDay = ConvertDate_IOS(fileDate);
        
        % pull data from each file
        if strcmp(dataType,'CBV') == true || strcmp(dataType,'CBV_HbT') == true
            if strcmp(dataType,'CBV') == true
                LH_AllData{o,1} = (ProcData.data.(dataType).LH - RestingBaselines.(baselineType).(dataType).LH.(AD_strDay))/RestingBaselines.(baselineType).(dataType).LH.(AD_strDay);
                RH_AllData{o,1} = (ProcData.data.(dataType).RH - RestingBaselines.(baselineType).(dataType).RH.(AD_strDay))/RestingBaselines.(baselineType).(dataType).RH.(AD_strDay);
            else
                LH_AllData{o,1} = ProcData.data.(dataType).LH;
                RH_AllData{o,1} = ProcData.data.(dataType).RH;
            end
        else
            LH_AllData{o,1} = (ProcData.data.cortical_LH.(dataType) - RestingBaselines.(baselineType).cortical_LH.(dataType).(AD_strDay))/RestingBaselines.(baselineType).cortical_LH.(dataType).(AD_strDay);
            RH_AllData{o,1} = (ProcData.data.cortical_RH.(dataType) - RestingBaselines.(baselineType).cortical_RH.(dataType).(AD_strDay))/RestingBaselines.(baselineType).cortical_RH.(dataType).(AD_strDay);
            Hip_AllData{o,1} = (ProcData.data.hippocampus.(dataType) - RestingBaselines.(baselineType).hippocampus.(dataType).(AD_strDay))/RestingBaselines.(baselineType).hippocampus.(dataType).(AD_strDay);
        end
    end
   
    % detend and lowpass filter each signal
    for p = 1:length(LH_AllData)
        LH_ProcAllData{p,1} = detrend(filtfilt(B,A,LH_AllData{p,1}),'constant');
        RH_ProcAllData{p,1} = detrend(filtfilt(B,A,RH_AllData{p,1}),'constant');
        if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
            Hip_ProcAllData{p,1} = detrend(filtfilt(B,A,Hip_AllData{p,1}),'constant');
        end
    end
    
    % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontally)
    LH_FinalAllData = zeros(length(LH_ProcAllData{1,1}),length(LH_ProcAllData));
    RH_FinalAllData = zeros(length(RH_ProcAllData{1,1}),length(RH_ProcAllData));
    if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
        Hip_FinalAllData = zeros(length(Hip_ProcAllData{1,1}),length(Hip_ProcAllData));
    end
    for q = 1:length(LH_ProcAllData)
        LH_FinalAllData(:,q) = LH_ProcAllData{q,1};
        RH_FinalAllData(:,q) = RH_ProcAllData{q,1};
        if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
            Hip_FinalAllData(:,q) = Hip_ProcAllData{q,1};
        end
    end
    
    % calculate the power spectra of the desired signals
    disp(['Analyzing the power spectrum of the LH all data ' dataType ' signal power...']); disp(' ')
    [LH_allData_S,LH_allData_f,LH_allData_sErr] = mtspectrumc_IOS(LH_FinalAllData,params);
    disp(['Analyzing the power spectrum of the RH all data ' dataType ' signal power...']); disp(' ')
    [RH_allData_S,RH_allData_f,RH_allData_sErr] = mtspectrumc_IOS(RH_FinalAllData,params);
    if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
        disp(['Analyzing the power spectrum of the Hippocampal all data ' dataType ' signal power...']); disp(' ')
        [Hip_allData_S,Hip_allData_f,Hip_allData_sErr] = mtspectrumc_IOS(Hip_FinalAllData,params);
    end
    
    % nboot = 1000;
    % LH_allDataCI = bootci(nboot,@mtspectrumc2,LH_FinalAllData',params);
    % RH_allDataCI = bootci(nboot,@mtspectrumc2,RH_FinalAllData',params);
    % if strcmp(dataType, 'CBV') == false && strcmp(dataType,'CBV_HbT') == false
    %     Hip_allDataCI = bootci(nboot,@mtspectrumc2,Hip_FinalAllData',params);
    % end
    
    % summary figures
    LH_allDataPower = figure;
    loglog(LH_allData_f,LH_allData_S,'k')
    hold on;
    loglog(LH_allData_f,LH_allData_sErr,'color',colors_IOS('battleship grey'))
    xlabel('Freq (Hz)');
    ylabel('Power');
    title([animalID  ' LH ' dataType ' Power during all data']);
    set(gca,'Ticklength',[0 0]);
    legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
    set(legend,'FontSize',6);
    xlim([0 1])
    axis square
    
    RH_allDataPower = figure;
    loglog(RH_allData_f,RH_allData_S,'k')
    hold on;
    loglog(RH_allData_f,RH_allData_sErr,'color',colors_IOS('battleship grey'))
    xlabel('Freq (Hz)');
    ylabel('Power');
    title([animalID  ' RH ' dataType ' Power during all data']);
    set(gca,'Ticklength',[0 0]);
    legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
    set(legend,'FontSize',6);
    xlim([0 1])
    axis square
    
    if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
        Hip_allDataPower = figure;
        loglog(Hip_allData_f,Hip_allData_S,'k')
        hold on;
        loglog(Hip_allData_f,Hip_allData_sErr,'color',colors_IOS('battleship grey'))
        xlabel('Freq (Hz)');
        ylabel('Power');
        title([animalID  ' Hippocampal ' dataType ' Power during all data']);
        set(gca,'Ticklength',[0 0]);
        legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
        set(legend,'FontSize',6);
        xlim([0 1])
        axis square
    end
    
    % save results
    AnalysisResults.PowerSpectra.AllData.(dataType).LH.S = LH_allData_S;
    AnalysisResults.PowerSpectra.AllData.(dataType).LH.f = LH_allData_f;
    AnalysisResults.PowerSpectra.AllData.(dataType).LH.sErr = LH_allData_sErr;
    AnalysisResults.PowerSpectra.AllData.(dataType).RH.S = RH_allData_S;
    AnalysisResults.PowerSpectra.AllData.(dataType).RH.f = RH_allData_f;
    AnalysisResults.PowerSpectra.AllData.(dataType).RH.sErr = RH_allData_sErr;
    if strcmp(dataType, 'CBV') == false && strcmp(dataType,'CBV_HbT') == false
        AnalysisResults.PowerSpectra.AllData.(dataType).Hip.S = Hip_allData_S;
        AnalysisResults.PowerSpectra.AllData.(dataType).Hip.f = Hip_allData_f;
        AnalysisResults.PowerSpectra.AllData.(dataType).Hip.sErr = Hip_allData_sErr;
    end
    
    % save figures
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Analysis Power Spectra/'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(LH_allDataPower, [dirpath animalID '_AllData_LH_' filterSet '_' dataType '_PowerSpectra']);
    savefig(RH_allDataPower, [dirpath animalID '_AllData_RH_' filterSet '_' dataType '_PowerSpectra']);
    if strcmp(dataType,'CBV') == false && strcmp(dataType,'CBV_HbT') == false
        savefig(Hip_allDataPower, [dirpath animalID '_AllData_Hippocampal_' filterSet '_' dataType '_PowerSpectra']);
    end
end

%% save results struct
save([animalID '_AnalysisResults.mat'], 'AnalysisResults');

end