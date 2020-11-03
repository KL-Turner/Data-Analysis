function [AnalysisResults] = AnalyzeCoherence_eLife2020(animalID,saveFigs,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the spectral coherence between bilateral hemodynamic [HbT] and neural signals (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
dataTypes = {'CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
modelType = 'Forest';
params.minTime.Rest = 10;
params.minTime.NREM = 30;
params.minTime.REM = 60;
%% only run analysis for valid animal IDs
if any(strcmp(animalIDs,animalID))
    dataLocation = [rootFolder '/' animalID '/Bilateral Imaging/'];
    cd(dataLocation)
    % character list of all ProcData file IDs
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
    % find and load RestData.mat struct
    restDataFileStruct = dir('*_RestData.mat');
    restDataFile = {restDataFileStruct.name}';
    restDataFileID = char(restDataFile);
    load(restDataFileID)
    % find and load manual baseline event information
    manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
    manualBaselineFile = {manualBaselineFileStruct.name}';
    manualBaselineFileID = char(manualBaselineFile);
    load(manualBaselineFileID)
    % find and load RestingBaselines.mat struct
    baselineDataFileStruct = dir('*_RestingBaselines.mat');
    baselineDataFile = {baselineDataFileStruct.name}';
    baselineDataFileID = char(baselineDataFile);
    load(baselineDataFileID)
    % find and load SleepData.mat struct
    sleepDataFileStruct = dir('*_SleepData.mat');
    sleepDataFile = {sleepDataFileStruct.name}';
    sleepDataFileID = char(sleepDataFile);
    load(sleepDataFileID)
    % find and load Forest_ScoringResults.mat struct
    forestScoringResultsFileID = 'Forest_ScoringResults.mat';
    load(forestScoringResultsFileID,'-mat')
    % lowpass filter
    samplingRate = RestData.CBV_HbT.adjLH.CBVCamSamplingRate;
    [z,p,k] = butter(4,1/(samplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    % criteria for resting
    RestCriteria.Fieldname = {'durations'};
    RestCriteria.Comparison = {'gt'};
    RestCriteria.Value = {params.minTime.Rest};
    RestPuffCriteria.Fieldname = {'puffDistances'};
    RestPuffCriteria.Comparison = {'gt'};
    RestPuffCriteria.Value = {5};
    % go through each valid data type for arousal-based coherence analysis
    for aa = 1:length(dataTypes)
        dataType = dataTypes{1,aa};
        %% analyze bilateral coherence during periods of rest
        % pull data from RestData.mat structure
        if strcmp(dataType,'CBV_HbT') == true
            [restLogical] = FilterEvents_IOS_eLife2020(RestData.(dataType).adjLH,RestCriteria);
            [puffLogical] = FilterEvents_IOS_eLife2020(RestData.(dataType).adjLH,RestPuffCriteria);
            combRestLogical = logical(restLogical.*puffLogical);
            restFileIDs = RestData.(dataType).adjLH.fileIDs(combRestLogical,:);
            restEventTimes = RestData.(dataType).adjLH.eventTimes(combRestLogical,:);
            restDurations = RestData.(dataType).adjLH.durations(combRestLogical,:);
            LH_unstimRestingData = RestData.(dataType).adjLH.data(combRestLogical,:);
            RH_unstimRestingData = RestData.(dataType).adjRH.data(combRestLogical,:);
        else
            [restLogical] = FilterEvents_IOS_eLife2020(RestData.cortical_LH.(dataType),RestCriteria);
            [puffLogical] = FilterEvents_IOS_eLife2020(RestData.cortical_LH.(dataType),RestPuffCriteria);
            combRestLogical = logical(restLogical.*puffLogical);
            restFileIDs = RestData.cortical_LH.(dataType).fileIDs(combRestLogical,:);
            restEventTimes = RestData.cortical_LH.(dataType).eventTimes(combRestLogical,:);
            restDurations = RestData.cortical_LH.(dataType).durations(combRestLogical,:);
            LH_unstimRestingData = RestData.cortical_LH.(dataType).NormData(combRestLogical,:);
            RH_unstimRestingData = RestData.cortical_RH.(dataType).NormData(combRestLogical,:);
        end
        % keep only the data that occurs within the manually-approved awake regions
        [LH_finalRestData,~,~,~] = RemoveInvalidData_IOS_eLife2020(LH_unstimRestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [RH_finalRestData,~,~,~] = RemoveInvalidData_IOS_eLife2020(RH_unstimRestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        clear LH_ProcRestData RH_ProcRestData
        % filter, detrend, and truncate data to minimum length to match events
        for bb = 1:length(LH_finalRestData)
            if length(LH_finalRestData{bb,1}) < params.minTime.Rest*samplingRate
                restChunkSampleDiff = params.minTime.Rest*samplingRate - length(LH_finalRestData{bb,1});
                LH_restPad = (ones(1,restChunkSampleDiff))*LH_finalRestData{bb,1}(end);
                RH_restPad = (ones(1,restChunkSampleDiff))*RH_finalRestData{bb,1}(end);
                LH_ProcRestData{bb,1} = horzcat(LH_finalRestData{bb,1},LH_restPad); %#ok<*AGROW>
                RH_ProcRestData{bb,1} = horzcat(RH_finalRestData{bb,1},RH_restPad);
                LH_ProcRestData{bb,1} = filtfilt(sos,g,detrend(LH_ProcRestData{bb,1},'constant'));
                RH_ProcRestData{bb,1} = filtfilt(sos,g,detrend(RH_ProcRestData{bb,1},'constant'));
            else
                LH_ProcRestData{bb,1} = filtfilt(sos,g,detrend(LH_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant'));
                RH_ProcRestData{bb,1} = filtfilt(sos,g,detrend(RH_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant'));
            end
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
        LH_restData = zeros(length(LH_ProcRestData{1,1}),length(LH_ProcRestData));
        RH_restData = zeros(length(RH_ProcRestData{1,1}),length(RH_ProcRestData));
        for cc = 1:length(LH_ProcRestData)
            LH_restData(:,cc) = LH_ProcRestData{cc,1};
            RH_restData(:,cc) = RH_ProcRestData{cc,1};
        end
        % parameters for coherencyc - information available in function
        params.tapers = [1,1];   % Tapers [n, 2n - 1]
        params.pad = 1;
        params.Fs = samplingRate;
        params.fpass = [0,0.5];   % Pass band [0, nyquist]
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the coherence between desired signals
        [C_RestData,~,~,~,~,f_RestData,confC_RestData,~,cErr_RestData] = coherencyc_eLife2020(LH_restData,RH_restData,params);
        % save results
        AnalysisResults.(animalID).Coherence.Rest.(dataType).C = C_RestData;
        AnalysisResults.(animalID).Coherence.Rest.(dataType).f = f_RestData;
        AnalysisResults.(animalID).Coherence.Rest.(dataType).confC = confC_RestData;
        AnalysisResults.(animalID).Coherence.Rest.(dataType).cErr = cErr_RestData;
        % save figures if desired
        if strcmp(saveFigs,'y') == true
            restCoherence = figure;
            semilogx(f_RestData,C_RestData,'k')
            hold on;
            semilogx(f_RestData,cErr_RestData,'color',colors_eLife2020('battleship grey'))
            xlabel('Freq (Hz)');
            ylabel('Coherence');
            title([animalID  ' ' dataType ' coherence for resting data']);
            set(gca,'Ticklength',[0,0]);
            legend('Coherence','Jackknife Lower','Jackknife Upper','Location','Southeast');
            set(legend,'FontSize',6);
            ylim([0,1])
            xlim([0,0.5])
            axis square
            set(gca,'box','off')
            [pathstr,~,~] = fileparts(cd);
            dirpath = [pathstr '/Figures/Coherence/'];
            if ~exist(dirpath,'dir')
                mkdir(dirpath);
            end
            savefig(restCoherence,[dirpath animalID '_Rest_' dataType '_Coherence']);
            close(restCoherence)
        end
        %% analyze bilateral coherence during periods of alert
        zz = 1;
        clear LH_AwakeData RH_AwakeData LH_ProcAwakeData RH_ProcAwakeData
        LH_AwakeData = [];
        for bb = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(bb,:);
            [~,allDataFileDate,allDataFileID] = GetFileInfo_IOS_eLife2020(procDataFileID);
            strDay = ConvertDate_IOS_eLife2020(allDataFileDate);
            scoringLabels = [];
            for cc = 1:length(ScoringResults.fileIDs)
                if strcmp(allDataFileID,ScoringResults.fileIDs{cc,1}) == true
                    scoringLabels = ScoringResults.labels{cc,1};
                end
            end
            % check labels to match arousal state
            if sum(strcmp(scoringLabels,'Not Sleep')) > 144   % 36 bins (180 total) or 3 minutes of sleep
                load(procDataFileID)
                puffs = ProcData.data.solenoids.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    if strcmp(dataType,'CBV_HbT') == true
                        LH_AwakeData{zz,1} = ProcData.data.(dataType).adjLH;
                        RH_AwakeData{zz,1} = ProcData.data.(dataType).adjRH;
                    else
                        LH_AwakeData{zz,1} = (ProcData.data.cortical_LH.(dataType) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay);
                        RH_AwakeData{zz,1} = (ProcData.data.cortical_RH.(dataType) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay);
                    end
                    zz = zz + 1;
                end
            end
        end
        % filter and detrend data
        if isempty(LH_AwakeData) == false
            for bb = 1:length(LH_AwakeData)
                LH_ProcAwakeData{bb,1} = filtfilt(sos,g,detrend(LH_AwakeData{bb,1},'constant'));
                RH_ProcAwakeData{bb,1} = filtfilt(sos,g,detrend(RH_AwakeData{bb,1},'constant'));
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            LH_awakeData = zeros(length(LH_ProcAwakeData{1,1}),length(LH_ProcAwakeData));
            RH_awakeData = zeros(length(RH_ProcAwakeData{1,1}),length(RH_ProcAwakeData));
            for cc = 1:length(LH_ProcAwakeData)
                LH_awakeData(:,cc) = LH_ProcAwakeData{cc,1};
                RH_awakeData(:,cc) = RH_ProcAwakeData{cc,1};
            end
            % parameters for coherencyc - information available in function
            params.tapers = [10,19];   % Tapers [n, 2n - 1]
            params.pad = 1;
            params.Fs = samplingRate;
            params.fpass = [0,0.5];   % Pass band [0, nyquist]
            params.trialave = 1;
            params.err = [2,0.05];
            % calculate the coherence between desired signals
            [C_AwakeData,~,~,~,~,f_AwakeData,confC_AwakeData,~,cErr_AwakeData] = coherencyc_eLife2020(LH_awakeData,RH_awakeData,params);
            % save results
            AnalysisResults.(animalID).Coherence.Awake.(dataType).C = C_AwakeData;
            AnalysisResults.(animalID).Coherence.Awake.(dataType).f = f_AwakeData;
            AnalysisResults.(animalID).Coherence.Awake.(dataType).confC = confC_AwakeData;
            AnalysisResults.(animalID).Coherence.Awake.(dataType).cErr = cErr_AwakeData;
            % save figures if desired
            if strcmp(saveFigs,'y') == true
                awakeCoherence = figure;
                semilogx(f_AwakeData,C_AwakeData,'k')
                hold on;
                semilogx(f_AwakeData,cErr_AwakeData,'color',colors_eLife2020('battleship grey'))
                xlabel('Freq (Hz)');
                ylabel('Coherence');
                title([animalID  ' ' dataType ' coherence for awake data']);
                set(gca,'Ticklength',[0,0]);
                legend('Coherence','Jackknife Lower','Jackknife Upper','Location','Southeast');
                set(legend,'FontSize',6);
                ylim([0,1])
                xlim([0,0.5])
                axis square
                set(gca,'box','off')
                [pathstr,~,~] = fileparts(cd);
                dirpath = [pathstr '/Figures/Coherence/'];
                if ~exist(dirpath,'dir')
                    mkdir(dirpath);
                end
                savefig(awakeCoherence,[dirpath animalID '_Awake_' dataType '_Coherence']);
                close(awakeCoherence)
            end
        else
            % save results
            AnalysisResults.(animalID).Coherence.Awake.(dataType).C = [];
            AnalysisResults.(animalID).Coherence.Awake.(dataType).f = [];
            AnalysisResults.(animalID).Coherence.Awake.(dataType).confC = [];
            AnalysisResults.(animalID).Coherence.Awake.(dataType).cErr = [];
        end
        %% analyze bilateral coherence during periods of asleep
        zz = 1;
        clear LH_SleepData RH_SleepData LH_ProcSleepData RH_ProcSleepData
        LH_SleepData = [];
        for bb = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(bb,:);
            [~,allDataFileDate,allDataFileID] = GetFileInfo_IOS_eLife2020(procDataFileID);
            strDay = ConvertDate_IOS_eLife2020(allDataFileDate);
            scoringLabels = [];
            for cc = 1:length(ScoringResults.fileIDs)
                if strcmp(allDataFileID,ScoringResults.fileIDs{cc,1}) == true
                    scoringLabels = ScoringResults.labels{cc,1};
                end
            end
            % check labels to match arousal state
            if sum(strcmp(scoringLabels,'Not Sleep')) < 36   % 36 bins (180 total) or 3 minutes of awake
                load(procDataFileID)
                puffs = ProcData.data.solenoids.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    if strcmp(dataType,'CBV_HbT') == true
                        LH_SleepData{zz,1} = ProcData.data.(dataType).adjLH;
                        RH_SleepData{zz,1} = ProcData.data.(dataType).adjRH;
                    else
                        LH_SleepData{zz,1} = (ProcData.data.cortical_LH.(dataType) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay);
                        RH_SleepData{zz,1} = (ProcData.data.cortical_RH.(dataType) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay);
                    end
                    zz = zz + 1;
                end
            end
        end
        % ilter and detrend data
        if isempty(LH_SleepData) == false
            for bb = 1:length(LH_SleepData)
                LH_ProcSleepData{bb,1} = filtfilt(sos,g,detrend(LH_SleepData{bb,1},'constant'));
                RH_ProcSleepData{bb,1} = filtfilt(sos,g,detrend(RH_SleepData{bb,1},'constant'));
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            LH_sleepData = zeros(length(LH_ProcSleepData{1,1}),length(LH_ProcSleepData));
            RH_sleepData = zeros(length(RH_ProcSleepData{1,1}),length(RH_ProcSleepData));
            for cc = 1:length(LH_ProcSleepData)
                LH_sleepData(:,cc) = LH_ProcSleepData{cc,1};
                RH_sleepData(:,cc) = RH_ProcSleepData{cc,1};
            end
            % parameters for coherencyc - information available in function
            params.tapers = [10,19];   % Tapers [n, 2n - 1]
            params.pad = 1;
            params.Fs = samplingRate;
            params.fpass = [0,0.5];   % Pass band [0, nyquist]
            params.trialave = 1;
            params.err = [2,0.05];
            % calculate the coherence between desired signals
            [C_SleepData,~,~,~,~,f_SleepData,confC_SleepData,~,cErr_SleepData] = coherencyc_eLife2020(LH_sleepData,RH_sleepData,params);
            % save results
            AnalysisResults.(animalID).Coherence.Sleep.(dataType).C = C_SleepData;
            AnalysisResults.(animalID).Coherence.Sleep.(dataType).f = f_SleepData;
            AnalysisResults.(animalID).Coherence.Sleep.(dataType).confC = confC_SleepData;
            AnalysisResults.(animalID).Coherence.Sleep.(dataType).cErr = cErr_SleepData;
            % save figures if desired
            if strcmp(saveFigs,'y') == true
                sleepCoherence = figure;
                semilogx(f_SleepData,C_SleepData,'k')
                hold on;
                semilogx(f_SleepData,cErr_SleepData,'color',colors_eLife2020('battleship grey'))
                xlabel('Freq (Hz)');
                ylabel('Coherence');
                title([animalID  ' ' dataType ' coherence for sleep data']);
                set(gca,'Ticklength',[0,0]);
                legend('Coherence','Jackknife Lower','Jackknife Upper','Location','Southeast');
                set(legend,'FontSize',6);
                ylim([0,1])
                xlim([0,0.5])
                axis square
                set(gca,'box','off')
                [pathstr,~,~] = fileparts(cd);
                dirpath = [pathstr '/Figures/Coherence/'];
                if ~exist(dirpath,'dir')
                    mkdir(dirpath);
                end
                savefig(sleepCoherence,[dirpath animalID '_Sleep_' dataType '_Coherence']);
                close(sleepCoherence)
            end
        else
            % save results
            AnalysisResults.(animalID).Coherence.Sleep.(dataType).C = [];
            AnalysisResults.(animalID).Coherence.Sleep.(dataType).f = [];
            AnalysisResults.(animalID).Coherence.Sleep.(dataType).confC = [];
            AnalysisResults.(animalID).Coherence.Sleep.(dataType).cErr = [];
        end
        %% analyze bilateral coherence during periods of all data
        zz = 1;
        clear LH_AllUnstimData RH_AllUnstimData LH_ProcAllUnstimData RH_ProcAllUnstimData
        LH_AllUnstimData = [];
        for bb = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(bb,:);
            [~,allUnstimDataFileDate,~] = GetFileInfo_IOS_eLife2020(procDataFileID);
            strDay = ConvertDate_IOS_eLife2020(allUnstimDataFileDate);
            load(procDataFileID)
            puffs = ProcData.data.solenoids.LPadSol;
            % don't include trials with stimulation
            if isempty(puffs) == true
                if strcmp(dataType,'CBV_HbT') == true
                    LH_AllUnstimData{zz,1} = ProcData.data.(dataType).adjLH;
                    RH_AllUnstimData{zz,1} = ProcData.data.(dataType).adjRH;
                else
                    LH_AllUnstimData{zz,1} = (ProcData.data.cortical_LH.(dataType) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay);
                    RH_AllUnstimData{zz,1} = (ProcData.data.cortical_RH.(dataType) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay);
                end
                zz = zz + 1;
            end
        end
        % filter and detrend data
        if isempty(LH_AllUnstimData) == false
            for bb = 1:length(LH_AllUnstimData)
                LH_ProcAllUnstimData{bb,1} = filtfilt(sos,g,detrend(LH_AllUnstimData{bb,1},'constant'));
                RH_ProcAllUnstimData{bb,1} = filtfilt(sos,g,detrend(RH_AllUnstimData{bb,1},'constant'));
            end
            % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            LH_allUnstimData = zeros(length(LH_ProcAllUnstimData{1,1}),length(LH_ProcAllUnstimData));
            RH_allUnstimData = zeros(length(RH_ProcAllUnstimData{1,1}),length(RH_ProcAllUnstimData));
            for cc = 1:length(LH_ProcAllUnstimData)
                LH_allUnstimData(:,cc) = LH_ProcAllUnstimData{cc,1};
                RH_allUnstimData(:,cc) = RH_ProcAllUnstimData{cc,1};
            end
            % parameters for coherencyc - information available in function
            params.tapers = [10,19];   % Tapers [n, 2n - 1]
            params.pad = 1;
            params.Fs = samplingRate;   % Sampling Rate
            params.fpass = [0,0.5];   % Pass band [0, nyquist]
            params.trialave = 1;
            params.err = [2,0.05];
            % calculate the coherence between desired signals
            [C_AllUnstimData,~,~,~,~,f_AllUnstimData,confC_AllUnstimData,~,cErr_AllUnstimData] = coherencyc_eLife2020(LH_allUnstimData,RH_allUnstimData,params);
            % save results
            AnalysisResults.(animalID).Coherence.All.(dataType).C = C_AllUnstimData;
            AnalysisResults.(animalID).Coherence.All.(dataType).f = f_AllUnstimData;
            AnalysisResults.(animalID).Coherence.All.(dataType).confC = confC_AllUnstimData;
            AnalysisResults.(animalID).Coherence.All.(dataType).cErr = cErr_AllUnstimData;
            % save figures if desired
            if strcmp(saveFigs,'y') == true
                allUnstimCoherence = figure;
                semilogx(f_AllUnstimData,C_AllUnstimData,'k')
                hold on;
                semilogx(f_AllUnstimData,cErr_AllUnstimData,'color',colors_eLife2020('battleship grey'))
                xlabel('Freq (Hz)');
                ylabel('Coherence');
                title([animalID  ' ' dataType ' coherence for all unstim data']);
                set(gca,'Ticklength',[0,0]);
                legend('Coherence','Jackknife Lower','Jackknife Upper','Location','Southeast');
                set(legend,'FontSize',6);
                ylim([0,1])
                xlim([0,0.5])
                axis square
                set(gca,'box','off')
                [pathstr,~,~] = fileparts(cd);
                dirpath = [pathstr '/Figures/Coherence/'];
                if ~exist(dirpath,'dir')
                    mkdir(dirpath);
                end
                savefig(allUnstimCoherence,[dirpath animalID '_AllUnstim_' dataType '_Coherence']);
                close(allUnstimCoherence)
            end
        end
        %% analyze bilateral coherence during periods of NREM sleep
        % pull data from SleepData.mat structure
        if strcmp(dataType,'CBV_HbT') == true
            [LH_nremData,~,~] = RemoveStimSleepData_IOS_eLife2020(animalID,SleepData.(modelType).NREM.data.(dataType).LH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            [RH_nremData,~,~] = RemoveStimSleepData_IOS_eLife2020(animalID,SleepData.(modelType).NREM.data.(dataType).RH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        else
            [LH_nremData,~,~] = RemoveStimSleepData_IOS_eLife2020(animalID,SleepData.(modelType).NREM.data.cortical_LH.(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            [RH_nremData,~,~] = RemoveStimSleepData_IOS_eLife2020(animalID,SleepData.(modelType).NREM.data.cortical_RH.(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        end
        % filter, detrend, and truncate data to minimum length to match events
        for ee = 1:length(LH_nremData)
            LH_nremData{ee,1} = filtfilt(sos,g,detrend(LH_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant'));
            RH_nremData{ee,1} = filtfilt(sos,g,detrend(RH_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant'));
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        LH_nrem = zeros(length(LH_nremData{1,1}),length(LH_nremData));
        RH_nrem = zeros(length(RH_nremData{1,1}),length(RH_nremData));
        for ff = 1:length(LH_nremData)
            LH_nrem(:,ff) = LH_nremData{ff,1};
            RH_nrem(:,ff) = RH_nremData{ff,1};
        end
        % parameters for coherencyc - information available in function
        params.tapers = [3,5];   % Tapers [n, 2n - 1]
        params.pad = 1;
        params.Fs = samplingRate;
        params.fpass = [0,0.5];   % Pass band [0, nyquist]
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the coherence between desired signals
        [C_nrem,~,~,~,~,f_nrem,confC_nrem,~,cErr_nrem] = coherencyc_eLife2020(LH_nrem,RH_nrem,params);
        % save results
        AnalysisResults.(animalID).Coherence.NREM.(dataType).C = C_nrem;
        AnalysisResults.(animalID).Coherence.NREM.(dataType).f = f_nrem;
        AnalysisResults.(animalID).Coherence.NREM.(dataType).confC = confC_nrem;
        AnalysisResults.(animalID).Coherence.NREM.(dataType).cErr = cErr_nrem;
        % save figures if desired
        if strcmp(saveFigs,'y') == true
            nremCoherence = figure;
            semilogx(f_nrem,C_nrem,'k')
            hold on;
            semilogx(f_nrem,cErr_nrem,'color',colors_eLife2020('battleship grey'))
            xlabel('Freq (Hz)');
            ylabel('Coherence');
            title([animalID  ' ' dataType ' coherence for ' modelType ' NREM data']);
            set(gca,'Ticklength',[0,0]);
            legend('Coherence','Jackknife Lower','Jackknife Upper','Location','Southeast');
            set(legend,'FontSize',6);
            ylim([0.1,0.5])
            xlim([0,0.5])
            axis square
            set(gca,'box','off')
            savefig(nremCoherence,[dirpath animalID '_' modelType '_NREM_' dataType '_Coherence']);
            close(nremCoherence)
        end
        %% analyze bilateral coherence during periods of REM sleep
        % pull data from SleepData.mat structure
        if strcmp(dataType,'CBV_HbT') == true
            [LH_remData,~,~] = RemoveStimSleepData_IOS_eLife2020(animalID,SleepData.(modelType).REM.data.(dataType).LH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            [RH_remData,~,~] = RemoveStimSleepData_IOS_eLife2020(animalID,SleepData.(modelType).REM.data.(dataType).RH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        else
            [LH_remData,~,~] = RemoveStimSleepData_IOS_eLife2020(animalID,SleepData.(modelType).REM.data.cortical_LH.(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            [RH_remData,~,~] = RemoveStimSleepData_IOS_eLife2020(animalID,SleepData.(modelType).REM.data.cortical_RH.(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        end
        % filter, detrend, and truncate data to minimum length to match events
        for gg = 1:length(LH_remData)
            LH_remData{gg,1} = filtfilt(sos,g,detrend(LH_remData{gg,1}(1:(params.minTime.REM*samplingRate)),'constant'));
            RH_remData{gg,1} = filtfilt(sos,g,detrend(RH_remData{gg,1}(1:(params.minTime.REM*samplingRate)),'constant'));
        end
        % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        LH_rem = zeros(length(LH_remData{1,1}),length(LH_remData));
        RH_rem = zeros(length(RH_remData{1,1}),length(RH_remData));
        for hh = 1:length(LH_remData)
            LH_rem(:,hh) = LH_remData{hh,1};
            RH_rem(:,hh) = RH_remData{hh,1};
        end
        % parameters for coherencyc - information available in function
        params.tapers = [5,9];   % Tapers [n, 2n - 1]
        params.pad = 1;
        params.Fs = samplingRate;
        params.fpass = [0,0.5];   % Pass band [0, nyquist]
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the coherence between desired signals
        [C_rem,~,~,~,~,f_rem,confC_rem,~,cErr_rem] = coherencyc_eLife2020(LH_rem,RH_rem,params);
        % save results
        AnalysisResults.(animalID).Coherence.REM.(dataType).C = C_rem;
        AnalysisResults.(animalID).Coherence.REM.(dataType).f = f_rem;
        AnalysisResults.(animalID).Coherence.REM.(dataType).confC = confC_rem;
        AnalysisResults.(animalID).Coherence.REM.(dataType).cErr = cErr_rem;
        % save figures if desired
        if strcmp(saveFigs,'y') == true
            remCoherence = figure;
            semilogx(f_rem,C_rem,'k')
            hold on;
            semilogx(f_rem,cErr_rem,'color',colors_eLife2020('battleship grey'))
            xlabel('Freq (Hz)');
            ylabel('Coherence');
            title([animalID  ' ' dataType ' coherence for ' modelType 'REM data']);
            set(gca,'Ticklength',[0,0]);
            legend('Coherence','Jackknife Lower','Jackknife Upper','Location','Southeast');
            set(legend,'FontSize',6);
            ylim([0.1,0.5])
            xlim([0,0.5])
            axis square
            set(gca,'box','off')
            savefig(remCoherence,[dirpath animalID '_' modelType '_REM_' dataType '_Coherence']);
            close(remCoherence)
        end
    end
    % save data
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end

end
