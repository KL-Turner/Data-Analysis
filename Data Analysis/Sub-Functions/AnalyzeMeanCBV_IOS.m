function [AnalysisResults] = AnalyzeMeanCBV_IOS(dataTypes,params,baselineType,AnalysisResults)
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

% list of unstim Procdata.mat files
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

% identify animal's ID and pull important infortmat
fileBreaks = strfind(restDataFileID, '_');
animalID = restDataFileID(1:fileBreaks(1)-1);
trialDuration_min = RestData.CBV.adjLH.trialDuration_sec/60;   % min
manualFileIDs = unique(RestingBaselines.manualSelection.baselineFileInfo.fileIDs);
fileTarget = params.targetMinutes/trialDuration_min;
filterSets = {'manualSelection','setDuration','entireDuration'};
samplingRate = RestData.CBV.adjLH.CBVCamSamplingRate;

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
        %% Analyze coherence during periods of rest
        % use the RestCriteria we specified earlier to find unstim resting events that are greater than the criteria
        if strcmp(dataType,'CBV') == true || strcmp(dataType,'CBV_HbT') == true
            [restLogical] = FilterEvents_IOS(RestData.(dataType).adjLH,RestCriteria);
            [puffLogical] = FilterEvents_IOS(RestData.(dataType).adjLH,PuffCriteria);
            combRestLogical = logical(restLogical.*puffLogical);
            restFiles = RestData.(dataType).adjLH.fileIDs(combRestLogical,:);
            if strcmp(dataType,'CBV') == true
                LH_RestingData = RestData.(dataType).adjLH.NormData(combRestLogical,:);
                RH_RestingData = RestData.(dataType).adjRH.NormData(combRestLogical,:);
            else
                LH_RestingData = RestData.(dataType).adjLH.data(combRestLogical,:);
                RH_RestingData = RestData.(dataType).adjRH.data(combRestLogical,:);
            end
        else
            break
        end
        
        % identify the unique days and the unique number of files from the list of unstim resting events
        restUniqueDays = GetUniqueDays_IOS(restFiles);
        restUniqueFiles = unique(restFiles);
        restNumberOfFiles = length(unique(restFiles));
        
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
        
        % extract unstim the resting events that correspond to the acceptable file list and the acceptable resting criteria
        clear restFileFilter
        filtRestFiles = restUniqueFiles(restFinalLogical,:);
        for f = 1:length(restFiles)
            restLogic = strcmp(restFiles{f},filtRestFiles);
            restLogicSum = sum(restLogic);
            if restLogicSum == 1
                restFileFilter(f,1) = 1;
            else
                restFileFilter(f,1) = 0;
            end
        end
        restFinalFileFilter = logical(restFileFilter);
        LH_finalRestData = LH_RestingData(restFinalFileFilter,:);
        RH_finalRestData = RH_RestingData(restFinalFileFilter,:);
        
        % only take the first 10 seconds of the epoch. occassionally a sample gets lost from rounding during the
        % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
        % lowpass filter and detrend each segment
        [B, A] = butter(3,1/(samplingRate/2),'low');
        clear LH_ProcRestData
        clear RH_ProcRestData
        for g = 1:length(LH_finalRestData)      
            LH_ProcRestData{g,1} = filtfilt(B,A,LH_finalRestData{g,1});
            RH_ProcRestData{g,1} = filtfilt(B,A,RH_finalRestData{g,1});
        end
        
        % analyze correlation coefficient between resting epochs
        for n = 1:length(LH_ProcRestData)
            LH_restCBVMean(n,1) = mean(LH_ProcRestData{n,1});
            RH_restCBVMean(n,1) = mean(RH_ProcRestData{n,1});
        end
        
        % save results
        AnalysisResults.MeanCBV.Rest.(dataType).(filterSet).LH = LH_restCBVMean;
        AnalysisResults.MeanCBV.Rest.(dataType).(filterSet).RH = RH_restCBVMean;
    end
    
    %% Analyze coherence during periods of NREM sleep
    % pull data from SleepData.mat structure
    if strcmp(dataType,'CBV') == true || strcmp(dataType,'CBV_HbT') == true
        LH_nremData = SleepData.NREM.data.(dataType).LH;
        RH_nremData = SleepData.NREM.data.(dataType).RH;
    else
        break
    end
  
    % analyze correlation coefficient between NREM epochs
    for n = 1:length(LH_nremData)
        LH_nremCBVMean(n,1) = mean(LH_nremData{n,1});
        RH_nremCBVMean(n,1) = mean(RH_nremData{n,1});
    end
    
    % save results
    AnalysisResults.MeanCBV.NREM.(dataType).LH = LH_nremCBVMean;
    AnalysisResults.MeanCBV.NREM.(dataType).RH = RH_nremCBVMean;
    
    %% Analyze coherence during periods of REM sleep
    % pull data from SleepData.mat structure
    if strcmp(dataType,'CBV') == true || strcmp(dataType,'CBV_HbT') == true
        LH_remData = SleepData.REM.data.(dataType).LH;
        RH_remData = SleepData.REM.data.(dataType).RH;
    else
        break
    end
    
    % analyze correlation coefficient between NREM epochs
    for n = 1:length(LH_remData)
        LH_remCBVMean(n,1) = mean(LH_remData{n,1});
        RH_remCBVMean(n,1) = mean(RH_remData{n,1});
    end

    % save results
    AnalysisResults.MeanCBV.REM.(dataType).LH = LH_remCBVMean;
    AnalysisResults.MeanCBV.REM.(dataType).RH = RH_remCBVMean;
    
    %% Analyze coherence during unstim unstimulated data
    for o = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(o,:);
        load(procDataFileID);
        if isempty(ProcData.data.solenoids.LPadSol) == true
            stimLogical(o,1) = 1;
        else
            stimLogical(o,1) = 0;
        end
    end
    stimLogical = logical(stimLogical);
    unstim_procDataFileIDs = procDataFileIDs(stimLogical,:);
    for p = 1:size(unstim_procDataFileIDs)
        unstim_procDataFileID = unstim_procDataFileIDs(p,:);
        load(unstim_procDataFileID)
        [~,fileDate,~] = GetFileInfo_IOS(unstim_procDataFileID);
        US_strDay = ConvertDate_IOS(fileDate);
        % pull data from each file
        if strcmp(dataType,'CBV') == true || strcmp(dataType,'CBV_HbT') == true
            if strcmp(dataType,'CBV') == true
                LH_UnstimData{p,1} = (ProcData.data.(dataType).adjLH - RestingBaselines.(baselineType).(dataType).adjLH.(US_strDay))/RestingBaselines.(baselineType).(dataType).adjLH.(US_strDay);
                RH_UnstimData{p,1} = (ProcData.data.(dataType).adjRH - RestingBaselines.(baselineType).(dataType).adjRH.(US_strDay))/RestingBaselines.(baselineType).(dataType).adjRH.(US_strDay);
            else
                LH_UnstimData{p,1} = ProcData.data.(dataType).adjLH;
                RH_UnstimData{p,1} = ProcData.data.(dataType).adjRH;
            end
            break
        end
    end
    
    % detend and lowpass filter each signal
    for q = 1:length(LH_UnstimData)
        LH_ProcUnstimData{q,1} = filtfilt(B,A,LH_UnstimData{q,1});
        RH_ProcUnstimData{q,1} = filtfilt(B,A,RH_UnstimData{q,1});
    end
    
    % analyze correlation coefficient between NREM epochs
    for n = 1:length(LH_ProcUnstimData)
        LH_unstimCBVMean(n,1) = mean(LH_ProcUnstimData{n,1});
        RH_unstimCBVMean(n,1) = mean(RH_ProcUnstimData{n,1});
    end

    % save results
    AnalysisResults.MeanCBV.Unstim.(dataType).LH = LH_unstimCBVMean;
    AnalysisResults.MeanCBV.Unstim.(dataType).RH = RH_unstimCBVMean;
    
    %% Analyze coherence during all data
    for p = 1:size(procDataFileIDs)
        procDataFileID = procDataFileIDs(p,:);
        load(procDataFileID)
        [~,fileDate,~] = GetFileInfo_IOS(procDataFileID);
        AD_strDay = ConvertDate_IOS(fileDate);
        % pull data from each file
        if strcmp(dataType,'CBV') == true || strcmp(dataType,'CBV_HbT') == true
            if strcmp(dataType,'CBV') == true
                LH_AllData{p,1} = (ProcData.data.(dataType).adjLH - RestingBaselines.(baselineType).(dataType).adjLH.(AD_strDay))/RestingBaselines.(baselineType).(dataType).adjLH.(AD_strDay);
                RH_AllData{p,1} = (ProcData.data.(dataType).adjRH - RestingBaselines.(baselineType).(dataType).adjRH.(AD_strDay))/RestingBaselines.(baselineType).(dataType).adjRH.(AD_strDay);
            else
                LH_AllData{p,1} = ProcData.data.(dataType).adjLH;
                RH_AllData{p,1} = ProcData.data.(dataType).adjRH;
            end
        else
            break
        end
    end
    
    % detend and lowpass filter each signal
    for q = 1:length(LH_AllData)
        LH_ProcAllData{q,1} = filtfilt(B,A,LH_AllData{q,1});
        RH_ProcAllData{q,1} = filtfilt(B,A,RH_AllData{q,1});
    end
    
    % analyze correlation coefficient between NREM epochs
    for n = 1:length(LH_ProcAllData)
        LH_allCBVMean(n,1) = mean(LH_ProcAllData{n,1});
        RH_allCBVMean(n,1) = mean(RH_ProcAllData{n,1});
    end

    % save results
    AnalysisResults.MeanCBV.All.(dataType).LH = LH_allCBVMean;
    AnalysisResults.MeanCBV.All.(dataType).RH = RH_allCBVMean;
end

%% save results strucure
save([animalID '_AnalysisResults.mat'], 'AnalysisResults');

end
