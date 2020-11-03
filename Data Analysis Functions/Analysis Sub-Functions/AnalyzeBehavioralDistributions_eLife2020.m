function [AnalysisResults] = AnalyzeBehavioralDistributions_eLife2020(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the arousal-state distribution of different behavioral measurements (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
%% only run analysis for valid animal IDs
if any(strcmp(animalIDs,animalID))
    dataLocation = [rootFolder '/' animalID '/Bilateral Imaging/'];
    cd(dataLocation)
    % character list of all ProcData files
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
    % find and load RestData.mat struct
    forestModelID = 'Forest_ScoringResults.mat';
    load(forestModelID,'-mat')
    % find and load RestingBaselines.mat strut
    baselineDataFileStruct = dir('*_RestingBaselines.mat');
    baselineDataFile = {baselineDataFileStruct.name}';
    baselineDataFileID = char(baselineDataFile);
    load(baselineDataFileID)
    % lowpass filter
    samplingRate = 30;
    [z,p,k] = butter(4,10/(samplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    %% go through each file ID and create determine the average EMG
    awakeEMG = []; nremEMG = []; remEMG = [];
    awakeHR = []; nremHR = []; remHR = [];
    awakeWhisk = []; nremWhisk = []; remWhisk = [];
    for aa = 1:length(procDataFileIDs)
        procDataFileID = procDataFileIDs(aa,:);
        [~,fileDate,fileID] = GetFileInfo_IOS_eLife2020(procDataFileID);
        strDay = ConvertDate_IOS_eLife2020(fileDate);
        load(procDataFileID,'-mat')
        for bb = 1:length(ScoringResults.fileIDs)
            if strcmp(fileID,ScoringResults.fileIDs{bb,1}) == true
                labels = ScoringResults.labels{bb,1};
            end
        end
        % process data for subdividing
        EMG = filtfilt(sos,g,ProcData.data.EMG.emg - RestingBaselines.manualSelection.EMG.emg.(strDay));
        whiskers = ProcData.data.whiskerAngle;
        HR = ProcData.data.heartRate;
        for cc = 1:length(labels)
            label = labels{cc,1};
            if cc == 1
                startPoint = 1;
                endPoint = 5*samplingRate;
                startPoint2 = 1;
                endPoint2 = 5;
            else
                startPoint = (cc-1)*5*samplingRate + 1;
                endPoint = cc*5*samplingRate;
                startPoint2 = (cc-1)*5 + 1;
                endPoint2 = cc*5;
            end
            if strcmp(label,'Not Sleep') == true
                awakeEMG = vertcat(awakeEMG,EMG(startPoint:endPoint)); %#ok<*AGROW>
                awakeWhisk = vertcat(awakeWhisk,whiskers(startPoint:endPoint));
                awakeHR = vertcat(awakeHR,HR(startPoint2:endPoint2));
            elseif strcmp(label,'NREM Sleep') == true
                nremEMG = vertcat(nremEMG,EMG(startPoint:endPoint));
                nremWhisk = vertcat(nremWhisk,whiskers(startPoint:endPoint));
                nremHR = vertcat(nremHR,HR(startPoint2:endPoint2));
            elseif strcmp(label,'REM Sleep') == true
                remEMG = vertcat(remEMG,EMG(startPoint:endPoint));
                remWhisk = vertcat(remWhisk,whiskers(startPoint:endPoint));
                remHR = vertcat(remHR,HR(startPoint2:endPoint2));
            end
        end
    end
    % save results
    AnalysisResults.(animalID).BehaviorDistributions = [];
    AnalysisResults.(animalID).BehaviorDistributions.Awake.EMG = mean(awakeEMG,2);
    AnalysisResults.(animalID).BehaviorDistributions.Awake.Whisk = var(awakeWhisk,0,2);
    AnalysisResults.(animalID).BehaviorDistributions.Awake.HR = mean(awakeHR,2);
    AnalysisResults.(animalID).BehaviorDistributions.NREM.EMG = mean(nremEMG,2);
    AnalysisResults.(animalID).BehaviorDistributions.NREM.Whisk = var(nremWhisk,0,2);
    AnalysisResults.(animalID).BehaviorDistributions.NREM.HR = mean(nremHR,2);
    AnalysisResults.(animalID).BehaviorDistributions.REM.EMG = mean(remEMG,2);
    AnalysisResults.(animalID).BehaviorDistributions.REM.Whisk = var(remWhisk,0,2);
    AnalysisResults.(animalID).BehaviorDistributions.REM.HR = mean(remHR,2);
end
% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults')

end

