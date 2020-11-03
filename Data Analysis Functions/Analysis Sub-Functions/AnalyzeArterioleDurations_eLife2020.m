function [AnalysisResults] = AnalyzeArterioleDurations_eLife2020(TwoP_animalIDs,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the time of each arousal-state data per artery (2PLSM)
%________________________________________________________________________________________________________________________

%% extract data from each animal's sleep scoring results
for aa = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,aa};
    dataLoc = [rootFolder '/' animalID '/2P Data/'];
    cd(dataLoc)
    % find and load RestingBaselines.mat strut
    baselineDataFileStruct = dir('*_RestingBaselines.mat');
    baselineDataFile = {baselineDataFileStruct.name}';
    baselineDataFileID = char(baselineDataFile);
    load(baselineDataFileID)
    % character list of all merged data file IDs
    mergedDirectory = dir('*_MergedData.mat');
    mergedDataFiles = {mergedDirectory.name}';
    mergedDataFileIDs = char(mergedDataFiles);
    % take out logicals from each vessel to determine how much data it has
    data.(animalID) = [];
    for bb = 1:size(mergedDataFileIDs,1)
        [~,~,~,~,~,vID] = GetFileInfo2_2P_eLife2020(mergedDataFileIDs(bb,:));
        if strcmp(vID(1),'V') == false
            load(mergedDataFileIDs(bb,:),'-mat')
            if isfield(data.(animalID),vID) == false
                data.(animalID).(vID).awakeData = sum(MergedData.sleep.logicals.Manual.awakeLogical);
                data.(animalID).(vID).nremData = sum(MergedData.sleep.logicals.Manual.nremLogical);
                data.(animalID).(vID).remData = sum(MergedData.sleep.logicals.Manual.remLogical);
            else
                data.(animalID).(vID).awakeData = data.(animalID).(vID).awakeData + sum(MergedData.sleep.logicals.Manual.awakeLogical);
                data.(animalID).(vID).nremData = data.(animalID).(vID).nremData + sum(MergedData.sleep.logicals.Manual.nremLogical);
                data.(animalID).(vID).remData =  data.(animalID).(vID).remData + sum(MergedData.sleep.logicals.Manual.remLogical);
            end
        end
    end
    % find baseline of each vessel ID
    uniqueVIDs = fieldnames(data.(animalID));
    for cc = 1:length(uniqueVIDs)
        diamBaseline = [];
        strDays = fieldnames(RestingBaselines.manualSelection.vesselDiameter.data.(uniqueVIDs{cc,1}));
        for dd = 1:length(strDays)
            diamBaseline(dd,1) = RestingBaselines.manualSelection.vesselDiameter.data.(uniqueVIDs{cc,1}).(strDays{dd,1}); %#ok<*AGROW>
        end
        data.(animalID).(uniqueVIDs{cc,1}).baseline = mean(diamBaseline);
    end
    cd(rootFolder)
end
% go through and put each animal/vessel ID into an table
gg = 1;
labelTime = 5;   % seconds
for ee = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,ee};
    uniqueVIDs = fieldnames(data.(animalID));
    TwoP_indTotalTimeAwake = [];
    TwoP_indTotalTimeNREM = [];
    TwoP_indTotalTimeREM = [];
    for ff = 1:length(uniqueVIDs)
        TwoP_animalIDList{gg,1} = [animalID '_' uniqueVIDs{ff,1}];
        TwoP_baselineDiams(gg,1) = round(data.(animalID).(uniqueVIDs{ff,1}).baseline,1);
        TwoP_totalTimeAwake(gg,1) = round(((data.(animalID).(uniqueVIDs{ff,1}).awakeData*labelTime)/60),1);   % sec -> min
        TwoP_totalTimeNREM(gg,1) = round(((data.(animalID).(uniqueVIDs{ff,1}).nremData*labelTime)/60),1);   % sec -> min -> hrs
        TwoP_totalTimeREM(gg,1) = round(((data.(animalID).(uniqueVIDs{ff,1}).remData*labelTime)/60),1);   % sec -> min -> hrs
        TwoP_totalTimeMins(gg,1) = round(TwoP_totalTimeAwake(gg,1) + TwoP_totalTimeNREM(gg,1) + TwoP_totalTimeREM(gg,1),1); 
        gg = gg + 1;
        % per animal
        TwoP_indTotalTimeAwake(ff,1) = ((data.(animalID).(uniqueVIDs{ff,1}).awakeData*labelTime)/60)/60;   % sec -> min -> hrs
        TwoP_indTotalTimeNREM(ff,1) = ((data.(animalID).(uniqueVIDs{ff,1}).nremData*labelTime)/60)/60;   % sec -> min -> hrs
        TwoP_indTotalTimeREM(ff,1) = ((data.(animalID).(uniqueVIDs{ff,1}).remData*labelTime)/60)/60;   % sec -> min -> hrs    
    end
    TwoP_totalTimePerAnimal(ee,1) = sum(TwoP_indTotalTimeAwake) + sum(TwoP_indTotalTimeNREM) + sum(TwoP_indTotalTimeREM);
end
TwoP_allTimeHours = round(sum(TwoP_totalTimeMins)/60);
TwoP_meanTimeHours = round(mean(TwoP_totalTimePerAnimal,1),1);
TwoP_stdTimeHours = round(std(TwoP_totalTimePerAnimal,0,1),1);
% save results
AnalysisResults.ArterioleDurations.TwoP_animalIDs = TwoP_animalIDList;
AnalysisResults.ArterioleDurations.TwoP_baselineDiams = TwoP_baselineDiams;
AnalysisResults.ArterioleDurations.TwoP_totalTimeAwake = TwoP_totalTimeAwake;
AnalysisResults.ArterioleDurations.TwoP_totalTimeNREM = TwoP_totalTimeNREM;
AnalysisResults.ArterioleDurations.TwoP_totalTimeREM = TwoP_totalTimeREM;
AnalysisResults.ArterioleDurations.TwoP_totalTimeMins = TwoP_totalTimeMins;
AnalysisResults.ArterioleDurations.TwoP_allTimeHours = TwoP_allTimeHours;
AnalysisResults.ArterioleDurations.TwoP_meanTimeHours = TwoP_meanTimeHours;
AnalysisResults.ArterioleDurations.TwoP_stdTimeHours = TwoP_stdTimeHours;
% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults')

end
