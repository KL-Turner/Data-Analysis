function [Results_VesselBaselineShift] = AnalyzeVesselBaselineShift(animalID,group,rootFolder,delim,Results_VesselBaselineShift)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the whisking-evoked arteriole D/D responses (2PLSM)
%________________________________________________________________________________________________________________________

%% only run analysis for valid animal IDs
baseLocation = [rootFolder delim group delim animalID delim 'Combined Imaging'];
cd(baseLocation)
% find and load RestingBaselines.mat strut
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID)
% cd to data location
dataLocation = [rootFolder delim group delim animalID delim 'Isoflurane Trials'];
cd(dataLocation)
% character list of all MergedData files
mergedDirectory = dir('*_MergedData.mat');
mergedDataFiles = {mergedDirectory.name}';
mergedDataFileIDs = char(mergedDataFiles);
for aa = 1:size(mergedDataFileIDs,1)
    mergedDataFileID = mergedDataFileIDs(aa,:);
    load(mergedDataFileID)
    [~,~,fileDate,~,~,vID] = GetFileInfo2_2P(mergedDataFileID);
    strDay = ConvertDate_2P(fileDate);
    % save results
    Results_VesselBaselineShift.(animalID).(vID).diameter = mean(MergedData.data.vesselDiameter.data(1:30*5));
    Results_VesselBaselineShift.(animalID).(vID).baseline = RestingBaselines.manualSelection.vesselDiameter.data.(vID).(strDay);
end
% save data
cd(rootFolder)
save('Results_VesselBaselineShift.mat','Results_VesselBaselineShift')

end
