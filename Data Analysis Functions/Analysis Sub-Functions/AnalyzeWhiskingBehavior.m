function [Results_WhiskBehav] = AnalyzeWhiskingBehavior(animalID,group,rootFolder,delim,Results_WhiskBehav)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: 
%________________________________________________________________________________________________________________________

%% only run analysis for valid animal IDs
dataLocation = [rootFolder delim group delim animalID delim 'Bilateral Imaging'];
cd(dataLocation)
% find and load EventData.mat struct
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID,'-mat')
%% analyze data during periods of whisking
% pull data from EventData.mat structure
whiskDurations = EventData.CBV_HbT.adjLH.whisk.duration;
fileIDs = unique(EventData.CBV_HbT.adjLH.whisk.fileIDs);
imagingDuration = 15*length(fileIDs);
% save results
Results_WhiskBehav.(animalID).whiskDurations = whiskDurations;
Results_WhiskBehav.(animalID).whiskDurationSec = sum(whiskDurations);
Results_WhiskBehav.(animalID).whiskDurationPerc = ((sum(whiskDurations)/60)/imagingDuration)*100;
% save data
cd(rootFolder)
save('Results_WhiskBehav.mat','Results_WhiskBehav')

end
