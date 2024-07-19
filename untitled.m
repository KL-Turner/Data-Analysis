zap;
procDataFileIDs = ls('*ProcData.mat');
rawDataFileIDs = ls('*RawData.mat');
for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    rawDataFileID = rawDataFileIDs(a,:);
    [animalID,fileDate,fileID] = GetFileInfo_IOS(procDataFileID);
    load(procDataFileID)
    load(rawDataFileID)
    newFileID = [animalID '_' fileID '_Data.mat'];
    Data.notes = ProcData.notes;
    Data.stimulations = ProcData.data.stimulations;
    Data.logEMG = ProcData.data.EMG.emg;
    Data.whiskerAngle = ProcData.data.whiskerAngle;
    try
        Data.EMG = RawData.data.EMG(1:18000000);
        Data.forceSensor = RawData.data.forceSensor(1:18000000);
        Data.flags = ProcData.flags;
        Data.isoTag = false;
    catch
        disp('iso?')
        Data.EMG = RawData.data.EMG(1:6000000);
        Data.forceSensor = RawData.data.forceSensor(1:6000000);
        Data.isoTag = true;
    end
    save(newFileID,'Data')
end

% RestingBaselines = RestingBaselines.manualSelection.baselineFileInfo;