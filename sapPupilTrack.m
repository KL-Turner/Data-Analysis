zap;
animalIDs = {'T142','T144','T159','T172','T182','T183','T184','T185','T198'};
curDir = cd;
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    dataDir = [curDir '/' animalID '/Imaging/'];
    cd(dataDir)
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
    [animalID,~,~] = GetFileInfo_IOS(procDataFileIDs(1,:));
    % track pupil area and blink detetction
    RunPupilTracker_IOS(procDataFileIDs)
    for bb = 1:size(procDataFileIDs,1)
        disp([num2str(aa) '/' num2str(length(animalIDs))  ' _ ' num2str(bb) '/' num2str(size(procDataFileIDs,1))]); disp(' ')
        PatchPupilArea_IOS(procDataFileIDs(bb,:))
        CheckPupilVideoFrames_IOS(procDataFileIDs(bb,:))
        CheckPupilBlinks_IOS(procDataFileIDs(bb,:))
        CheckPupilDiameter_IOS(procDataFileIDs(bb,:))
        ConvertPupilAreaToDiameter_IOS(procDataFileIDs(bb,:))
    end
    % add pupil area to RestData.mat
    dataTypes = {'pupilArea','diameter','mmArea','mmDiameter','patchCentroidX','patchCentroidY'};
    ExtractPupilRestingData_IOS(procDataFileIDs,dataTypes);
    % add pupil baseline to Restingbaselines.mat
    [RestingBaselines] = AddPupilRestingBaseline_IOS();
    % zScore pupil data
    for ff = 1:size(procDataFileIDs,1)
        disp(['Z-scoring pupil data of file ' num2str(ff) '/' num2str(size(procDataFileIDs,1))]); disp(' ')
        procDataFileID = procDataFileIDs(ff,:);
        zScorePupilData_IOS(procDataFileID,RestingBaselines)
    end
    % add pupil area to RestData.mat
    dataTypes = {'pupilArea','diameter','mmArea','mmDiameter','zArea','zDiameter','patchCentroidX','patchCentroidY','LH_HbT','RH_HbT','LH_gammaBandPower','RH_gammaBandPower'};
    [RestData] = ExtractPupilRestingData_IOS(procDataFileIDs,dataTypes);
    % add pupil area to EventData.mat
    [EventData] = ExtractPupilEventTriggeredData_IOS(procDataFileIDs,dataTypes);
    % normalize Rest/Event data structures
    [RestData] = NormRestDataStruct_IOS(RestData,RestingBaselines,'manualSelection');
    [EventData] = NormEventDataStruct_IOS(EventData,RestingBaselines,'manualSelection');
    % add pupil data to SleepData.mat
    AddPupilSleepParameters_IOS(procDataFileIDs)
    UpdatePupilSleepData_IOS(procDataFileIDs)
    cd(curDir)
end