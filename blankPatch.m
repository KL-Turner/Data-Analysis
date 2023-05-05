zap;
animalIDs = {'T224','T228','T233','T263','T264','T265','T267'};
curDir = cd;
for aa = 1:length(animalIDs)
    disp(num2str(aa))
    animalID = animalIDs{1,aa};
    dataDir = [curDir '/' animalID '/Imaging/'];
    cd(dataDir)
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);

    baselineFile = ls('*RestingBaselines.mat');
    load(baselineFile)

    [RestData] = ExtractRestingData_IOS(procDataFileIDs,2);
    [RestData] = NormRestDataStruct_IOS(animalID,RestData,RestingBaselines,'manualSelection');

    [EventData] = ExtractEventTriggeredData_IOS(procDataFileIDs);
    [EventData] = NormEventDataStruct_IOS(animalID,EventData,RestingBaselines,'manualSelection');

    cd(curDir)
end