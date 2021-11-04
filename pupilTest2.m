zap;
curDir = cd;
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123','T141','T155','T156','T157','T186','T187','T189','T220'};
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    cd([animalID '/Bilateral Imaging/']);
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
    AddPupilSleepParameters_IOS(procDataFileIDs)
    UpdatePupilSleepData_IOS(procDataFileIDs)
    cd(curDir)
end