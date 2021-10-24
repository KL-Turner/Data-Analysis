zap;
curDir = cd;
animalIDs = {'T110','T111','T119','T120','T121','T122','T123','T141','T155','T156','T157','T186','T187','T189','T220'};
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    cd([animalID '/Bilateral Imaging/']);
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
    TrackPupilDiameter_IOS(procDataFileIDs)
    cd(curDir)
end