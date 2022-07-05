zap;
curDir = cd;
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123','T141','T155','T156','T157','T186','T187','T189','T220'};
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    cd([animalID '/Bilateral Imaging/']);
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        load(procDataFileID);
        pupilArea = ProcData.data.Pupil.pupilArea;
        diameter = sqrt(pupilArea./pi)*2;
        ProcData.data.Pupil.pixelDiameter = diameter;
        ProcData.data.Pupil.pupilMinor = ProcData.data.Pupil.pupiMinor;
        ProcData.data.Pupil = rmfield(ProcData.data.Pupil,'pupiMinor');
        ProcData.data.Pupil.mmPerPixel = 0.018;
        ProcData.data.Pupil.mmDiameter = diameter.*ProcData.data.Pupil.mmPerPixel;
        ProcData.data.Pupil.mmArea = pupilArea.*(ProcData.data.Pupil.mmPerPixel^2);
        save(procDataFileID,'ProcData')
    end
    StageThreeProcessing_IOS()
    % find and load RestingBaselines.mat struct
    baselineDataFileStruct = dir('*_RestingBaselines.mat');
    baselineDataFile = {baselineDataFileStruct.name}';
    baselineDataFileID = char(baselineDataFile);
    load(baselineDataFileID,'-mat')
    for cc = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(cc,:);
        load(procDataFileID)
        ProcData.data.Pupil.zDiameter = (ProcData.data.Pupil.mmDiameter - RestingBaselines.manualSelection);
        ProcData.data.Pupil.zArea = (ProcData.data.Pupil.mmDiameter - RestingBaselines.manualSelection);
    end
    cd(curDir)
end