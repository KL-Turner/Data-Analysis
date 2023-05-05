zap;
animalIDs = {'T141','T155','T156','T157','T186','T187','T189','T220'};
curDir = cd;
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    dataDir = [curDir '/' animalID '/Imaging/'];
    cd(dataDir)
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);

    for bb = 1:size(procDataFileIDs,1)
        load(procDataFileIDs(bb,:));

        if isfield(ProcData.notes,'motionArtifact') == false
            disp([num2str(aa) ' ' num2str(bb)])
            ProcData.notes.motionArtifact = 0;
            save(procDataFileIDs(bb,:),'ProcData')
        end

    end

    cd(curDir)
end