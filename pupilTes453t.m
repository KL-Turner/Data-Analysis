clear; clc;
curDir = cd;
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123','T141','T155','T156','T157','T186','T187','T189','T220'};
goodFiles = []; xx = 1;
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    cd([animalID '/Bilateral Imaging/']);
    procDataFileIDs = ls('*_ProcData.mat');
    for bb = 1:size(procDataFileIDs,1)
        load(procDataFileIDs(bb,:));
        if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
            goodFiles(xx,1) = 1;
        elseif strcmp(ProcData.data.Pupil.diameterCheck,'n') == true
            goodFiles(xx,1) = 0;
        end
        xx = xx + 1;
    end
    cd(curDir)
end
test = sum(goodFiles)/length(goodFiles)
totalHours = (length(goodFiles)*15)/60
testHours = (sum(goodFiles)*15)/60