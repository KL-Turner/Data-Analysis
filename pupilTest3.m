clear; clc;
curDir = cd;
animalIDs = {'T120','T121','T122','T123','T141','T155'}; %'T156','T157','T186','T187','T189','T220'};
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
%     cd([animalID '/Bilateral Imaging/']);
%     StageThreeProcessing_IOS('bilateral','single','lime')
%     cd(curDir)
%     cd(animalID)
%     SleepScoreMainScript_IOS()
%     cd(curDir)
    cd([animalID '/Bilateral Imaging/']);
    Pupil_CoreProcessing()
    cd(curDir)
end