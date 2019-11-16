%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%
%   Last Revised: Oct 1st, 2019
%________________________________________________________________________________________________________________________
clear
clc

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111'};
driveLetters = {'E','E','E','F','F','F','D','D','D'};

for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\' animalID '\Combined Imaging\'];
    disp(['Re-running data analysis for animal: ' animalID]); disp(' ')
    cd(dataPath)
    DataAnalysis_IOS
end

%% Average figure generation
AvgCoherence_IOS
AvgPowerSpectra_IOS
AvgXCorr_IOS
AvgStim_IOS
AvgWhisk_IOS
AvgCorrCoeff_IOS
AvgCBVandHeartRate_IOS
AvgResponseFunctionPredictions_IOS
