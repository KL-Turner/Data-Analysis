%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Calculate the average coherence of different behavioral states
%________________________________________________________________________________________________________________________
%
%   Inputs: none
%
%   Outputs: Generates summary figures saved to C: drive Documents folder
%
%   Last Revised: Oct 1st, 2019
%________________________________________________________________________________________________________________________
clear
clc

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110'};
driveLetters = {'E','E','E','F','F','F','D','D'};
hemDataTypes = {'LH','RH'};
behavFields = {'Rest','NREM','REM','Unstim','All'};
baselineTypes = {'manualSelection','setDuration','entireDuration'};
coherr_dataTypes = {'gammaBandPower','muaPower'};
colorbrewer_setA_colorA = [0.520000 0.520000 0.510000];
colorbrewer_setA_colorB = [(31/256) (120/256) (180/256)];
colorbrewer_setA_colorC = [(255/256) (0/256) (115/256)];
colorbrewer_setA_colorD = [(51/256) (160/256) (44/256)];
colorbrewer_setA_colorE = [(255/256) (140/256) (0/256)];

%% cd through each animal's directory and extract the appropriate analysis results
figure;
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\' animalID '\Combined Imaging\'];
    cd(dataPath)
    load([animalID '_AnalysisResults.mat']);
    for b = 1:length(hemDataTypes)
        hemDataType = hemDataTypes{1,b};
        for c = 1:length(coherr_dataTypes)
            coherr_dataType = coherr_dataTypes{1,c};
            plot(AnalysisResults.HRFs.(hemDataType).(coherr_dataType).kernelT,AnalysisResults.HRFs.(hemDataType).(coherr_dataType).bestKernel)
            hold on
        end
    end
end