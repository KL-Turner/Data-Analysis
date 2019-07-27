function [] = CreateTrainingDataSet_SVM(procDataFileIDs, RestingBaselines)
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
%   Last Revised: July 26th, 2019
%________________________________________________________________________________________________________________________

for a = 1:size(procDataFileIDs, 1)
    procDataFileID = procDataFileIDs(a,:);
    load(procDataFileID)
    
    variableNames = {'maxLH_CBV', 'maxRH_CBV', 'maxLH_Delta', 'maxRH_Delta', 'maxLH_Theta', 'maxRH_Theta',...
        'maxLH_Gamma', 'maxRH_Gamma', 'numWhiskEvents', 'numForceEvents', 'numEMGEvents'};
    
    [figHandle] = GenerateSingleFigures_SVM(procDataFileID, RestingBaselines);
    trialDuration = ProcData.notes.trialDuration_sec;
    numBins = trialDuration/5;
    
    for b = 1:numBins
        subplot(6,1,3)
    
    keyboard
   
    
%     T = table(

end

end
