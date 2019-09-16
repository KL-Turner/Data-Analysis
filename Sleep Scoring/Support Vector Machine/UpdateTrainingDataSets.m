function [] = UpdateTrainingDataSets(procDataFileIDs)
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

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    modelDataSetID = [procDataFileID(1:end-12) 'ModelData.mat'];
    trainingDataSetID = [procDataFileID(1:end-12) 'TrainingData.mat'];
    if exist(trainingDataSetID) > 1
        load(modelDataSetID)
        load(trainingDataSetID)
        disp(['Updating training data set for ' trainingDataSetID '...' ]); disp(' ')
        paramsTable.behavState = trainingTable.behavState;
        trainingTable = paramsTable;
        save(trainingDataSetID, 'trainingTable')
    end
end