function [SVMResults] = PredictBehaviorEvents_SVM(modelDataFileIDs, SVMModel)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs: List of processed data file IDs (char) and the loaded RestingBaselines structure
%
%   Outputs:
%
%   Last Revised: July 26th, 2019
%________________________________________________________________________________________________________________________

[animalID,~,~] = GetFileInfo_IOS(modelDataFileIDs(1,:));
disp('Predicting behavior events using SVM model'); disp(' ')

for a = 1:size(modelDataFileIDs,1)
    modelDataFileID = modelDataFileIDs(a,:);
    if a == 1
        load(modelDataFileID)
        joinedTable = paramsTable;
        joinedFileList = cell(size(paramsTable,1),1);
        joinedFileList(:) = {modelDataFileID};
    else
        load(modelDataFileID)
        fileIDCells = cell(size(paramsTable,1),1);
        fileIDCells(:) = {modelDataFileID};
        joinedTable = vertcat(joinedTable, paramsTable);
        joinedFileList = vertcat(joinedFileList, fileIDCells);
    end
end

%% Obtain the label/score from the model
scoringTable = joinedTable;
[label,score] = predict(SVMModel,scoringTable);
SVMResults.fileIDs = joinedFileList;
SVMResults.inputData = scoringTable;
SVMResults.labels = label;
SVMResults.scores = score;
saveID = [animalID '_SVM_SleepScoringResults'];
save(saveID, 'SVMResults')

end