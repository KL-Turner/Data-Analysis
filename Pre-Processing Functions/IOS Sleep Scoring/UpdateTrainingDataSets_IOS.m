function [] = UpdateTrainingDataSets_IOS(procDataFileIDs,TrainingFiles)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Update training data file with most recent predictors
%________________________________________________________________________________________________________________________

cc = 1;
% reduce file list to those with the training dates
for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    [~,fileDate,~] = GetFileInfo_IOS(procDataFileID);
    if strcmp(fileDate,TrainingFiles.day1) == true || strcmp(fileDate,TrainingFiles.day2) == true
        trainingFileList(cc,:) = procDataFileID;
        cc = cc + 1;
    end
end
% go through each training file
for bb = 1:size(trainingFileList,1)
    disp(['Updating training table from file (' num2str(bb) '/' num2str(size(trainingFileList,1)) ')']); disp(' ')
    procDataFileID = trainingFileList(bb,:);
    modelDataSetID = [procDataFileID(1:end - 12) 'ModelData.mat'];
    trainingDataSetID = [procDataFileID(1:end - 12) 'TrainingData.mat'];
    load(modelDataSetID)
    load(trainingDataSetID)
    % update training data decisions with most recent predictor table
    paramsTable.behavState = trainingTable.behavState;
    trainingTable = paramsTable;
    save(trainingDataSetID,'trainingTable')
end

end
