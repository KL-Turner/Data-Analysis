function [SVMModel] = TrainModel_SVM(animalIDs)
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
%   Last Revised: July 27th, 2019
%________________________________________________________________________________________________________________________

for a = 1:length(animalIDs)
    startingDirectory = cd;
    trainingDirectory = [animalIDs{1,a} '\SVM Training Set\'];
    cd(trainingDirectory)
    % character list of all training files
    trainingDataFileStruct = dir('*_TrainingData.mat');
    trainingDataFiles = {trainingDataFileStruct.name}';
    trainingDataFileIDs = char(trainingDataFiles);
    % Load each updated training set and concatenate the data into table
    for b = 1:size(trainingDataFileIDs,1)
        trainingTableFileID = trainingDataFileIDs(b,:);
        if a == 1 && b == 1
           load(trainingTableFileID)
           joinedTable = trainingTable;
       else
           load(trainingTableFileID)
           joinedTable = vertcat(joinedTable, trainingTable); %#ok<AGROW>
       end
    end
   cd(startingDirectory)
end

X = joinedTable(:,1:end-1);
Y = joinedTable(:,end);
t = templateSVM('Standardize',true,'KernelFunction','gaussian');
disp('Training Support Vector Machine...'); disp(' ')
SVMModel = fitcecoc(X,Y,'Learners',t,'FitPosterior',true,'ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'},'Verbose',2);
disp('Cross-validating (10-fold) the support vector machine classifier...'); disp(' ')
CVSVMModel = crossval(SVMModel);
loss = kfoldLoss(CVSVMModel);
disp(['k-fold loss classification error: ' num2str(loss*100) '%']); disp(' ')
saveLoc = startingDirectory;
save([saveLoc '\IOS_SVM_SleepScoringModel.mat\'],'SVMModel')
save([saveLoc '\IOS_SVM_ModelCrossValidation.mat\'],'CVSVMModel')

end
