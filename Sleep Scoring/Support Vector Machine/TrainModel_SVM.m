function [SVMModel] = TrainModel_SVM(animalIDs, driveLetters)
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
   fileLoc = [driveLetters{1,a} ':\' animalIDs{1,a} '\Combined Imaging\SVM Training Set'];
   cd(fileLoc)
   
   % Character list of all '*_ProcData.mat' files
   procDataFileStruct = dir('*_ProcData.mat');
   procDataFiles = {procDataFileStruct.name}';
   procDataFileIDs = char(procDataFiles);
   
   cd ..
   % Resting baseline file structure
   baselinesFileStruct = dir('*_RestingBaselines.mat');
   baselinesFile = {baselinesFileStruct.name}';
   baselinesFileID = char(baselinesFile);
   load(baselinesFileID)
   cd(fileLoc)
   
   % Update parameters for the training data set
   AddSleepParameters_SVM(procDataFileIDs, RestingBaselines)
   CreateModelDataSet_SVM(procDataFileIDs)
   UpdateTrainingDataSets(procDataFileIDs)

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

saveLoc = 'C:\Users\klt8\Documents\';
save([saveLoc 'SVM_SleepScoringModel.mat'], 'SVMModel')
save([saveLoc 'SVM_ModelCrossValidation.mat'], 'CVSVMModel')

end

