function [SVMModel] = TrainModel_SVM(trainingDataFileIDs)
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

for a = 1:size(trainingDataFileIDs,1)
   trainingTableFileID = trainingDataFileIDs(a,:); 
   if a == 1
       load(trainingTableFileID)
       joinedTable = T;
   else
       load(trainingTableFileID)
       joinedTable = vertcat(joinedTable, T);
   end
end

[animalID,~,~] = GetFileInfo_IOS(trainingTableFileID);
X = joinedTable(:,1:12);
Y = joinedTable(:,13);
t = templateSVM('Standardize',true,'KernelFunction','gaussian');

disp('Training Support Vector Machine...'); disp(' ')
SVMModel = fitcecoc(X,Y,'Learners',t,'FitPosterior',true,'ClassNames',{'Awake','Transition','NREM','REM'},'Verbose',2);

disp('Cross-validating (10-fold) the support vector machine classifier...'); disp(' ')
CVSVMModel = crossval(SVMModel);

loss = kfoldLoss(CVSVMModel);
disp(['k-fold loss classification error: ' num2str(loss*100) '%']); disp(' ')

save([animalID '_SVM_SleepScoringModel.mat'], 'SVMModel')
save([animalID '_SVM_ModelCrossValidation.mat'], 'CVSVMModel')

end

