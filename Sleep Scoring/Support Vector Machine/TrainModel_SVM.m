function [SVMModel] = TrainModel_SVM(trainingTableFileIDs)
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

for a = 1:size(trainingTableFileIDs,1)
   trainingTableFileID = trainingTableFileIDs(a,:); 
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
SVMModel = fitcecoc(X, Y,'Learners',t,'FitPosterior',true,'ClassNames',{'Awake','Transition', 'NREM','REM'}, 'Verbose',2);
save([animalID '_SVM_SleepScoringModel.mat'], 'SVMModel')

end

