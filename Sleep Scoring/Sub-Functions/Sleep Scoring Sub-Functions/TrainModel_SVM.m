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
       joinedTable = join(joinedTable, T);
   end
end

SVMModel = fitcecoc(joinedTable,'behaviorState');

end

