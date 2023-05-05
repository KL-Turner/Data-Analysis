function [Results_ModelAccuracy_GCaMP] = AnalyzeModelAccuracy_GCaMP(animalID,group,set,rootFolder,delim,Results_ModelAccuracy_GCaMP)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataLocation = [rootFolder delim 'Data' delim group delim set delim animalID delim 'Imaging'];
cd(dataLocation)
% load the random forest model for evaluating the cross-validation
modelName = [animalID '_IOS_RF_SleepScoringModel.mat'];
load(modelName,'-mat')
iterations = 100;
X = RF_MDL.X;
Y = RF_MDL.Y;
% determine the misclassification probability (for classification trees) for out-of-bag observations in the training data
Results_ModelAccuracy_GCaMP.(group).(animalID).ModelAccuracy.oobErr = oobError(RF_MDL,'Mode','Ensemble')*100;
% re-create the model 100 times with shuffled data and determine the oobError distribution
for aa = 1:iterations
    shuffYIdx = randperm(numel(Y));
    shuffY = Y(shuffYIdx);
    numTrees = 128;
    shuffRF_MDL = TreeBagger(numTrees,X,shuffY,'Method','Classification','Surrogate','all','OOBPrediction','on','ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
    Results_ModelAccuracy_GCaMP.(group).(animalID).ModelAccuracy.shuff_oobErr(aa,1) = oobError(shuffRF_MDL,'Mode','Ensemble')*100;
end
% save results
cd([rootFolder delim 'Results_Turner'])
save('Results_ModelAccuracy_GCaMP.mat','Results_ModelAccuracy_GCaMP')
cd([rootFolder delim 'Data'])