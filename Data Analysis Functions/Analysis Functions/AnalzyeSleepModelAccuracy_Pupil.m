function [Results_SleepModel] = AnalzyeSleepModelAccuracy_Pupil(animalID,rootFolder,delim,Results_SleepModel)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________

%% function parameters
dataLocation = [rootFolder delim 'Data' delim animalID delim 'Training Data'];
cd(dataLocation)
% ProcData file IDs
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% prepare training data by updating parameters
AddPupilSleepParameters_IOS(procDataFileIDs)
CreatePupilModelDataSet_IOS(procDataFileIDs)
CreatePupilTrainingDataSet_IOS(procDataFileIDs)
% training data file IDs
trainingDataFileStruct = dir('*_TrainingData.mat');
trainingDataFiles = {trainingDataFileStruct.name}';
trainingDataFileIDs = char(trainingDataFiles);
% Load each updated training set and concatenate the data into table
for bb = 1:size(trainingDataFileIDs,1)
    trainingTableFileID = trainingDataFileIDs(bb,:);
    if bb == 1
        load(trainingTableFileID)
        dataLength = size(trainingTable,1);
        joinedTableOdd = trainingTable;
    elseif bb == 2
        load(trainingTableFileID)
        joinedTableEven = trainingTable;
    elseif rem(bb,2) == 1
        load(trainingTableFileID)
        joinedTableOdd = vertcat(joinedTableOdd,trainingTable); %#ok<*AGROW>
    elseif rem(bb,2) == 0
        load(trainingTableFileID)
        joinedTableEven = vertcat(joinedTableEven,trainingTable);
    end
end
% train on odd data
Xodd = joinedTableOdd(:,1:end - 1);
Yodd = joinedTableOdd(:,end);
% test on even data
Xeven = joinedTableEven(:,1:end - 1);
Yeven = joinedTableEven(:,end);
%% Train Support Vector Machine (SVM) classifier
t = templateSVM('Standardize',true,'KernelFunction','gaussian');
disp('Training Support Vector Machine...'); disp(' ')
SVM_MDL = fitcecoc(Xodd,Yodd,'Learners',t,'FitPosterior',true,'ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'},'Verbose',2);
% save model in desired location
% save([dirpath animalID '_IOS_SVM_PupilSleepScoringModel.mat'],'SVM_MDL')
% determine k-fold loss of the model
disp('Cross-validating (3-fold) the support vector machine classifier...'); disp(' ')
CV_SVM_MDL = crossval(SVM_MDL,'kfold',3);
loss = kfoldLoss(CV_SVM_MDL);
disp(['k-fold loss classification error: ' num2str(loss*100) '%']); disp(' ')
% use the model to generate a set of scores for the even set of data
[XoddLabels,~] = predict(SVM_MDL,Xodd);
[XevenLabels,~] = predict(SVM_MDL,Xeven);
% apply a logical patch on the REM events
oddREMindex = strcmp(XoddLabels,'REM Sleep');
evenREMindex = strcmp(XevenLabels,'REM Sleep');
oddNumFiles = length(XoddLabels)/dataLength;
evenNumFiles = length(XevenLabels)/dataLength;
oddReshapedREMindex = reshape(oddREMindex,dataLength,oddNumFiles);
evenReshapedREMindex = reshape(evenREMindex,dataLength,evenNumFiles);
oddPatchedREMindex = [];
evenPatchedREMindex = [];
% training data - patch missing REM indeces due to theta band falling off
for ii = 1:size(oddReshapedREMindex,2)
    oddREMArray = oddReshapedREMindex(:,ii);
    oddPatchedREMarray = LinkBinaryEvents_IOS(oddREMArray',[5,0]);
    oddPatchedREMindex = vertcat(oddPatchedREMindex,oddPatchedREMarray');
end
% testing data - patch missing REM indeces due to theta band falling off
for ii = 1:size(evenReshapedREMindex,2)
    evenREMArray = evenReshapedREMindex(:,ii);
    evenPatchedREMarray = LinkBinaryEvents_IOS(evenREMArray',[5,0]);
    evenPatchedREMindex = vertcat(evenPatchedREMindex,evenPatchedREMarray');
end
% training data - change labels for each event
for jj = 1:length(XoddLabels)
    if oddPatchedREMindex(jj,1) == 1
        XoddLabels{jj,1} = 'REM Sleep';
    end
end
% testing data - change labels for each event
for jj = 1:length(XevenLabels)
    if evenPatchedREMindex(jj,1) == 1
        XevenLabels{jj,1} = 'REM Sleep';
    end
end
% save labels for later confusion matrix
Results_SleepModel.(animalID).SVM.trainYlabels = Yodd.behavState;
Results_SleepModel.(animalID).SVM.trainXlabels = XoddLabels;
Results_SleepModel.(animalID).SVM.testYlabels = Yeven.behavState;
Results_SleepModel.(animalID).SVM.testXlabels = XevenLabels;
%% Ensemble classification - AdaBoostM2, Subspace, Bag, LPBoost,RUSBoost, TotalBoost
disp('Training Ensemble Classifier...'); disp(' ')
t = templateTree('Reproducible',true);
EC_MDL = fitcensemble(Xodd,Yodd,'OptimizeHyperparameters','auto','Learners',t,'HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'),'ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
% save model in desired location
% save([dirpath animalID '_IOS_EC_PupilSleepScoringModel.mat'],'EC_MDL')
% determine k-fold loss of the model
disp('Cross-validating (3-fold) the ensemble classifier...'); disp(' ')
CV_EC_MDL = crossval(EC_MDL,'kfold',3);
loss = kfoldLoss(CV_EC_MDL);
disp(['k-fold loss classification error: ' num2str(loss*100) '%']); disp(' ')
% use the model to generate a set of scores for the even set of data
[XoddLabels,~] = predict(EC_MDL,Xodd);
[XevenLabels,~] = predict(EC_MDL,Xeven);
% apply a logical patch on the REM events
oddREMindex = strcmp(XoddLabels,'REM Sleep');
evenREMindex = strcmp(XevenLabels,'REM Sleep');
oddNumFiles = length(XoddLabels)/dataLength;
evenNumFiles = length(XevenLabels)/dataLength;
oddReshapedREMindex = reshape(oddREMindex,dataLength,oddNumFiles);
evenReshapedREMindex = reshape(evenREMindex,dataLength,evenNumFiles);
oddPatchedREMindex = [];
evenPatchedREMindex = [];
% training data - patch missing REM indeces due to theta band falling off
for ii = 1:size(oddReshapedREMindex,2)
    oddREMArray = oddReshapedREMindex(:,ii);
    oddPatchedREMarray = LinkBinaryEvents_IOS(oddREMArray',[5,0]);
    oddPatchedREMindex = vertcat(oddPatchedREMindex,oddPatchedREMarray');
end
% testing data - patch missing REM indeces due to theta band falling off
for ii = 1:size(evenReshapedREMindex,2)
    evenREMArray = evenReshapedREMindex(:,ii);
    evenPatchedREMarray = LinkBinaryEvents_IOS(evenREMArray',[5,0]);
    evenPatchedREMindex = vertcat(evenPatchedREMindex,evenPatchedREMarray');
end
% training data - change labels for each event
for jj = 1:length(XoddLabels)
    if oddPatchedREMindex(jj,1) == 1
        XoddLabels{jj,1} = 'REM Sleep';
    end
end
% testing data - change labels for each event
for jj = 1:length(XevenLabels)
    if evenPatchedREMindex(jj,1) == 1
        XevenLabels{jj,1} = 'REM Sleep';
    end
end
% save labels for later confusion matrix
Results_SleepModel.(animalID).EC.trainYlabels= Yodd.behavState;
Results_SleepModel.(animalID).EC.trainXlabels = XoddLabels;
Results_SleepModel.(animalID).EC.testYlabels = Yeven.behavState;
Results_SleepModel.(animalID).EC.testXlabels = XevenLabels;
%% Decision Tree classification
disp('Training Decision Tree Classifier...'); disp(' ')
DT_MDL = fitctree(Xodd,Yodd,'ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
% save model in desired location
% save([dirpath animalID '_IOS_DT_PupilSleepScoringModel.mat'],'DT_MDL')
% use the model to generate a set of scores for the even set of data
[XoddLabels,~] = predict(DT_MDL,Xodd);
[XevenLabels,~] = predict(DT_MDL,Xeven);
% apply a logical patch on the REM events
oddREMindex = strcmp(XoddLabels,'REM Sleep');
evenREMindex = strcmp(XevenLabels,'REM Sleep');
oddNumFiles = length(XoddLabels)/dataLength;
evenNumFiles = length(XevenLabels)/dataLength;
oddReshapedREMindex = reshape(oddREMindex,dataLength,oddNumFiles);
evenReshapedREMindex = reshape(evenREMindex,dataLength,evenNumFiles);
oddPatchedREMindex = [];
evenPatchedREMindex = [];
% training data - patch missing REM indeces due to theta band falling off
for ii = 1:size(oddReshapedREMindex,2)
    oddREMArray = oddReshapedREMindex(:,ii);
    oddPatchedREMarray = LinkBinaryEvents_IOS(oddREMArray',[5,0]);
    oddPatchedREMindex = vertcat(oddPatchedREMindex,oddPatchedREMarray');
end
% testing data - patch missing REM indeces due to theta band falling off
for ii = 1:size(evenReshapedREMindex,2)
    evenREMArray = evenReshapedREMindex(:,ii);
    evenPatchedREMarray = LinkBinaryEvents_IOS(evenREMArray',[5,0]);
    evenPatchedREMindex = vertcat(evenPatchedREMindex,evenPatchedREMarray');
end
% training data - change labels for each event
for jj = 1:length(XoddLabels)
    if oddPatchedREMindex(jj,1) == 1
        XoddLabels{jj,1} = 'REM Sleep';
    end
end
% testing data - change labels for each event
for jj = 1:length(XevenLabels)
    if evenPatchedREMindex(jj,1) == 1
        XevenLabels{jj,1} = 'REM Sleep';
    end
end
% save labels for later confusion matrix
Results_SleepModel.(animalID).DT.trainYlabels = Yodd.behavState;
Results_SleepModel.(animalID).DT.trainXlabels = XoddLabels;
Results_SleepModel.(animalID).DT.testYlabels = Yeven.behavState;
Results_SleepModel.(animalID).DT.testXlabels = XevenLabels;
%% Random forest
disp('Training Random Forest Classifier...'); disp(' ')
numTrees = 128;
RF_MDL = TreeBagger(numTrees,Xodd,Yodd,'Method','Classification','Surrogate','all','OOBPrediction','on','ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
% save model in desired location
% save([dirpath animalID '_IOS_RF_PupilSleepScoringModel.mat'],'RF_MDL')
% determine the misclassification probability (for classification trees) for out-of-bag observations in the training data
RF_OOBerror = oobError(RF_MDL,'Mode','Ensemble');
disp(['Random Forest out-of-bag error: ' num2str(RF_OOBerror*100) '%']); disp(' ')
% use the model to generate a set of scores for the even set of data
[XoddLabels,~] = predict(RF_MDL,Xodd);
[XevenLabels,~] = predict(RF_MDL,Xeven);
% apply a logical patch on the REM events
oddREMindex = strcmp(XoddLabels,'REM Sleep');
evenREMindex = strcmp(XevenLabels,'REM Sleep');
oddNumFiles = length(XoddLabels)/dataLength;
evenNumFiles = length(XevenLabels)/dataLength;
oddReshapedREMindex = reshape(oddREMindex,dataLength,oddNumFiles);
evenReshapedREMindex = reshape(evenREMindex,dataLength,evenNumFiles);
oddPatchedREMindex = [];
evenPatchedREMindex = [];
% training data - patch missing REM indeces due to theta band falling off
for ii = 1:size(oddReshapedREMindex,2)
    oddREMArray = oddReshapedREMindex(:,ii);
    oddPatchedREMarray = LinkBinaryEvents_IOS(oddREMArray',[5,0]);
    oddPatchedREMindex = vertcat(oddPatchedREMindex,oddPatchedREMarray');
end
% testing data - patch missing REM indeces due to theta band falling off
for ii = 1:size(evenReshapedREMindex,2)
    evenREMArray = evenReshapedREMindex(:,ii);
    evenPatchedREMarray = LinkBinaryEvents_IOS(evenREMArray',[5,0]);
    evenPatchedREMindex = vertcat(evenPatchedREMindex,evenPatchedREMarray');
end
% training data - change labels for each event
for jj = 1:length(XoddLabels)
    if oddPatchedREMindex(jj,1) == 1
        XoddLabels{jj,1} = 'REM Sleep';
    end
end
% testing data - change labels for each event
for jj = 1:length(XevenLabels)
    if evenPatchedREMindex(jj,1) == 1
        XevenLabels{jj,1} = 'REM Sleep';
    end
end
% save labels for later confusion matrix
Results_SleepModel.(animalID).RF.trainYlabels = Yodd.behavState;
Results_SleepModel.(animalID).RF.trainXlabels = XoddLabels;
Results_SleepModel.(animalID).RF.testYlabels = Yeven.behavState;
Results_SleepModel.(animalID).RF.testXlabels = XevenLabels;
%% k-nearest neighbor classifier
disp('Training k-nearest neighbor Classifier...'); disp(' ')
t = templateKNN('NumNeighbors',5,'Standardize',1);
KNN_MDL = fitcecoc(Xodd,Yodd,'Learners',t);
% use the model to generate a set of scores for the even set of data
[XoddLabels,~] = predict(KNN_MDL,Xodd);
[XevenLabels,~] = predict(KNN_MDL,Xeven);
% apply a logical patch on the REM events
oddREMindex = strcmp(XoddLabels,'REM Sleep');
evenREMindex = strcmp(XevenLabels,'REM Sleep');
oddNumFiles = length(XoddLabels)/dataLength;
evenNumFiles = length(XevenLabels)/dataLength;
oddReshapedREMindex = reshape(oddREMindex,dataLength,oddNumFiles);
evenReshapedREMindex = reshape(evenREMindex,dataLength,evenNumFiles);
oddPatchedREMindex = [];
evenPatchedREMindex = [];
% training data - patch missing REM indeces due to theta band falling off
for ii = 1:size(oddReshapedREMindex,2)
    oddREMArray = oddReshapedREMindex(:,ii);
    oddPatchedREMarray = LinkBinaryEvents_IOS(oddREMArray',[5,0]);
    oddPatchedREMindex = vertcat(oddPatchedREMindex,oddPatchedREMarray');
end
% testing data - patch missing REM indeces due to theta band falling off
for ii = 1:size(evenReshapedREMindex,2)
    evenREMArray = evenReshapedREMindex(:,ii);
    evenPatchedREMarray = LinkBinaryEvents_IOS(evenREMArray',[5,0]);
    evenPatchedREMindex = vertcat(evenPatchedREMindex,evenPatchedREMarray');
end
% training data - change labels for each event
for jj = 1:length(XoddLabels)
    if oddPatchedREMindex(jj,1) == 1
        XoddLabels{jj,1} = 'REM Sleep';
    end
end
% testing data - change labels for each event
for jj = 1:length(XevenLabels)
    if evenPatchedREMindex(jj,1) == 1
        XevenLabels{jj,1} = 'REM Sleep';
    end
end
% save labels for later confusion matrix
Results_SleepModel.(animalID).KNN.trainYlabels = Yodd.behavState;
Results_SleepModel.(animalID).KNN.trainXlabels = XoddLabels;
Results_SleepModel.(animalID).KNN.testYlabels = Yeven.behavState;
Results_SleepModel.(animalID).KNN.testXlabels = XevenLabels;
%% Naive Bayes classifier
disp('Training naive Bayes Classifier...'); disp(' ')
NB_MDL = fitcnb(Xodd,Yodd,'ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
% use the model to generate a set of scores for the even set of data
[XoddLabels,~] = predict(NB_MDL,Xodd);
[XevenLabels,~] = predict(NB_MDL,Xeven);
% apply a logical patch on the REM events
oddREMindex = strcmp(XoddLabels,'REM Sleep');
evenREMindex = strcmp(XevenLabels,'REM Sleep');
oddNumFiles = length(XoddLabels)/dataLength;
evenNumFiles = length(XevenLabels)/dataLength;
oddReshapedREMindex = reshape(oddREMindex,dataLength,oddNumFiles);
evenReshapedREMindex = reshape(evenREMindex,dataLength,evenNumFiles);
oddPatchedREMindex = [];
evenPatchedREMindex = [];
% training data - patch missing REM indeces due to theta band falling off
for ii = 1:size(oddReshapedREMindex,2)
    oddREMArray = oddReshapedREMindex(:,ii);
    oddPatchedREMarray = LinkBinaryEvents_IOS(oddREMArray',[5,0]);
    oddPatchedREMindex = vertcat(oddPatchedREMindex,oddPatchedREMarray');
end
% testing data - patch missing REM indeces due to theta band falling off
for ii = 1:size(evenReshapedREMindex,2)
    evenREMArray = evenReshapedREMindex(:,ii);
    evenPatchedREMarray = LinkBinaryEvents_IOS(evenREMArray',[5,0]);
    evenPatchedREMindex = vertcat(evenPatchedREMindex,evenPatchedREMarray');
end
% training data - change labels for each event
for jj = 1:length(XoddLabels)
    if oddPatchedREMindex(jj,1) == 1
        XoddLabels{jj,1} = 'REM Sleep';
    end
end
% testing data - change labels for each event
for jj = 1:length(XevenLabels)
    if evenPatchedREMindex(jj,1) == 1
        XevenLabels{jj,1} = 'REM Sleep';
    end
end
% save labels for later confusion matrix
Results_SleepModel.(animalID).NB.trainYlabels = Yodd.behavState;
Results_SleepModel.(animalID).NB.trainXlabels = XoddLabels;
Results_SleepModel.(animalID).NB.testYlabels = Yeven.behavState;
Results_SleepModel.(animalID).NB.testXlabels = XevenLabels;
% save datae
cd([rootFolder delim])
save('Results_SleepModel.mat','Results_SleepModel')

end
