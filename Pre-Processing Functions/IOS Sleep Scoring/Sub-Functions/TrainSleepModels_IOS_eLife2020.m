function [] = TrainSleepModels_IOS_eLife2020(animalIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Train several machine learning techniques on manually scored sleep data, and evaluate each model's accuracy
%________________________________________________________________________________________________________________________

%% load in all the data to create a table of values
for aa = 1:length(animalIDs)
    clearvars -except animalIDs aa ConfusionData
    startingDirectory = cd;
    trainingDirectory = [animalIDs{1,aa} '\Training Data\'];
    cd(trainingDirectory)
    % character list of all training files
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
    % directory path for saving data
    dirpath = [startingDirectory '\' animalIDs{1,aa} '\Figures\Sleep Models\'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    
    %% Train Support Vector Machine (SVM) classifier
    t = templateSVM('Standardize',true,'KernelFunction','gaussian');
    disp('Training Support Vector Machine...'); disp(' ')
    SVM_MDL = fitcecoc(Xodd,Yodd,'Learners',t,'FitPosterior',true,'ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'},'Verbose',2);
    % save model in desired location
    save([dirpath animalIDs{1,aa} '_IOS_SVM_SleepScoringModel.mat'],'SVM_MDL')
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
        oddPatchedREMarray = LinkBinaryEvents_IOS_eLife2020(oddREMArray',[5,0]);
        oddPatchedREMindex = vertcat(oddPatchedREMindex,oddPatchedREMarray');
    end
    % testing data - patch missing REM indeces due to theta band falling off
    for ii = 1:size(evenReshapedREMindex,2)
        evenREMArray = evenReshapedREMindex(:,ii);
        evenPatchedREMarray = LinkBinaryEvents_IOS_eLife2020(evenREMArray',[5,0]);
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
    ConfusionData.SVM.trainYlabels{aa,1} = Yodd.behavState;
    ConfusionData.SVM.trainXlabels{aa,1} = XoddLabels;
    ConfusionData.SVM.testYlabels{aa,1} = Yeven.behavState;
    ConfusionData.SVM.testXlabels{aa,1} = XevenLabels;
    % confusion matrix
    SVM_confMat = figure;
    sgtitle('Support Vector Machine Classifier Confusion Matrix')
    % training data confusion chart
    subplot(1,2,1)
    oddCM = confusionchart(Yodd.behavState,XoddLabels);
    oddCM.ColumnSummary = 'column-normalized';
    oddCM.RowSummary = 'row-normalized';
    oddCM.Title = 'Training Data';
    oddConfVals = oddCM.NormalizedValues;
    oddTotalScores = sum(oddConfVals(:));
    oddSVM_accuracy = (sum(oddConfVals([1,5,9])/oddTotalScores))*100;
    disp(['Support Vector Machine model prediction accuracy (training): ' num2str(oddSVM_accuracy) '%']); disp(' ')
    % testing data confusion chart
    subplot(1,2,2)
    evenCM = confusionchart(Yeven.behavState,XevenLabels);
    evenCM.ColumnSummary = 'column-normalized';
    evenCM.RowSummary = 'row-normalized';
    evenCM.Title = 'Testing Data';
    evenConfVals = evenCM.NormalizedValues;
    evenTotalScores = sum(evenConfVals(:));
    evenSVM_accuracy = (sum(evenConfVals([1,5,9])/evenTotalScores))*100;
    disp(['Support Vector Machine model prediction accuracy (testing): ' num2str(evenSVM_accuracy) '%']); disp(' ')
    % save model and figure
    savefig(SVM_confMat,[dirpath animalIDs{1,aa} '_IOS_SVM_ConfusionMatrix']);
    close(SVM_confMat)
    
    %% Ensemble classification - AdaBoostM2, Subspace, Bag, LPBoost,RUSBoost, TotalBoost
    disp('Training Ensemble Classifier...'); disp(' ')
    t = templateTree('Reproducible',true);
    EC_MDL = fitcensemble(Xodd,Yodd,'OptimizeHyperparameters','auto','Learners',t,'HyperparameterOptimizationOptions',...
        struct('AcquisitionFunctionName','expected-improvement-plus'),'ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
    % save model in desired location
    save([dirpath animalIDs{1,aa} '_IOS_EC_SleepScoringModel.mat'],'EC_MDL')
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
        oddPatchedREMarray = LinkBinaryEvents_IOS_eLife2020(oddREMArray',[5,0]);
        oddPatchedREMindex = vertcat(oddPatchedREMindex,oddPatchedREMarray');
    end
    % testing data - patch missing REM indeces due to theta band falling off
    for ii = 1:size(evenReshapedREMindex,2)
        evenREMArray = evenReshapedREMindex(:,ii);
        evenPatchedREMarray = LinkBinaryEvents_IOS_eLife2020(evenREMArray',[5,0]);
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
    ConfusionData.EC.trainYlabels{aa,1} = Yodd.behavState;
    ConfusionData.EC.trainXlabels{aa,1} = XoddLabels;
    ConfusionData.EC.testYlabels{aa,1} = Yeven.behavState;
    ConfusionData.EC.testXlabels{aa,1} = XevenLabels;
    % confusion matrix
    EC_confMat = figure;
    sgtitle('Support Vector Machine Classifier Confusion Matrix')
    % training data confusion chart
    subplot(1,2,1)
    oddCM = confusionchart(Yodd.behavState,XoddLabels);
    oddCM.ColumnSummary = 'column-normalized';
    oddCM.RowSummary = 'row-normalized';
    oddCM.Title = 'Training Data';
    oddConfVals = oddCM.NormalizedValues;
    oddTotalScores = sum(oddConfVals(:));
    oddEC_accuracy = (sum(oddConfVals([1,5,9])/oddTotalScores))*100;
    disp(['Ensemble model prediction accuracy (training): ' num2str(oddEC_accuracy) '%']); disp(' ')
    % testing data confusion chart
    subplot(1,2,2)
    evenCM = confusionchart(Yeven.behavState,XevenLabels);
    evenCM.ColumnSummary = 'column-normalized';
    evenCM.RowSummary = 'row-normalized';
    evenCM.Title = 'Testing Data';
    evenConfVals = evenCM.NormalizedValues;
    evenTotalScores = sum(evenConfVals(:));
    evenEC_accuracy = (sum(evenConfVals([1,5,9])/evenTotalScores))*100;
    disp(['Ensemble model prediction accuracy (testing): ' num2str(evenEC_accuracy) '%']); disp(' ')
    % save model and figure
    savefig(EC_confMat,[dirpath animalIDs{1,aa} '_IOS_EC_ConfusionMatrix']);
    close(EC_confMat)
    
    %% Decision Tree classification
    disp('Training Decision Tree Classifier...'); disp(' ')
    DT_MDL = fitctree(Xodd,Yodd,'ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
    % save model in desired location
    save([dirpath animalIDs{1,aa} '_IOS_DT_SleepScoringModel.mat'],'DT_MDL')
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
        oddPatchedREMarray = LinkBinaryEvents_IOS_eLife2020(oddREMArray',[5,0]);
        oddPatchedREMindex = vertcat(oddPatchedREMindex,oddPatchedREMarray');
    end
    % testing data - patch missing REM indeces due to theta band falling off
    for ii = 1:size(evenReshapedREMindex,2)
        evenREMArray = evenReshapedREMindex(:,ii);
        evenPatchedREMarray = LinkBinaryEvents_IOS_eLife2020(evenREMArray',[5,0]);
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
    ConfusionData.DT.trainYlabels{aa,1} = Yodd.behavState;
    ConfusionData.DT.trainXlabels{aa,1} = XoddLabels;
    ConfusionData.DT.testYlabels{aa,1} = Yeven.behavState;
    ConfusionData.DT.testXlabels{aa,1} = XevenLabels;
    % confusion matrix
    DT_confMat = figure;
    sgtitle('Decision Tree Classifier Confusion Matrix')
    % training data confusion chart
    subplot(1,2,1)
    oddCM = confusionchart(Yodd.behavState,XoddLabels);
    oddCM.ColumnSummary = 'column-normalized';
    oddCM.RowSummary = 'row-normalized';
    oddCM.Title = 'Training Data';
    oddConfVals = oddCM.NormalizedValues;
    oddTotalScores = sum(oddConfVals(:));
    oddDT_accuracy = (sum(oddConfVals([1,5,9])/oddTotalScores))*100;
    disp(['Support Vector Machine model prediction accuracy (training): ' num2str(oddDT_accuracy) '%']); disp(' ')
    % testing data confusion chart
    subplot(1,2,2)
    evenCM = confusionchart(Yeven.behavState,XevenLabels);
    evenCM.ColumnSummary = 'column-normalized';
    evenCM.RowSummary = 'row-normalized';
    evenCM.Title = 'Testing Data';
    evenConfVals = evenCM.NormalizedValues;
    evenTotalScores = sum(evenConfVals(:));
    evenDT_accuracy = (sum(evenConfVals([1,5,9])/evenTotalScores))*100;
    disp(['Decision Tree model prediction accuracy (testing): ' num2str(evenDT_accuracy) '%']); disp(' ')
    % save model and figure
    savefig(DT_confMat,[dirpath animalIDs{1,aa} '_IOS_DT_ConfusionMatrix']);
    close(DT_confMat)
    
    %% Random forest
    disp('Training Random Forest Classifier...'); disp(' ')
    numTrees = 128;
    RF_MDL = TreeBagger(numTrees,Xodd,Yodd,'Method','Classification','Surrogate','all','OOBPrediction','on','ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
    % save model in desired location
    save([dirpath animalIDs{1,aa} '_IOS_RF_SleepScoringModel.mat'],'RF_MDL')
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
        oddPatchedREMarray = LinkBinaryEvents_IOS_eLife2020(oddREMArray',[5,0]);
        oddPatchedREMindex = vertcat(oddPatchedREMindex,oddPatchedREMarray');
    end
    % testing data - patch missing REM indeces due to theta band falling off
    for ii = 1:size(evenReshapedREMindex,2)
        evenREMArray = evenReshapedREMindex(:,ii);
        evenPatchedREMarray = LinkBinaryEvents_IOS_eLife2020(evenREMArray',[5,0]);
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
    ConfusionData.RF.trainYlabels{aa,1} = Yodd.behavState;
    ConfusionData.RF.trainXlabels{aa,1} = XoddLabels;
    ConfusionData.RF.testYlabels{aa,1} = Yeven.behavState;
    ConfusionData.RF.testXlabels{aa,1} = XevenLabels;
    % confusion matrix
    RF_confMat = figure;
    sgtitle('Random Forest Classifier Confusion Matrix')
    % training data confusion chart
    subplot(1,2,1)
    oddCM = confusionchart(Yodd.behavState,XoddLabels);
    oddCM.ColumnSummary = 'column-normalized';
    oddCM.RowSummary = 'row-normalized';
    oddCM.Title = 'Training Data';
    oddConfVals = oddCM.NormalizedValues;
    oddTotalScores = sum(oddConfVals(:));
    oddRF_accuracy = (sum(oddConfVals([1,5,9])/oddTotalScores))*100;
    disp(['Random Forest model prediction accuracy (training): ' num2str(oddRF_accuracy) '%']); disp(' ')
    % testing data confusion chart
    subplot(1,2,2)
    evenCM = confusionchart(Yeven.behavState,XevenLabels);
    evenCM.ColumnSummary = 'column-normalized';
    evenCM.RowSummary = 'row-normalized';
    evenCM.Title = 'Testing Data';
    evenConfVals = evenCM.NormalizedValues;
    evenTotalScores = sum(evenConfVals(:));
    evenRF_accuracy = (sum(evenConfVals([1,5,9])/evenTotalScores))*100;
    disp(['Random Forest model prediction accuracy (testing): ' num2str(evenRF_accuracy) '%']); disp(' ')
    % save model and figure
    savefig(RF_confMat,[dirpath animalIDs{1,aa} '_IOS_RF_ConfusionMatrix']);
    close(RF_confMat)
    
    %% k-nearest neighbor classifier
    disp('Training k-nearest neighbor Classifier...'); disp(' ')
    t = templateKNN('NumNeighbors',5,'Standardize',1);
    KNN_MDL = fitcecoc(Xodd,Yodd,'Learners',t);
    % save model in desired location
    save([dirpath animalIDs{1,aa} '_IOS_KNN_SleepScoringModel.mat'],'KNN_MDL')
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
        oddPatchedREMarray = LinkBinaryEvents_IOS_eLife2020(oddREMArray',[5,0]);
        oddPatchedREMindex = vertcat(oddPatchedREMindex,oddPatchedREMarray');
    end
    % testing data - patch missing REM indeces due to theta band falling off
    for ii = 1:size(evenReshapedREMindex,2)
        evenREMArray = evenReshapedREMindex(:,ii);
        evenPatchedREMarray = LinkBinaryEvents_IOS_eLife2020(evenREMArray',[5,0]);
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
    ConfusionData.KNN.trainYlabels{aa,1} = Yodd.behavState;
    ConfusionData.KNN.trainXlabels{aa,1} = XoddLabels;
    ConfusionData.KNN.testYlabels{aa,1} = Yeven.behavState;
    ConfusionData.KNN.testXlabels{aa,1} = XevenLabels;
    % confusion matrix
    KNN_confMat = figure;
    sgtitle('Support Vector Machine Classifier Confusion Matrix')
    % training data confusion chart
    subplot(1,2,1)
    oddCM = confusionchart(Yodd.behavState,XoddLabels);
    oddCM.ColumnSummary = 'column-normalized';
    oddCM.RowSummary = 'row-normalized';
    oddCM.Title = 'Training Data';
    oddConfVals = oddCM.NormalizedValues;
    oddTotalScores = sum(oddConfVals(:));
    oddKNN_accuracy = (sum(oddConfVals([1,5,9])/oddTotalScores))*100;
    disp(['k-Nearest Neighbor model prediction accuracy (training): ' num2str(oddKNN_accuracy) '%']); disp(' ')
    % testing data confusion chart
    subplot(1,2,2)
    evenCM = confusionchart(Yeven.behavState,XevenLabels);
    evenCM.ColumnSummary = 'column-normalized';
    evenCM.RowSummary = 'row-normalized';
    evenCM.Title = 'Testing Data';
    evenConfVals = evenCM.NormalizedValues;
    evenTotalScores = sum(evenConfVals(:));
    evenKNN_accuracy = (sum(evenConfVals([1,5,9])/evenTotalScores))*100;
    disp(['k-Nearest Neightbor model prediction accuracy (testing): ' num2str(evenKNN_accuracy) '%']); disp(' ')
    % save model and figure
    savefig(KNN_confMat,[dirpath animalIDs{1,aa} '_IOS_KNN_ConfusionMatrix']);
    close(KNN_confMat)
    
    %% Naive Bayes classifier
    disp('Training naive Bayes Classifier...'); disp(' ')
    NB_MDL = fitcnb(Xodd,Yodd,'ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
    % save model in desired location
    save([dirpath animalIDs{1,aa} '_IOS_NB_SleepScoringModel.mat'],'NB_MDL')
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
        oddPatchedREMarray = LinkBinaryEvents_IOS_eLife2020(oddREMArray',[5,0]);
        oddPatchedREMindex = vertcat(oddPatchedREMindex,oddPatchedREMarray');
    end
    % testing data - patch missing REM indeces due to theta band falling off
    for ii = 1:size(evenReshapedREMindex,2)
        evenREMArray = evenReshapedREMindex(:,ii);
        evenPatchedREMarray = LinkBinaryEvents_IOS_eLife2020(evenREMArray',[5,0]);
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
    ConfusionData.NB.trainYlabels{aa,1} = Yodd.behavState;
    ConfusionData.NB.trainXlabels{aa,1} = XoddLabels;
    ConfusionData.NB.testYlabels{aa,1} = Yeven.behavState;
    ConfusionData.NB.testXlabels{aa,1} = XevenLabels;
    % confusion matrix
    NB_confMat = figure;
    sgtitle('Naive Bayes Classifier Confusion Matrix')
    % training data confusion chart
    subplot(1,2,1)
    oddCM = confusionchart(Yodd.behavState,XoddLabels);
    oddCM.ColumnSummary = 'column-normalized';
    oddCM.RowSummary = 'row-normalized';
    oddCM.Title = 'Training Data';
    oddConfVals = oddCM.NormalizedValues;
    oddTotalScores = sum(oddConfVals(:));
    oddNB_accuracy = (sum(oddConfVals([1,5,9])/oddTotalScores))*100;
    disp(['Naive Bayes model prediction accuracy (training): ' num2str(oddNB_accuracy) '%']); disp(' ')
    % testing data confusion chart
    subplot(1,2,2)
    evenCM = confusionchart(Yeven.behavState,XevenLabels);
    evenCM.ColumnSummary = 'column-normalized';
    evenCM.RowSummary = 'row-normalized';
    evenCM.Title = 'Testing Data';
    evenConfVals = evenCM.NormalizedValues;
    evenTotalScores = sum(evenConfVals(:));
    evenNB_accuracy = (sum(evenConfVals([1,5,9])/evenTotalScores))*100;
    disp(['Naive Bayes model prediction accuracy (testing): ' num2str(evenNB_accuracy) '%']); disp(' ')
    % save model and figure
    savefig(NB_confMat,[dirpath animalIDs{1,aa} '_IOS_NB_ConfusionMatrix']);
    close(NB_confMat)
    cd(startingDirectory)
end
% save confusion matrix results
dirpath = [startingDirectory '\Summary Figures and Structures\Confusion Matricies\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
saveLoc = [startingDirectory '\Summary Figures and Structures\Confusion Matricies\'];
save([saveLoc 'ConfusionData.mat'],'ConfusionData')

end
