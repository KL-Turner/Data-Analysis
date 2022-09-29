function [] = PupilSleepModelAccuracy2_Pupil(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

resultsStructA = 'Results_PhysioSleepModel';
load(resultsStructA);
resultsStructB = 'Results_PupilSleepModel';
load(resultsStructB);
animalIDs = fieldnames(Results_PhysioSleepModel);
modelNames = {'SVM','EC','DT','RF','KNN','NB'};
holdXlabels.physio = []; holdYlabels.physio = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    for bb = 1:length(modelNames)
        modelName = modelNames{1,bb};
        if isfield(holdXlabels.physio,modelName) == false
            holdXlabels.physio.(modelName) = {};
            holdYlabels.physio.(modelName) = {};
        end
        for cc = 1:length(Results_PhysioSleepModel.(animalID).(modelName).testXlabels)
            holdXlabels.physio.(modelName) = vertcat(holdXlabels.physio.(modelName),Results_PhysioSleepModel.(animalID).(modelName).testXlabels{cc,1});
            holdYlabels.physio.(modelName) = vertcat(holdYlabels.physio.(modelName),Results_PhysioSleepModel.(animalID).(modelName).testYlabels{cc,1});
        end
    end
end
%
holdXlabels.pupil = []; holdYlabels.pupil = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    for bb = 1:length(modelNames)
        modelName = modelNames{1,bb};
        if isfield(holdXlabels.pupil,modelName) == false
            holdXlabels.pupil.(modelName) = {};
            holdYlabels.pupil.(modelName) = {};
        end
        for cc = 1:length(Results_PupilSleepModel.(animalID).(modelName).testXlabels)
            holdXlabels.pupil.(modelName) = vertcat(holdXlabels.pupil.(modelName),Results_PupilSleepModel.(animalID).(modelName).testXlabels{cc,1});
            holdYlabels.pupil.(modelName) = vertcat(holdYlabels.pupil.(modelName),Results_PupilSleepModel.(animalID).(modelName).testYlabels{cc,1});
        end
    end
end
%% confusion matrix
summaryFigure = figure;
sgtitle('Standard physiological sleep scoring classification model')
subplot(2,3,1)
cm = confusionchart(holdYlabels.physio.SVM,holdXlabels.physio.SVM);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,4])/totalScores))*100,1);
cm.Title = {'Support Vector Machine',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
subplot(2,3,2)
cm = confusionchart(holdYlabels.physio.EC,holdXlabels.physio.EC);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,4])/totalScores))*100,1);
cm.Title = {'Ensemble Classifier',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
subplot(2,3,3)
cm = confusionchart(holdYlabels.physio.DT,holdXlabels.physio.DT);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,4])/totalScores))*100,1);
cm.Title = {'Decision Tree',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
subplot(2,3,4)
cm = confusionchart(holdYlabels.physio.RF,holdXlabels.physio.RF);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,4])/totalScores))*100,1);
cm.Title = {'Bootstrap Aggregate Random Forest',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
subplot(2,3,5)
cm = confusionchart(holdYlabels.physio.KNN,holdXlabels.physio.KNN);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,4])/totalScores))*100,1);
cm.Title = {'K-Nearest Neighbor',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
subplot(2,3,6)
cm = confusionchart(holdYlabels.physio.NB,holdXlabels.physio.NB);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,4])/totalScores))*100,1);
cm.Title = {'Naive Bayes',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Sleep Model Accuracy' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Physio_SleepModel']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Physio_SleepModel'])
end
%% confusion matrix
summaryFigure2 = figure;
sgtitle('Pupil/Whisker sleep scoring classification model')
subplot(2,3,1)
cm = confusionchart(holdYlabels.pupil.SVM,holdXlabels.pupil.SVM);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,4])/totalScores))*100,1);
cm.Title = {'Support Vector Machine',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
subplot(2,3,2)
cm = confusionchart(holdYlabels.pupil.EC,holdXlabels.pupil.EC);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,4])/totalScores))*100,1);
cm.Title = {'Ensemble Classifier',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
subplot(2,3,3)
cm = confusionchart(holdYlabels.pupil.DT,holdXlabels.pupil.DT);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,4])/totalScores))*100,1);
cm.Title = {'DT unseen data confusion matrix',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
subplot(2,3,4)
cm = confusionchart(holdYlabels.pupil.RF,holdXlabels.pupil.RF);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,4])/totalScores))*100,1);
cm.Title = {'Bootstrap Aggregate Random Forest',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
subplot(2,3,5)
cm = confusionchart(holdYlabels.pupil.KNN,holdXlabels.pupil.KNN);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,4])/totalScores))*100,1);
cm.Title = {'K-Nearest Neighbor',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
subplot(2,3,6)
cm = confusionchart(holdYlabels.pupil.NB,holdXlabels.pupil.NB);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,4])/totalScores))*100,1);
cm.Title = {'Naive Bayes',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Sleep Model Accuracy' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure2,[dirpath 'Pupil_SleepModel']);
    set(summaryFigure2,'PaperPositionMode','auto');
    print('-vector','-dpdf','-bestfit',[dirpath 'Pupil_SleepModel'])
end

end
