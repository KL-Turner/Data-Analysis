function [] = PupilSleepModelAccuracy_Pupil(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

resultsStructA = 'Results_SleepModel';
load(resultsStructA);
resultsStructB = 'Results_PupilSleepModel';
load(resultsStructB);
animalIDs = fieldnames(Results_SleepModel);
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
        for cc = 1:length(Results_SleepModel.(animalID).(modelName).testXlabels)
            holdXlabels.physio.(modelName) = vertcat(holdXlabels.physio.(modelName),Results_SleepModel.(animalID).(modelName).testXlabels{cc,1});
            holdYlabels.physio.(modelName) = vertcat(holdYlabels.physio.(modelName),Results_SleepModel.(animalID).(modelName).testYlabels{cc,1});
        end
    end
end
for dd = 1:length(modelNames)
    modelName = modelNames{1,dd};
    holdXlabels.pupil.(modelName) = Results_PupilSleepModel.(modelName).testXlabels;
    holdYlabels.pupil.(modelName) = Results_PupilSleepModel.(modelName).testYlabels;
end
% confusion matrix
figure;
sgtitle('Standard physiological sleep scoring classification model')
subplot(2,3,1)
cm = confusionchart(holdYlabels.physio.SVM,holdXlabels.physio.SVM);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
cm.Title = {'Support Vector Machine',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
subplot(2,3,2)
cm = confusionchart(holdYlabels.physio.EC,holdXlabels.physio.EC);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
cm.Title = {'Ensemble Classifier',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
subplot(2,3,3)
cm = confusionchart(holdYlabels.physio.DT,holdXlabels.physio.DT);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
cm.Title = {'Decision Tree',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
subplot(2,3,4)
cm = confusionchart(holdYlabels.physio.RF,holdXlabels.physio.RF);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
cm.Title = {'Bootstrap Aggregate Random Forest',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
subplot(2,3,5)
cm = confusionchart(holdYlabels.physio.KNN,holdXlabels.physio.KNN);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
cm.Title = {'K-Nearest Neighbor',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
subplot(2,3,6)
cm = confusionchart(holdYlabels.physio.NB,holdXlabels.physio.NB);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
cm.Title = {'Naive Bayes',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
%% confusion matrix
figure;
sgtitle('Pupil/Whisker sleep scoring classification model')
subplot(2,3,1)
cm = confusionchart(holdYlabels.pupil.SVM,holdXlabels.pupil.SVM);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
cm.Title = {'Support Vector Machine',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
subplot(2,3,2)
cm = confusionchart(holdYlabels.pupil.EC,holdXlabels.pupil.EC);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
cm.Title = {'Ensemble Classifier',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
subplot(2,3,3)
cm = confusionchart(holdYlabels.pupil.DT,holdXlabels.pupil.DT);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
cm.Title = {'DT unseen data confusion matrix',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
subplot(2,3,4)
cm = confusionchart(holdYlabels.pupil.RF,holdXlabels.pupil.RF);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
cm.Title = {'Bootstrap Aggregate Random Forest',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
subplot(2,3,5)
cm = confusionchart(holdYlabels.pupil.KNN,holdXlabels.pupil.KNN);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
cm.Title = {'K-Nearest Neighbor',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
subplot(2,3,6)
cm = confusionchart(holdYlabels.pupil.NB,holdXlabels.pupil.NB);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
cm.Title = {'Naive Bayes',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
%%
% % save figure(s)
% if saveFigs == true
%     dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Stimulus-evoked Pupil Area' delim];
%     if ~exist(dirpath,'dir')
%         mkdir(dirpath);
%     end
%     savefig(summaryFigure2,[dirpath 'Stimulus_PupilArea']);
%     set(summaryFigure2,'PaperPositionMode','auto');
%     print('-painters','-dpdf','-bestfit',[dirpath 'Stimulus_PupilArea'])
% end

end
