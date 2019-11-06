function [] = VerifyModelPredictions_SVM(animalIDs,driveLetters,saveFigs,baselineType)
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

modelLocation = 'C:\Users\klt8\Documents\';
cd(modelLocation)
load('SVM_SleepScoringModel.mat')

for a = 1:length(animalIDs)
    baseLoc = [driveLetters{1,a} ':\' animalIDs{1,a} '\Combined Imaging\'];
    cd(baseLoc)
    
    % Character name for the '*_RestingBaselines.mat' structure
    baselinesFileStruct = dir('*_RestingBaselines.mat');
    baselinesFile = {baselinesFileStruct.name}';
    baselinesFileID = char(baselinesFile);
    load(baselinesFileID)
    
    dataLoc = [driveLetters{1,a} ':\' animalIDs{1,a} '\SVM Validation Set'];
    cd(dataLoc)
    
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
    
    AddSleepParameters_SVM(procDataFileIDs, RestingBaselines, baselineType)
    CreateModelDataSet_SVM(procDataFileIDs)
    modelDataFileStruct = dir('*_ModelData.mat');
    modelDataFiles = {modelDataFileStruct.name}';
    modelDataFileIDs = char(modelDataFiles);
    
    PredictBehaviorEvents_SVM(modelDataFileIDs, SVMModel)
    
    scoringResultsDataFileStruct = dir('*_SleepScoringResults.mat');
    scoringResultsDataFiles = {scoringResultsDataFileStruct.name}';
    scoringResultsDataFileID = char(scoringResultsDataFiles);
    load(scoringResultsDataFileID)
    
    trainingDataFileStruct = dir('*_TrainingData.mat');
    trainingDataFiles = {trainingDataFileStruct.name}';
    trainingDataFileIDs = char(trainingDataFiles);
    
    for b = 1:size(trainingDataFileIDs,1)
%         disp(['Generating SVM verification figure ' num2str(b) ' of ' num2str(size(trainingDataFileIDs,1)) '...']); disp(' ')
        trainingDataFileID = trainingDataFileIDs(b,:);
        load(trainingDataFileID)
        procDataFileID = procDataFileIDs(b,:);
        modelDataFileID = [procDataFileID(1:end-12) 'ModelData.mat'];
        c = 1;
        clear fileSVMScores
        for d = 1:length(SVMResults.fileIDs)
            if strcmp(modelDataFileID, SVMResults.fileIDs{d,1}) == true
                fileSVMScores{c,1} = SVMResults.labels{d,1};
                c = c + 1;
            end
        end
        fileManualScores = trainingTable.behavState;
        
        if b == 1
            allSVMScores = fileSVMScores;
            allManualScores = fileManualScores;
        else
            allSVMScores = vertcat(allSVMScores, fileSVMScores);
            allManualScores = vertcat(allManualScores, fileManualScores);
        end
        
        if strcmp(saveFigs, 'y') == true
            svmNotSleepInds = double(strcmp(fileSVMScores, 'Not Sleep'));
            svmNotSleepInds(svmNotSleepInds~=1) = NaN;
            
            svmNREMSleepInds = double(strcmp(fileSVMScores, 'NREM Sleep'));
            svmNREMSleepInds(svmNREMSleepInds~=1) = NaN;
            
            svmREMSleepInds = double(strcmp(fileSVMScores, 'REM Sleep'));
            svmREMSleepInds(svmREMSleepInds~=1) = NaN;
            
            manualNotSleepInds = double(strcmp(fileManualScores, 'Not Sleep'));
            manualNotSleepInds(manualNotSleepInds~=1) = NaN;
            
            manualNREMSleepInds = double(strcmp(fileManualScores, 'NREM Sleep'));
            manualNREMSleepInds(manualNREMSleepInds~=1) = NaN;
            
            manualREMSleepInds = double(strcmp(fileManualScores, 'REM Sleep'));
            manualREMSleepInds(manualREMSleepInds~=1) = NaN;
            
            xInds = 1:5:900;
            
            baselineType = 'manualSelection';
            [figHandle] = GenerateSingleFigures_SVM(procDataFileID, RestingBaselines, baselineType);
            
            subplot(6,1,1)
            yyaxis left
            ylimits1 = ylim;
            yMax1 = ylimits1(2);
            yInds_svmNotSleep = svmNotSleepInds*yMax1*1.2;
            yInds_svmNREMSleep = svmNREMSleepInds*yMax1*1.2;
            yInds_svmREMSleep = svmREMSleepInds*yMax1*1.2;
            yInds_manualNotSleep = manualNotSleepInds*yMax1*1.4;
            yInds_manualNREMSleep = manualNREMSleepInds*yMax1*1.4;
            yInds_manualREMSleep = manualREMSleepInds*yMax1*1.4;
            hold on
            scatter(xInds,yInds_svmNotSleep, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
            scatter(xInds,yInds_svmNREMSleep, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
            scatter(xInds,yInds_svmREMSleep, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
            scatter(xInds,yInds_manualNotSleep, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
            scatter(xInds,yInds_manualNREMSleep, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
            scatter(xInds,yInds_manualREMSleep, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
            
            subplot(6,1,2)
            yyaxis left
            ylimits2 = ylim;
            yMax2 = ylimits2(2);
            yInds_svmNotSleep = svmNotSleepInds*yMax2*1.1;
            yInds_svmNREMSleep = svmNREMSleepInds*yMax2*1.1;
            yInds_svmREMSleep = svmREMSleepInds*yMax2*1.1;
            yInds_manualNotSleep = manualNotSleepInds*yMax2*1.2;
            yInds_manualNREMSleep = manualNREMSleepInds*yMax2*1.2;
            yInds_manualREMSleep = manualREMSleepInds*yMax2*1.2;
            hold on
            scatter(xInds,yInds_svmNotSleep, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
            scatter(xInds,yInds_svmNREMSleep, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
            scatter(xInds,yInds_svmREMSleep, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
            scatter(xInds,yInds_manualNotSleep, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
            scatter(xInds,yInds_manualNREMSleep, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
            scatter(xInds,yInds_manualREMSleep, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
            ylim([-20 80])
            yyaxis right
            ylim([5 15])
            
            subplot(6,1,3)
            yyaxis left
            ylimits3 = ylim;
            yMax3 = ylimits3(2);
            yInds_svmNotSleep = svmNotSleepInds*yMax3*1.2;
            yInds_svmNREMSleep = svmNREMSleepInds*yMax3*1.2;
            yInds_svmREMSleep = svmREMSleepInds*yMax3*1.2;
            yInds_manualNotSleep = manualNotSleepInds*yMax3*1.4;
            yInds_manualNREMSleep = manualNREMSleepInds*yMax3*1.4;
            yInds_manualREMSleep = manualREMSleepInds*yMax3*1.4;
            hold on
            scatter(xInds,yInds_svmNotSleep,'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
            scatter(xInds,yInds_svmNREMSleep, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
            scatter(xInds,yInds_svmREMSleep, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
            scatter(xInds,yInds_manualNotSleep, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
            scatter(xInds,yInds_manualNREMSleep, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
            scatter(xInds,yInds_manualREMSleep, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
            
            
            subplot(6,1,4)
            yyaxis left
            yMax4 = 100;
            yInds_svmNotSleep = svmNotSleepInds*yMax4*1.2;
            yInds_svmNREMSleep = svmNREMSleepInds*yMax4*1.2;
            yInds_svmREMSleep = svmREMSleepInds*yMax4*1.2;
            yInds_manualNotSleep = manualNotSleepInds*yMax4*1.6;
            yInds_manualNREMSleep = manualNREMSleepInds*yMax4*1.6;
            yInds_manualREMSleep = manualREMSleepInds*yMax4*1.6;
            hold on
            scatter(xInds,yInds_svmNotSleep, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
            scatter(xInds,yInds_svmNREMSleep, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
            scatter(xInds,yInds_svmREMSleep, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
            scatter(xInds,yInds_manualNotSleep, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
            scatter(xInds,yInds_manualNREMSleep, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
            scatter(xInds,yInds_manualREMSleep, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
            
            
            subplot(6,1,5)
            yyaxis left
            yMax5 = 100;
            yInds_svmNotSleep = svmNotSleepInds*yMax5*1.2;
            yInds_svmNREMSleep = svmNREMSleepInds*yMax5*1.2;
            yInds_svmREMSleep = svmREMSleepInds*yMax5*1.2;
            yInds_manualNotSleep = manualNotSleepInds*yMax5*1.6;
            yInds_manualNREMSleep = manualNREMSleepInds*yMax5*1.6;
            yInds_manualREMSleep = manualREMSleepInds*yMax5*1.6;
            hold on
            scatter(xInds,yInds_svmNotSleep, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
            scatter(xInds,yInds_svmNREMSleep, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
            scatter(xInds,yInds_svmREMSleep, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
            scatter(xInds,yInds_manualNotSleep, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
            scatter(xInds,yInds_manualNREMSleep, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
            scatter(xInds,yInds_manualREMSleep, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
            
            subplot(6,1,6)
            yyaxis left
            yMax6 = 100;
            yInds_svmNotSleep = svmNotSleepInds*yMax6*1.2;
            yInds_svmNREMSleep = svmNREMSleepInds*yMax6*1.2;
            yInds_svmREMSleep = svmREMSleepInds*yMax6*1.2;
            yInds_manualNotSleep = manualNotSleepInds*yMax6*1.6;
            yInds_manualNREMSleep = manualNREMSleepInds*yMax6*1.6;
            yInds_manualREMSleep = manualREMSleepInds*yMax6*1.6;
            hold on
            scatter(xInds,yInds_svmNotSleep, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
            scatter(xInds,yInds_svmNREMSleep, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
            scatter(xInds,yInds_svmREMSleep, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
            scatter(xInds,yInds_manualNotSleep, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
            scatter(xInds,yInds_manualNREMSleep, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
            scatter(xInds,yInds_manualREMSleep, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
            
            dirpath = [driveLetters{1,a} ':\' animalIDs{1,a} '\SVM Model Validation\Figures\'];
            
            if ~exist(dirpath, 'dir')
                mkdir(dirpath);
            end
            
            savefig(figHandle, [dirpath procDataFileID(1:end-12) 'SVM_ModelValidation']);
            close(figHandle)
        end
    end
    
    confMat = figure;
    cm = confusionchart(allManualScores,allSVMScores);
    cm.ColumnSummary = 'column-normalized';
    cm.RowSummary = 'row-normalized';
    cm.Title = 'T109 SVM Confusion Matrix';
    ConfMat.manualScores = allManualScores;
    ConfMat.svmScore = allSVMScores;
    dirpathA = [driveLetters{1,a} ':\' animalIDs{1,a} '\SVM Validation Set\Figures\'];
    if ~exist(dirpathA, 'dir')
        mkdir(dirpathA);
    end
    save([dirpathA animalIDs{1,a} '_SVM_ConfusionInputs'],'ConfMat')
    dirpathB = [driveLetters{1,a} ':\' animalIDs{1,a} '\SVM Validation Set\'];
    if ~exist(dirpathB, 'dir')
        mkdir(dirpathB);
    end
    savefig(confMat, [dirpathB animalIDs{1,a} '_SVM_ConfusionMatrix']);
    close(confMat)
    
    totalScores = length(allSVMScores);
    positiveScores = sum(strcmp(allSVMScores, allManualScores));
    disp([animalIDs{1,a} ' SVM model prediction success: ' num2str((positiveScores/totalScores)*100) '%']); disp(' ')
end