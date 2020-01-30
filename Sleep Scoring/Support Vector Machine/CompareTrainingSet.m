function [] = CompareTrainingSet(modelDataFileIDs, RestingBaselines)
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

[animalID,~,~] = GetFileInfo_IOS(modelDataFileIDs(1,:));
load([animalID '_SVM_SleepScoringModel'])
load([animalID '_SVM_SleepScoringResults'])

for a = 1:size(modelDataFileIDs,1)
    disp(['Generating SVM comparison figure ' num2str(a) ' of ' num2str(size(modelDataFileIDs,1)) '...']); disp(' ')
    modelDataFileID = modelDataFileIDs(a,:);
    procDataFileID = [modelDataFileID(1:end-14) '_ProcData.mat'];
    b = 1;
    for c = 1:length(SVMResults.fileIDs)
        if strcmp(modelDataFileID, SVMResults.fileIDs{c,1}) == true
            fileManualScores{b,1} = SVMModel.Y{c,1};
            fileSVMScores{b,1} = SVMResults.labels{c,1};
            b = b + 1;
        end
    end
    
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
    
    [figHandle] = GenerateSingleFigures_SVM(procDataFileID, RestingBaselines);

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
    
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Figures/SVM Training Comparison/'];
    
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    
    savefig(figHandle, [dirpath modelDataFileID(1:end-14) '_SVMcomparison']);
    close(figHandle)
end