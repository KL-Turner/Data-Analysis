%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Calculate the average correlation coefficient of the different behavioral states
%________________________________________________________________________________________________________________________
%
%   Inputs: none
%
%   Outputs: Generates summary figures saved to C: drive Documents folder
%
%   Last Revised: Oct 1st, 2019
%________________________________________________________________________________________________________________________
clear
clc

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111'};
driveLetters = {'E','E','E','F','F','F','D','D','D'};
behavFields = {'Whisk','Rest','NREM','REM','Unstim','AllData'};
neuralBands = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','muaPower'};
baselineTypes = {'manualSelection','setDuration','entireDuration'};
fileSets = {'fileSetA','fileSetB'};
hemDataTypes = {'LH','RH'};
cbv_dataTypes = {'CBV','CBV_HbT'};
colorbrewer_setA_colorA = [0.520000 0.520000 0.510000];
colorbrewer_setA_colorB = [(31/256) (120/256) (180/256)];
colorbrewer_setA_colorC = [(255/256) (0/256) (115/256)];
colorbrewer_setA_colorD = [(51/256) (160/256) (44/256)];
colorbrewer_setA_colorE = [(255/256) (140/256) (0/256)];
colorbrewer_setA_colorF = [0.750000 0.000000 1.000000];

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\' animalID '\Combined Imaging\'];
    cd(dataPath)
    load([animalID '_AnalysisResults.mat']);
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        for x = 1:length(fileSets)
            fileSet = fileSets{1,x};
            for c = 1:length(cbv_dataTypes)
                cbv_dataType = cbv_dataTypes{1,c};
                for y = 1:length(hemDataTypes)
                    hemDataType = hemDataTypes{1,y};
                    for z = 1:length(neuralBands)
                        neuralBand = neuralBands{1,z};
                        if strcmp(behavField,'Rest') == true || strcmp(behavField,'Whisk') == true
                            for d = 1:length(baselineTypes)
                                baselineType = baselineTypes{1,d};
                                data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).(baselineType).R(a,1) = mean(AnalysisResults.HRF_Predictions.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).(baselineType).R);
                                data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).(baselineType).allR{a,1} = AnalysisResults.HRF_Predictions.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).(baselineType).R;
                                data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).(baselineType).R2(a,1) = mean(AnalysisResults.HRF_Predictions.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).(baselineType).R2);
                                data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).(baselineType).allR2{a,1} = AnalysisResults.HRF_Predictions.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).(baselineType).R2;
                            end
                        else
                            data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).R(a,1) = mean(AnalysisResults.HRF_Predictions.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).entireDuration.R);
                            data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).allR{a,1} = AnalysisResults.HRF_Predictions.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).entireDuration.R;
                            data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).R2(a,1) = mean(AnalysisResults.HRF_Predictions.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).entireDuration.R2);
                            data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).allR2{a,1} = AnalysisResults.HRF_Predictions.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).entireDuration.R2;
                        end
                    end
                end
            end
        end
    end
end

for b = 1:length(behavFields)
    behavField = behavFields{1,b};
    for x = 1:length(fileSets)
        fileSet = fileSets{1,x};
        for c = 1:length(cbv_dataTypes)
            cbv_dataType = cbv_dataTypes{1,c};
            for y = 1:length(hemDataTypes)
                hemDataType = hemDataTypes{1,y};
                for z = 1:length(neuralBands)
                    neuralBand = neuralBands{1,z};
                    if strcmp(behavField,'Rest') == true || strcmp(behavField,'Whisk') == true
                        for d = 1:length(baselineTypes)
                            baselineType = baselineTypes{1,d};
                            data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).(baselineType).meanR = mean(data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).(baselineType).R);
                            data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).(baselineType).stdR = std(data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).(baselineType).R,0,2);
                            data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).(baselineType).meanR2 = mean(data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).(baselineType).R2);
                            data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).(baselineType).stdR2 = std(data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).(baselineType).R2,0,2);
                        end
                    else
                        data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).meanR = mean(data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).R);
                        data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).stdR = std(data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).R,0,2);
                        data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).meanR2 = mean(data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).R2);
                        data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).stdR2 = std(data.(behavField).(fileSet).(cbv_dataType).(hemDataType).(neuralBand).R2,0,2);
                    end
                end
            end
        end
    end
end

