function [] = CreateModelDataSet_IOS(procDataFileIDs,imagingType)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Arrange data into a table of most-relevant parameters for model training/classification
%________________________________________________________________________________________________________________________

for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    disp(['Creating model data set for ProcData file (' num2str(aa) '/' num2str(size(procDataFileIDs,1)) ')']); disp(' ')
    modelDataSetID = [procDataFileID(1:end - 12) 'ModelData.mat'];
    load(procDataFileID)
    % pre-allocation
    maxCortDeltaColumn = zeros(180,1);
    maxCortBetaColumn = zeros(180,1);
    maxCortGammaColumn = zeros(180,1);
    maxHippThetaColumn = zeros(180,1);
    numWhiskEventsColumn = zeros(180,1);
    medEMGColumn = zeros(180,1);
    if strcmpi(imagingType,'bilateral') == true
        avgHeartRateColumn = zeros(180,1);
    end
    % extract relevant parameters from each epoch
    for bb = 1:length(maxCortDeltaColumn)
        % cortical delta
        maxLHcortDelta = mean(cell2mat(ProcData.sleep.parameters.cortical_LH.specDeltaBandPower{bb,1}));
        maxRHcortDelta = mean(cell2mat(ProcData.sleep.parameters.cortical_RH.specDeltaBandPower{bb,1}));
        if strcmpi(imagingType,'GCaMP') == true
            maxCortDeltaColumn(bb,1) = maxRHcortDelta;
        else
            if maxLHcortDelta >= maxRHcortDelta
                maxCortDeltaColumn(bb,1) = maxLHcortDelta;
            else
                maxCortDeltaColumn(bb,1) = maxRHcortDelta;
            end
        end
        % cortical beta
        maxLHcortBeta = mean(cell2mat(ProcData.sleep.parameters.cortical_LH.specBetaBandPower{bb,1}));
        maxRHcortBeta = mean(cell2mat(ProcData.sleep.parameters.cortical_RH.specBetaBandPower{bb,1}));
        if strcmpi(imagingType,'GCaMP') == true
            maxCortBetaColumn(bb,1) = maxRHcortBeta;
        else
            if maxLHcortBeta >= maxRHcortBeta
                maxCortBetaColumn(bb,1) = maxLHcortBeta;
            else
                maxCortBetaColumn(bb,1) = maxRHcortBeta;
            end
        end
        % cortical gamma
        maxLHcortGamma = mean(cell2mat(ProcData.sleep.parameters.cortical_LH.specGammaBandPower{bb,1}));
        maxRHcortGamma = mean(cell2mat(ProcData.sleep.parameters.cortical_RH.specGammaBandPower{bb,1}));
        if strcmpi(imagingType,'GCaMP') == true
            maxCortGammaColumn(bb,1) = maxRHcortGamma;
        else
            if maxLHcortGamma >= maxRHcortGamma
                maxCortGammaColumn(bb,1) = maxLHcortGamma;
            else
                maxCortGammaColumn(bb,1) = maxRHcortGamma;
            end
        end
        % hippocampal theta
        maxHippThetaColumn(bb,1) = mean(cell2mat(ProcData.sleep.parameters.hippocampus.specThetaBandPower{bb,1}));
        % number of binarized whisking events
        numWhiskEventsColumn(bb,1) = sum(ProcData.sleep.parameters.binWhiskerAngle{bb,1});
        % average of the log of the EMG profile
        EMG = ProcData.sleep.parameters.EMG{bb,1};
        medEMGColumn(bb,1) = median(EMG);
        % average heart rate
        if strcmpi(imagingType,'GCaMP') == false
            avgHeartRateColumn(bb,1) = round(mean(ProcData.sleep.parameters.heartRate{bb,1}),1);
        end
    end
    % create table
    if strcmpi(imagingType,'GCaMP') == true
        variableNames = {'maxCortDelta','maxCortBeta','maxCortGamma','maxHippTheta','numWhiskEvents','avgEMG'};
        paramsTable = table(maxCortDeltaColumn,maxCortBetaColumn,maxCortGammaColumn,maxHippThetaColumn,numWhiskEventsColumn,medEMGColumn,'VariableNames',variableNames);
    elseif strcmpi(imagingType,'bilateral') == true
        variableNames = {'maxCortDelta','maxCortBeta','maxCortGamma','maxHippTheta','numWhiskEvents','avgEMG','avgHeartRate'};
        paramsTable = table(maxCortDeltaColumn,maxCortBetaColumn,maxCortGammaColumn,maxHippThetaColumn,numWhiskEventsColumn,medEMGColumn,avgHeartRateColumn,'VariableNames',variableNames);
    end
    save(modelDataSetID,'paramsTable')
end

end
