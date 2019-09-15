function [] = CreateModelDataSet_SVM(procDataFileIDs)
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
%   Last Revised: July 26th, 2019
%________________________________________________________________________________________________________________________

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    modelDataSetID = [procDataFileID(1:end-12) 'ModelData.mat'];
        load(procDataFileID)
        disp(['Creating model data set for ' procDataFileID '... (' num2str(a) '/' num2str(size(procDataFileIDs,1)) ')' ]); disp(' ')
        %% Create table to send into SVM model
        variableNames = {'minRefl', 'maxCortDelta', 'maxCortGamma', 'maxHippTheta', 'numWhiskEvents', 'numForceEvents', 'avgEMG', 'avgHeartRate'};
        % pre-allocation
        minRefl_column = zeros(180,1);
        maxCortDelta_column = zeros(180,1);
        maxCortGamma_column = zeros(180,1);
        maxHippTheta_column = zeros(180,1);
        numWhiskEvents_column = zeros(180,1);
        numForceEvents_column = zeros(180,1);
        avgEMG_column = zeros(180,1);
        avgHeartRate_column = zeros(180,1);
        % extract relevant parameters from each epoch
        for b = 1:length(minRefl_column)
            % CBV Reflectance
            minLHrefl = round(min(ProcData.sleep.parameters.CBV.LH{b,1})*100,1);
            minRHrefl = round(min(ProcData.sleep.parameters.CBV.RH{b,1})*100,1);
            if minLHrefl <= minRHrefl 
                minRefl_column(b,1) = minLHrefl;
            else
                minRefl_column(b,1) = minRHrefl;
            end
            % Cortical delta
            maxLHcortDelta = round(max(ProcData.sleep.parameters.cortical_LH.deltaBandPower{b,1}),1);
            maxRHcortDelta = round(max(ProcData.sleep.parameters.cortical_RH.deltaBandPower{b,1}),1);
            if maxLHcortDelta >= maxRHcortDelta
                maxCortDelta_column(b,1) = maxLHcortDelta;
            else
                maxCortDelta_column(b,1) = maxRHcortDelta;
            end
            % Cortical gamma
            maxLHcortGamma = round(max(ProcData.sleep.parameters.cortical_LH.gammaBandPower{b,1}),1);
            maxRHcortGamma = round(max(ProcData.sleep.parameters.cortical_RH.gammaBandPower{b,1}),1);
            if maxLHcortGamma >= maxRHcortGamma
                maxCortGamma_column(b,1) = maxLHcortGamma;
            else
                maxCortGamma_column(b,1) = maxRHcortGamma;
            end
            % Hippocampal theta
            maxHippTheta_column(b,1) = round(max(ProcData.sleep.parameters.hippocampus.thetaBandPower{b,1}),1);
            % number of binarized whisking events
            numWhiskEvents_column(b,1) = sum(ProcData.sleep.parameters.binWhiskerAngle{b,1});
            % number of binarized force sensor events
            numForceEvents_column(b,1) = sum(ProcData.sleep.parameters.binForceSensor{b,1});
            % average of the log of the EMG profile
            TF = ~isinf(ProcData.sleep.parameters.EMG{b,1});
            EMG2 = ProcData.sleep.parameters.EMG{b,1}(TF);
            avgEMG_column(b,1) = round(mean(EMG2),1);
            % average heart rate
            avgHeartRate_column(b,1) = round(mean(ProcData.sleep.parameters.heartRate{b,1}),1);
        end
        % create table
        paramsTable = table(minRefl_column, maxCortDelta_column, maxCortGamma_column, maxHippTheta_column, ...
            numWhiskEvents_column, numForceEvents_column, avgEMG_column, avgHeartRate_column,...
            'VariableNames', variableNames);
        save(modelDataSetID, 'paramsTable')
end

end
