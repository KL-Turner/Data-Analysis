function [] = PredictBehaviorEvents_SVM(procDataFileIDs, SVMModel)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs: List of processed data file IDs (char) and the loaded RestingBaselines structure
%
%   Outputs:
%
%   Last Revised: July 26th, 2019
%________________________________________________________________________________________________________________________

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    disp(['Sleep scoring ' procDataFileID ' through support vector machine model...']); disp(' ')
    load(procDataFileID)
    if ~isfield(ProcData.sleep, 'SVM')
        %% Create table to send into SVM model
        variableNames = {'maxLH_CBV', 'maxRH_CBV', 'maxLH_Delta', 'maxRH_Delta', 'maxLH_Theta', 'maxRH_Theta',...
            'maxLH_Gamma', 'maxRH_Gamma', 'numWhiskEvents', 'numForceEvents', 'avgEMG', 'avgHeartRate'};
        % pre-allocation
        maxLH_CBV_column = cell(180,1);
        maxRH_CBV_column = cell(180,1);
        maxLH_Delta_column = cell(180,1);
        maxRH_Delta_column = cell(180,1);
        maxLH_Theta_column = cell(180,1);
        maxRH_Theta_column = cell(180,1);
        maxLH_Gamma_column = cell(180,1);
        maxRH_Gamma_column = cell(180,1);
        numWhiskEvents_column = cell(180,1);
        numForceEvents_column = cell(180,1);
        avgEMG_column = cell(180,1);
        avgHeartRate_column = cell(180,1);
        % extract relevant parameters from each epoch
        for b = 1:length(maxLH_CBV_column)
            maxLH_CBV_column{b,1} = min(ProcData.sleep.parameters.CBV.LH{b,1});
            maxRH_CBV_column{b,1} = min(ProcData.sleep.parameters.CBV.RH{b,1});
            maxLH_Delta_column{b,1} = max(ProcData.sleep.parameters.deltaBandPower.LH{b,1});
            maxRH_Delta_column{b,1} = max(ProcData.sleep.parameters.deltaBandPower.RH{b,1});
            maxLH_Theta_column{b,1} = max(ProcData.sleep.parameters.thetaBandPower.LH{b,1});
            maxRH_Theta_column{b,1} = max(ProcData.sleep.parameters.thetaBandPower.RH{b,1});
            maxLH_Gamma_column{b,1} = max(ProcData.sleep.parameters.gammaBandPower.LH{b,1});
            maxRH_Gamma_column{b,1} = max(ProcData.sleep.parameters.gammaBandPower.RH{b,1});
            numWhiskEvents_column{b,1} = sum(ProcData.sleep.parameters.binWhiskerAngle{b,1});
            numForceEvents_column{b,1} = sum(ProcData.sleep.parameters.binForceSensor{b,1});
            TF = ~isinf(ProcData.sleep.parameters.EMG{b,1});
            EMG2 = ProcData.sleep.parameters.EMG{b,1}(TF);
            avgEMG_column{b,1} = mean(EMG2);
            avgHeartRate_column{b,1} = mean(ProcData.sleep.parameters.heartRate{b,1});
        end
        % create table
        paramsTable = table(maxLH_CBV_column, maxRH_CBV_column, maxLH_Delta_column, maxRH_Delta_column,...
            maxLH_Theta_column, maxRH_Theta_column, maxLH_Gamma_column, maxRH_Gamma_column,...
            numWhiskEvents_column, numForceEvents_column, avgEMG_column, avgHeartRate_column,...
            'VariableNames', variableNames);
        
        %% Obtain the label/score from the model
        [label,score] = predict(SVMModel, paramsTable);
        ProcData.sleep.SVM.label = label;
        ProcData.sleep.SVM.score = score;
        save(procDataFileID, 'ProcData')
    else
        disp([procDataFileID ' already evaluated. Continuing...']); disp(' ')
    end
end
