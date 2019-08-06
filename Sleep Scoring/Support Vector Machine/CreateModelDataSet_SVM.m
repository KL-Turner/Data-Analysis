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

for a = 1:size(procDataFileIDs,1)    procDataFileID = procDataFileIDs(a,:);
    modelDataSetID = [procDataFileID(1:end-12) 'ModelData.mat'];
    if ~exist(modelDataSetID)
        load(procDataFileID)
        disp(['Creating model data set for ' procDataFileID '(' num2str(a) '/' num2str(size(procDataFileIDs,1)) ')' ]); disp(' ')
        %% Create table to send into SVM model
        variableNames = {'maxLH_CBV', 'maxRH_CBV', 'maxLH_Delta', 'maxRH_Delta', 'maxLH_Theta', 'maxRH_Theta',...
            'maxLH_Gamma', 'maxRH_Gamma', 'numWhiskEvents', 'numForceEvents', 'avgEMG', 'avgHeartRate'};
        % pre-allocation
        maxLH_CBV_column = zeros(180,1);
        maxRH_CBV_column = zeros(180,1);
        maxLH_Delta_column = zeros(180,1);
        maxRH_Delta_column = zeros(180,1);
        maxLH_Theta_column = zeros(180,1);
        maxRH_Theta_column = zeros(180,1);
        maxLH_Gamma_column = zeros(180,1);
        maxRH_Gamma_column = zeros(180,1);
        numWhiskEvents_column = zeros(180,1);
        numForceEvents_column = zeros(180,1);
        avgEMG_column = zeros(180,1);
        avgHeartRate_column = zeros(180,1);
        % extract relevant parameters from each epoch
        for b = 1:length(maxLH_CBV_column)
            maxLH_CBV_column(b,1) = round(min(ProcData.sleep.parameters.CBV.LH{b,1})*100,1);
            maxRH_CBV_column(b,1) = round(min(ProcData.sleep.parameters.CBV.RH{b,1})*100,1);
            maxLH_Delta_column(b,1) = round(max(ProcData.sleep.parameters.deltaBandPower.LH{b,1}),1);
            maxRH_Delta_column(b,1) = round(max(ProcData.sleep.parameters.deltaBandPower.RH{b,1}),1);
            maxLH_Theta_column(b,1) = round(max(ProcData.sleep.parameters.thetaBandPower.LH{b,1}),1);
            maxRH_Theta_column(b,1) = round(max(ProcData.sleep.parameters.thetaBandPower.RH{b,1}),1);
            maxLH_Gamma_column(b,1) = round(max(ProcData.sleep.parameters.gammaBandPower.LH{b,1}),1);
            maxRH_Gamma_column(b,1) = round(max(ProcData.sleep.parameters.gammaBandPower.RH{b,1}),1);
            numWhiskEvents_column(b,1) = sum(ProcData.sleep.parameters.binWhiskerAngle{b,1});
            numForceEvents_column(b,1) = sum(ProcData.sleep.parameters.binForceSensor{b,1});
            TF = ~isinf(ProcData.sleep.parameters.EMG{b,1});
            EMG2 = ProcData.sleep.parameters.EMG{b,1}(TF);
            avgEMG_column(b,1) = round(mean(EMG2),1);
            avgHeartRate_column(b,1) = round(mean(ProcData.sleep.parameters.heartRate{b,1}),1);
        end
        % create table
        paramsTable = table(maxLH_CBV_column, maxRH_CBV_column, maxLH_Delta_column, maxRH_Delta_column,...
            maxLH_Theta_column, maxRH_Theta_column, maxLH_Gamma_column, maxRH_Gamma_column,...
            numWhiskEvents_column, numForceEvents_column, avgEMG_column, avgHeartRate_column,...
            'VariableNames', variableNames);
        save(modelDataSetID, 'paramsTable')
    end
end

end
