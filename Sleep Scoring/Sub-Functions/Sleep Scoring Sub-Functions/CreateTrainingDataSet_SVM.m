function [] = CreateTrainingDataSet_SVM(procDataFileIDs, RestingBaselines)
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

for a = 1:size(procDataFileIDs, 1)
    procDataFileID = procDataFileIDs(a,:);
    tableFileID = [procDataFileID(1:end-12) 'TrainingTable.mat'];
    if ~exist(tableFileID, 'file')
        disp(['Loading ' procDataFileID ' for manual sleep scoring.' ]); disp(' ')
        load(procDataFileID)
        
        [figHandle] = GenerateSingleFigures_SVM(procDataFileID, RestingBaselines);
        trialDuration = ProcData.notes.trialDuration_sec;
        numBins = trialDuration/5;
        
        behavioralState = cell(180,1);
        for b = 1:numBins
            global buttonState %#ok<TLEV>
            buttonState = 0;
            subplot(6,1,3)
            ylimits = ylim;
            yMax = ylimits(2);
            yInds = ones(1,5)*yMax*1.2;
            xStartVal = (b*5)-4;
            xEndVal = b*5;
            xInds = xStartVal:1:xEndVal;
            h1 = scatter(xInds,yInds);
            if b <= 60
                xlim([1 300])
            elseif b >= 61 && b <= 120
                xlim([300 600])
            elseif b >= 121 && b <= 180
                xlim([600 900])
            end
            
            [updatedGUI] = BehaviorStateButtonGUI_SVM;
            while buttonState == 0
                drawnow()
                if buttonState == 1
                    guiResults = guidata(updatedGUI);
                    if guiResults.pushbutton1.Value == true
                        behavioralState{b,1} = 'Awake';
                    elseif guiResults.pushbutton2.Value == true
                        behavioralState{b,1} = 'Transition';
                    elseif guiResults.pushbutton3.Value == true
                        behavioralState{b,1} = 'NREM';
                    elseif guiResults.pushbutton4.Value == true
                        behavioralState{b,1} = 'REM';
                    else
                        disp('No button pressed'); disp(' ')
                        keyboard
                    end
                    close(updatedGUI)
                    break;
                end
                ...
            end
        delete(h1)
        end
        close(figHandle)
        
        variableNames = {'maxLH_CBV', 'maxRH_CBV', 'maxLH_Delta', 'maxRH_Delta', 'maxLH_Theta', 'maxRH_Theta',...
            'maxLH_Gamma', 'maxRH_Gamma', 'numWhiskEvents', 'numForceEvents', 'avgEMG', 'avgHeartRate', 'behavState'};
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
        
        for c = 1:length(maxLH_CBV_column)
            maxLH_CBV_column(c,1) = min(ProcData.sleep.parameters.CBV.LH{c,1});
            maxRH_CBV_column(c,1) = min(ProcData.sleep.parameters.CBV.RH{c,1});
            maxLH_Delta_column(c,1) = max(ProcData.sleep.parameters.deltaBandPower.LH{c,1});
            maxRH_Delta_column(c,1) = max(ProcData.sleep.parameters.deltaBandPower.RH{c,1});
            maxLH_Theta_column(c,1) = max(ProcData.sleep.parameters.thetaBandPower.LH{c,1});
            maxRH_Theta_column(c,1) = max(ProcData.sleep.parameters.thetaBandPower.RH{c,1});
            maxLH_Gamma_column(c,1) = max(ProcData.sleep.parameters.gammaBandPower.LH{c,1});
            maxRH_Gamma_column(c,1) = max(ProcData.sleep.parameters.gammaBandPower.RH{c,1});
            numWhiskEvents_column(c,1) = sum(ProcData.sleep.parameters.binWhiskerAngle{c,1});
            numForceEvents_column(c,1) = sum(ProcData.sleep.parameters.binForceSensor{c,1});
            TF = ~isinf(ProcData.sleep.parameters.EMG{c,1});
            EMG2 = ProcData.sleep.parameters.EMG{c,1}(TF);
            avgEMG_column(c,1) = mean(EMG2);
            avgHeartRate_column(c,1) = mean(ProcData.sleep.parameters.heartRate{c,1});
        end
        
        T = table(maxLH_CBV_column, maxRH_CBV_column, maxLH_Delta_column, maxRH_Delta_column,...
            maxLH_Theta_column, maxRH_Theta_column, maxLH_Gamma_column, maxRH_Gamma_column,...
            numWhiskEvents_column, numForceEvents_column, avgEMG_column, avgHeartRate_column,...
            behavioralState, 'VariableNames', variableNames);
        save(tableFileID, 'T')
    end
    disp([tableFileID ' already exists. Continuing...']); disp(' ')
end

end
