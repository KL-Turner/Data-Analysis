function [] = CreateTrainingDataSet_SVM(RestingBaselines)
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

procDataFileIDs = uigetfile('*_ProcData.mat','Multiselect','on');
for a = 1:length(procDataFileIDs)
    procDataFileID = procDataFileIDs{1,a};
    modelDataFileID = [procDataFileID(1:end-12) 'ModelData.mat'];
    trainingDataFileID = [procDataFileID(1:end-12) 'TrainingData.mat'];
    if ~exist(trainingDataFileID, 'file')
        disp(['Loading ' procDataFileID ' for manual sleep scoring.' ]); disp(' ')
        load(procDataFileID)
        load(modelDataFileID)
        [figHandle] = GenerateSingleFigures_SVM(procDataFileID, RestingBaselines);
        trialDuration = ProcData.notes.trialDuration_sec;
        numBins = trialDuration/5;
        
        behavioralState = cell(180,1);
        for b = 1:numBins
            global buttonState %#ok<TLEV>
            buttonState = 0;
            subplot(7,1,3)
            ylimits = ylim;
            yMax = ylimits(2);
            yInds = ones(1,5)*yMax*1.2;
            xStartVal = (b*5)-4;
            xEndVal = b*5;
            xInds = xStartVal:1:xEndVal;
            h1 = scatter(xInds,yInds);
%             if b <= 60
%                 xlim([1 300])
%             elseif b >= 61 && b <= 120
%                 xlim([300 600])
%             elseif b >= 121 && b <= 180
%                 xlim([600 900])
%             end
            
            [updatedGUI] = SelectBehavioralStateGUI_SVM;
            while buttonState == 0
                drawnow()
                if buttonState == 1
                    guiResults = guidata(updatedGUI);
                    if guiResults.togglebutton1.Value == true
                        behavioralState{b,1} = 'Not Sleep';
                    elseif guiResults.togglebutton2.Value == true
                        behavioralState{b,1} = 'NREM Sleep';
                    elseif guiResults.togglebutton3.Value == true
                        behavioralState{b,1} = 'REM Sleep';
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
        paramsTable.behavState = behavioralState;
        trainingTable = paramsTable;
        save(trainingDataFileID, 'trainingTable')
    else
        disp([trainingDataFileID ' already exists. Continuing...']); disp(' ')
    end
end

end
