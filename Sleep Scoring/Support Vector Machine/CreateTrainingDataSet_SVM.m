function [] = CreateTrainingDataSet_SVM(procDataFileIDs,RestingBaselines)
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

for a = 1:length(procDataFileIDs)
    procDataFileID = procDataFileIDs(a,:);
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
            
            xStartVal = (b*5)-4;
            xEndVal = b*5;
            xInds = xStartVal:1:xEndVal;
                 
%             subplot(6,1,1)
%             yyaxis left
%             ylimits1 = ylim;
%             yMax1 = ylimits1(2);
%             yInds1 = ones(1,5)*yMax1*1.2;
%             hold on
%             h1 = scatter(xInds,yInds1);
%             
%             subplot(6,1,2)
%             yyaxis left
%             ylimits2 = ylim;
%             yMax2 = ylimits2(2);
%             yInds2 = ones(1,5)*yMax2*1.2;
%             ylim([-20 yInds2(1)])
%             hold on
%             h2 = scatter(xInds,yInds2);
%             yyaxis right
%             ylim([5 15])
            
            subplot(6,1,3)
            yyaxis left
            ylimits3 = ylim;
            yMax3 = ylimits3(2);
            yInds3 = ones(1,5)*yMax3*1.2;
            hold on
            h3 = scatter(xInds,yInds3);
%             
%             subplot(6,1,4)
%             yyaxis left
%             ylimits4 = ylim;
%             yMax4 = ylimits4(2);
%             yInds4 = ones(1,5)*110;
%             hold on
%             h4 = scatter(xInds,yInds4);
%             
%             subplot(6,1,5)
%             yyaxis left
%             ylimits5 = ylim;
%             yMax5 = ylimits5(2);
%             yInds5 = ones(1,5)*110;
%             hold on
%             h5 = scatter(xInds,yInds5);
%             
%             subplot(6,1,6)
%             yyaxis left
%             ylimits6 = ylim;
%             yMax6 = ylimits6(2);
%             yInds6 = ones(1,5)*110;
%             hold on
%             h6 = scatter(xInds,yInds6);
            
            if b <= 60
                xlim([1 300])
            elseif b >= 61 && b <= 120
                xlim([300 600])
            elseif b >= 121 && b <= 180
                xlim([600 900])
            end
            
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
%         delete(h1)
%         delete(h2)
        delete(h3)
%         delete(h4)
%         delete(h5)
%         delete(h6)
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
