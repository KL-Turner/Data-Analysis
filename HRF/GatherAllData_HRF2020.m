function [NeuralData,HemoData] = GatherAllData_HRF2020(neuralBand,hemisphere,behavior,RestingBaselines,ScoringResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner - adapted from code written by Aaron T. Winder
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: calculate the hemodynamic response function from neural data
%________________________________________________________________________________________________________________________

% Character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
NeuralData = {};
HemoData = {};
if strcmp(behavior,'All') == true
    % load each file and put processed data into each structure
    for aa = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(aa,:);
        load(procDataFileID,'-mat')
        [~,fileDate,~] = GetFileInfo_HRF2020(procDataFileID);
        strDay = ConvertDate_HRF2020(fileDate);
        if strcmp(hemisphere(1:4),'cort') == true
            NeuralData{aa,1} = (ProcData.data.(['cortical_' hemisphere(end - 1:end)]).(neuralBand) - RestingBaselines.manualSelection.(['cortical_' hemisphere(end - 1:end)]).(neuralBand).(strDay))./RestingBaselines.manualSelection.(['cortical_' hemisphere(end - 1:end)]).(neuralBand).(strDay); %#ok<*AGROW>
        elseif strcmp(hemisphere(1:4),'hipp') == true
            NeuralData{aa,1} = (ProcData.data.hippocampus.(neuralBand) - RestingBaselines.manualSelection.hippocampus.(neuralBand).(strDay))./RestingBaselines.manualSelection.hippocampus.(neuralBand).(strDay);
        end
        HemoData{aa,1} = ProcData.data.CBV_HbT.(['adj' hemisphere(end - 1:end)]);
    end
elseif strcmp(behavior,'Alert') == true
    dd = 1;
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,alertDataFileDate,alertDataFileID] = GetFileInfo_HRF2020(procDataFileID);
        strDay = ConvertDate_HRF2020(alertDataFileDate);
        scoringLabels = [];
        for cc = 1:length(ScoringResults.fileIDs)
            if strcmp(alertDataFileID,ScoringResults.fileIDs{cc,1}) == true
                scoringLabels = ScoringResults.labels{cc,1};
            end
        end
        % check labels to match arousal state
        if sum(strcmp(scoringLabels,'Not Sleep')) > 144   % 36 bins (180 total) or 3 minutes of sleep
            load(procDataFileID)
            % puffs = ProcData.data.solenoids.LPadSol;
            % % don't include trials with stimulation
            % if isempty(puffs) == true
            if strcmp(hemisphere(1:4),'cort') == true
                NeuralData{dd,1} = (ProcData.data.(['cortical_' hemisphere(end - 1:end)]).(neuralBand) - RestingBaselines.manualSelection.(['cortical_' hemisphere(end - 1:end)]).(neuralBand).(strDay))./RestingBaselines.manualSelection.(['cortical_' hemisphere(end - 1:end)]).(neuralBand).(strDay); %#ok<*AGROW>
            elseif strcmp(hemisphere(1:4),'hipp') == true
                NeuralData{dd,1} = (ProcData.data.hippocampus.(neuralBand) - RestingBaselines.manualSelection.hippocampus.(neuralBand).(strDay))./RestingBaselines.manualSelection.hippocampus.(neuralBand).(strDay);
            end
            HemoData{dd,1} = ProcData.data.CBV_HbT.(['adj' hemisphere(end - 1:end)]);
            dd = dd + 1;
            % end
        end
    end
elseif strcmp(behavior,'Asleep') == true
    gg = 1;
    for ee = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(ee,:);
        [~,alertDataFileDate,alertDataFileID] = GetFileInfo_HRF2020(procDataFileID);
        strDay = ConvertDate_HRF2020(alertDataFileDate);
        scoringLabels = [];
        for ff = 1:length(ScoringResults.fileIDs)
            if strcmp(alertDataFileID,ScoringResults.fileIDs{ff,1}) == true
                scoringLabels = ScoringResults.labels{ff,1};
            end
        end
        % check labels to match arousal state
        if sum(strcmp(scoringLabels,'Not Sleep')) < 36   % 36 bins (180 total) or 3 minutes of awake
            load(procDataFileID)
            % puffs = ProcData.data.solenoids.LPadSol;
            % % don't include trials with stimulation
            % if isempty(puffs) == true
            if strcmp(hemisphere(1:4),'cort') == true
                NeuralData{gg,1} = (ProcData.data.(['cortical_' hemisphere(end - 1:end)]).(neuralBand) - RestingBaselines.manualSelection.(['cortical_' hemisphere(end - 1:end)]).(neuralBand).(strDay))./RestingBaselines.manualSelection.(['cortical_' hemisphere(end - 1:end)]).(neuralBand).(strDay); %#ok<*AGROW>
            elseif strcmp(hemisphere(1:4),'hipp') == true
                NeuralData{gg,1} = (ProcData.data.hippocampus.(neuralBand) - RestingBaselines.manualSelection.hippocampus.(neuralBand).(strDay))./RestingBaselines.manualSelection.hippocampus.(neuralBand).(strDay);
            end
            HemoData{gg,1} = ProcData.data.CBV_HbT.(['adj' hemisphere(end - 1:end)]);
            gg = gg + 1;
            % end
        end
    end
end

end
