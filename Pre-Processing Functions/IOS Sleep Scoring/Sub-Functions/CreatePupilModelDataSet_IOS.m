function [] = CreatePupilModelDataSet_IOS(procDataFileIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Arrange data into a table of most-relevant parameters for model training/classification
%________________________________________________________________________________________________________________________

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    pupilModelDataSetID = [procDataFileID(1:end-12) 'PupilModelData.mat'];
    load(procDataFileID)
    if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
        %% Create table to send into model
        variableNames = {'pupilArea','numWhiskEvents'};
        % pre-allocation
        numWhiskEvents_column = zeros(180,1);
        avgPupilArea_column = zeros(180,1);
        % extract relevant parameters from each epoch
        for b = 1:length(numWhiskEvents_column)
            % number of binarized whisking events
            numWhiskEvents_column(b,1) = sum(ProcData.sleep.parameters.binWhiskerAngle{b,1});
            % average pupil area
            avgPupilArea_column(b,1) = round(nanmean(ProcData.sleep.parameters.Pupil.pupilArea{b,1}),1);
        end
        % create table
        pupilParamsTable = table(avgPupilArea_column,numWhiskEvents_column,'VariableNames',variableNames);
        save(pupilModelDataSetID,'pupilParamsTable')
    end
end

end
