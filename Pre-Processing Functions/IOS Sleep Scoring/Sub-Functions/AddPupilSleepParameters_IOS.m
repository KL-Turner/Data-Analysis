function [] = AddPupilSleepParameters_IOS(procDataFileIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Organize data into appropriate bins for sleep scoring characterization
%________________________________________________________________________________________________________________________

% restingBaselineDataFileStruct = dir('*_RestingBaselines.mat');
% restingBaselineDataFiles = {restingBaselineDataFileStruct.name}';
% restingBaselineDataFileID = char(restingBaselineDataFiles);
% load(restingBaselineDataFileID)
% baselineType = 'manualSelection';
for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    disp(['Adding sleep scoring parameters to ' procDataFileID '... (' num2str(a) '/' num2str(size(procDataFileIDs,1)) ')']); disp(' ')
%     [~,fileDate,~] = GetFileInfo_IOS(procDataFileID);
%     strDay = ConvertDate_IOS(fileDate);
    load(procDataFileID)
    if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
        %% BLOCK PURPOSE: Create folder for the Neural data of each electrode
        % cortical delta
        pupilArea = ProcData.data.Pupil.pupilArea;
        % pupilBaseline = RestingBaselines.(baselineType).Pupil.pupilArea.(strDay);
        % pupilArea = (pupilArea - pupilBaseline)/pupilBaseline;
        % Divide the neural signals into five second bins and put them in a cell array
        pupilStruct = cell(180,1);
        % loop through all samples across the 15 minutes in 5 second bins (180 total)
        for b = 1:180
            if b == 1
                pupilStruct(b,1) = {pupilArea(b:150)};
            elseif b == 180
                pupilStruct(b,1) = {pupilArea((((150*(b - 1)) + 1)):end)};
            else
                pupilStruct(b,1) = {pupilArea((((150*(b - 1)) + 1)):(150*b))};
            end
        end
        % save data under ProcData file
        ProcData.sleep.parameters.Pupil.pupilArea = pupilStruct;
        save(procDataFileID,'ProcData');
    end
end

end
