function [] = GenerateTrialFigures_IOS(procDataFileIDs,RestingBaselines)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Create a summary figure for a single n minute IOS trial
%________________________________________________________________________________________________________________________

for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    % load file and gather information
    load(procDataFileID)
    [animalID,~,fileID] = GetFileInfo_IOS(procDataFileID);
    % imaging type
    imagingType = ProcData.notes.imagingType;
    if strcmpi(imagingType,'Single ROI (SI)') == true
        [figHandle] = GenerateTrialFigures_SI_IOS(procDataFileID,RestingBaselines);
    elseif strcmpi(imagingType,'Single ROI (SSS)') == true
        [figHandle] = GenerateTrialFigures_SSS_IOS(procDataFileID,RestingBaselines);
    elseif strcmpi(imagingType,'Bilateral ROI (SI)') == true
        [figHandle] = GenerateTrialFigures_bilatSI_IOS(procDataFileID,RestingBaselines);
    elseif strcmpi(imagingType,'Bilateral ROI (SI,FC)') == true
        [figHandle] = GenerateTrialFigures_bilatSIFC_IOS(procDataFileID,RestingBaselines);
    end
    % save the file to directory.
    [pathstr,~,~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Single Trial Figures/'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(figHandle,[dirpath animalID '_' fileID '_SingleTrialFig']);
    close(figHandle)
end
