function [] = ProcessIntrinsicData_IOS(procDataFileIDs)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
load(procDataFileIDs(1,:));
[animalID,~,~] = GetFileInfo_IOS(procDataFileIDs(1,:));
imagingType = ProcData.notes.imagingType;
imagingWavelengths = ProcData.notes.imagingWavelengths;
% use imaging type to determine ROI names and typical lens magnification
if strcmpi(imagingType,'Single ROI (SI)') == true
    ROInames = {'barrels'};
elseif strcmpi(imagingType,'Single ROI (SSS)') == true
    ROInames = {'SSS','lSSS','rSSS'};
elseif strcmpi(imagingType,'Bilateral ROI (SI)') == true
    ROInames = {'LH','RH'};
elseif strcmpi(imagingType,'Bilateral ROI (SI,FC)') == true
    ROInames = {'LH','RH','fLH','fRH'};
end
lensMag = ProcData.notes.lensMag;
% create/load pre-existing ROI file with the coordinates
ROIFileDir = dir('*_ROIs.mat');
if isempty(ROIFileDir) == true
    ROIs = [];
else
    ROIFileName = {ROIFileDir.name}';
    ROIFileID = char(ROIFileName);
    load(ROIFileID);
end
% check whether or not each ROI already exists
[ROIs] = CheckROIDates_IOS(animalID,ROIs,ROInames,lensMag,imagingType,imagingWavelengths);
% extract CBV data from each ROI for each RawData file in the directory that hasn't been processed yet.
if any(strcmp(imagingWavelengths,{'Red, Green, & Blue','Red, Lime, & Blue'})) == true
    ExtractTriWavelengthData_IOS(ROIs,ROInames,procDataFileIDs)
elseif any(strcmp(imagingWavelengths,{'Green & Blue','Lime & Blue'})) == true
    ExtractDualWavelengthData_IOS(ROIs,ROInames,procDataFileIDs)
else
    ExtractSingleWavelengthData_IOS(ROIs,ROInames,procDataFileIDs)
end