function [] = ProcessIntrinsicData_IOS(rawDataFileIDs,procDataFileIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Base function to run the functions necessary for IOS data extraction from drawn ROIs over the images
%________________________________________________________________________________________________________________________

load(rawDataFileIDs(1,:));
[animalID,~,~] = GetFileInfo_IOS(rawDataFileIDs(1,:));
imagingType = RawData.notes.imagingType;
imagingWavelengths = RawData.notes.imagingWavelengths;
% use imaging type to determine ROI names and typical lens magnification
if strcmpi(imagingType,'Single ROI (SI)') == true
    lensMag = '2.0';
    ROInames = {'barrels'};
elseif strcmpi(imagingType,'Single ROI (SSS)') == true
    lensMag = '3.0X';
    ROInames = {'SSS','lSSS','rSSS'};
elseif strcmpi(imagingType,'Bilateral ROI (SI)') == true
    lensMag = '1.5X';
    ROInames = {'LH','RH'};
elseif strcmpi(imagingType,'Bilateral ROI (SI,FC)') == true
    lensMag = '1.5X';
    ROInames = {'LH','RH','fLH','fRH'};
end
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
    ExtractTriWavelengthData_IOS(ROIs,ROInames,rawDataFileIDs,procDataFileIDs)
elseif any(strcmp(imagingWavelengths,{'Green & Blue','Lime & Blue'})) == true
    ExtractDualWavelengthData_IOS(ROIs,ROInames,rawDataFileIDs,procDataFileIDs)
else
    ExtractSingleWavelengthData_IOS(ROIs,ROInames,rawDataFileIDs)
end
% go through each ProcData file and add the pixel data to each
for aa = 1:size(procDataFileIDs,1)
    disp(['Adding IOS CBV data to ProcData file (' num2str(aa) '/' num2str(size(procDataFileIDs,1)) ')']); disp(' ')
    procDataFileID = procDataFileIDs(aa,:);
    load(procDataFileID)
    rawDataFileID = rawDataFileIDs(aa,:);
    load(rawDataFileID)
    [~,fileDate,~] = GetFileInfo_IOS(rawDataFileID);
    strDay = ConvertDate_IOS(fileDate);
    for bb = 1:length(ROInames)
        if any(strcmp(imagingWavelengths,{'Red, Green, & Blue','Red, Lime, & Blue'})) == true
            ProcData.data.CBV.(ROInames{1,bb}) = RawData.data.CBV.([ROInames{1,bb} '_' strDay]);
            ProcData.data.GCaMP7s.(ROInames{1,bb}) = RawData.data.GCaMP7s.([ROInames{1,bb} '_' strDay]);
            ProcData.data.Deoxy.(ROInames{1,bb}) = RawData.data.Deoxy.([ROInames{1,bb} '_' strDay]);
            ProcData.notes.CBVCamSamplingRate = RawData.notes.CBVCamSamplingRate/3;
        elseif any(strcmp(imagingWavelengths,{'Green & Blue','Lime & Blue'})) == true
            ProcData.data.CBV.(ROInames{1,bb}) = RawData.data.CBV.([ROInames{1,bb} '_' strDay]);
            ProcData.data.GCaMP7s.(ROInames{1,bb}) = RawData.data.GCaMP7s.([ROInames{1,bb} '_' strDay]);
            ProcData.notes.CBVCamSamplingRate = RawData.notes.CBVCamSamplingRate/2;
        else
            ProcData.data.CBV.(ROInames{1,bb}) = RawData.data.CBV.([ROInames{1,bb} '_' strDay])(1:end - 1);
        end
        save(procDataFileID,'ProcData')
    end
end

end

