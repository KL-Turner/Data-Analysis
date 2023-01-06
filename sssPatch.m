zap;

% select imaging type
imagingOptions = {'Single ROI (SI)','Single ROI (SSS)','Bilateral ROI (SI)','Bilateral ROI (SI,FC)'};
imagingType = SelectImagingType_IOS(imagingOptions);
% select imaging type
wavelengthOptions = {'Green','Lime','Blue','Green & Blue','Lime & Blue','Red, Green, & Blue','Red, Lime, & Blue'};
imagingWavelengths = SelectImagingType_IOS(wavelengthOptions);

rawDataFileStruct = dir('*_RawData.mat');
rawDataFiles = {rawDataFileStruct.name}';
rawDataFileIDs = char(rawDataFiles);
for aa = 1:size(rawDataFileIDs)
    rawDataFileID = rawDataFileIDs(aa,:);
    load(rawDataFileID)
    RawData.notes.imagingType = imagingType;
    RawData.notes.imagingWavelengths = imagingWavelengths;
    save(rawDataFileID,'RawData')
end

procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
for aa = 1:size(procDataFileIDs)
    procDataFileID = procDataFileIDs(aa,:);
    load(procDataFileID)
    ProcData.notes.imagingType = imagingType;
    ProcData.notes.imagingWavelengths = imagingWavelengths;
    save(procDataFileID,'ProcData')
end