zap
% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
imagingOptions = {'Single ROI (SI)','Single ROI (SSS)','Bilateral ROI (SI)','Bilateral ROI (SI,FC)'};
imagingType = SelectImagingType_IOS(imagingOptions);
% select imaging type
wavelengthOptions = {'Green','Lime','Blue','Green & Blue','Lime & Blue','Red, Green, & Blue','Red, Lime, & Blue'};
imagingWavelengths = SelectImagingType_IOS(wavelengthOptions);
% edit fieldname
for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    load(procDataFileID);
    disp(num2str(aa))
    ProcData.notes.imagingType = imagingType;
    ProcData.notes.imagingWavelengths = imagingWavelengths;
    save(procDataFileID,'ProcData')
end
