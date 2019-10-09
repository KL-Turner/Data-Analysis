function ExtractCBVData_IOS(ROIs, ROInames, rawDataFileIDs)
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
%   Last Revised: February 29th, 2019
%________________________________________________________________________________________________________________________

for a = 1:size(rawDataFileIDs,1)
    rawDataFile = rawDataFileIDs(a,:);
    disp(['Analyzing RawData file ' num2str(a) ' of ' num2str(size(rawDataFileIDs, 1)) '...']); disp(' ')
    [~, fileDate, fileID] = GetFileInfo_IOS(rawDataFile);
    strDay = ConvertDate_IOS(fileDate);
    load(rawDataFile)
    RawData.data.CBV = [];
%     if ~isfield(RawData.data,'CBV')
        [frames] = ReadDalsaBinary_IOS([fileID '_WindowCam.bin'],RawData.notes.CBVCamPixelHeight,RawData.notes.CBVCamPixelWidth);
        for b = 1:length(ROInames)
            ROIname = [ROInames{1,b} '_' strDay];
            disp(['Extracting ' ROIname ' ROI CBV data from ' rawDataFile '...']); disp(' ')
            xi = ROIs.(ROIname).xi;
            yi = ROIs.(ROIname).yi;
            mask = roipoly(frames{1},xi,yi);
            meanIntensity = BinToIntensity_IOS([fileID '_WindowCam.bin'],mask,frames);
            RawData.data.CBV.(ROIname) = meanIntensity;
        end
        save(rawDataFile, 'RawData')
%     else
%         disp([rawDataFile ' CBV already extracted. Continuing...']); disp(' ')
%     end
end

end

