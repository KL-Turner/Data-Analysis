function [] = ExtractCBVData_IOS(ROIs, ROInames, rawDataFileIDs)
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
    [animalID,fileDate,fileID] = GetFileInfo_IOS(rawDataFile);
    strDay = ConvertDate_IOS(fileDate);
    load(rawDataFile)
    for b = 1:length(ROInames)
        ROIshortName = ROInames{1,b}; 
        ROIname = [ROInames{1,b} '_' strDay];
        if ~isfield(RawData.data,'CBV')
            RawData.data.CBV = [];
        end
%         if ~isfield(RawData.data.CBV,ROIname)
            [frames] = ReadDalsaBinary_IOS(animalID,[fileID '_WindowCam.bin']);
            disp(['Extracting ' ROIname ' ROI CBV data from ' rawDataFile '...']); disp(' ')
            if strcmp(ROIshortName,'LH') == true || strcmp(ROIshortName,'RH') == true
                circROI = drawcircle('Center',ROIs.(ROIname).circPosition,'Radius',ROIs.(ROIname).circRadius);
                mask = createMask(circROI,frames{1});
            else
                mask = roipoly(frames{1},ROIs.(ROIname).xi,ROIs.(ROIname).yi);
            end
            meanIntensity = BinToIntensity_IOS([fileID '_WindowCam.bin'],mask,frames);
            RawData.data.CBV.(ROIname) = meanIntensity;
            save(rawDataFile,'RawData')
%         else
%             disp([ROIname ' ROI CBV data from ' rawDataFile ' already extracted. Continuing...']); disp(' ')
%         end
    end
end

end

