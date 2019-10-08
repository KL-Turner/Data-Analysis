function [] = CreateIntrinsicTiffStacks_IOS(rawDataFileIDs)
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
%   Last Revised: 
%________________________________________________________________________________________________________________________

for a = 1:length(rawDataFileIDs,1)
    disp(['Processing WindowCam file (' num2str(a) '/' num2str(size(rawDataFileIDs,1)) ')']); disp(' ')
    rawDataFileID = rawDataFileIDs(a,:);
    load(rawDataFileID)
    [~,~,fileID] = GetFileInfo_IOS(rawDataFileID);
    % extract camera frames from WindowCam binary file
    [frames] = ReadDalsaBinary_IOS([fileID '_WindowCam.bin'],RawData.notes.CBVCamPixelHeight,RawData.notes.CBVCamPixelWidth);
    % write each frame to a Tiff stack with an identical name
    tiffStackFileID = [fileID '_WindowCam.tif'];
    for b = 1:length(frames)
        frame = frames{1,b};
        disp(['Writing frame (' num2str(b) '/' num2str(length(frames)) ') to TIFF stack']); disp(' ')
        uintImage = imcomplement(im2uint16(frame));
        imwrite(uintImage,tiffStackFileID,'WriteMode','append')
    end
end

end

