function [] = DrawVesselROIs_IOS(procDataFileIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________

ROIFileDir = dir('*_ROIs.mat');
ROIFileName = {ROIFileDir.name}';
ROIFileID = char(ROIFileName);
load(ROIFileID);

procDataFileID = procDataFileIDs(1,:);
[~,~,fileID] = GetFileInfo_IOS(procDataFileID);
windowCamFileID = [fileID '_WindowCam.bin'];
load(procDataFileID)

imageHeight = ProcData.notes.CBVCamPixelHeight;
imageWidth = ProcData.notes.CBVCamPixelWidth;
samplingRate = ProcData.notes.CBVCamSamplingRate;
trialDuration_sec = ProcData.notes.trialDuration_sec;

frameInds = 1:samplingRate*trialDuration_sec;
[imageStack] =  GetCBVFrameSubset_IOS(windowCamFileID,imageWidth,imageHeight,frameInds);

drawVesselLines = 'y';
b = 1;
while strcmp(drawVesselLines,'y') == true
    fig = figure;
    imagesc(imageStack(:,:,10));
    axis image;
    colormap gray;
    hold on;
    [x1,y1] = ginput(1);
    h1 = plot(gca,x1,y1,'bo',x1,y1,'r+','markersize',10);
    [x2,y2] = ginput(1);
    h2 = plot(gca,x2,y2,'bo',x1,y1,'r+','markersize',10);
    delete([h1,h2]);
    imline(gca,[x1,y1;x2,y2]);
    hold off;
    ROIs.vesselIDs{b,1} = input('Input the vessel name and type (A#,V#): ','s'); disp(' ')
    ROIs.x_endpoints{b,1} = [x1,x2];
    ROIs.y_endpoints{b,1} = [y1,y2];
    b = b + 1;
    drawVesselLines = input('Continue on to another vessel? (y/n): ','s'); disp(' ')
    if strcmp(drawVesselLines,'y') == true || strcmp(drawVesselLines,'n') == true
        continue
    else
        drawVesselLines = 'n';
    end
end
save(ROIFileID,'ROIs')

end
