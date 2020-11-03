function [MscanData] = DiamCalcSurfaceVessel_2P(tempData,imageID)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab
%________________________________________________________________________________________________________________________
%
%   Purpose: Draw an ROI for the surface vessel and fill in notes information.
%________________________________________________________________________________________________________________________

%% Pull information from graphics file
movieInfo = imfinfo([imageID '.TIF']);
% take file info and extrac magnification and frame rate
tempData.notes.header.fileName = movieInfo(1).Filename;
tempData.notes.header.frameWidth = num2str(movieInfo(1).Width);
tempData.notes.header.frameHeight = num2str(movieInfo(1).Height);
tempData.notes.header.numberOfFrames = length(movieInfo);
tempData.notes.xSize = str2double(tempData.notes.header.frameWidth);
tempData.notes.ySize = str2double(tempData.notes.header.frameHeight);
% Read header and take further action based on header information
textHold = strread(movieInfo(1).ImageDescription, '%s', 'delimiter', '\n'); %#ok<*DSTRRD>
magStart = strfind(textHold{20},': ');
tempData.notes.header.magnification = textHold{20}(magStart + 2:end);
rotationStart = strfind(textHold{19},': ');
tempData.notes.header.rotation = textHold{19}(rotationStart + 2:end);
frameRateStart = strfind(textHold{24},': ');
tempData.notes.header.frameRate = (textHold{24}(frameRateStart + 2:end - 3));
tempData.notes.frameRate = 1/str2num(tempData.notes.header.frameRate); %#ok<*ST2NM>
tempData.notes.startframe = 1;
tempData.notes.endframe = tempData.notes.header.numberOfFrames;
if tempData.notes.objectiveID == 1       %10X
    micronsPerPixel = 1.2953;
elseif tempData.notes.objectiveID == 2   % Small 20X
    micronsPerPixel = 0.5595; 
elseif tempData.notes.objectiveID == 3   % Big 20X
    micronsPerPixel = 0.64;  
elseif tempData.notes.objectiveID == 4   % 40X
    micronsPerPixel = 0.3619; 
elseif tempData.notes.objectiveID == 5   % 16X
    micronsPerPixel = 0.825;
end
tempData.notes.micronsPerPixel = micronsPerPixel;
tempData.notes.header.timePerLine = 1/(tempData.notes.frameRate*str2num(tempData.notes.header.frameHeight));
xFactor = micronsPerPixel/(str2num(tempData.notes.header.magnification(1:end - 1)));
% generate fig
image = imread(imageID,'TIFF','Index',1);
vesROI = figure;
imagesc(double(image))
title([tempData.notes.animalID '_' tempData.notes.date '_' tempData.notes.imageID])
colormap('gray');
axis image
xlabel('pixels')
ylabel('pixels')
% check input
yString = 'y';
theInput = 'n';
xSize = size(image,2);
ySize = size(image,1);
area = impoly(gca,[1,1;1,20;20,20;20,1]); %#ok<*IMPOLY>    
while strcmp(yString,theInput) ~= 1
    theInput = input('Is the diameter of the box ok? (y/n): ', 's'); disp(' ')
end
if strcmp(yString,theInput)
    get_API = iptgetapi(area);
    tempData.notes.vesselROI.boxPosition.xy = get_API.getPosition();
    tempData.notes.vesselROI.xSize = xSize;
    tempData.notes.vesselROI.ySize = ySize;
    theInput = 'n';
end
diamAxis = imline(gca,round(xSize*[.25,.75]),round(ySize*[.25,.75])); %#ok<*IMLINE>
while strcmp(yString,theInput) ~= 1
    theInput = input('Is the line along the diameter axis ok? (y/n): ','s'); disp(' ')
end
if strcmp(yString,theInput)
    get_API = iptgetapi(diamAxis);
    tempData.notes.vesselROI.vesselLine.position.xy = get_API.getPosition();
end
tempData.notes.xFactor = xFactor;
close(vesROI);
MscanData = tempData;

end
