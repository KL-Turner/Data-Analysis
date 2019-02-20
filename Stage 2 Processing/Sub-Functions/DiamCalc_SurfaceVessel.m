function MscanData = DiamCalc_SurfaceVessel(tempData, imageID)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%
% Originally written by Patrick J. Drew and Yurong Gao
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%
%   Last Revised: February 19th, 2019    
%________________________________________________________________________________________________________________________

clear MscanData
movieInfo = imfinfo([imageID '.TIF']);   % Pull information from graphics file

%take file info and extrac magnification and frame rate
tempData.Notes.Header.Filename = movieInfo(1).Filename;
tempData.Notes.Header.frameWidth = num2str(movieInfo(1).Width);
tempData.Notes.Header.frameHeight = num2str(movieInfo(1).Height);
tempData.Notes.Header.numberOfFrames = length(movieInfo);
tempData.Notes.Header.frameCount = tempData.Notes.Header.numberOfFrames;
tempData.Notes.xSize = str2double(tempData.Notes.Header.frameWidth);
tempData.Notes.ySize = str2double(tempData.Notes.Header.frameHeight);

% Read header and take further action based on header information
textHold = strread(movieInfo(1).ImageDescription, '%s', 'delimiter', '\n'); %#ok<*DSTRRD>
magStart = strfind(textHold{20}, ': ');
tempData.Notes.Header.magnification = textHold{20}(magStart + 2:end);
rotationStart = strfind(textHold{19}, ': ');
tempData.Notes.Header.rotation = textHold{19}(rotationStart + 2:end);
frameRateStart = strfind(textHold{24}, ': ');
tempData.Notes.Header.frameRate=(textHold{24}(frameRateStart + 2:end - 3));
tempData.Notes.frameRate = 1/str2num(tempData.Notes.Header.frameRate); %#ok<*ST2NM>
tempData.Notes.startframe = 1;
tempData.Notes.endframe = tempData.Notes.Header.numberOfFrames;

if tempData.Notes.objectiveID == 1   %10X
    micronsPerPixel = 1.2953;
    
elseif tempData.Notes.objectiveID == 2   % Small 20X
    micronsPerPixel = 0.5595;
    
elseif tempData.Notes.objectiveID == 3   % Big 20X
    micronsPerPixel = 0.64;
    
elseif tempData.Notes.objectiveID == 4   %40X
    micronsPerPixel = 0.3619;
    
elseif tempData.Notes.objectiveID == 5   %16X
    micronsPerPixel = 0.825;
end

tempData.Notes.micronsPerPixel = micronsPerPixel;
tempData.Notes.Header.timePerLine = 1/(tempData.Notes.frameRate*str2num(tempData.Notes.Header.frameHeight));
xFactor = micronsPerPixel/(str2num(tempData.Notes.Header.magnification(1:end - 1)));

image = imread(imageID, 'TIFF', 'Index', 1);
vesselROI = figure;
imagesc(double(image))
title([tempData.Notes.animalID '_' tempData.Notes.date '_' tempData.Notes.imageID])
colormap('gray');
axis image
xlabel('pixels')
ylabel('pixels')

yString = 'y';
theInput = 'n';
xSize = size(image, 2);
ySize = size(image, 1);
area = impoly(gca, [1 1; 1 20; 20 20; 20 1]); %#ok<*IMPOLY>
        
while strcmp(yString, theInput) ~= 1
    theInput = input('Is the diameter of the box ok? (y/n): ', 's');
end
disp(' ')

if strcmp(yString, theInput)
    get_API = iptgetapi(area);
    tempData.Notes.vessel.boxPosition.xy = get_API.getPosition();
    tempData.Notes.vessel.xSize = xSize;
    tempData.Notes.vessel.ySize = ySize;
    theInput = 'n';
end

diamAxis = imline(gca, round(xSize*[.25 .75]), round(ySize*[.25 .75])); %#ok<*IMLINE>
while strcmp(yString, theInput) ~= 1
    theInput = input('Is the line along the diameter axis ok? (y/n): ', 's');
end
disp(' ')

if strcmp(yString, theInput)
    get_API = iptgetapi(diamAxis);
    tempData.Notes.vessel.vesselLine.position.xy = get_API.getPosition();
end

tempData.Notes.xFactor = xFactor;
tempData.Notes.vesselType = input('What is the type of this vessel? (A, V, PA, AV): ', 's'); disp(' ');
tempData.Notes.vessel.order = input('what is the order of this vessel? (number): ', 's'); disp(' ')

MscanData = tempData;

%% Save the file to directory.
[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Vessel ROIs/'];

if ~exist(dirpath, 'dir') 
    mkdir(dirpath); 
end

savefig(vesselROI, [dirpath tempData.Notes.animalID '_' tempData.Notes.date '_' tempData.Notes.imageID '_VesselROI']);

end
