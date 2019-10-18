function [] = ExtractTiffData_IOS(rawDataFileIDs)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab and Yurong Gao
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________

% loop through each raw data file and save the analyzed vessels diameter to it
for a = 1:size(rawDataFileIDs)
    rawDataFileID = rawDataFileIDs(a,:);
    load(rawDataFileID)
    [~,~,fileID] = GetFileInfo_IOS(rawDataFileID);
    tiffStackFileID = [fileID '_WindowCam.tif'];
    structFields = fieldnames(RawData.notes);
    % identify how many individual vessels to analyze for this file
    d = 1;
    b = 0;
    vesselIDs = [];
    while b == 0
        for c = 1:length(structFields)
            structFieldName = structFields{c,1};
            vesselFieldName = ['V' num2str(d)];
            checkD = d;
            if strcmp(structFieldName,vesselFieldName) == true
                vesselIDs{d,1} = vesselFieldName; %#ok<*AGROW>
                d = d + 1;
                break
            end
        end
        if d == checkD
            b = 1;
        end
    end
    % loop through each vessel and perform diameter analysis
    for e = 1:length(vesselIDs)
        vID = vesselIDs{e,1};
        %% Takes a tiff file (movie) and an analog ascii file, and extracts the diameters
        disp('Analyzing vessel projections from defined polygons...'); disp(' ');
        [RawData] = GetDiameterFromMovie_IOS(RawData,vID,tiffStackFileID);
        try
            [RawData] = FWHM_MovieProjection_IOS(RawData,vID,[RawData.notes.startframe RawData.notes.endframe]);
        catch
            disp('FWHM calculation failed.'); disp(' ')
            keyboard
        end
        try
            % 1 dural/vein, >40% changes spline, artery: >60% spline
            % 2 dural/vein, >30% changes interpolate, artery: >50% interpolate
            if strcmp(RawData.notes.(vID).vesselType,'D') == true || strcmp(RawData.notes.(vID).vesselType,'V') == true
                RawData.data.(vID).vesselDiameter = RemoveMotion_IOS(RawData.data.(vID).tempVesselDiameter,RawData.notes.(vID).vesselROI.modalFixedDiameter,2,0.3);
            else
                RawData.data.(vID).vesselDiameter = RemoveMotion_IOS(RawData.data.(vID).tempVesselDiameter,RawData.notes.(vID).vesselROI.modalFixedDiameter,2,0.5);
            end
            [diamPerc,powerS,powerf] = DiamPercPower_IOS(RawData.data.(vID).vesselDiameter,RawData.notes.(vID).vesselROI.modalFixedDiameter,RawData.notes.CBVCamSamplingRate);
            RawData.notes.(vID).vessel.diamPerc = diamPerc;
            RawData.notes.(vID).vessel.powerf = powerf;
            RawData.notes.(vID).vessel.powerS = powerS;
        catch
            disp('Diameter percentage analysis failed.'); disp(' ')
            keyboard
        end
    end
    save(rawDataFileID,'RawData')
end

end

%% Opens the tiff file and gets the  vessel projections from the defined polygons
function [RawData] = GetDiameterFromMovie_IOS(RawData,vID,tiffStackFileID)
RawData.notes.(vID).firstFrame = imread(tiffStackFileID,'TIFF','Index',1);
fftFirstFrame = fft2(double(RawData.notes.(vID).firstFrame));
X = repmat(1:RawData.notes.xSize, RawData.notes.ySize, 1);
Y = repmat((1:RawData.notes.ySize)',1,RawData.notes.xSize);
RawData.notes.(vID).vesselROI.projectionAngle = atand(diff(RawData.notes.(vID).vesselROI.vesselLine.position.xy(:,1))/diff(RawData.notes.(vID).vesselROI.vesselLine.position.xy(:,2)));
atand(diff(RawData.notes.(vID).vesselROI.vesselLine.position.xy(:,1))/diff(RawData.notes.(vID).vesselROI.vesselLine.position.xy(:,2)));
% Find the projection for each frame based on the Radon transform
for f = RawData.notes.startframe:RawData.notes.endframe
    rawFrame = imread(tiffStackFileID,'TIFF','Index',f);
    fftRawFrame = fft2(double(rawFrame));
    [RawData.notes.(vID).pixelShift(:,f),~] = DftRegistration_IOS(fftFirstFrame,fftRawFrame,1);
    inpolyFrame = inpolygon(X + RawData.notes.(vID).pixelShift(3,f),Y + RawData.notes.(vID).pixelShift(4,f),RawData.notes.(vID).vesselROI.boxPosition.xy(:,1),RawData.notes.(vID).vesselROI.boxPosition.xy(:,2));
    boundedrawFrame = rawFrame.*uint16(inpolyFrame);
    RawData.notes.(vID).vesselROI.projection(f,:) = radon(boundedrawFrame,RawData.notes.(vID).vesselROI.projectionAngle);
end

end

%% Calculate diameter using FWHM and get the baseline diameter
function [RawData] = FWHM_MovieProjection_IOS(RawData,vID,theFrames)
for g = min(theFrames):max(theFrames)
    % Add in a 5 pixel median filter
    RawData.data.(vID).rawVesselDiameter(g) = CalcFWHM_IOS(medfilt1(RawData.notes.(vID).vesselROI.projection(g,:),5));
end
RawData.data.(vID).tempVesselDiameter = RawData.data.(vID).rawVesselDiameter*RawData.notes.xFactor;
[holdHist,d] = hist(RawData.data.(vID).tempVesselDiameter,0:0.25:100);
[~,maxD] = max(holdHist);
RawData.notes.(vID).vesselROI.modalFixedDiameter = d(maxD);

end

