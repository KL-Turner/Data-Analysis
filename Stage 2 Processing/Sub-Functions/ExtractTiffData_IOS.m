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
%
%   Inputs:
%
%   Outputs:
%
%   Last Revised:
%________________________________________________________________________________________________________________________

for a = 1:size(rawDataFileIDs)
    rawDataFileID = rawDataFileIDs(a,:);
    load(rawDataFileID)
    [~,~,fileID] = GetFileInfo_IOS(rawDataFileID);
    tiffStackFileID = [fileID '_WindowCam.tif'];
    
    %% Takes a tiff file (movie) and an analog ascii file, and extracts the diameters
    MScan_analogData = [fileID '.TXT'];
    disp(['Loading MScan file: ' MScan_analogData '...']); disp(' ');
    analogData = load(MScan_analogData, '-ascii');
    RawData.data.EMG = analogData(:, 2);
    RawData.data.forceSensor = analogData(:, 3);
    RawData.data.rawNeuralData = analogData(:, 4);
    RawData.notes.analogSamplingRate = 20000;
    
    disp('Analyzing vessel projections from defined polygons...'); disp(' ');
    [RawData] = GetDiameterFromMovie_2P(RawData, fileID);
    
    try
        [RawData] = FWHM_MovieProjection_2P(RawData, [RawData.notes.startframe RawData.notes.endframe]);
    catch
        disp([RawData.notes.imageID ' FWHM calculation failed!']); disp(' ')
    end
    
    try
        % 1 dural/vein, >40% changes spline, artery: >60% spline
        % 2 dural/vein, >30% changes interpolate, artery: >50% interpolate
        if strcmp(RawData.notes.vesselType, 'D') || strcmp(RawData.notes.vesselType, 'V')
            RawData.data.vesselDiameter = RemoveMotion_2P(RawData.data.tempVesselDiameter, RawData.notes.vesselROI.modalFixedDiameter, 2, 0.3);
        else
            RawData.data.vesselDiameter = RemoveMotion_2P(RawData.data.tempVesselDiameter, RawData.notes.vesselROI.modalFixedDiameter, 2, 0.5);
        end
        [diamPerc, S, f] = DiamPercPower_2P(RawData.data.vesselDiameter, RawData.notes.vesselROI.modalFixedDiameter, RawData.notes.frameRate);
        RawData.notes.vessel.diamPerc = diamPerc;
        RawData.notes.vessel.power_f = f;
        RawData.notes.vessel.power_S = S;
    catch
        disp([RawData.notes.imageID ' Diameter percentage analysis failed!']); disp(' ')
    end
    
end

%% Opens the tiff file and gets the  vessel projections from the defined polygons
    function [MScanData] = GetDiameterFromMovie_2P(MScanData, fileID)
        MScanData.notes.firstFrame = imread(fileID, 'TIFF', 'Index', 1);
        fftFirstFrame = fft2(double(MScanData.notes.firstFrame));
        X = repmat(1:MScanData.notes.xSize, MScanData.notes.ySize, 1);
        Y = repmat((1:MScanData.notes.ySize)', 1, MScanData.notes.xSize);
        MScanData.notes.vesselROI.projectionAngle = atand(diff(MScanData.notes.vesselROI.vesselLine.position.xy(:, 1))/diff(MScanData.notes.vesselROI.vesselLine.position.xy(:, 2)));
        atand(diff(MScanData.notes.vesselROI.vesselLine.position.xy(:, 1))/diff(MScanData.notes.vesselROI.vesselLine.position.xy(:, 2)));
        
        for theFrame = MScanData.notes.startframe:MScanData.notes.endframe
            rawFrame = imread(fileID, 'TIFF', 'Index', theFrame);
            fftRawFrame = fft2(double(rawFrame));
            
            [MScanData.notes.pixelShift(:, theFrame), ~] = DftRegistration_2P(fftFirstFrame, fftRawFrame, 1);
            
            inpolyFrame = inpolygon(X + MScanData.notes.pixelShift(3, theFrame), Y + MScanData.notes.pixelShift(4, theFrame), MScanData.notes.vesselROI.boxPosition.xy(:, 1), MScanData.notes.vesselROI.boxPosition.xy(:, 2));
            boundedrawFrame = rawFrame.*uint16(inpolyFrame);
            MScanData.notes.vesselROI.projection(theFrame, :) = radon(boundedrawFrame, MScanData.notes.vesselROI.projectionAngle);
        end
        
    end

%% Calculate diameter using FWHM and get the baseline diameter
    function [MScanData] = FWHM_MovieProjection_2P(MScanData, theFrames)
        
        for f = min(theFrames):max(theFrames)
            % Add in a 5 pixel median filter
            MScanData.data.rawVesselDiameter(f) = CalcFWHM_2P(medfilt1(MScanData.notes.vesselROI.projection(f, :), 5));
        end
        
        MScanData.data.tempVesselDiameter = MScanData.data.rawVesselDiameter*MScanData.notes.xFactor;
        [holdHist, d] = hist(MScanData.data.tempVesselDiameter, 0:.25:100);
        [~, maxD] = max(holdHist);
        MScanData.notes.vesselROI.modalFixedDiameter = d(maxD);
        
    end

end

