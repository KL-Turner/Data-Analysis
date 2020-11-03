function [MScanData] = ExtractTiffAnalogData_2P_eLife2020(MScanData, fileID)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab and Yurong Gao
%________________________________________________________________________________________________________________________
%
%   Purpose: Uses the TIFF image stack and txt files from MScan to extract the vessel diameter and the analog data.
%________________________________________________________________________________________________________________________

%% Takes a tiff file (movie) and an analog ascii file, and extracts the diameters
MScan_analogData = [fileID '.TXT'];
disp(['Loading MScan file: ' MScan_analogData '...']); disp(' ');
analogData = load(MScan_analogData,'-ascii');
MScanData.data.corticalNeural = analogData(:,2);
MScanData.data.EMG = analogData(:,3);
MScanData.data.hippocampalNeural = analogData(:,4);
MScanData.data.forceSensor = analogData(:,5);
MScanData.notes.analogSamplingRate = 20000;
% 
% figure; plot(MScanData.data.corticalNeural)
% figure; plot(MScanData.data.EMG)
% figure; plot(MScanData.data.hippocampalNeural)
% figure; plot(MScanData.data.forceSensor)

disp('Analyzing vessel projections from defined polygons...'); disp(' ');
[MScanData] = GetDiameterFromMovie_2P_eLife2020(MScanData,fileID);

try
    [MScanData] = FWHM_MovieProjection_2P_eLife2020(MScanData,[MScanData.notes.startframe,MScanData.notes.endframe]);
catch
    disp([MScanData.notes.imageID ' FWHM calculation failed!']); disp(' ')
end

try
    % 1 dural/vein, >40% changes spline, artery: >60% spline
    % 2 dural/vein, >30% changes interpolate, artery: >50% interpolate
    if strcmp(MScanData.notes.vesselType,'D') || strcmp(MScanData.notes.vesselType,'V')
        MScanData.data.vesselDiameter = RemoveMotion_2P_eLife2020(MScanData.data.tempVesselDiameter,MScanData.notes.vesselROI.modalFixedDiameter,2,0.3);
    else
        MScanData.data.vesselDiameter = RemoveMotion_2P_eLife2020(MScanData.data.tempVesselDiameter,MScanData.notes.vesselROI.modalFixedDiameter,2,0.5);
    end
    [diamPerc,S,f] = DiamPercPower_2P_eLife2020(MScanData.data.vesselDiameter,MScanData.notes.vesselROI.modalFixedDiameter,MScanData.notes.frameRate);
    MScanData.notes.vessel.diamPerc = diamPerc;
    MScanData.notes.vessel.power_f = f;
    MScanData.notes.vessel.power_S = S;
catch
    disp([MScanData.notes.imageID ' Diameter percentage analysis failed!']); disp(' ')
    keyboard
end

end

%% Opens the tiff file and gets the  vessel projections from the defined polygons
function [MScanData] = GetDiameterFromMovie_2P_eLife2020(MScanData,fileID)
MScanData.notes.firstFrame = imread(fileID,'TIFF','Index',1);
fftFirstFrame = fft2(double(MScanData.notes.firstFrame));
X = repmat(1:MScanData.notes.xSize,MScanData.notes.ySize,1);
Y = repmat((1:MScanData.notes.ySize)',1,MScanData.notes.xSize);
MScanData.notes.vesselROI.projectionAngle = atand(diff(MScanData.notes.vesselROI.vesselLine.position.xy(:,1))/diff(MScanData.notes.vesselROI.vesselLine.position.xy(:,2)));
for theFrame = MScanData.notes.startframe:MScanData.notes.endframe
    rawFrame = imread(fileID,'TIFF','Index',theFrame);
    fftRawFrame = fft2(double(rawFrame));
    [MScanData.notes.pixelShift(:,theFrame),~] = DftRegistration_2P_eLife2020(fftFirstFrame,fftRawFrame,1);
    inpolyFrame = inpolygon(X + MScanData.notes.pixelShift(3,theFrame),Y + MScanData.notes.pixelShift(4,theFrame),MScanData.notes.vesselROI.boxPosition.xy(:,1),MScanData.notes.vesselROI.boxPosition.xy(:,2));
    boundedrawFrame = rawFrame.*uint16(inpolyFrame);
    MScanData.notes.vesselROI.projection(theFrame,:) = radon(boundedrawFrame,MScanData.notes.vesselROI.projectionAngle);
end

end

%% Calculate diameter using FWHM and get the baseline diameter
function [MScanData] = FWHM_MovieProjection_2P_eLife2020(MScanData, theFrames)
for f = min(theFrames):max(theFrames)
    % Add in a 5 pixel median filter
    MScanData.data.rawVesselDiameter(f) = CalcFWHM_2P_eLife2020(medfilt1(MScanData.notes.vesselROI.projection(f,:),5));
end
MScanData.data.tempVesselDiameter = MScanData.data.rawVesselDiameter*MScanData.notes.xFactor;
[holdHist,d] = hist(MScanData.data.tempVesselDiameter,0:.25:100); %#ok<HIST>
[~,maxD] = max(holdHist);
MScanData.notes.vesselROI.modalFixedDiameter = d(maxD);

end
