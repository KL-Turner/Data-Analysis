function [OrgData] = Extract_Tiff_Analog_Data(OrgData, fileID)
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

%% Takes a tiff file (movie) and an analog ascii file, and extracts the diameters
MScan_analogData = [fileID '.TXT'];
disp(['Loading MScan file: ' MScan_analogData '...']); disp(' ');
analogData = load(MScan_analogData);
OrgData.Data.MScan_Force_Sensor = analogData(:, 2);
OrgData.Data.MScan_Neural_Data = analogData(:, 3);
OrgData.Notes.MScan_analogSamplingRate = 20000;

disp('Analyzing vessel projections from defined polygons...'); disp(' ');
[OrgData] = GetDiametersFromMovie(OrgData, fileID);

try
    [OrgData] = FWHM_MovieProjection(OrgData, [OrgData.Notes.startframe OrgData.Notes.endframe]);
catch
    disp([OrgData.Notes.imageID ' FWHM calculation failed!'])
end

try
    % 1 dural/vein, >40% changes spline, artery: >60% spline
    % 2 dural/vein, >30% changes interpolate, artery: >50% interpolate
    if strcmp(OrgData.Notes.vesselType, 'D') || strcmp(OrgData.Notes.vesselType, 'V')
        OrgData.Data.Vessel_Diameter = Remove_Motion(OrgData.Data.Vessel_TempDiameter, OrgData.Notes.vessel.modalFixedDiameter, 2, 0.3);
    else
        OrgData.Data.Vessel_Diameter = Remove_Motion(OrgData.Data.Vessel_TempDiameter, OrgData.Notes.vessel.modalFixedDiameter, 2, 0.5);
    end
    [diamPerc, S, f] = diamPerc_Power(OrgData.Data.Vessel_Diameter, OrgData.Notes.vessel.modalFixedDiameter, OrgData.Notes.frameRate);
    OrgData.Notes.vessel.diamPerc = diamPerc;
    OrgData.Notes.vessel.power_f = f;
    OrgData.Notes.vessel.power_S = S;
catch
    disp([OrgData.Notes.imageID ' Diameter percentage analysis failed!'])
end

end

%% Opens the tiff file and gets the  vessel projections from the defined polygons
function [OrgData] = GetDiametersFromMovie(OrgData, fileID)
OrgData.Notes.firstFrame = imread(fileID, 'TIFF', 'Index', 1);
fft_firstFrame = fft2(double(OrgData.Notes.firstFrame));
X = repmat(1:OrgData.Notes.xSize, OrgData.Notes.ySize, 1);
Y = repmat((1:OrgData.Notes.ySize)', 1, OrgData.Notes.xSize);
OrgData.Notes.vessel.projectionAngle = atand(diff(OrgData.Notes.vessel.vesselLine.position.xy(:, 1))/diff(OrgData.Notes.vessel.vesselLine.position.xy(:, 2)));
atand(diff(OrgData.Notes.vessel.vesselLine.position.xy(:, 1))/diff(OrgData.Notes.vessel.vesselLine.position.xy(:, 2)));

for theFrame = OrgData.Notes.startframe:OrgData.Notes.endframe
    rawFrame = imread(fileID, 'TIFF', 'Index', theFrame);
    fft_rawFrame = fft2(double(rawFrame));
    
    [OrgData.Notes.pixelShift(:, theFrame), ~] = dftregistration(fft_firstFrame, fft_rawFrame, 1);
    
    inpolyFrame = inpolygon(X + OrgData.Notes.pixelShift(3, theFrame), Y + OrgData.Notes.pixelShift(4, theFrame), OrgData.Notes.vessel.boxPosition.xy(:, 1), OrgData.Notes.vessel.boxPosition.xy(:, 2));
    bounded_rawFrame = rawFrame.*uint16(inpolyFrame);
    OrgData.Notes.vessel.projection(theFrame, :) = radon(bounded_rawFrame, OrgData.Notes.vessel.projectionAngle);
end

end

%% Calculate diameter using FWHM and get the baseline diameter
function [OrgData] = FWHM_MovieProjection(OrgData, theFrames)
for f = min(theFrames):max(theFrames)
    %a Add in a 5 pixel median filter
    OrgData.Data.Vessel_RawDiameter(f) = calcFWHM(medfilt1(OrgData.Notes.vessel.projection(f, :), 5));
end

OrgData.Data.Vessel_TempDiameter = OrgData.Data.Vessel_RawDiameter*OrgData.Notes.xFactor;
[holdHist, d] = hist(OrgData.Data.Vessel_TempDiameter, 0:.25:100);
[~, maxD] = max(holdHist);
OrgData.Notes.vessel.modalFixedDiameter = d(maxD);
end
