function [MscanData] = ExtractTiffAnalogData(MscanData, fileID)
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
MscanData.Data.MScan_Force_Sensor = analogData(:, 2);
MscanData.Data.MScan_Neural_Data = analogData(:, 3);
MscanData.Notes.MScan_analogSamplingRate = 20000;

disp('Analyzing vessel projections from defined polygons...'); disp(' ');
[MscanData] = GetDiametersFromMovie(MscanData, fileID);

try
    [MscanData] = FWHM_MovieProjection(MscanData, [MscanData.Notes.startframe MscanData.Notes.endframe]);
catch
    disp([MscanData.Notes.imageID ' FWHM calculation failed!'])
end

try
    % 1 dural/vein, >40% changes spline, artery: >60% spline
    % 2 dural/vein, >30% changes interpolate, artery: >50% interpolate
    if strcmp(MscanData.Notes.vesselType, 'D') || strcmp(MscanData.Notes.vesselType, 'V')
        MscanData.Data.Vessel_Diameter = RemoveMotion(MscanData.Data.Vessel_TempDiameter, MscanData.Notes.vessel.modalFixedDiameter, 2, 0.3);
    else
        MscanData.Data.Vessel_Diameter = RemoveMotion(MscanData.Data.Vessel_TempDiameter, MscanData.Notes.vessel.modalFixedDiameter, 2, 0.5);
    end
    [diamPerc, S, f] = DiamPercPower(MscanData.Data.Vessel_Diameter, MscanData.Notes.vessel.modalFixedDiameter, MscanData.Notes.frameRate);
    MscanData.Notes.vessel.diamPerc = diamPerc;
    MscanData.Notes.vessel.power_f = f;
    MscanData.Notes.vessel.power_S = S;
catch
    disp([MscanData.Notes.imageID ' Diameter percentage analysis failed!'])
end

end

%% Opens the tiff file and gets the  vessel projections from the defined polygons
function [MscanData] = GetDiametersFromMovie(MscanData, fileID)
MscanData.Notes.firstFrame = imread(fileID, 'TIFF', 'Index', 1);
fft_firstFrame = fft2(double(MscanData.Notes.firstFrame));
X = repmat(1:MscanData.Notes.xSize, MscanData.Notes.ySize, 1);
Y = repmat((1:MscanData.Notes.ySize)', 1, MscanData.Notes.xSize);
MscanData.Notes.vessel.projectionAngle = atand(diff(MscanData.Notes.vessel.vesselLine.position.xy(:, 1))/diff(MscanData.Notes.vessel.vesselLine.position.xy(:, 2)));
atand(diff(MscanData.Notes.vessel.vesselLine.position.xy(:, 1))/diff(MscanData.Notes.vessel.vesselLine.position.xy(:, 2)));

for theFrame = MscanData.Notes.startframe:MscanData.Notes.endframe
    rawFrame = imread(fileID, 'TIFF', 'Index', theFrame);
    fft_rawFrame = fft2(double(rawFrame));
    
    [MscanData.Notes.pixelShift(:, theFrame), ~] = DftRegistration(fft_firstFrame, fft_rawFrame, 1);
    
    inpolyFrame = inpolygon(X + MscanData.Notes.pixelShift(3, theFrame), Y + MscanData.Notes.pixelShift(4, theFrame), MscanData.Notes.vessel.boxPosition.xy(:, 1), MscanData.Notes.vessel.boxPosition.xy(:, 2));
    bounded_rawFrame = rawFrame.*uint16(inpolyFrame);
    MscanData.Notes.vessel.projection(theFrame, :) = radon(bounded_rawFrame, MscanData.Notes.vessel.projectionAngle);
end

end

%% Calculate diameter using FWHM and get the baseline diameter
function [MscanData] = FWHM_MovieProjection(MscanData, theFrames)
for f = min(theFrames):max(theFrames)
    %a Add in a 5 pixel median filter
    MscanData.Data.Vessel_RawDiameter(f) = CalcFWHM(medfilt1(MscanData.Notes.vessel.projection(f, :), 5));
end

MscanData.Data.Vessel_TempDiameter = MscanData.Data.Vessel_RawDiameter*MscanData.Notes.xFactor;
[holdHist, d] = hist(MscanData.Data.Vessel_TempDiameter, 0:.25:100);
[~, maxD] = max(holdHist);
MscanData.Notes.vessel.modalFixedDiameter = d(maxD);
end
