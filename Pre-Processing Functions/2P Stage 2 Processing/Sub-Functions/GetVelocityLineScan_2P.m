function [MScanData] = GetVelocityLineScan_2P(MScanData,fileID)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________

% takes a tiff file (movie) and an analog ascii file, and extracts the line scan velocity
MScan_analogData = [fileID '.TXT'];
disp(['Loading MScan file: ' MScan_analogData '...']); disp(' ');
analogData = load(MScan_analogData,'-ascii');
MScanData.data.corticalNeural = analogData(:,2);
MScanData.data.EMG = analogData(:,3);
MScanData.data.hippocampalNeural = analogData(:,4);
MScanData.data.forceSensor = analogData(:,5);
MScanData.notes.analogSamplingRate = 20000;
maxVelocity = 10000; % maximum physiological velocity in micrometers/sec
angleSpan = 15; % number of degrees (+/-) around previous angle to look
angleResolution = .1; % how accurate to determine the angle
scanMirrors = 1; % we only use 6215 - but 1) 6210(fast) or 2)6215(slow)
if scanMirrors == 1
    MScanData.notes.scanMirrors = 6210;
elseif scanMirrors==2
    MScanData.notes.scanMirrors = 6215;
end
disp(['Processing ' MScanData.notes.date '_' MScanData.notes.imageID '...']); disp(' ')
currentFilePath = mfilename('fullpath');
MScanData.notes.code = GetFunctionCode_2P(currentFilePath);
[MScanData] = GetVelocityRadon_2P(fileID,angleSpan,angleResolution,maxVelocity,scanMirrors,MScanData);
try
    MScanData.notes.radon_code = GetFunctionCode_2P(currentFilePath);
catch
    disp('code read fail!')
    keyboard
end
[MScanData,xcfig]=linescan_xcov_velocity_TIFF_2P(MScanData);%calculate the velocity using the cross correlation method of Kim et al, 2012

end
