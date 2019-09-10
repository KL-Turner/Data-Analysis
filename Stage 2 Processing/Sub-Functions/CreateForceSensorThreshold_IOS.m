function [thresh] = CreateForceSensorThreshold_IOS(PSWF)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs: 
%
%   Last Revised: February 29th, 2019
%________________________________________________________________________________________________________________________

y = hilbert(diff(PSWF));
force = abs(y);
forceFig = figure;
isok = 'n';

while strcmp(isok,'y') == 0
    plot(force, 'k');
    thresh = input('No Threshold to binarize pressure sensor found. Please enter a threshold: '); disp(' ')
    binPSWF = BinarizeForceSensor_IOS(PSWF,thresh);
    binInds = find(binPSWF);
    subplot(211)
    plot(PSWF, 'k') 
    axis tight
    hold on
    scatter(binInds, max(PSWF)*ones(size(binInds)),'r');
    subplot(212) 
    plot(force, 'k')
    axis tight
    hold on
    scatter(binInds, max(force)*ones(size(binInds)),'r');
    isok = input('Is this threshold okay? (y/n) ','s'); disp(' ')
    hold off
end
close(forceFig)

end
