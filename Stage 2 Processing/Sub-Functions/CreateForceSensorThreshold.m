function [thresh] = CreateForceSensorThreshold(PSWF)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%
% Originally written by Aaron T. Winder
%
%   Last Revised: August 4th, 2018
%________________________________________________________________________________________________________________________
%
%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, Dec 2013
%   Version 1
%
%   SUMMARY:
%_______________________________________________________________
%   INPUTS:
%_______________________________________________________________
%   OUTPUTS:
%_______________________________________________________________

y = hilbert(diff(PSWF));
force = abs(y);
figure;
isok = 'n';

while strcmp(isok,'y') == 0
    plot(force, 'k');
    thresh = input('No Threshold to binarize pressure sensor found. Please enter a threshold: '); disp(' ')
    binPSWF = BinarizeForceSensor(PSWF,thresh);
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

end
