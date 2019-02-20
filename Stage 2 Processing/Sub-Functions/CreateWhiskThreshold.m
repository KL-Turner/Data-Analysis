function [thresh1,thresh2] = CreateWhiskThreshold(angl, fs)
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

isok = 'n';
dd_wwf = abs((diff(angl, 2)))*fs^2;

while strcmp(isok,'y') == 0
    close all; 
    plot(dd_wwf,'k');
    thresh2 = input('No Threshold for volitional whisks found. Please enter a threshold: '); disp(' ')
    thresh1 = input('No Threshold for resting behavior found. Please enter a threshold: '); disp(' ')
    bin_wwf = BinarizeWhiskers(angl, fs, thresh1, thresh2);
    ax1 = subplot(311); 
    plot(angl, 'k'); 
    axis tight; 
    ylabel('Angle')
    ax2 = subplot(312);
    plot(abs(diff(angl, 2))*fs^2, 'k'); 
    axis tight; 
    ylabel('Second Derivative')
    ax3 = subplot(313); 
    plot(bin_wwf, 'k'); 
    axis tight;
    ylabel('Binarization')
    linkaxes([ax1, ax2, ax3], 'x');
    isok = input('Is this threshold okay? (y/n) ','s'); disp(' ')
end

end

