function [thresh1,thresh2] = CreateWhiskThreshold_IOS(angl, fs)
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

isok = 'n';
dd_wwf = abs((diff(angl, 2)))*fs^2;
whiskFig = figure;

while strcmp(isok,'y') == 0
    plot(dd_wwf,'k');
    thresh2 = input('No Threshold for volitional whisks found. Please enter a threshold: '); disp(' ')
    thresh1 = input('No Threshold for resting behavior found. Please enter a threshold: '); disp(' ')
    bin_wwf = BinarizeWhiskers_IOS(angl, fs, thresh1, thresh2);
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
close(whiskFig)

end

