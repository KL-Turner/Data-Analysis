function [imagingColors] = SelectWavelengthType_IOS(wavelengthOptions)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Select option from menu
%________________________________________________________________________________________________________________________

tf = 0;
while tf == 0
    disp('Select imaging wavelengths'); disp(' ')
    [indx,tf] = listdlg('PromptString',{'Select imaging wavelength(s)',''},'SelectionMode','single','ListString',wavelengthOptions);
    if tf == 0
        disp('Please select wavelength type'); disp(' ');
    else
        disp(['Wavelength(s): ' num2str(wavelengthOptions{1,indx}) ' selected']); disp(' ')
    end
end
wavelengthType = wavelengthOptions{1,indx};
% use wavelength type to determine color of imaging/spectroscopy 
if strcmpi(wavelengthType,'1 wavelength (G/L)') == true
    imagingColors = 'G';
elseif strcmpi(wavelengthType,'1 wavelength (B)') == true
    imagingColors = 'B';
elseif strcmpi(wavelengthType,'2 wavelenths (G/L & B)') == true
    imagingColors = 'GB';
elseif strcmpi(wavelengthType,'3 wavelengths (R, G/L, & B)') == true
    imagingColors = 'RGB';
end

end