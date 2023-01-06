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


end