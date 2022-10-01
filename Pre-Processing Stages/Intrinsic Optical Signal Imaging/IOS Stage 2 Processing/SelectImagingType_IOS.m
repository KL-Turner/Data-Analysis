function [imagingType] = SelectImagingType_IOS(imagingOptions)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Select option from menu
%________________________________________________________________________________________________________________________

tf = 0;
while tf == 0
    disp('Select imaging type'); disp(' ')
    [indx,tf] = listdlg('PromptString',{'Select an ROI option',''},'SelectionMode','single','ListString',imagingOptions);
    if tf == 0
        disp('Please select an imaging type'); disp(' ');
    else
        disp(['Imaging type: ' num2str(imagingOptions{1,indx}) ' selected']); disp(' ')
    end
end
imagingType = imagingOptions{1,indx};

end