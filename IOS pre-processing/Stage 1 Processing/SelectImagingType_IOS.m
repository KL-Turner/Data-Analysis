function [imagingType] = SelectImagingType_IOS(imagingOptions)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
tf = 0;
while tf == 0
    disp('Select imaging type'); disp(' ')
    [indx,tf] = listdlg('PromptString',{'Select an option',''},'SelectionMode','single','ListString',imagingOptions);
    if tf ~= 0
        disp(['Option selected: ' num2str(imagingOptions{1,indx}) ' selected']); disp(' ')
    end
end
imagingType = imagingOptions{1,indx};