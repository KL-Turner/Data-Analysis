function [] = CheckForNaNs_IOS(ProcData,imagingType)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: 
%________________________________________________________________________________________________________________________
%
%   Inputs: 
%
%   Outputs:     
%
%   Last Revised: June 27th, 2019    
%________________________________________________________________________________________________________________________

if strcmp(imagingType,'bilateral') == true
    ROInames = {'LH','RH'};
elseif strcmp(imagingType,'single') == true
    ROInames = {'Barrels'};
end

for b = 1:length(ROInames)
    nanCheck{b,1} = sum(isnan(ProcData.data.CBV.(ROInames{1,b}))); %#ok<AGROW>
end

for b = 1:length(nanCheck)
    if nanCheck{b,1} ~= 0 
        disp('WARNING - NaNs found in CBV array'); disp(' ')
        keyboard
    end
end

end