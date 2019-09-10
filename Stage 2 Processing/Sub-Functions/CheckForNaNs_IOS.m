function [] = CheckForNaNs_IOS(ProcData)
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

CBVfields = fieldnames(ProcData.data.CBV);
for b = 1:length(CBVfields)
    nanCheck{b,1} = sum(isnan(ProcData.data.CBV.(CBVfields{b})));
end

for b = 1:length(nanCheck)
    if nanCheck{b,1} ~= 0
        disp('WARNING - NaNs found in CBV array'); disp(' ')
%         keyboard
    end
end

end