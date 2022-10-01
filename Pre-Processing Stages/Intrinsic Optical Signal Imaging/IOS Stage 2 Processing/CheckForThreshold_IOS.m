function [ok] = CheckForThreshold_IOS(sfield,animal)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Check to see if a threshold has already been created for a given field in a given day.
%________________________________________________________________________________________________________________________

% begin Check
ok = 0;
if exist([animal '_Thresholds.mat'],'file') == 2
    load([animal '_Thresholds.mat']);
    if isfield(Thresholds,sfield)
        ok = 1;
    end
end

end
