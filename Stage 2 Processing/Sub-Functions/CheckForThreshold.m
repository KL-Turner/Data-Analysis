function [ok] = CheckForThreshold(sfield, animal)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%
% Originally written by Aaron T. Winder
%
%   Last Revised: August 8th, 2018
%________________________________________________________________________________________________________________________
%
%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, Nov 2013
%   Version 1
%
%   SUMMARY: This function checks shared variable for baselines that can 
%               be used to normalize data.
%_______________________________________________________________
%   INPUTS:
%                           sfield - a string giving the subfield name
%                           (second level) of the shared variables
%                           structure
%                                       BinWWF - for binarizing whisker
%                                       angle measurements
%                                       BinPSWF - for binarizing pressure
%                                       sensor measurements
%                                       Spike - for detecting a spike
%
%                           animal - animal name as a string
%
%                           hem - recorded hemisphere as a string
%                                   [RH/LH]
%_______________________________________________________________
%   OUTPUTS:
%                           ok - output of 1 or 0 to indicate whether
%                           hemo baseline exist for the given day
%
%_______________________________________________________________
%   REQUIRED SCRIPTS:
%                           [Q] = DetectMachine(prevpath)
%_______________________________________________________________
%   CALLED BY:
%                           HemoBaseline.m
%_______________________________________________________________
%   FUTURE VERSIONS: Make the code generalizable to other shared variables
%_______________________________________________________________
%   CHANGES FROM PREV VERS:
%_______________________________________________________________

% Navigate to Shared Variables folder
disp(['CheckForThreshold.m: Checking for Threshold field: ' sfield '...']); disp(' ')
% Begin Check
ok = 0;
if exist([animal '_Thresholds.mat'],'file') == 2
    load([animal '_Thresholds.mat']);
    if isfield(Thresholds, sfield)
        ok = 1;
        disp(['CheckForThreshold.m: Threshold: ' sfield ' found.']); disp(' ')
    else
        disp(['CheckForThreshold.m: Threshold: ' sfield ' not found.']); disp(' ')
    end
end

end
