function [linkedWF] = LinkBinaryEvents(binWF, dCrit)
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
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Takes a binary waveform and links the peaks which
%   occur within a close period of time of each other creating a single 
%   peak.   
%_______________________________________________________________
%   PARAMETERS:             
%                   bin_wf - [Array] The binary waveform
%
%                   dCrit - [1x2 Array] The distances, in samples, between 
%                   the falling edge of the previous event and the rising 
%                   edge of the current waveform. 
%                       Input should be given as a 2D array: 
%                           [dCrit for breaks ,dCrit for peaks]. 
%                       Any events with breaks between them less than dCrit
%                           for breaks will be merged. 
%                       Any events that have duration less than dCrit for 
%                           peaks will be removed. 
%                       Either dCrit value can be zero to avoid merging or 
%                           erasing events.                    
%_______________________________________________________________
%   RETURN:                     
%                   linked_wf - [Array] The new, linked binary waveform as 
%                   an array.
%_______________________________________________________________

%% Identify Edges, control for trial start/stop
dBinWF = diff(gt(binWF, 0));
upInd = find(dBinWF == 1);
downInd = find(dBinWF == -1);
if binWF(end) > 0
    downInd = [downInd length(binWF)];
end
if binWF(1) > 0
    upInd = [1 upInd];
end

%% Link periods of bin_wf==0 together if less than dCrit(1)
% Calculate time between events
brkTimes = upInd(2:length(upInd)) - downInd(1:(length(downInd) - 1));
% Identify times less than user-defined period
sub_dCritDowns = find(lt(brkTimes, dCrit(1)));

% Link any identified breaks together
if isempty(sub_dCritDowns) == 0
    for d = 1:length(sub_dCritDowns)
        start = downInd(sub_dCritDowns(d));
        stop = upInd(sub_dCritDowns(d) + 1);
        binWF(start:stop) = 1;
    end
end

%% Link periods of bin_wf==1 together if less than dCrit(2)
hitimes = downInd - upInd;
blips = find(lt(hitimes, dCrit(2)) == 1);
if isempty(blips) == 0
    for b = 1:length(blips)
        start = upInd(blips(b));
        stop = downInd(blips(b));
        binWF(start:stop) = 0;
    end
end

linkedWF = binWF;

end
