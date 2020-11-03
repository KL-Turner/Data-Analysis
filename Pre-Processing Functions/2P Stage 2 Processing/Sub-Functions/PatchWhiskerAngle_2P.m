function [patchedWhiskerAngle,sampleDiff] = PatchWhiskerAngle_2P(whiskerAngle,fs,expectedDuration_Sec,droppedFrameIndex)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: The whisker camera occasionally drops packets of frames. We can calculate the difference in the number
%            of expected frames as well as the indeces that LabVIEW found the packets lost. This is a rough fix, as
%            we are not sure the exact number of frames at each index, only the total number.
%________________________________________________________________________________________________________________________

expectedSamples = expectedDuration_Sec*fs;
droppedFrameIndex = str2num(droppedFrameIndex); %#ok<ST2NM>
sampleDiff = expectedSamples - length(whiskerAngle);
framesPerIndex = ceil(sampleDiff/length(droppedFrameIndex));
% loop through each index, linear interpolate the values between the index and the right edge, then shift the samples.
% take into account that the old dropped frame indeces will no longer correspond to the new length of the array.
if ~isempty(droppedFrameIndex)
    % each dropped index
    for x = 1:length(droppedFrameIndex)
        % for the first event, it's okay to start at the actual index
        if x == 1
            leftEdge = (droppedFrameIndex(1,x));
        else
        % for all other dropped frames after the first, we need to correct for the fact that index is shifted right.
            leftEdge = (droppedFrameIndex(1,x)) + ((x - 1)*framesPerIndex);
        end
        % set the edges for the interpolation points. we want n number of samples between the two points,vthe left and
        % right edge values. This equates to having a 1/(dropped frames + 1) step size between the edges.
        rightEdge = leftEdge + 1;
        patchFrameInds = leftEdge:(1/(framesPerIndex + 1)):rightEdge;
        % concatenate the original whisker angle for the first index, then the new patched angle for all subsequent
        % indeces. Take the values from 1:left edge, add in the new frames, then right edge to end.
        if x == 1
            patchFrameVals = interp1(1:length(whiskerAngle),whiskerAngle,patchFrameInds);   % linear interp
            snipPatchFrameVals = patchFrameVals(2:end - 1);
            patchedWhiskerAngle = horzcat(whiskerAngle(1:leftEdge),snipPatchFrameVals,whiskerAngle(rightEdge:end));
        else
            patchFrameVals = interp1(1:length(patchedWhiskerAngle),patchedWhiskerAngle,patchFrameInds);   % linear interp
            snipPatchFrameVals = patchFrameVals(2:end - 1);
            patchedWhiskerAngle = horzcat(patchedWhiskerAngle(1:leftEdge),snipPatchFrameVals,patchedWhiskerAngle(rightEdge:end));
        end
    end
    patchedWhiskerAngle = patchedWhiskerAngle(1:expectedSamples);
else
    patchedWhiskerAngle = whiskerAngle(1:expectedSamples);
end
% due to rounding up on the number of dropped frames per index, we have a few extra frames. Snip them off.

end
