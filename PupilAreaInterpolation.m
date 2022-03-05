clear; clc; close all;
% load example file
exampleProcDataFileID = 'T141_201105_12_05_20_ProcData.mat';
load(exampleProcDataFileID,'-mat')
% pupil area and derivative
pupilArea = ProcData.data.Pupil.pupilArea;
diffArea = abs(diff(pupilArea));
% threshold for interpolation
threshold = 250;
diffIndex = diffArea > threshold;
% derivative figure
figure; 
p1 = plot(pupilArea);
hold on
p2 = plot(diffArea);
p3 = yline(threshold);
axis tight
legend([p1,p2,p3],'original','abs(1st deriv)','manual threshold')
% derivative index figure
figure
p1 = plot(pupilArea);
hold on
p2 = plot(diffIndex*max(pupilArea));
axis tight
legend([p1,p2],'original','diff index')
% link adjacent indeces
[linkedDiffIndex] = LinkBinaryEvents(gt(diffIndex,0),[30,0]);
figure
p1 = plot(pupilArea);
hold on
p2 = plot(linkedDiffIndex*max(pupilArea));
axis tight
legend([p1,p2],'original','linked diff index')
% obtain edge pairs for interpolation
edgeFound = false;
xx = 1;
for aa = 1:length(linkedDiffIndex)
    if edgeFound == false
        if linkedDiffIndex(1,aa) == 1
            startEdge(xx,1) = aa;
            edgeFound = true;
        end
    elseif edgeFound == true
        if linkedDiffIndex(1,aa) == 0
            endEdge(xx,1) = aa;
            xx = xx + 1;
            edgeFound = false;
        end
    end
end
% edge figure
figure
p1 = plot(pupilArea,'k');
hold on
for aa = 1:length(startEdge)
    x1 = xline(startEdge(aa,1),'r');
end
for aa = 1:length(endEdge)
    x2 = xline(endEdge(aa,1),'b');
end
axis tight
legend([p1,x1,x2],'original','start edge','end edge')
% set time between edges as NaN
nanPupilArea = pupilArea;
for aa = 1:length(startEdge)
    startTime = startEdge(aa,1);
    endTime = endEdge(aa,1);
    nanPupilArea(startTime:endTime) = NaN;
    patchedPupilArea = fillmissing(nanPupilArea,'spline');
end
% comparison figure
figure
p1 = plot(pupilArea);
hold on
p2 = plot(patchedPupilArea);
axis tight
legend([p1,p2],'original','patched')
% check any nans in final array
disp(['Number of NaNs in final array: ' num2str(sum(isnan(patchedPupilArea)))]); disp(' ')

function [linkedWF] = LinkBinaryEvents(binWF,dCrit)

% Identify Edges, control for trial start/stop
dBinWF = diff(gt(binWF,0));
upInd = find(dBinWF == 1);
downInd = find(dBinWF == -1);
if binWF(end) > 0
    downInd = [downInd,length(binWF)];
end
if binWF(1) > 0
    upInd = [1,upInd];
end
% Link periods of bin_wf==0 together if less than dCrit(1). Calculate time between events
brkTimes = upInd(2:length(upInd)) - downInd(1:(length(downInd) - 1));
% Identify times less than user-defined period
sub_dCritDowns = find(lt(brkTimes,dCrit(1)));
% Link any identified breaks together
if isempty(sub_dCritDowns) == 0
    for d = 1:length(sub_dCritDowns)
        start = downInd(sub_dCritDowns(d));
        stop = upInd(sub_dCritDowns(d) + 1);
        binWF(start:stop) = 1;
    end
end
% Link periods of bin_wf==1 together if less than dCrit(2)
hitimes = downInd - upInd;
blips = find(lt(hitimes,dCrit(2)) == 1);
if isempty(blips) == 0
    for b = 1:length(blips)
        start = upInd(blips(b));
        stop = downInd(blips(b));
        binWF(start:stop) = 0;
    end
end
linkedWF = binWF;

end
