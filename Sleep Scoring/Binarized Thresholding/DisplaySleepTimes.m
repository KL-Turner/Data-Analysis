function [] = DisplaySleepTimes(animal, SleepData)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: 
%________________________________________________________________________________________________________________________
%
%   Inputs: 
%   
%   Outputs: 
%________________________________________________________________________________________________________________________

%% NREM
fileList = ls('*_WindowCam.bin');
uniqueDays = GetUniqueDays(fileList);
numberOfDays = length(uniqueDays);
numberOfFiles = length(fileList);
filesPerDay = numberOfFiles / numberOfDays;
NREMsleepFileList = SleepData.NREM.FileIDs;
NREMsleepLogical = zeros(1, length(fileList));

numberOfNREMSleepEvents = length(SleepData.NREM.BinTimes);
totalNREMSleepTime = 0;
for ii = 1:length(SleepData.NREM.BinTimes)
    startTime = SleepData.NREM.BinTimes{ii, 1}(1);
    endTime = SleepData.NREM.BinTimes{ii, 1}(end);
    sleepTime = endTime - startTime;
    totalNREMSleepTime = totalNREMSleepTime + sleepTime;
end
totalNREMSleepTime = totalNREMSleepTime / 60; % Minutes

for ii = 1:size(fileList, 1)
    NREMfileID = fileList(ii, :);
    NREMsleepVal = NaN;
    for iii = 1:length(NREMsleepFileList)
    NREMsleepFileID = NREMsleepFileList{iii, 1};
        if strcmp(NREMfileID(1:15), NREMsleepFileID)
            NREMsleepVal = 1;
        end
    end
    NREMsleepLogical(ii) = NREMsleepVal;
end

NREMdayIndex(1) = 0;
for x = 1:length(uniqueDays)
    NREMdayIndex(x + 1) = x*filesPerDay;
end

%% REM
if ~isempty(SleepData.REM)
    REMsleepFileList = SleepData.REM.FileIDs;
    REMsleepLogical = zeros(1, length(fileList));
    numberOfREMSleepEvents = length(SleepData.REM.BinTimes);
    
    totalREMSleepTime = 0;
    for ii = 1:length(SleepData.REM.BinTimes)
        startTime = SleepData.REM.BinTimes{ii, 1}(1);
        endTime = SleepData.REM.BinTimes{ii, 1}(end);
        sleepTime = endTime - startTime;
        totalREMSleepTime = totalREMSleepTime + sleepTime;
    end
    totalREMSleepTime = totalREMSleepTime / 60; % Minutes
    
    for ii = 1:length(fileList)
        REMfileID = fileList(ii, :);
        REMsleepVal = NaN;
        for iii = 1:length(REMsleepFileList)
            REMsleepFileID = REMsleepFileList{iii, 1};
            if strcmp(REMfileID(1:15), REMsleepFileID)
                REMsleepVal = 1.2;
            end
        end
        REMsleepLogical(ii) = REMsleepVal;
    end
    
    REMdayIndex(1) = 0;
    for x = 1:length(uniqueDays)
        REMdayIndex(x + 1) = x*filesPerDay;
    end
end

%% Figure
SleepTimes = figure;
scatter(1:length(NREMsleepLogical), NREMsleepLogical, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c')
hold on
if ~isempty(SleepData.REM)
    scatter(1:length(REMsleepLogical), REMsleepLogical, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors('deep carrot orange'))
    ylim([0 2])
end

y = get(gca,'ylim');
for ii = 1:length(NREMdayIndex)
    plot([NREMdayIndex(ii) NREMdayIndex(ii)], y, 'k')
end

title([animal ' Sleeping Events Across Imaging Days'])
if ~isempty(SleepData.REM)
    xlabel({'File number (Vertical lines denote a new imaging day)'; [num2str(totalNREMSleepTime) ' minutes of NREM Sleep found across ' num2str(numberOfNREMSleepEvents) ' unique events.']; ...
        [num2str(totalREMSleepTime) ' minutes of REM Sleep found across ' num2str(numberOfREMSleepEvents) ' unique events.']});
else
    xlabel({'File number (Vertical lines denote a new imaging day)'; [num2str(totalNREMSleepTime) ' minutes of NREM Sleep found across ' num2str(numberOfNREMSleepEvents) ' unique events.']});
end
ylabel('Sleeping Events')
legend('NREM Sleep Events', 'REM Sleep Events')
xlim([0 numberOfFiles])
ylim([0 2])
set(gca, 'Ticklength', [0 0])

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Sleep Times/'];

if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end

savefig(SleepTimes, [dirpath animal '_SleepTimes']);

end