function [ComparisonData] = AnalyzeCBV_STD(animal, dataType, params, procDataFiles, RestingBaselines, RestData, SleepData, ComparisonData)
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
%   Outputs: //
%________________________________________________________________________________________________________________________

%% Load in relevant data from the RestData struct:
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};

PuffCriteria.Fieldname = {'puffDistances'};
PuffCriteria.Comparison = {'gt'};
PuffCriteria.Value = {5};

% Use the RestCriteria we specified earlier to find all resting events that are greater than the criteria
[restLogical] = FilterEvents(RestData.CBV.(dataType), RestCriteria);   % Output is a logical
[puffLogical] = FilterEvents(RestData.CBV.(dataType), PuffCriteria);   % Output is a logical
combRestLogical = logical(restLogical.*puffLogical);
allRestFiles = RestData.CBV.(dataType).fileIDs(combRestLogical, :);   % Overall logical for all resting file names that meet criteria

allRestingCBVData = RestData.CBV.(dataType).NormData(combRestLogical, :);   % Pull out data from all those resting files that meet criteria

uniqueDays = GetUniqueDays(RestData.CBV.(dataType).fileIDs);   % Find the unique days of imaging
uniqueFiles = unique(RestData.CBV.(dataType).fileIDs);   % Find the unique files from the filelist. This removes duplicates
% since most files have more than one resting event
numberOfFiles = length(unique(RestData.CBV.(dataType).fileIDs));   % Find the number of unique files
fileTarget = params.targetMinutes / 5;   % Divide that number of unique files by 5 (minutes) to get the number of files that
% corresponds to the desired targetMinutes

% Loop through each unique day in order to create a logical to filter the file list so that it only includes the first
% x number of files that fall within the targetMinutes requirement
for uD = 1:length(uniqueDays)
    day = uniqueDays(uD);
    x = 1;
    for nOF = 1:numberOfFiles
        file = uniqueFiles(nOF);
        fileID = file{1}(1:6);
        if strcmp(params.Infusion, 'y')
            if strcmp(day, fileID) && x > fileTarget
                filtLogical{uD, 1}(nOF, 1) = 1;
            else
                filtLogical{uD, 1}(nOF, 1) = 0;
                x = x + 1;
            end
        else
            if strcmp(day, fileID) && x <= fileTarget
                filtLogical{uD, 1}(nOF, 1) = 1;
                x = x + 1;
            else
                filtLogical{uD, 1}(nOF, 1) = 0;
            end
        end
    end
end

% Combine the 3 logicals so that it reflects the first "x" number of files from each day
finalLogical = any(sum(cell2mat(filtLogical'), 2), 2);

% Now that the appropriate files from each day are identified, loop through each file name with respect to the original
% list of all relevant files, only keeping the ones that fall within the first targetMinutes of each day.
filtRestFiles = uniqueFiles(finalLogical, :);
for rF = 1:length(allRestFiles)
    logic = strcmp(allRestFiles{rF}, filtRestFiles);
    logicSum = sum(logic);
    if logicSum == 1
        fileFilter(rF, 1) = 1;
    else
        fileFilter(rF, 1) = 0;
    end
end

finalFileFilter = logical(fileFilter);
finalRestCBVData = allRestingCBVData(finalFileFilter, :);

[B, A] = butter(4, 2 / (30 / 2), 'low');
for ii = 1:length(finalRestCBVData)
    finalRestCBV_means(ii, 1) = mean(filtfilt(B, A, finalRestCBVData{ii, 1}));
    finalRestCBV_STDs(ii, 1) = std(filtfilt(B, A, finalRestCBVData{ii, 1}));
end

nboot = 1000;
restCBV_CI = bootci(nboot, @mean, finalRestCBVData);
restCBV_mean = mean(finalRestCBV_means);
restCBV_STD = mean(finalRestCBV_STDs);

disp([(dataType) ' Resting CBV mean is  ' num2str(restCBV_mean) ' with a STD of: ' num2str(restCBV_STD)]); disp(' ')

ComparisonData.AvgCBV_STD.Rest.(dataType).mean = restCBV_mean;
ComparisonData.AvgCBV_STD.Rest.(dataType).std = restCBV_STD;

%% Load in relevant data from SleepData
if ~isempty(SleepData)
    nremCBV = SleepData.NREM.Data.CBV.(dataType);
    
    for n = 1:length(nremCBV)
        nremCBV_means(:, n) = mean(nremCBV{n, 1});
        nremCBV_STDs(:, n) = std(nremCBV{n, 1});
    end
    
    nremCBV_CI = bootci(nboot, @mean, nremCBV);
    nremCBV_mean = mean(nremCBV_means);
    nremCBV_STD = mean(nremCBV_STDs);
    
    disp([(dataType) ' NREM CBV mean is  ' num2str(nremCBV_mean) ' with a STD of: ' num2str(nremCBV_STD)]); disp(' ')
    
    ComparisonData.AvgCBV_STD.NREM.(dataType).mean = nremCBV_mean;
    ComparisonData.AvgCBV_STD.NREM.(dataType).std = nremCBV_STD;
    
    if ~isempty(SleepData.REM)
        remCBV = SleepData.REM.Data.CBV.(dataType);
        
        for n = 1:length(remCBV)
            remCBV_means(:, n) = mean(remCBV{n, 1});
            remCBV_STDs(:, n) = std(remCBV{n, 1});
        end
        
        remCBV_mean = mean(remCBV_means);
        remCBV_STD = mean(remCBV_STDs);
        
        disp([(dataType) ' REM CBV mean is  ' num2str(remCBV_mean) ' with a STD of: ' num2str(remCBV_STD)]); disp(' ')
        
        ComparisonData.AvgCBV_STD.REM.(dataType).mean = remCBV_mean;
        ComparisonData.AvgCBV_STD.REM.(dataType).std = remCBV_STD;
    end
end

%% Load in relevant data from all ProcDataFiles
for pDF = 1:size(procDataFiles, 1)
    procDataFile = procDataFiles(pDF, :);
    load(procDataFile);
    [~, ~, fileDate, fileID] = GetFileInfo_IOS(procDataFile);
    strDay = ConvertDate(fileDate);
    CBV{pDF, 1} = (ProcData.Data.CBV.(dataType) - RestingBaselines.CBV.(dataType).(strDay)) / RestingBaselines.CBV.(dataType).(strDay);
end

for ii = 1:length(CBV)
    allCBV_means(ii, 1) = mean(filtfilt(B, A, CBV{ii, 1}));
    allCBV_STDs(ii, 1) = std(filtfilt(B, A, CBV{ii, 1}));
end

allCBV_mean = mean(allCBV_means);
allCBV_STD = mean(allCBV_STDs);

disp([(dataType) ' all data CBV mean is  ' num2str(allCBV_mean) ' with a STD of: ' num2str(allCBV_STD)]); disp(' ')

ComparisonData.AvgCBV_STD.AllData.(dataType).mean = allCBV_mean;
ComparisonData.AvgCBV_STD.AllData.(dataType).std = allCBV_STD;

save([animal '_ComparisonData.mat'], 'ComparisonData');

end