function [RestingBaselines] = CalculateSpectrogramBaselines_2P(animal, trialDuration_Sec, specDataFiles, RestingBaselines)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: Uses the resting time indeces to extract the average resting power in each frequency bin during periods of
%            rest to normalize the spectrogram data.
%________________________________________________________________________________________________________________________
%
%   Inputs: animal ID, trialDuration of the session, list of the SpecData.mat files, and the RestingBaselines.mat struct.
%
%   Outputs: Updates to the RestingBaselines.mat structure containing a resting frequency-dependent power for each day.
%________________________________________________________________________________________________________________________

dsFs = 30;
restFileList = unique(RestingBaselines.baselineFileInfo.fileIDs);      % Obtain the list of unique fileIDs
restS1 = cell(size(restFileList,1), 1);
restS5 = cell(size(restFileList,1), 1);
% Obtain the spectrogram information from all the resting files
for a = 1:length(restFileList)
    fileID = restFileList{a, 1};   % FileID of currently loaded file
    % Load in neural data from current file
    for b = 1:size(specDataFiles, 1)
        [~, ~, specDataFile, ~, ~] = GetFileInfo2_2P(specDataFiles(b,:));
        if strcmp(fileID, specDataFile)
            load(specDataFiles(b,:), '-mat')
            S1 = SpecData.oneSec.S;
            S5 = SpecData.fiveSec.S;
            break
        end
    end
    restS1{a,1} = S1;
    restS5{a,1} = S5;
end

for c = 1:length(restFileList)
    fileID = restFileList{c,1};
    strDay = ConvertDate_2P(fileID(1:6));
    S1_data = restS1{c,1};
    S5_data = restS5{c,1};
    s1Length = size(S1_data,2);
    s5Length = size(S5_data,2);                               
    binSize1 = ceil(s1Length/trialDuration_Sec); 
    binSize5 = ceil(s5Length/trialDuration_Sec);
    samplingDiff1 = dsFs/binSize1;
    samplingDiff5 = dsFs/binSize5;  
    S1_trialRest = [];
    S5_trialRest = [];
    for d = 1:length(RestingBaselines.baselineFileInfo.fileIDs)
        restFileID = RestingBaselines.baselineFileInfo.fileIDs{d, 1};
        if strcmp(fileID, restFileID)
            restDuration1 = floor(floor(RestingBaselines.baselineFileInfo.durations(d, 1)*dsFs) / samplingDiff1);
            restDuration5 = floor(floor(RestingBaselines.baselineFileInfo.durations(d, 1)*dsFs) / samplingDiff5);
            startTime1 = floor(floor(RestingBaselines.baselineFileInfo.eventTimes(d, 1)*dsFs) / samplingDiff1);
            startTime5 = floor(floor(RestingBaselines.baselineFileInfo.eventTimes(d, 1)*dsFs) / samplingDiff5);
            try
                S1_single_rest = S1_data(:, (startTime1:(startTime1 + restDuration1)));
                S5_single_rest = S5_data(:, (startTime5:(startTime5 + restDuration5)));
            catch
                S1_single_rest = S1_data(:, end - restDuration1:end);
                S5_single_rest = S5_data(:, end - restDuration5:end);
            end
            S1_trialRest = [S1_single_rest, S1_trialRest];
            S5_trialRest = [S5_single_rest, S5_trialRest];
        end
    end
    S_trialAvg1 = mean(S1_trialRest, 2);
    S_trialAvg5 = mean(S5_trialRest, 2);
    trialRestData.([strDay '_' fileID]).oneSec.S_avg = S_trialAvg1;
    trialRestData.([strDay '_' fileID]).fiveSec.S_avg = S_trialAvg5;
end

fields = fieldnames(trialRestData);
uniqueDays = GetUniqueDays_2P(RestingBaselines.baselineFileInfo.fileIDs);

for e = 1:length(uniqueDays)
    f = 1;
    for field = 1:length(fields)
        if strcmp(fields{field}(7:12), uniqueDays{e})
            stringDay = ConvertDate_2P(uniqueDays{e});
            S_avgs.oneSec.(stringDay){f, 1} = trialRestData.(fields{field}).oneSec.S_avg;
            S_avgs.fiveSec.(stringDay){f, 1} = trialRestData.(fields{field}).fiveSec.S_avg;
            f = f + 1;
        end
    end
end

dayFields = fieldnames(S_avgs.oneSec);
for g = 1:length(dayFields)
    dayVals1 = [];
    dayVals5 = [];
    for f = 1:length(S_avgs.oneSec.(dayFields{g}))
        dayVals1 = [dayVals1, S_avgs.oneSec.(dayFields{g}){f, 1}];
        dayVals5 = [dayVals5, S_avgs.fiveSec.(dayFields{g}){f, 1}];
    end
    disp(['Adding spectrogram baseline to baseline file for ' dayFields{g} '...']); disp(' ')
    RestingBaselines.Spectrograms.oneSec.(dayFields{g}) = mean(dayVals1, 2);
    RestingBaselines.Spectrograms.fiveSec.(dayFields{g}) = mean(dayVals5, 2);
end

save([animal '_RestingBaselines.mat'], 'RestingBaselines');

end
