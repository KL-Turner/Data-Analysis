function [RestingBaselines] = CalculateSpectrogramBaselines_IOS(animal, neuralDataTypes, trialDuration_sec, specDataFiles, RestingBaselines, baselineType)
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

for a = 1:length(neuralDataTypes)
    neuralDataType = neuralDataTypes{1,a};
    dsFs = 30;
    restFileList = unique(RestingBaselines.(baselineType).baselineFileInfo.fileIDs);      % Obtain the list of unique fileIDs
    restS1 = cell(size(restFileList,1), 1);
    restS5 = cell(size(restFileList,1), 1);
    % Obtain the spectrogram information from all the resting files
    for b = 1:length(restFileList)
        fileID = restFileList{b, 1};   % FileID of currently loaded file
        % Load in neural data from current file
        for c = 1:size(specDataFiles, 1)
            [~, ~, specDataFile] = GetFileInfo_IOS(specDataFiles(c,:));
            if strcmp(fileID, specDataFile)
                load(specDataFiles(c,:), '-mat')
                S1 = SpecData.(neuralDataType).oneSec.S;
                S5 = SpecData.(neuralDataType).fiveSec.S;
                break
            end
        end
        restS1{b,1} = S1;
        restS5{b,1} = S5;
    end
    
    for d = 1:length(restFileList)
        fileID = restFileList{d,1};
        strDay = ConvertDate_IOS(fileID(1:6));
        S1_data = restS1{d,1};
        S5_data = restS5{d,1};
        s1Length = size(S1_data,2);
        s5Length = size(S5_data,2);
        binSize1 = ceil(s1Length/trialDuration_sec);
        binSize5 = ceil(s5Length/trialDuration_sec);
        samplingDiff1 = dsFs/binSize1;
        samplingDiff5 = dsFs/binSize5;
        S1_trialRest = [];
        S5_trialRest = [];
        for e = 1:length(RestingBaselines.(baselineType).baselineFileInfo.fileIDs)
            restFileID = RestingBaselines.(baselineType).baselineFileInfo.fileIDs{e, 1};
            if strcmp(fileID, restFileID)
                restDuration1 = floor(floor(RestingBaselines.(baselineType).baselineFileInfo.durations(e, 1)*dsFs) / samplingDiff1);
                restDuration5 = floor(floor(RestingBaselines.(baselineType).baselineFileInfo.durations(e, 1)*dsFs) / samplingDiff5);
                startTime1 = floor(floor(RestingBaselines.(baselineType).baselineFileInfo.eventTimes(e, 1)*dsFs) / samplingDiff1);
                startTime5 = floor(floor(RestingBaselines.(baselineType).baselineFileInfo.eventTimes(e, 1)*dsFs) / samplingDiff5);
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
    uniqueDays = GetUniqueDays_IOS(RestingBaselines.(baselineType).baselineFileInfo.fileIDs);
    
    for f = 1:length(uniqueDays)
        g = 1;
        for field = 1:length(fields)
            if strcmp(fields{field}(7:12), uniqueDays{f})
                stringDay = ConvertDate_IOS(uniqueDays{f});
                S_avgs.oneSec.(stringDay){g, 1} = trialRestData.(fields{field}).oneSec.S_avg;
                S_avgs.fiveSec.(stringDay){g, 1} = trialRestData.(fields{field}).fiveSec.S_avg;
                g = g + 1;
            end
        end
    end
    
    dayFields = fieldnames(S_avgs.oneSec);
    for h = 1:length(dayFields)
        dayVals1 = [];
        dayVals5 = [];
        for j = 1:length(S_avgs.oneSec.(dayFields{h}))
            dayVals1 = [dayVals1, S_avgs.oneSec.(dayFields{h}){j, 1}];
            dayVals5 = [dayVals5, S_avgs.fiveSec.(dayFields{h}){j, 1}];
        end
        disp(['Adding spectrogram baseline to baseline file for ' neuralDataType ' on ' dayFields{h} '...']); disp(' ')
        RestingBaselines.Spectrograms.(neuralDataType).oneSec.(dayFields{h}) = mean(dayVals1, 2);
        RestingBaselines.Spectrograms.(neuralDataType).fiveSec.(dayFields{h}) = mean(dayVals5, 2);
    end
end

save([animal '_RestingBaselines.mat'], 'RestingBaselines');

end
