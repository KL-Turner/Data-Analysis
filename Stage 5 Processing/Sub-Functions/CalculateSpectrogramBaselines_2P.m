function [RestingBaselines] = CalculateSpectrogramBaselines_2P(animal, RestingBaselines, SpectrogramData)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: //
%________________________________________________________________________________________________________________________
%
%   Inputs: //
%
%   Outputs: //
%________________________________________________________________________________________________________________________


restFileList = unique(RestingBaselines.baselineFileInfo.fileIDs);      % Obtain the list of unique fileIDs

% Obtain the spectrogram information from all the resting files
for ii = 1:length(restFileList)
    fileID = restFileList{ii, 1};   % FileID of currently loaded file
    % Load in neural data from current file
    for iii = 1:length(SpectrogramData.FileIDs)
        if strcmp(fileID, SpectrogramData.FileIDs{iii, 1})
            S1 = SpectrogramData.OneSec.S{iii, 1};
            T1 = SpectrogramData.OneSec.T{iii, 1};
            F1 = SpectrogramData.OneSec.F{iii, 1};
            S5 = SpectrogramData.FiveSec.S{iii, 1};
            T5 = SpectrogramData.FiveSec.T{iii, 1};
            F5 = SpectrogramData.FiveSec.F{iii, 1};
        end
    end
    Neural_restS1{ii, 1} = S1;
    Neural_restT1{ii, 1} = T1;
    Neural_restF1{ii, 1} = F1;
    Neural_restS5{ii, 1} = S5;
    Neural_restT5{ii, 1} = T5;
    Neural_restF5{ii, 1} = F5;
end

for ii = 1:length(restFileList)
    fileID = restFileList{ii, 1};
    strDay = ConvertDate(fileID(1:6));
    S1_data = Neural_restS1{ii, 1};
    S5_data = Neural_restS5{ii, 1};
    s1Length = size(S1_data, 2);
    s5Length = size(S5_data, 2);                                  % Length of the data across time (number of samples)
    binSize1 = ceil(s1Length / 300);                              % Find the number of bins needed to divide this into 290 seconds
    binSize5 = ceil(s5Length / 300);                              % Find the number of bins needed to divide this into 290 seconds
    samplingRate = 28;
    samplingDiff1 = samplingRate / binSize1;
    samplingDiff5 = samplingRate / binSize5;  
    S1_trialRest = [];
    S5_trialRest = [];
    for iii = 1:length(RestingBaselines.baselineFileInfo.fileIDs)
        restFileID = RestingBaselines.baselineFileInfo.fileIDs{iii, 1};
        if strcmp(fileID, restFileID)
            restDuration1 = floor(floor(RestingBaselines.baselineFileInfo.durations(iii, 1)*samplingRate) / samplingDiff1);
            restDuration5 = floor(floor(RestingBaselines.baselineFileInfo.durations(iii, 1)*samplingRate) / samplingDiff5);
            startTime1 = floor(floor(RestingBaselines.baselineFileInfo.eventTimes(iii, 1)*samplingRate) / samplingDiff1);
            startTime5 = floor(floor(RestingBaselines.baselineFileInfo.eventTimes(iii, 1)*samplingRate) / samplingDiff5);
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
    trialRestData.([strDay '_' fileID]).OneSec.S_avg = S_trialAvg1;
    trialRestData.([strDay '_' fileID]).FiveSec.S_avg = S_trialAvg5;
end

fields = fieldnames(trialRestData);
uniqueDays = GetUniqueDays(RestingBaselines.baselineFileInfo.fileIDs);

for day = 1:length(uniqueDays)
    x = 1;
    for field = 1:length(fields)
        if strcmp(fields{field}(7:12), uniqueDays{day})
            stringDay = ConvertDate(uniqueDays{day});
            S_avgs.OneSec.(stringDay){x, 1} = trialRestData.(fields{field}).OneSec.S_avg;
            S_avgs.FiveSec.(stringDay){x, 1} = trialRestData.(fields{field}).FiveSec.S_avg;
            x = x + 1;
        end
    end
end

dayFields = fieldnames(S_avgs.OneSec);
for day = 1:length(dayFields)
    dayVals1 = [];
    dayVals5 = [];
    for x = 1:length(S_avgs.OneSec.(dayFields{day}))
        dayVals1 = [dayVals1, S_avgs.OneSec.(dayFields{day}){x, 1}];
        dayVals5 = [dayVals5, S_avgs.FiveSec.(dayFields{day}){x, 1}];
    end
    disp(['Adding spectrogram baseline to baseline file for ' dayFields{day} '...']); disp(' ')
    RestingBaselines.Spectrograms.OneSec.(dayFields{day}) = mean(dayVals1, 2);
    RestingBaselines.Spectrograms.FiveSec.(dayFields{day}) = mean(dayVals5, 2);
end

save([animal '_RestingBaselines.mat'], 'RestingBaselines');

end
