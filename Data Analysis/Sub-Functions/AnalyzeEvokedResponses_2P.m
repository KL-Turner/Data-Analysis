function [ComparisonData] = AnalyzeEvokedResponses_2P(animalID, RestingBaselines, EventData, SpectrogramData)
%___________________________________________________________________________________________________
% Written by Kevin L. Turner, Jr.
% Adapted from codes credited to Dr. Patrick J. Drew and Aaron T. Winder
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%___________________________________________________________________________________________________
%
%   Purpose:
%___________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%___________________________________________________________________________________________________

p2Fs = 20;

%%
whiskCriteria.Fieldname{1,1} = {'duration', 'duration', 'puffDistance'};
whiskCriteria.Comparison{1,1} = {'gt','lt','gt'};
whiskCriteria.Value{1,1} = {0.5, 2, 5};

whiskCriteria.Fieldname{2,1} = {'duration', 'duration', 'puffDistance'};
whiskCriteria.Comparison{2,1} = {'gt','lt','gt'};
whiskCriteria.Value{2,1} = {2, 5, 5};

whiskCriteria.Fieldname{3,1} = {'duration', 'duration', 'puffDistance'};
whiskCriteria.Comparison{3,1} = {'gt', 'lt', 'gt'};
whiskCriteria.Value{3,1} = {5, 10, 5};

whiskData = cell(length(whiskCriteria.Fieldname), 1);
whiskVesselIDs = cell(length(whiskCriteria.Fieldname), 1);
whiskEventTimes = cell(length(whiskCriteria.Fieldname), 1);
whiskFileIDs = cell(length(whiskCriteria.Fieldname), 1);
for x = 1:length(whiskCriteria.Fieldname)
    criteria.Fieldname = whiskCriteria.Fieldname{x,1};
    criteria.Comparison = whiskCriteria.Comparison{x,1};
    criteria.Value = whiskCriteria.Value{x,1};
    whiskFilter = FilterEvents(EventData.Vessel_Diameter.whisk, criteria);
    [tempWhiskData] = EventData.Vessel_Diameter.whisk.data(whiskFilter, :);
    whiskData{x,1} = tempWhiskData;
    [tempWhiskVesselIDs] = EventData.Vessel_Diameter.whisk.vesselIDs(whiskFilter, :);
    whiskVesselIDs{x,1} = tempWhiskVesselIDs;
    [tempWhiskEventTimes] = EventData.Vessel_Diameter.whisk.eventTime(whiskFilter, :);
    whiskEventTimes{x,1} = tempWhiskEventTimes;
    [tempWhiskFileIDs] = EventData.Vessel_Diameter.whisk.fileIDs(whiskFilter, :);
    whiskFileIDs{x,1} = tempWhiskFileIDs;
end

%%
[B, A] = butter(4, 2/(p2Fs/2), 'low');
processedWhiskData.data = cell(length(whiskData), 1);
for x = 1:length(whiskData)
    uniqueVesselIDs = unique(whiskVesselIDs{x,1});
    for y = 1:length(uniqueVesselIDs)
        uniqueVesselID = uniqueVesselIDs{y,1};
        w = 1;
        for z = 1:length(whiskVesselIDs{x,1})
            vesselID = whiskVesselIDs{x,1}{z,1};
            if strcmp(uniqueVesselID, vesselID) == 1
                fileID = whiskFileIDs{x,1}{z,1};
                strDay = ConvertDate(fileID(1:6));
                vesselDiam = whiskData{x,1}(z,:);
                normVesselDiam = (vesselDiam - RestingBaselines.(uniqueVesselID).(strDay).Vessel_Diameter.baseLine)./(RestingBaselines.(vesselID).(strDay).Vessel_Diameter.baseLine);
                filtVesselDiam = sgolayfilt(normVesselDiam, 10, 101)*100;
%                 filtVesselDiam = (filtfilt(B, A, normVesselDiam))*100;
                processedWhiskData.data{x,1}{y,1}(w,:) = filtVesselDiam;
                processedWhiskData.vesselIDs{x,1}{y,1}{w} = vesselID;
                w = w + 1;
            end
        end
    end
end

whiskCritMeans.data = cell(length(processedWhiskData.data),1);
whiskCritSTD = cell(length(processedWhiskData.data),1);
for x = 1:length(processedWhiskData.data)
    for y = 1:length(processedWhiskData.data{x,1})
        whiskCritMeans.data{x,1}{y,1} = mean(processedWhiskData.data{x,1}{y,1},1);
        whiskCritMeans.vesselIDs{x,1}{y,1} = unique(processedWhiskData.vesselIDs{x,1}{y,1});
        whiskCritSTD{x,1}{y,1} = std(processedWhiskData.data{x,1}{y,1},1, 1);
    end
end

for x = 1:length(whiskCritMeans.data)
    legendIDs = [];
    figure('NumberTitle', 'off', 'Name', ['Whisker Criteria ' num2str(x)]);
    for y = 1:length(whiskCritMeans.data{x,1})
        plot(((1:length(whiskCritMeans.data{x,1}{y,1}))/p2Fs)-2, whiskCritMeans.data{x,1}{y,1} - mean(whiskCritMeans.data{x,1}{y,1}(1:40)));
        vID = join([string(animalID) string(whiskCritMeans.vesselIDs{x,1}{y,1})]);
        vID = strrep(vID, ' ', '');
        legendIDs = [legendIDs vID];
        hold on
    end
    title('Whisking evoked vessel diameter')
    xlabel('Peri-whisk time (sec)')
    ylabel('Diameter change {%}')
    legend(legendIDs)
end

%%
whiskZhold = [];
sFiles = whiskFileIDs{3,1};
sEventTimes = whiskEventTimes{3,1};
for x = 1:length(sFiles)   % Loop through each non-unique file
    whiskFileID = sFiles{x, 1}; 
    % Load in Neural Data from rest period
    for s = 1:length(SpectrogramData.FileIDs)
        if strcmp(whiskFileID, SpectrogramData.FileIDs{s, 1})
            whiskS_Data = SpectrogramData.OneSec.S_Norm{s, 1};  % S data for this specific file
        end
    end
    whiskSLength = size(whiskS_Data, 2);                                  % Length of the data across time (number of samples)
    whiskBinSize = ceil(whiskSLength / 280);                              % Find the number of bins needed to divide this into 300 seconds
    whiskSamplingDiff = p2Fs/whiskBinSize;                   % Number to divide by for data to be at 5 Hz
    
    % Find the start time and duration
    whiskDuration = whiskBinSize*12;
    startTime = floor(floor(sEventTimes(x,1)*p2Fs)/whiskSamplingDiff);
    if startTime == 0
        startTime = 1;
    end
    
    % Take the S_data from the start time throughout the duration
    try
        whiskS_Vals = whiskS_Data(:, (startTime - (2*whiskBinSize)):(startTime + ((12 - 2)*whiskBinSize)));
    catch
        whiskS_Vals = whiskS_Data(:, end - (12*whiskBinSize):end);
    end
    
    whiskZhold = cat(3, whiskZhold, whiskS_Vals);
end

T = 1:whiskDuration;
timevec = (T/whiskBinSize)-2;
F = SpectrogramData.OneSec.F{1, 1};
whiskS = mean(whiskZhold, 3);

figure;
imagesc(timevec,F,whiskS);
axis xy
colorbar
caxis([-1 1])
title('Hippocampal LFP during extended whisking')
ylabel('Frequency (Hz)')
xlabel('Peri-whisk time (sec)')

ComparisonData.Whisk.data = whiskCritMeans.data;
ComparisonData.Whisk.vesselIDs = whiskCritMeans.vesselIDs;
ComparisonData.Whisk.std = whiskCritSTD;
ComparisonData.Whisk.LFP.S = whiskS;
ComparisonData.Whisk.LFP.T = timevec;
ComparisonData.Whisk.LFP.F = F;

save([animalID '_ComparisonData.mat'], 'ComparisonData');

end
