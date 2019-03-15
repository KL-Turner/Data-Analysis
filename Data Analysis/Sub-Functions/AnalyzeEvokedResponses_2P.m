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
processedWhiskData = cell(length(whiskData), 1);
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
                filtVesselDiam = (filtfilt(B, A, normVesselDiam))*100;
                processedWhiskData{x,1}{y,1}(w,:) = filtVesselDiam;
                w = w + 1;
            end
        end
    end
end

whiskCritMeans = cell(length(processedWhiskData),1);
whiskCritSTD = cell(length(processedWhiskData),1);
for x = 1:length(processedWhiskData)
    for y = 1:length(processedWhiskData{x,1})
        whiskCritMeans{x,1}{y,1} = mean(processedWhiskData{x,1}{y,1},1);
        whiskCritSTD{x,1}{y,1} = std(processedWhiskData{x,1}{y,1},1, 1);
    end
end

for x = 1:length(whiskCritMeans)
    figure('NumberTitle', 'off', 'Name', ['Whisker Criteria ' num2str(x)]);
    for y = 1:length(whiskCritMeans{x,1})
        plot(((1:length(whiskCritMeans{x,1}{y,1}))/p2Fs)-2, whiskCritMeans{x,1}{y,1} - mean(whiskCritMeans{x,1}{y,1}(1:40)));
        hold on
    end
    title('Whisking evoked vessel diameter')
    xlabel('Peri-whisk time (sec)')
    ylabel('Diameter change {%}')
    legend('T76A1', 'T76A2', 'T76A3')
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
        whiskS_Vals = whiskS_Data(:, (startTime - (2*whiskBinSize)):(startTime + ((whiskDuration - 2)*whiskBinSize)));
    catch
        whiskS_Vals = whiskS_Data(:, end - (whiskDuration*whiskBinSize):end);
    end
    
    % Mean subtract each row with detrend
    transpWhiskS_Vals = whiskS_Vals';   % Transpose since detrend goes down columns
    dTWhiskS_Vals = detrend(transpWhiskS_Vals);   % detrend
    whiskZhold = cat(3, whiskZhold, dTWhiskS_Vals');   % transpose back to original orientation and save into Data.S
end
% 
T = 1:whiskDuration;
timevec = (T/whiskBinSize)-2;
F = SpectrogramData.OneSec.F{1, 1};
whiskS = mean(whiskZhold, 3);

figure;
imagesc(timevec,F,whiskS);
axis xy
caxis([-1 2])
ylabel(' ')

S = SpectrogramData.OneSec.S_Norm{3,1};
T = SpectrogramData.OneSec.T{1,1};
F = SpectrogramData.OneSec.F{1,1};
figure;
imagesc(T,F,S)
axis xy
caxis([-1 2])

end

%%
% whiskSpecs = cell(length(whiskData), 1);
% for x = 1:length(whiskData)
%     uniqueVesselIDs = unique(whiskVesselIDs{x,1});
%     for y = 1:length(uniqueVesselIDs)
%         uniqueVesselID = uniqueVesselIDs{y,1};
%         whiskHold = [];
%         for z = 1:length(whiskVesselIDs{x,1})
%             vesselID = whiskVesselIDs{x,1}{z,1};
%             if strcmp(uniqueVesselID, vesselID) == 1
%                 fileID = whiskFileIDs{x,1}{z,1};
%                 for v = 1:length(SpectrogramData.FileIDs)
%                     specFileID = (SpectrogramData.FileIDs{v,1});
%                     if strcmp(specFileID, fileID)
%                         whiskSData = SpectrogramData.OneSec.S_Norm{v,1};
%                     end
%                 end
%                 
%                 whiskSLength = size(whiskSData, 2);
%                 whiskBinSize = ceil(whiskSLength/280);
%                 whiskSamplingDiff = p2Fs/whiskBinSize;
%                 whiskDuration = (floor(whiskSLength/280))*12;
%                 startTime = floor(floor(whiskEventTimes{x,1}(z,1)*p2Fs)/whiskSamplingDiff);
%                 if startTime == 0
%                     startTime = 1;
%                 end
%              
%                 try
%                     whiskSVals = whiskSData(:, (startTime - (2*whiskBinSize)):(startTime + ((whiskDuration - 2)*whiskBinSize)));
%                 catch
%                     whiskSVals = whiskSData(:, end - (whiskDuration*whiskBinSize):end);
%                 end
%                 transpWhiskSVals = whiskSVals';
%                 dTWhiskSVals = detrend(transpWhiskSVals);
%                 whiskHold = cat(3, whiskHold, dTWhiskSVals');
%             end
%         end
%         whiskSpecs{x,1}{y,1} = mean(whiskHold, 3);
%     end
% end
% 
% T = 1:whiskDuration;
% F = SpectrogramData.OneSec.F{1, 1};
% 
% % ax3 =subplot(3,1,3);
% % plot(timeVector, meanWhiskMUAData)
% % axis tight
% % xlabel('Peristimulus time (sec)')
% % ylabel('Normalized Power')
% % 
% % ComparisonData.Evoked.Whisk.(dataType).CBV.data = meanWhiskCBVData;
% % ComparisonData.Evoked.Whisk.(dataType).CBV.timeVector = timeVector;
% % ComparisonData.Evoked.Whisk.(dataType).MUA.data = meanWhiskMUAData;
% % ComparisonData.Evoked.Whisk.(dataType).MUA.timeVector = timeVector;
% % ComparisonData.Evoked.Whisk.(dataType).LFP.data = whiskS;
% % ComparisonData.Evoked.Whisk.(dataType).LFP.T = T;
% % ComparisonData.Evoked.Whisk.(dataType).LFP.F = F;
% 
% save([animal '_ComparisonData.mat'], 'ComparisonData');
% 
% [pathstr, ~, ~] = fileparts(cd);
% dirpath = [pathstr '/Figures/Evoked Averages/'];
% 
% if ~exist(dirpath, 'dir')
%     mkdir(dirpath);
% end
% 
% savefig(whiskEvoked, [dirpath animal '_' dataType '_WhiskEvokedAverages']);
% 
% end
% 
% 
% 
