function [ComparisonData] = AnalyzeXCorr_2P(animalID, mergedDataFiles, ComparisonData)
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

p2Fs = 20;
downSampledFs = 30;

%%
vesselIDs = {};
for a = 1:size(mergedDataFiles, 1)
    mergedDataFile = mergedDataFiles(a,:);
    [~,~,~, vID] = GetFileInfo_2P(mergedDataFile);
    vesselIDs{a,1} = vID;
end

uniqueVesselIDs = unique(vesselIDs);
[B, A] = butter(4, 2/(p2Fs/2), 'low');
for b = 1:length(uniqueVesselIDs)
    uniqueVesselID = string(uniqueVesselIDs{b,1});
    d = 1;
    for c = 1:size(mergedDataFiles, 1)
        mergedDataFile = mergedDataFiles(c,:);
        [~,~,~, mdID] = GetFileInfo_2P(mergedDataFile);
        if strcmp(uniqueVesselID, mdID) == true
            load(mergedDataFile);
            uniqueVesselData{b,1}(:,d) = detrend(filtfilt(B, A, MergedData.Data.Vessel_Diameter(2:end - 1)), 'constant');
            uniqueWhiskerData{b,1}(:,d) = detrend(abs(diff(resample(MergedData.Data.Whisker_Angle, p2Fs, downSampledFs), 2)), 'constant');
            d = d + 1;
        end
    end
end

for k = 1:length(uniqueVesselIDs)
    uniqueVesselIDs{k,1} = [animalID uniqueVesselIDs{k,1}];
end

%%
z_hold = [];
lagTime = 25;       % Seconds
frequency = 20;     % Hz
maxLag = lagTime*frequency;    % Number of points
for x = 1:length(uniqueVesselIDs)
    z_hold = [];
    for y = 1:size(uniqueVesselData{x, 1}, 2)
        vesselArray = uniqueVesselData{x,1}(:,y);
        whiskArray = uniqueWhiskerData{x,1}(:,y);
        [XC_Vals(y, :), lags] = xcorr(vesselArray, whiskArray, maxLag, 'coeff');
    end
    XC_means{x,1} = mean(XC_Vals, 1);
end

%%
lags = lags/frequency;
figure;
for f = 1:length(XC_means)
    plot(lags, XC_means{f,1});
    hold on
end
title('Cross Correlation between whisker acceleration and vessel diameter')
xlabel('Lags (sec)')
ylabel('Correlation')
legend(uniqueVesselIDs)

ComparisonData.WhiskVessel_XCorr.XC_means = XC_means;
ComparisonData.WhiskVessel_XCorr.lags = lags;
ComparisonData.WhiskVessel_XCorr.vesselIDs = uniqueVesselIDs;

save([animalID '_ComparisonData.mat'], 'ComparisonData');

end