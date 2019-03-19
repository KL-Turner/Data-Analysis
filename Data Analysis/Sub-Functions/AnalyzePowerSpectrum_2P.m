function [ComparisonData] = AnalyzePowerSpectrum_2P(animalID, mergedDataFiles, ComparisonData)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs: 
%
%   Last Revised: March 18th, 2019
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
t = 1;
for b = 1:length(uniqueVesselIDs)
    uniqueVesselID = string(uniqueVesselIDs{b,1});
    d = 1;
    for c = 1:size(mergedDataFiles, 1)
        mergedDataFile = mergedDataFiles(c,:);
        [~,~,~, mdID] = GetFileInfo_2P(mergedDataFile);
        if strcmp(uniqueVesselID, mdID) == true
            load(mergedDataFile);
            uniqueVesselData{b,1}(:,d) = detrend(filtfilt(B, A, MergedData.Data.Vessel_Diameter), 'constant');
            whiskerData(:,t) = detrend(abs(diff(resample(MergedData.Data.Whisker_Angle, p2Fs, downSampledFs), 2)), 'constant');
            d = d + 1;
            t = t + 1;
        end
    end
end

for k = 1:length(uniqueVesselIDs)
    uniqueVesselIDs{k,1} = [animalID uniqueVesselIDs{k,1}];
end

params.tapers = [3 5];
params.pad = 1;
params.Fs = p2Fs;
params.fpass = [0 1]; 
params.trialave = 1;
params.err = [2 0.05];

%%
for e = 1:length(uniqueVesselData)
    [S, f, sErr] = mtspectrumc(uniqueVesselData{e,1}, params);
    allS{e,1} = S;
    allf{e,1} = f;
    allsErr{e,1} = sErr;
end

[wS, wf, wsErr] = mtspectrumc(whiskerData, params);

figure;
for f = 1:length(allS)
    plot(allf{f,1}, allS{f,1});
    hold on
end
title('Power Spectrum of vessel diameter')
xlabel('Frequency (Hz)')
ylabel('Power')
legend(uniqueVesselIDs)

figure;
plot(wf, wS, 'k')
title('Power Spectrum of abs(whiskerAccel)')
xlabel('Frequency (Hz)')
ylabel('Power')

ComparisonData.Vessel_PowerSpec.S = allS;
ComparisonData.Vessel_PowerSpec.f = allf;
ComparisonData.Vessel_PowerSpec.sErr = allsErr;
ComparisonData.Vessel_PowerSpec.vesselIDs = uniqueVesselIDs;
ComparisonData.Whisk_PowerSpec.S = wS;
ComparisonData.Whisk_PowerSpec.f = wf;
ComparisonData.Whisk_PowerSpec.sErr = wsErr;
% ComparisonData.Whisk_PowerSpec.tblVals = tblVals;

save([animalID '_ComparisonData.mat'], 'ComparisonData');

end