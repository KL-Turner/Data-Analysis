function [ComparisonData] = AnalyzeCoherence_2P(animalID, mergedDataFiles, ComparisonData)
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
for b = 1:length(uniqueVesselIDs)
    uniqueVesselID = string(uniqueVesselIDs{b,1});
    d = 1;
    for c = 1:size(mergedDataFiles, 1)
        mergedDataFile = mergedDataFiles(c,:);
        [~,~,~, mdID] = GetFileInfo_2P(mergedDataFile);
        if strcmp(uniqueVesselID, mdID) == true
            load(mergedDataFile);
            uniqueVesselData{b,1}(:,d) = detrend(filtfilt(B, A, MergedData.Data.Vessel_Diameter(1:end - 2)), 'constant');
            uniqueWhiskerData{b,1}(:,d) = detrend(abs(diff(resample(MergedData.Data.Whisker_Angle, p2Fs, downSampledFs), 2)), 'constant');
            d = d + 1;
        end
    end
end

for k = 1:length(uniqueVesselIDs)
    uniqueVesselIDs{k,1} = [animalID uniqueVesselIDs{k,1}];
end

%%
params.tapers = [3 5];
params.pad = 1;
params.Fs = p2Fs; 
params.fpass = [0 1]; 
params.trialave = 1;
params.err = [2 0.05];

for e = 1:length(uniqueVesselData)
    [C, ~, ~, ~, ~, f, ~, ~, cErr] = coherencyc(uniqueVesselData{e,1}, uniqueWhiskerData{e,1}, params);
    allC{e,1} = C;
    allf{e,1} = f;
    allcErr{e,1} = cErr;
end

%%
figure;
for f = 1:length(allC)
    plot(allf{f,1}, allC{f,1});
    hold on
end
title('Coherence bewteen whisking and vessel diameter')
xlabel('Frequency (Hz)')
ylabel('Coherence')
legend(uniqueVesselIDs)

ComparisonData.WhiskVessel_Coherence.C = allC;
ComparisonData.WhiskVessel_Coherence.f = allf;
ComparisonData.WhiskVessel_Coherence.cErr = allcErr;
ComparisonData.WhiskVessel_Coherence.vesselIDs = uniqueVesselIDs;

save([animalID '_ComparisonData.mat'], 'ComparisonData');

end
