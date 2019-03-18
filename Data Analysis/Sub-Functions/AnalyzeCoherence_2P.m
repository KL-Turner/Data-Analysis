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

%%
vesselIDs = {};
for a = 1:length(mergedDataFiles)
    mergedDataFile = mergedDataFiles(a,:);
    [~,~,~, vID] = GetFileInfo_2P(mergedDataFile);
    vesselIDs = vID{a,1};
end

uniqueVesselIDs = unique(vesselIDs);
[B, A] = butter(4, 2/(20/2), 'low');
for b = 1:length(uniqueVesselIDs)
    uniqueVesselID = string(uniqueVesselIDs{b,1});
    d = 1;
    for c = 1:length(mergedDataFiles)
        [~,~,~, mdID] = GetFileInfo_2P(mergedDataFile(c,:));
        if strcmp(uniqueVesselID, mdID) == true
            load(mergedDataFile);
            uniqueVesselData{b,1}(d,:) = detrend(filtfilt(B, A, MergedData.Data.Vessel_Diameter), 'constant');
            uniqueWhiskerData{b,1}(d,:) = abs(diff(MergedData.Data.Whisker_Angle, 2));
            d = d + 1;
        end
    end
end

%%
params.tapers = [3 5];
params.pad = 1;
params.Fs = 20; 
params.fpass = [0 2]; 
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
