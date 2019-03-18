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

%%
vesselIDs = {};
for a = 1:length(mergedDataFiles)
    mergedDataFile = mergedDataFiles(a,:);
    [~,~,~, vID] = GetFileInfo_2P(mergedDataFile);
    vesselIDs = vID{a,1};
end

uniqueVesselIDs = unique(vesselIDs);
[B, A] = butter(4, 2/(p2Fs/2), 'low');
for b = 1:length(uniqueVesselIDs)
    uniqueVesselID = string(uniqueVesselIDs{b,1});
    d = 1;
    for c = 1:length(mergedDataFiles)
        [~,~,~, mdID] = GetFileInfo_2P(mergedDataFile(c,:));
        if strcmp(uniqueVesselID, mdID) == true
            load(mergedDataFile);
            uniqueVesselData{b,1}(d,:) = detrend(filtfilt(B, A, MergedData.Data.Vessel_Diameter), 'constant');
            d = d + 1;
        end
    end
end

params.tapers = [1 1];
params.pad = 1;
params.Fs = MergedData.Notes.p2Fs;
params.fpass = [0 2]; 
params.trialave = 1;
params.err = [2 0.05];

%%
for e = 1:length(uniqueVesselData)
    [S, f, sErr] = mtspectrumc(uniqueVesselData{e,1}, params);
    allS{e,1} = S;
    allf{e,1} = f;
    allsErr{e,1} = sErr;
end

figure;
for f = 1:length(allS)
    plot(allf{f,1}, allS{f,1});
    hold on
end
title('Power Spectrum of vessel diameter')
xlabel('Frequency (Hz)')
ylabel('Power')
legend(uniqueVesselIDs)

ComparisonData.Vessel_PowerSpec.S = allS;
ComparisonData.Vessel_PowerSpec.f = allf;
ComparisonData.Vessel_PowerSpec.sErr = allsErr;
ComparisonData.Vessel_PowerSpec.vesselIDs = uniqueVesselIDs;

save([animalID '_ComparisonData.mat'], 'ComparisonData');

end