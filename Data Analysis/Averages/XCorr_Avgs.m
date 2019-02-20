%_______________________________________________________________________________________________
% Written by Kevin L. Turner, Jr. 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%_______________________________________________________________________________________________
%
%   Purpose: 
%_______________________________________________________________________________________________
%
%   Inputs: 
%
%   Outputs: 
%________________________________________________________________________________________________

%% Rest and All data
restAnimalIDs = {'T48', 'T49', 'T52', 'T61', 'T62', 'T64', 'T65', 'T66'};
for n = 1:length(restAnimalIDs)
    animalID = restAnimalIDs{n};
    if strcmp(animalID, 'T48' || 'T49')
        cd(['G:\' animalID '\Combined Imaging\']);
    else
        cd(['I:\' animalID '\Combined Imaging\']);
    end 
        load([animalID '_AnalysisAverages.mat']);
        LH_restXC_Vals(:, :, n) = ComparisonData.Xcorr.Rest.LH.XC_Vals;
        LH_restLags(:, :, n) = ComparisonData.XCorr.Rest.LH.lags;
        LH_restF(:, :, n) = ComparisonData.XCorr.Rest.LH.F;

        RH_restXC_Vals(:, :, n) = ComparisonData.Xcorr.Rest.RH.XC_Vals;
        RH_restLags(:, :, n) = ComparisonData.XCorr.Rest.RH.lags;
        RH_restF(:, :, n) = ComparisonData.XCorr.Rest.RH.F;
        
        LH_allXC_Vals(:, :, n) = ComparisonData.XCorr.AllData.LH.XC_Vals;
        LH_allLags(:, :, n) = ComparisonData.XCorr.AllData.LH.lags;
        LH_allF(:, :, n) = ComparisonData.XCorr.AllData.LH.F;
        
        RH_allXC_Vals(:, :, n) = ComparisonData.XCorr.AllData.RH.XC_Vals;
        RH_allLags(:, :, n) = ComparisonData.XCorr.AllData.RH.lags;
        RH_allF(:, :, n) = ComparisonData.XCorr.AllData.RH.F;
end

restXC_Vals = concat(LH_restXC_Vals,RH_restXC_Vals);
restLags = concat(LH_restLags, RH_restLags);
restF = concat(LH_restF, RH_restF);

allXC_Vals = concat(LH_allXC_Vals,RH_allXC_Vals);
allLags = concat(LH_allLags, RH_allLags);
allF = concat(LH_allF, RH_allF);

meanRestXC_Vals = mean(restXC_Vals, 3);
meanRestLags = mean(restLags, 2);
meanRestF = mean(restF, 2);

meanAllXC_Vals = mean(allXC_Vals, 3);
meanAllLags = mean(allLags, 2);
meanAllF = mean(allF, 2);







































































figure('name', 'First ten seconds of sleep');
subplot(1,2,1)
imagesc(meanRestLags, meanRestF, meanRestMeanXC_Vals);
colormap parula
colorbar
title('Rest XCorr (Combined)')
ylabel('Frequency (Hz)')
xlabel('Time (sec)')
xticks([-maxLag2 -maxLag2/2 0 maxLag2/2 maxLag2])
xticklabels({'-5', '-2.5', '0', '2.5' '5'})
axis square
axis xy
ylim([1 100])

subplot(1,2,2)
imagesc(meanSleepLags, meanSleepF, meanSleepMeanXC_Vals);
colormap parula
colorbar
title('Sleep XCorr (Combined)')
ylabel('Frequency (Hz)')
xlabel('Time (sec)')
xticks([-maxLag -maxLag/2 0 maxLag/2 maxLag])
xticklabels({'-10', '-5', '0', '5' '10'})
axis square
axis xy
ylim([1 100])
