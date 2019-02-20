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

clear
clc

%% Rest and All data
restAnimalIDs = {'T48', 'T49', 'T52', 'T61', 'T62', 'T64', 'T65', 'T66'};
for n = 1:length(restAnimalIDs)
    animalID = restAnimalIDs{n};
    if strcmp(animalID, 'T48') || strcmp(animalID, 'T49')
        cd(['G:\' animalID '\Combined Imaging\']);
    else
        cd(['I:\' animalID '\Combined Imaging\']);
    end 
        load([animalID '_ComparisonData.mat']);
        restCBV_mean(1,1) = ComparisonData.AvgCBV_STD.Rest.LH.mean;
        restCBV_std(1,1) = ComparisonData.AvgCBV_STD.Rest.LH.std;
        restCBV_mean(1,2) = ComparisonData.AvgCBV_STD.Rest.RH.mean;
        restCBV_std(1,2) = ComparisonData.AvgCBV_STD.Rest.RH.std;    
        restCBV_means(n, 1) = mean(restCBV_mean);
        restCBV_stds(n, 1) = mean(restCBV_std);
        
        allCBV_mean(1,1) = ComparisonData.AvgCBV_STD.AllData.LH.mean;
        allCBV_std(1,1) = ComparisonData.AvgCBV_STD.AllData.LH.std;
        allCBV_mean(1,2) = ComparisonData.AvgCBV_STD.AllData.RH.mean;
        allCBV_std(1,2) = ComparisonData.AvgCBV_STD.AllData.RH.std;
        allCBV_means(n, 1) = mean(allCBV_mean);
        allCBV_stds(n, 1) = mean(allCBV_std);
end

figure;
errorbar(1:length(restCBV_means), restCBV_means*100, restCBV_stds*100, 'o', 'LineStyle', 'none')
title('Rest CBV')
ylim([-15 5])
axis square

figure;
errorbar(1:length(allCBV_means), allCBV_means*100, allCBV_stds*100, 's', 'LineStyle', 'none')
title('All Data CBV')
ylim([-15 5])
axis square

%% NREM data
nremAnimalIDs = {'T48', 'T49', 'T52', 'T61', 'T62', 'T64'};
for n = 1:length(nremAnimalIDs)
    animalID = nremAnimalIDs{n};
    if strcmp(animalID, 'T48') || strcmp(animalID, 'T49')
        cd(['G:\' animalID '\Combined Imaging\']);
    else
        cd(['I:\' animalID '\Combined Imaging\']);
    end 
        load([animalID '_ComparisonData.mat']);
        nremCBV_mean(1,1) = ComparisonData.AvgCBV_STD.NREM.LH.mean;
        nremCBV_std(1,1) = ComparisonData.AvgCBV_STD.NREM.LH.std;
        nremCBV_mean(1,2) = ComparisonData.AvgCBV_STD.NREM.RH.mean;
        nremCBV_std(1,2) = ComparisonData.AvgCBV_STD.NREM.RH.std;    
        nremCBV_means(n, 1) = mean(nremCBV_mean);
        nremCBV_stds(n, 1) = mean(nremCBV_std);
end

figure;
errorbar(1:length(nremCBV_means), nremCBV_means*100, nremCBV_stds*100, 'h', 'LineStyle', 'none')
title('NREM CBV')
ylim([-15 5])
axis square

%% REM data
remAnimalIDs = {'T48', 'T52', 'T61', 'T62', 'T64'};
for n = 1:length(remAnimalIDs)
    animalID = remAnimalIDs{n};
    if strcmp(animalID, 'T48') || strcmp(animalID, 'T49')
        cd(['G:\' animalID '\Combined Imaging\']);
    else
        cd(['I:\' animalID '\Combined Imaging\']);
    end 
        load([animalID '_ComparisonData.mat']);
        remCBV_mean(1,1) = ComparisonData.AvgCBV_STD.REM.LH.mean;
        remCBV_std(1,1) = ComparisonData.AvgCBV_STD.REM.LH.std;
        remCBV_mean(1,2) = ComparisonData.AvgCBV_STD.REM.RH.mean;
        remCBV_std(1,2) = ComparisonData.AvgCBV_STD.REM.RH.std;    
        remCBV_means(n, 1) = mean(remCBV_mean);
        remCBV_stds(n, 1) = mean(remCBV_std);
end

figure;
errorbar(1:length(remCBV_means), remCBV_means*100, remCBV_stds*100, 'd', 'LineStyle', 'none')
title('REM CBV')
ylim([-15 5])
axis square

[a b] = find(restCBV_means==0);
restCBV_means(a,b) = NaN
anova1(restCBV_means)

