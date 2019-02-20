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
        restCBV_C(:, n) = ComparisonData.Coherence.Rest.CBV.C;
        restCBV_f(:, n) = ComparisonData.Coherence.Rest.CBV.f;
        restGAM_C(:, n) = ComparisonData.Coherence.Rest.GAM.C;
        restGAM_f(:, n) = ComparisonData.Coherence.Rest.GAM.f;
        
        allCBV_C(:, n) = ComparisonData.Coherence.AllData.CBV.C;
        allCBV_f(:, n) = ComparisonData.Coherence.AllData.CBV.f;
        allGAM_C(:, n) = ComparisonData.Coherence.AllData.GAM.C;
        allGAM_f(:, n) = ComparisonData.Coherence.AllData.GAM.f;
end

meanRestCBV_C = mean(restCBV_C, 2);
meanRestCBV_f = mean(restCBV_f, 2);
meanRestGAM_C = mean(restGAM_C, 2);
meanRestGAM_f = mean(restGAM_f, 2);

stdRestCBV_C = (std(restCBV_C, 0, 2));
stdRestGAM_C = (std(restGAM_C, 0, 2));

meanAllCBV_C = mean(allCBV_C, 2);
meanAllCBV_f = mean(allCBV_f, 2);
meanAllGAM_C = mean(allGAM_C, 2);
meanAllGAM_f = mean(allGAM_f, 2);

stdAllCBV_C = (std(allCBV_C, 0, 2));
stdAllGAM_C = (std(allGAM_C, 0, 2));

%% NREM Sleep
nremAnimalIDs = {'T48', 'T49', 'T52', 'T61', 'T62', 'T64'};
for n = 1:length(nremAnimalIDs)
    animalID = nremAnimalIDs{n};
    if strcmp(animalID, 'T48') || strcmp(animalID, 'T49')
        cd(['G:\' animalID '\Combined Imaging\']);
    else
        cd(['I:\' animalID '\Combined Imaging\']);
    end 
        load([animalID '_ComparisonData.mat']);
        nremCBV_C(:, n) = ComparisonData.Coherence.NREM.CBV.C;
        nremCBV_f(:, n) = ComparisonData.Coherence.NREM.CBV.f; 
        nremGAM_C(:, n) = ComparisonData.Coherence.NREM.GAM.C;
        nremGAM_f(:, n) = ComparisonData.Coherence.NREM.GAM.f;
end

meanNREMCBV_C = mean(nremCBV_C, 2);
meanNREMCBV_f = mean(nremCBV_f, 2);
meanNREMGAM_C = mean(nremGAM_C, 2);
meanNREMGAM_f = mean(nremGAM_f, 2);

stdNREMCBV_C = (std(nremCBV_C, 0, 2));
stdNREMGAM_C = (std(nremGAM_C, 0, 2));

%% REM Sleep
remAnimalIDs = {'T48', 'T52', 'T61', 'T62', 'T64'};
for n = 1:length(remAnimalIDs)
    animalID = remAnimalIDs{n};
    if strcmp(animalID, 'T48') || strcmp(animalID, 'T49')
        cd(['G:\' animalID '\Combined Imaging\']);
    else
        cd(['I:\' animalID '\Combined Imaging\']);
    end 
        load([animalID '_ComparisonData.mat']);
        remCBV_C(:, n) = ComparisonData.Coherence.REM.CBV.C;
        remCBV_f(:, n) = ComparisonData.Coherence.REM.CBV.f; 
        remGAM_C(:, n) = ComparisonData.Coherence.REM.GAM.C;
        remGAM_f(:, n) = ComparisonData.Coherence.REM.GAM.f;
end

meanREMCBV_C = mean(remCBV_C, 2);
meanREMCBV_f = mean(remCBV_f, 2);
meanREMGAM_C = mean(remGAM_C, 2);
meanREMGAM_f = mean(remGAM_f, 2);

stdREMCBV_C = (std(remCBV_C, 0, 2));
stdREMGAM_C = (std(remGAM_C, 0, 2));

%% Summary Fig
figure;
plot(meanRestCBV_f, meanRestCBV_C, 'r');
hold on
plot(meanRestCBV_f, meanRestCBV_C + stdRestCBV_C, '--r');
plot(meanRestCBV_f, meanRestCBV_C - stdRestCBV_C, '--r');

plot(meanAllCBV_f, meanAllCBV_C, 'k');
plot(meanAllCBV_f, meanAllCBV_C + stdAllCBV_C, '--k');
plot(meanAllCBV_f, meanAllCBV_C - stdAllCBV_C, '--k');

plot(meanNREMCBV_f, meanNREMCBV_C, 'g');
plot(meanNREMCBV_f, meanNREMCBV_C + stdNREMCBV_C, '--g');
plot(meanNREMCBV_f, meanNREMCBV_C - stdNREMCBV_C, '--g');

plot(meanREMCBV_f, meanREMCBV_C, 'b');
plot(meanREMCBV_f, meanREMCBV_C + stdREMCBV_C, '--b');
plot(meanREMCBV_f, meanREMCBV_C - stdREMCBV_C, '--b');

title('L/R CBV Coherence')
ylabel('Magnitude of Coherence')
xlabel('Frequency (Hz)')
legend('Rest Mean', 'upper Rest std', 'lower Rest std', 'All std', 'upper All std', 'lower All std', 'NREM Mean',...
    'upper NREM std', 'lower NREM std', 'REM Mean', 'upper REM std', 'lower REM std')
xlim([0 1])
ylim([0 1])
axis square

figure;
plot(meanRestGAM_f, meanRestGAM_C, 'r');
hold on
plot(meanRestGAM_f, meanRestGAM_C + stdRestGAM_C, '--r');
plot(meanRestGAM_f, meanRestGAM_C - stdRestGAM_C, '--r');

plot(meanAllGAM_f, meanAllGAM_C, 'k');
plot(meanAllGAM_f, meanAllGAM_C + stdAllGAM_C, '--k');
plot(meanAllGAM_f, meanAllGAM_C - stdAllGAM_C, '--k');

plot(meanNREMGAM_f, meanNREMGAM_C, 'g');
plot(meanNREMGAM_f, meanNREMGAM_C + stdNREMGAM_C, '--g');
plot(meanNREMGAM_f, meanNREMGAM_C - stdNREMGAM_C, '--g');

plot(meanREMGAM_f, meanREMGAM_C, 'b');
plot(meanREMGAM_f, meanREMGAM_C + stdREMGAM_C, '--b');
plot(meanREMGAM_f, meanREMGAM_C - stdREMGAM_C, '--b');

title('L/R GAM Coherence')
ylabel('Magnitude of Coherence')
xlabel('Frequency (Hz)')
legend('Rest Mean', 'upper Rest std', 'lower Rest std', 'All std', 'upper All std', 'lower All std', 'NREM Mean',...
    'upper NREM std', 'lower NREM std', 'REM Mean', 'upper REM std', 'lower REM std')
xlim([0 1])
ylim([0 1])
axis square
