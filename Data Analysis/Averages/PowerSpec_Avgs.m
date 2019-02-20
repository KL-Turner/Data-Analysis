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
        LH_restCBV_S(:, n) = ComparisonData.PowerSpectrum.Rest.CBV.LH.S;
        LH_restCBV_f(:, n) = ComparisonData.PowerSpectrum.Rest.CBV.LH.f;
        LH_restGAM_S(:, n) = ComparisonData.PowerSpectrum.Rest.GAM.LH.S;
        LH_restGAM_f(:, n) = ComparisonData.PowerSpectrum.Rest.GAM.LH.f;
        
        RH_restCBV_S(:, n) = ComparisonData.PowerSpectrum.Rest.CBV.RH.S;
        RH_restCBV_f(:, n) = ComparisonData.PowerSpectrum.Rest.CBV.RH.f;
        RH_restGAM_S(:, n) = ComparisonData.PowerSpectrum.Rest.GAM.RH.S;
        RH_restGAM_f(:, n) = ComparisonData.PowerSpectrum.Rest.GAM.RH.f;
        
        LH_allCBV_S(:, n) = ComparisonData.PowerSpectrum.AllData.CBV.LH.S;
        LH_allCBV_f(:, n) = ComparisonData.PowerSpectrum.AllData.CBV.LH.f;
        LH_allGAM_S(:, n) = ComparisonData.PowerSpectrum.AllData.GAM.LH.S;
        LH_allGAM_f(:, n) = ComparisonData.PowerSpectrum.AllData.GAM.LH.f;
          
        RH_allCBV_S(:, n) = ComparisonData.PowerSpectrum.AllData.CBV.RH.S;
        RH_allCBV_f(:, n) = ComparisonData.PowerSpectrum.AllData.CBV.RH.f;
        RH_allGAM_S(:, n) = ComparisonData.PowerSpectrum.AllData.GAM.RH.S;
        RH_allGAM_f(:, n) = ComparisonData.PowerSpectrum.AllData.GAM.RH.f;
end

restCBV_S = cat(2, LH_restCBV_S, RH_restCBV_S);
restCBV_f = cat(2, LH_restCBV_f, RH_restCBV_f);
restGAM_S = cat(2, LH_restGAM_S, RH_restGAM_S);
restGAM_f = cat(2, LH_restGAM_f, RH_restGAM_f);

allCBV_S = cat(2, LH_allCBV_S, RH_allCBV_S);
allCBV_f = cat(2, LH_allCBV_f, RH_allCBV_f);
allGAM_S = cat(2, LH_allGAM_S, RH_allGAM_S);
allGAM_f = cat(2, LH_allGAM_f, RH_allGAM_f);

meanRestCBV_S = mean(restCBV_S, 2);
meanRestCBV_f = mean(restCBV_f, 2);
meanRestGAM_S = mean(restGAM_S, 2);
meanRestGAM_f = mean(restGAM_f, 2);

stdRestCBV_S = (std(restCBV_S, 0, 2));
stdRestGAM_S = (std(restGAM_S, 0, 2));

meanAllCBV_S = mean(allCBV_S, 2);
meanAllCBV_f = mean(allCBV_f, 2);
meanAllGAM_S = mean(allGAM_S, 2);
meanAllGAM_f = mean(allGAM_f, 2);

stdAllCBV_S = (std(allCBV_S, 0, 2));
stdAllGAM_S = (std(allGAM_S, 0, 2));

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
        LH_nremCBV_S(:, n) = ComparisonData.PowerSpectrum.NREM.CBV.LH.S;
        LH_nremCBV_f(:, n) = ComparisonData.PowerSpectrum.NREM.CBV.LH.f; 
        LH_nremGAM_S(:, n) = ComparisonData.PowerSpectrum.NREM.GAM.LH.S;
        LH_nremGAM_f(:, n) = ComparisonData.PowerSpectrum.NREM.GAM.LH.f;
        
        RH_nremCBV_S(:, n) = ComparisonData.PowerSpectrum.NREM.CBV.RH.S;
        RH_nremCBV_f(:, n) = ComparisonData.PowerSpectrum.NREM.CBV.RH.f; 
        RH_nremGAM_S(:, n) = ComparisonData.PowerSpectrum.NREM.GAM.RH.S;
        RH_nremGAM_f(:, n) = ComparisonData.PowerSpectrum.NREM.GAM.RH.f;
end

nremCBV_S = cat(2, LH_nremCBV_S, RH_nremCBV_S);
nremCBV_f = cat(2, LH_nremCBV_f, RH_nremCBV_f);
nremGAM_S = cat(2, LH_nremGAM_S, RH_nremGAM_S);
nremGAM_f = cat(2, LH_nremGAM_f, RH_nremGAM_f);

meanNREMCBV_S = mean(nremCBV_S, 2);
meanNREMCBV_f = mean(nremCBV_f, 2);
meanNREMGAM_S = mean(nremGAM_S, 2);
meanNREMGAM_f = mean(nremGAM_f, 2);

stdNREMCBV_S = (std(nremCBV_S, 0, 2));
stdNREMGAM_S = (std(nremGAM_S, 0, 2));

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
        LH_remCBV_S(:, n) = ComparisonData.PowerSpectrum.REM.CBV.LH.S;
        LH_remCBV_f(:, n) = ComparisonData.PowerSpectrum.REM.CBV.LH.f; 
        LH_remGAM_S(:, n) = ComparisonData.PowerSpectrum.REM.GAM.LH.S;
        LH_remGAM_f(:, n) = ComparisonData.PowerSpectrum.REM.GAM.LH.f;
        
        RH_remCBV_S(:, n) = ComparisonData.PowerSpectrum.REM.CBV.RH.S;
        RH_remCBV_f(:, n) = ComparisonData.PowerSpectrum.REM.CBV.RH.f; 
        RH_remGAM_S(:, n) = ComparisonData.PowerSpectrum.REM.GAM.RH.S;
        RH_remGAM_f(:, n) = ComparisonData.PowerSpectrum.REM.GAM.RH.f;
end

remCBV_S = cat(2, LH_remCBV_S, RH_remCBV_S);
remCBV_f = cat(2, LH_remCBV_f, RH_remCBV_f);
remGAM_S = cat(2, LH_remGAM_S, RH_remGAM_S);
remGAM_f = cat(2, LH_remGAM_f, RH_remGAM_f);

meanREMCBV_S = mean(remCBV_S, 2);
meanREMCBV_f = mean(remCBV_f, 2);
meanREMGAM_S = mean(remGAM_S, 2);
meanREMGAM_f = mean(remGAM_f, 2);

stdREMCBV_S = (std(remCBV_S, 0, 2));
stdREMGAM_S = (std(remGAM_S, 0, 2));

%% Summary Fig
figure;
semilogx(meanRestCBV_f, meanRestCBV_S, 'r');
hold on
loglog(meanRestCBV_f, meanRestCBV_S + stdRestCBV_S, '--r');
loglog(meanRestCBV_f, meanRestCBV_S - stdRestCBV_S, '--r');

semilogx(meanAllCBV_f, meanAllCBV_S, 'k');
loglog(meanAllCBV_f, meanAllCBV_S + stdAllCBV_S, '--k');
loglog(meanAllCBV_f, meanAllCBV_S - stdAllCBV_S, '--k');

semilogx(meanNREMCBV_f, meanNREMCBV_S, 'g');
loglog(meanNREMCBV_f, meanNREMCBV_S + stdNREMCBV_S, '--g');
loglog(meanNREMCBV_f, meanNREMCBV_S - stdNREMCBV_S, '--g');

semilogx(meanREMCBV_f, meanREMCBV_S, 'b');
loglog(meanREMCBV_f, meanREMCBV_S + stdREMCBV_S, '--b');
loglog(meanREMCBV_f, meanREMCBV_S - stdREMCBV_S, '--b');

title('L/R CBV Power Spec')
ylabel('Magnitude of PowerSpectrum')
xlabel('Frequency (Hz)')
legend('Rest Mean', 'upper Rest std', 'lower Rest std', 'All Mean', 'upper All std', 'lower All std', 'NREM Mean',...
    'upper NREM std', 'lower NREM std', 'REM Mean', 'upper REM std', 'lower REM std')
xlim([0.0293 1])
axis square

figure;
loglog(meanRestGAM_f, meanRestGAM_S, 'r');
hold on
loglog(meanRestGAM_f, meanRestGAM_S + stdRestGAM_S, '--r');
loglog(meanRestGAM_f, meanRestGAM_S - stdRestGAM_S, '--r');

loglog(meanAllGAM_f, meanAllGAM_S, 'k');
loglog(meanAllGAM_f, meanAllGAM_S + stdAllGAM_S, '--k');
loglog(meanAllGAM_f, meanAllGAM_S - stdAllGAM_S, '--k');

loglog(meanNREMGAM_f, meanNREMGAM_S, 'g');
loglog(meanNREMGAM_f, meanNREMGAM_S + stdNREMGAM_S, '--g');
loglog(meanNREMGAM_f, meanNREMGAM_S - stdNREMGAM_S, '--g');

loglog(meanREMGAM_f, meanREMGAM_S, 'b');
loglog(meanREMGAM_f, meanREMGAM_S + stdREMGAM_S, '--b');
loglog(meanREMGAM_f, meanREMGAM_S - stdREMGAM_S, '--b');

title('L/R GAM Power Spec')
ylabel('Magnitude of PowerSpectrum')
xlabel('Frequency (Hz)')
legend('Rest Mean', 'upper Rest std', 'lower Rest std', 'All std', 'upper All std', 'lower All std', 'NREM Mean',...
    'upper NREM std', 'lower NREM std', 'REM Mean', 'upper REM std', 'lower REM std')
xlim([0.0293 1])
axis square
