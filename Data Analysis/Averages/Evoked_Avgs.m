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

%% Whisk Data
whiskAnimalIDs = {'T48', 'T49', 'T52', 'T61', 'T62', 'T64', 'T65', 'T66'};
for n = 1:length(whiskAnimalIDs)
    animalID = whiskAnimalIDs{n};
    if strcmp(animalID, 'T48') || strcmp(animalID, 'T49')
        cd(['G:\' animalID '\Combined Imaging\']);
    else
        cd(['I:\' animalID '\Combined Imaging\']);
    end 
        load([animalID '_ComparisonData.mat']);
        LH_whiskCBV(:, n) = ComparisonData.Evoked.Whisk.LH.CBV.data;
        LH_whiskMUA(:, n) = ComparisonData.Evoked.Whisk.LH.MUA.data;
        LH_whiskLFP(:, :, n) = ComparisonData.Evoked.Whisk.LH.LFP.data;
        
        RH_whiskCBV(:, n) = ComparisonData.Evoked.Whisk.RH.CBV.data;
        RH_whiskMUA(:, n) = ComparisonData.Evoked.Whisk.RH.MUA.data;
        RH_whiskLFP(:, :, n) = ComparisonData.Evoked.Whisk.RH.LFP.data;

        T = 1:61;
        F = ComparisonData.Evoked.Whisk.LH.LFP.F;
        timeVector = ComparisonData.Evoked.Whisk.LH.CBV.timeVector;
end

whiskCBV = cat(2, LH_whiskCBV, RH_whiskCBV);
whiskMUA = cat(2, LH_whiskMUA, RH_whiskMUA);
whiskLFP = cat(3, LH_whiskLFP, RH_whiskLFP);

meanWhiskCBV = mean(whiskCBV, 2);
meanWhiskMUA = mean(whiskMUA, 2);
meanWhiskLFP = mean(whiskLFP, 3);

stdWhiskCBV = std(whiskCBV, 0, 2);
stdWhiskMUA = std(whiskMUA, 0, 2);

%% Fig
figure;
ax1 = subplot(3,1,1);
plot(timeVector, meanWhiskMUA)
hold on
plot(timeVector, meanWhiskMUA + stdWhiskMUA)
plot(timeVector, meanWhiskMUA - stdWhiskMUA)
axis tight
ylabel('Normalized Power')
axis square

ax2 = subplot(3,1,2);
imagesc(T, F, meanWhiskLFP)
set(gca, 'Ticklength', [0 0])
ylim([1 100])
ylabel('Freq (Hz)')
colorbar
% caxis([-0.5 1])
axis xy
axis square

ax3 = subplot(3,1,3);
plot(timeVector, meanWhiskCBV*100)
hold on
plot(timeVector, meanWhiskCBV*100 + stdWhiskCBV*100)
plot(timeVector, meanWhiskCBV*100 - stdWhiskCBV*100)
ylabel('Reflectance (%)')
axis tight
title('Stimulus-evoked averages')
xlabel('Peristimulus time (sec)')
axis square


%% Stim Data
stimAnimalIDs = {'T48', 'T49', 'T52', 'T65', 'T66'};
for n = 1:length(stimAnimalIDs)
    animalID = stimAnimalIDs{n};
    if strcmp(animalID, 'T48') || strcmp(animalID, 'T49')
        cd(['G:\' animalID '\Combined Imaging\']);
    else
        cd(['I:\' animalID '\Combined Imaging\']);
    end 
        load([animalID '_ComparisonData.mat']);
        LH_stimCBV(:, n) = ComparisonData.Evoked.Stim.LH.CBV.data;
        LH_stimMUA(:, n) = ComparisonData.Evoked.Stim.LH.MUA.data;
        LH_stimLFP(:, :, n) = ComparisonData.Evoked.Stim.LH.LFP.data;
        
        RH_stimCBV(:, n) = ComparisonData.Evoked.Stim.RH.CBV.data;
        RH_stimMUA(:, n) = ComparisonData.Evoked.Stim.RH.MUA.data;
        RH_stimLFP(:, :, n) = ComparisonData.Evoked.Stim.RH.LFP.data;
end

stimCBV = cat(2, LH_stimCBV, RH_stimCBV);
stimMUA = cat(2, LH_stimMUA, RH_stimMUA);
stimLFP = cat(3, LH_stimLFP, RH_stimLFP);

meanStimCBV = mean(stimCBV, 2);
meanStimMUA = mean(stimMUA, 2);
meanStimLFP = mean(stimLFP, 3);

stdStimCBV = std(stimCBV, 0, 2);
stdStimMUA = std(stimMUA, 0, 2);

%% Fig
figure;
ax1 = subplot(3,1,1);
plot(timeVector, meanStimMUA)
hold on
plot(timeVector, meanStimMUA + stdStimMUA)
plot(timeVector, meanStimMUA - stdStimMUA)
axis tight
ylabel('Normalized Power')
axis square

ax2 = subplot(3,1,2);
imagesc(T, F, meanStimLFP)
set(gca, 'Ticklength', [0 0])
ylim([1 100])
ylabel('Freq (Hz)')
colorbar
caxis([-1 1])
axis xy
axis square

ax3 = subplot(3,1,3);
plot(timeVector, meanStimCBV*100)
hold on
plot(timeVector, meanStimCBV*100 + stdStimCBV*100)
plot(timeVector, meanStimCBV*100 - stdStimCBV*100)
ylabel('Reflectance (%)')
axis tight
title('Stimulus-evoked averages')
xlabel('Peristimulus time (sec)')
axis square

