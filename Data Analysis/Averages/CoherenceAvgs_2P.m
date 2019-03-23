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
%   Last Revised: March 22nd, 2019
%________________________________________________________________________________________________________________________

%%
whiskAnimalIDs = {'T72', 'T73', 'T74', 'T75', 'T76'};
x = 1;
for a = 1:length(whiskAnimalIDs)
    animalID = whiskAnimalIDs{a};
    cd(['I:\' animalID '\Combined Imaging\']);
    load([animalID '_ComparisonData.mat']);
  
    for b = 1:length(ComparisonData.WhiskVessel_Coherence.C)
        coherenceData(x, :) = ComparisonData.WhiskVessel_Coherence.C{b,1};
        vIDs{x,1} =  ComparisonData.Vessel_PowerSpec.vesselIDs{b,1};
        x = x + 1;
    end 
end

f = ComparisonData.WhiskVessel_Coherence.f{1,1};
coherenceMean = mean(coherenceData, 1);
coherenceSTD = std(coherenceData, 1, 1);

%%
cohAvgs = figure;
ax1 = subplot(1,2,1);
plot(f, coherenceMean, 'k')
hold on
plot(f, coherenceMean + coherenceSTD)
plot(f, coherenceMean - coherenceSTD)
title('Mean coherence Abs(whiskAccel) vs. vessel diameter')
xlabel('Frequency (Hz)')
ylabel('Coherence')
xlim([0 0.5])

ax2 = subplot(1,2,2);
for c = 1:size(coherenceData, 1)
    plot(f, coherenceData(c,:));
    hold on
end
title('Ind coherence Abs(whiskAccel) vs. vessel diameter')
xlabel('Frequency (Hz)')
ylabel('Coherence')
legend(vIDs)
xlim([0 0.5])
linkaxes([ax1 ax2], 'xy')
