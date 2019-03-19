%___________________________________________________________________________________________________
% Written by Kevin L. Turner, Jr.
% Adapted from codes credited to Dr. Patrick J. Drew and Aaron T. Winder
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%___________________________________________________________________________________________________
%
%   Purpose:
%___________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%___________________________________________________________________________________________________

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
title('Abs(whiskAccel) vs. vessel diameter')
xlabel('Frequency (Hz)')
ylabel('Coherence')

ax2 = subplot(1,2,2);
for c = 1:size(coherenceData, 1)
    plot(f, coherenceData(c,:));
    hold on
end
title('Abs(whiskAccel) vs. vessel diameter')
xlabel('Frequency (Hz)')
ylabel('Coherence')
legend(vIDs)
linkaxes([ax1 ax2], 'xy')
