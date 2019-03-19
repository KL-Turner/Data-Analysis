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
  
    for b = 1:length(ComparisonData.Vessel_PowerSpec.S)
        powerspecVesselData(x, :) = ComparisonData.Vessel_PowerSpec.S{b,1};
        powerspecWhiskData(x, :) = ComparisonData.Whisk_PowerSpec.S;
        vIDs{x,1} =  ComparisonData.Vessel_PowerSpec.vesselIDs{b,1};
        x = x + 1;
    end 
    powerspecWhiskData(a, :) = ComparisonData.Whisk_PowerSpec.S;

end

vf = ComparisonData.Vessel_PowerSpec.f{1,1};
powerspecVesselMean = mean(powerspecVesselData, 1);
powerspecVesselSTD = std(powerspecVesselData, 1, 1);

wf = ComparisonData.Whisk_PowerSpec.f;
powerspecWhiskMean = mean(powerspecWhiskData, 1);
powerspecWhiskSTD = std(powerspecWhiskData, 1, 1);

%%
specAvgs = figure;
ax1 = subplot(2,2,1);
plot(vf, powerspecVesselMean, 'k')
hold on
plot(vf, powerspecVesselMean + powerspecVesselSTD)
plot(vf, powerspecVesselMean - powerspecVesselSTD)
title('Power spec vessel diameter')
xlabel('Frequency (Hz)')
ylabel('Power')

ax2 = subplot(2,2,2);
for c = 1:size(powerspecVesselData, 1)
    plot(vf, powerspecVesselData(c,:));
    hold on
end
title('Power spec vessel diameter')
xlabel('Frequency (Hz)')
ylabel('Power')
legend(vIDs)
linkaxes([ax1 ax2], 'xy')

ax3 = subplot(2,2,3);
plot(wf, powerspecWhiskMean, 'k')
hold on
plot(wf, powerspecWhiskMean + powerspecWhiskSTD)
plot(wf, powerspecWhiskMean - powerspecWhiskSTD)
title('Power spec abs(whiskerAccel)')
xlabel('Frequency (Hz)')
ylabel('Power')

ax4 = subplot(2,2,4);
for c = 1:size(powerspecWhiskData, 1)
    plot(wf, powerspecWhiskData(c,:));
    hold on
end
title('Power spec abs(whiskerAccel)')
xlabel('Frequency (Hz)')
ylabel('Power')
legend(whiskAnimalIDs)
linkaxes([ax3 ax4], 'xy')
