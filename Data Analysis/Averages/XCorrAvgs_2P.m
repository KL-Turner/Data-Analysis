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
        XCorrData(x, :) = ComparisonData.WhiskVessel_XCorr.XC_means{b,1};
        vIDs{x,1} =  ComparisonData.WhiskVessel_XCorr.vesselIDs{b,1};
        x = x + 1;
    end 
end

lags = ComparisonData.WhiskVessel_XCorr.lags;
xcorrMean = mean(XCorrData, 1);
xcorrSTD = std(XCorrData, 1, 1);

%%
xcorrAvgs = figure;
ax1 = subplot(1,2,1);
plot(lags, xcorrMean, 'k')
hold on
plot(lags, xcorrMean + xcorrSTD)
plot(lags, xcorrMean - xcorrSTD)
title('Mean XCorr Abs(whiskAccel) vs. vessel diameter')
ylabel('Correlation')
xlabel('Lags (sec)')

ax2 = subplot(1,2,2);
for c = 1:size(XCorrData, 1)
    plot(lags, XCorrData(c,:));
    hold on
end
title('Ind XCorr Abs(whiskAccel) vs. vessel diameter')
ylabel('Correlation')
xlabel('Lags (sec)')
legend(vIDs)
linkaxes([ax1 ax2], 'xy')
