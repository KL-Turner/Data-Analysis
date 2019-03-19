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
y = 1;
z = 1;
for t = 1:length(whiskAnimalIDs)
    animalID = whiskAnimalIDs{t};
    cd(['I:\' animalID '\Combined Imaging\']);
    load([animalID '_ComparisonData.mat']);
    
    for a = 1:length(ComparisonData.Whisk.data{1,1})
        whiskData1(x, :) = ComparisonData.Whisk.data{1,1}{a,1};
        vesselID = join([string(animalID) string(ComparisonData.Whisk.vesselIDs{1,1}{a,1})]);
        vIDs1{x,1} = strrep(vesselID, ' ', '');
        x = x + 1;
    end
    
    for b = 1:length(ComparisonData.Whisk.data{2,1})
        whiskData2(y, :) = ComparisonData.Whisk.data{2,1}{b,1};
        vesselID = join([string(animalID) string(ComparisonData.Whisk.vesselIDs{2,1}{b,1})]);
        vIDs2{y,1} = strrep(vesselID, ' ', '');
        y = y + 1;
    end
    
    for c = 1:length(ComparisonData.Whisk.data{3,1})
        whiskData3(z,:) = ComparisonData.Whisk.data{3,1}{c,1};
        vesselID = join([string(animalID) string(ComparisonData.Whisk.vesselIDs{3,1}{c,1})]);
        vIDs3{z,1} = strrep(vesselID, ' ', '');
        z = z + 1;
    end
    
    whiskLFP1(:, :, t) = ComparisonData.Whisk.LFP.S{1,1};
    whiskLFP2(:, :, t) = ComparisonData.Whisk.LFP.S{2,1};
    whiskLFP3(:, :, t) = ComparisonData.Whisk.LFP.S{3,1};
    
end
T = ComparisonData.Whisk.LFP.T;
F = ComparisonData.Whisk.LFP.F;

whiskData_C1 = mean(whiskData1, 1);
whiskData_C2 = mean(whiskData2, 1);
whiskData_C3 = mean(whiskData3, 1);
whiskLFP_C1 = mean(whiskLFP1, 3);
whiskLFP_C2 = mean(whiskLFP2, 3);
whiskLFP_C3 = mean(whiskLFP3, 3);

whiskSTD_C1 = std(whiskData1, 1, 1);
whiskSTD_C2 = std(whiskData2, 1, 1);
whiskSTD_C3 = std(whiskData3, 1, 1);

T = ComparisonData.Whisk.LFP.T;
F = ComparisonData.Whisk.LFP.F;

%%
timeVec = ((1:length(whiskData_C1))/20) - 2;
evokedAvgs = figure;
ax1 = subplot(3,3,1);
legendIDs = [];
for x = 1:size(whiskData1, 1)
    plot(timeVec, whiskData1(x,:));
    vID = vIDs1{x,1};
    hold on
end
title('0.5 to 2 seconds')
xlabel('Peri-whisk time (sec)')
ylabel('\Delta Diameter (%)')

ax2 = subplot(3,3,2);
legendIDs = [];
for y = 1:size(whiskData2, 1)
    plot(timeVec, whiskData2(y,:));
    vID = vIDs1{y,1};
    hold on
end
title('2 to 5 seconds')
xlabel('Peri-whisk time (sec)')
ylabel('\Delta Diameter (%)')

ax3 = subplot(3,3,3);
legendIDs = [];
for z = 1:size(whiskData3, 1)
    plot(timeVec, whiskData3(z,:));
    vID = vIDs1{z,1};
    legendIDs = [legendIDs vID];
    hold on
end
title('5 to 10')
xlabel('Peri-whisk time (sec)')
ylabel('\Delta Diameter (%)')
legend(legendIDs)
linkaxes([ax1 ax2 ax3], 'xy')

%%
ax4 = subplot(3,3,4);
plot(timeVec, whiskData_C1, 'k')
hold on
plot(timeVec, whiskData_C1 + whiskSTD_C1)
plot(timeVec, whiskData_C1 - whiskSTD_C1)
xlabel('Peri-whisk time (sec)')
ylabel('\Delta Diameter (%)')

ax5 = subplot(3,3,5);
plot(timeVec, whiskData_C2, 'k')
hold on
plot(timeVec, whiskData_C2 + whiskSTD_C2)
plot(timeVec, whiskData_C2 - whiskSTD_C2)
xlabel('Peri-whisk time (sec)')
ylabel('\Delta Diameter (%)')

ax6 = subplot(3,3,6);
plot(timeVec, whiskData_C3, 'k')
hold on
plot(timeVec, whiskData_C3 + whiskSTD_C3)
plot(timeVec, whiskData_C3 - whiskSTD_C3)
xlabel('Peri-whisk time (sec)')
ylabel('\Delta Diameter (%)')
linkaxes([ax4 ax5 ax6], 'xy')

%%
subplot(3,3,7);
imagesc(T,F,whiskLFP_C1);
axis xy
caxis([-1 1])
ylabel('Frequency (Hz)')
xlabel('Peri-whisk time (sec)')

subplot(3,3,8);
imagesc(T,F,whiskLFP_C2);
axis xy
caxis([-1 1])
ylabel('Frequency (Hz)')
xlabel('Peri-whisk time (sec)')

subplot(3,3,9);
imagesc(T,F,whiskLFP_C3);
axis xy
caxis([-1 1])
ylabel('Frequency (Hz)')
xlabel('Peri-whisk time (sec)')
