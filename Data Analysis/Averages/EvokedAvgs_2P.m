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
whiskAnimalIDs = {'72', 'T73', 'T74', 'T75', 'T76'};
x = 1;
y = 1;
z = 1;
for t = 1:length(whiskAnimalIDs)
    animalID = whiskAnimalIDs{t};
    cd(['I:\' animalID '\Combined Imaging\']);
    load([animalID '_ComparisonData.mat']);
    
    for a = 1:length(ComparisonData.Whisk.data{1,1})
        whiskData1(x, :) = ComparisonData.Whisk.data{1,1}{a,1};
        vesselID = join([string(animalID) string(ComparisonData.vesselIDs{1,1}{a,1})]);
        vIDs1{x,1} = strrep(vesselID, ' ', '');
        x = 1 + 1;
    end
    
    for b = 1:length(ComparisonData.Whisk.data{2,1})
        whiskData2(y, :) = ComparisonData.Whisk.data{2,1}{b,1};
        vesselID = join([string(animalID) string(ComparisonData.vesselIDs{2,1}{b,1})]);
        vIDs2{y,1} = strrep(vesselID, ' ', '');
        y = 1 + 1;
    end
    
    for c = 1:length(ComparisonData.Whisk.data{3,1})
        whiskData3(z,:) = ComparisonData.Whisk.data{3,1}{c,1};
        vesselID = join([string(animalID) string(ComparisonData.vesselIDs{3,1}{c,1})]);
        vIDs3{z,1} = strrep(vesselID, ' ', '');
        z = 1 + 1;
    end
    
    whiskLFP(:, :, t) = ComparisonData.Whisk.LFP.S;
end
T = ComparisonData.Whisk.LFP.T;
F = ComparisonData.Whisk.LFP.F;

whiskData_C1 = mean(whiskData1, 1);
whiskData_C2 = mean(whiskData2, 1);
whiskData_C3 = mean(whiskData3, 1);
whiskLFP_All = mean(whiskLFP, 3);

whiskSTD_C1 = std(whiskData1, 1, 1);
whiskSTD_C2 = std(whiskData2, 1, 1);
whiskSTD_C3 = std(whiskData3, 1, 1);

%%
timeVec = (1:length(whiskData_C1) - 2)/20;
allAvgs = figure;
subplot(1,3,1)
plot(timeVec, whiskData_C1, 'k')
hold on
plot(timeVec, whiskData_C1 + whiskSTD_C1)
plot(timeVec, whiskData_C1 - whiskSTD_C1)
title('Mean vessel diameter response to whisking events 0.5 to 2 seconds long')
xlabel('Peri-whisk time (sec)')
ylabel('\Delta Diameter (%)')

subplot(1,3,2)
plot(timeVec, whiskData_C2, 'k')
hold on
plot(timeVec, whiskData_C2 + whiskSTD_C2)
plot(timeVec, whiskData_C2 - whiskSTD_C2)
title('Mean vessel diameter response to whisking events 2 to 5 seconds long')
xlabel('Peri-whisk time (sec)')
ylabel('\Delta Diameter (%)')

subplot(1,3,3)
plot(timeVec, whiskData_C3, 'k')
hold on
plot(timeVec, whiskData_C3 + whiskSTD_C3)
plot(timeVec, whiskData_C3 - whiskSTD_C3)
title('Mean vessel diameter response to whisking events 5 to 10 seconds long')
xlabel('Peri-whisk time (sec)')
ylabel('\Delta Diameter (%)')

%%
indAvgs = figure;
subplot(1,3,1)
legendIDs = [];
for x = 1:length(whiskData1)
    plot(timeVec, whiskData1(x,:));
    vID = vIDs1{x,1};
    legendIDs = [legendIDs vID];
    hold on
end
title('Mean vessel diameter response to whisking events 0.5 to 2 seconds long')
xlabel('Peri-whisk time (sec)')
ylabel('\Delta Diameter (%)')
legend(legendIDs)

subplot(1,3,2)
legendIDs = [];
for y = 1:length(whiskData2)
    plot(timeVec, whiskData2(y,:));
    vID = vIDs1{y,1};
    legendIDs = [legendIDs vID];
    hold on
end
title('Mean vessel diameter response to whisking events 2 to 5 seconds long')
xlabel('Peri-whisk time (sec)')
ylabel('\Delta Diameter (%)')
legend(legendIDs)

subplot(1,3,3)
legendIDs = [];
for z = 1:length(whiskData3)
    plot(timeVec, whiskData3(z,:));
    vID = vIDs1{z,1};
    legendIDs = [legendIDs vID];
    hold on
end
title('Mean vessel diameter response to whisking events 5 to 10 seconds long')
xlabel('Peri-whisk time (sec)')
ylabel('\Delta Diameter (%)')
legend(legendIDs)

%%
lfpAvg = figure;
imagesc(T,F,S);
axis xy
colorbar
caxis([-1 1])
title('Mean hippocampal LFP during extended whisking')
ylabel('Frequency (Hz)')
xlabel('Peri-whisk time (sec)')
