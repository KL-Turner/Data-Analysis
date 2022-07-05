function [AnalysisResults] = PupilHbTRelationship(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Generate figure panel 8 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

% colorBlack = [(0/256),(0/256),(0/256)];
% colorGrey = [(209/256),(211/256),(212/256)];
% colorRfcAwake = [(0/256),(64/256),(64/256)];
% colorRfcNREM = [(0/256),(174/256),(239/256)];
% colorRfcREM = [(190/256),(30/256),(45/256)];
colorRest = [(0/256),(166/256),(81/256)];
% colorWhisk = [(31/256),(120/256),(179/256)];
% colorStim = [(255/256),(28/256),(206/256)];
colorNREM = [(191/256),(0/256),(255/256)];
colorREM = [(254/256),(139/256),(0/256)];
colorAlert = [(255/256),(191/256),(0/256)];
colorAsleep = [(0/256),(128/256),(255/256)];
colorAll = [(183/256),(115/256),(51/256)];
% colorIso = [(0/256),(256/256),(256/256)];
%% variables for loops
resultsStruct = 'Results_PupilHbTRelationship';
load(resultsStruct);
animalIDs = fieldnames(Results_PupilHbTRelationship);
behavFields = {'Awake','NREM','REM'};
%% take data from each animal corresponding to the CBV-gamma relationship
catHbT = []; catPupil = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        if isfield(catHbT,behavField) == false
            catHbT.(behavField) = []; 
            catPupil.(behavField) = [];
        end
        catHbT.(behavField) = cat(1,catHbT.(behavField),Results_PupilHbTRelationship.(animalID).(behavField).HbT);
        catPupil.(behavField) = cat(1,catPupil.(behavField),Results_PupilHbTRelationship.(animalID).(behavField).Pupil);
    end
end
%% Fig. 8
% summaryFigure = figure('Name','Fig8 (a-d)');
% sgtitle('Figure 8 - Turner et al. 2020')
%% [8d] histogram images for HbT vs. gamma-band power during different arousal-states
% histogram for Awake
awakeHist = figure;
h1 = histogram2(catPupil.Awake,catHbT.Awake,'DisplayStyle','tile','ShowEmptyBins','on','XBinedges',-5:0.025:3,'YBinedges',-25:2.5:125,'Normalization','probability');
h1Vals = h1.Values;
% RGB image for Awake
awakeRGB = figure;
s = pcolor(-4.975:0.025:3,-22.5:2.5:125,h1Vals');
s.FaceColor = 'interp';
set(s,'EdgeColor','none');
n = 50;         
R = linspace(1,0,n);
B = linspace(1,0,n);
G = linspace(1,0,n); 
colormap(flipud([R(:),G(:),B(:)]));
h1Frame = getframe(gcf);
h1Img = frame2im(h1Frame);
close(awakeHist)
close(awakeRGB)
summaryFigure = figure;
subplot(1,3,1)
imshow(h1Img)
axis off
title('Awake Pupil-[HbT]')
xlabel('Diameter(z-units)')
ylabel('\Delta[HbT]')
% histogram for NREM
nremHist = figure;
h2 = histogram2(catPupil.NREM,catHbT.NREM,'DisplayStyle','tile','ShowEmptyBins','on','XBinedges',-5:0.025:3,'YBinedges',-25:2.5:125,'Normalization','probability');
h2Vals = h2.Values;
% RGB image for NREM
nremRGB = figure;
s = pcolor(-4.975:0.025:3,-22.5:2.5:125,h2Vals');
s.FaceColor = 'interp';
set(s,'EdgeColor','none');
n = 50;         
R = linspace(0,0,n);
B = linspace(1,0,n);
G = linspace(1,0,n); 
colormap(flipud([R(:),G(:),B(:)]));
h2Frame = getframe(gcf);
h2Img = frame2im(h2Frame);
close(nremHist)
close(nremRGB)
subplot(1,3,2)
imshow(h2Img)
axis off
title('NREM Pupil-[HbT]')
xlabel('Diameter(z-units)')
ylabel('\Delta[HbT]')
% histogram for REM
remHist = figure;
h3 = histogram2(catPupil.REM,catHbT.REM,'DisplayStyle','tile','ShowEmptyBins','on','XBinedges',-5:0.025:3,'YBinedges',-25:2.5:125,'Normalization','probability');
h3Vals = h3.Values;
% RGB image for REM
remRGB = figure;
s = pcolor(-4.975:0.025:3,-22.5:2.5:125,h3Vals');
s.FaceColor = 'interp';
set(s,'EdgeColor','none');
n = 50;         
R = linspace(1,0,n);
B = linspace(0,0,n);
G = linspace(0,0,n); 
colormap(flipud([R(:),G(:),B(:)]));
h3Frame = getframe(gcf);
h3Img = frame2im(h3Frame);
close(remHist)
close(remRGB)
subplot(1,3,3)
imshow(h3Img)
axis off
title('REM Pupil-[HbT]')
xlabel('Diameter(z-units)')
ylabel('\Delta[HbT]')
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Pupil-HbT Relationship' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    imwrite(h1Img,[dirpath 'Fig8d_Awake.png'])
    imwrite(h2Img,[dirpath 'Fig8d_NREM.png'])
    imwrite(h3Img,[dirpath 'Fig8d_REM.png'])
    savefig(summaryFigure,[dirpath 'PupilHbT_Relationship']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'PupilHbT_Relationship'])
end

end
