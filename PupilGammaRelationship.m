function [AnalysisResults] = PupilGammaRelationship(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Generate figure panel 8 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% variables for loops
resultsStruct = 'Results_PupilGammaRelationship';
load(resultsStruct);
animalIDs = fieldnames(Results_PupilGammaRelationship);
behavFields = {'Awake','NREM','REM'};
%% take data from each animal corresponding to the CBV-gamma relationship
catGamma = []; catPupil = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        if isfield(catGamma,behavField) == false
            catGamma.(behavField) = []; 
            catPupil.(behavField) = [];
        end
        catGamma.(behavField) = cat(1,catGamma.(behavField),Results_PupilGammaRelationship.(animalID).(behavField).Gamma);
        catPupil.(behavField) = cat(1,catPupil.(behavField),Results_PupilGammaRelationship.(animalID).(behavField).Pupil);
    end
end
%% Fig. 8
% summaryFigure = figure('Name','Fig8 (a-d)');
% sgtitle('Figure 8 - Turner et al. 2020')
%% [8d] histogram images for Gamma vs. gamma-band power during different arousal-states
% histogram for Awake
awakeHist = figure;
h1 = histogram2(catPupil.Awake,catGamma.Awake,'DisplayStyle','tile','ShowEmptyBins','on','XBinedges',-5:0.025:3,'YBinedges',-0.25:0.025:1,'Normalization','probability');
h1Vals = h1.Values;
% RGB image for Awake
awakeRGB = figure;
s = pcolor(-4.975:0.025:3,-0.225:0.025:1,h1Vals');
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
title('Awake Pupil-Gamma')
xlabel('Diameter(z-units)')
ylabel('Gamma-band power (%)') 
% histogram for NREM
nremHist = figure;
h2 = histogram2(catPupil.NREM,catGamma.NREM,'DisplayStyle','tile','ShowEmptyBins','on','XBinedges',-5:0.025:3,'YBinedges',-0.25:0.025:1,'Normalization','probability');
h2Vals = h2.Values;
% RGB image for NREM
nremRGB = figure;
s = pcolor(-4.975:0.025:3,-0.225:0.025:1,h2Vals');
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
title('NREM Pupil-Gamma')
xlabel('Diameter(z-units)')
ylabel('Gamma-band power (%)') 
% histogram for REM
remHist = figure;
h3 = histogram2(catPupil.REM,catGamma.REM,'DisplayStyle','tile','ShowEmptyBins','on','XBinedges',-5:0.025:3,'YBinedges',-0.25:0.025:1,'Normalization','probability');
h3Vals = h3.Values;
% RGB image for REM
remRGB = figure;
s = pcolor(-4.975:0.025:3,-0.225:0.025:1,h3Vals');
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
title('REM Pupil-Gamma')
xlabel('Diameter(z-units)')
ylabel('Gamma-band power (%)') 
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Pupil-Gamma Relationship' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    imwrite(h1Img,[dirpath 'Fig8d_Awake.png'])
    imwrite(h2Img,[dirpath 'Fig8d_NREM.png'])
    imwrite(h3Img,[dirpath 'Fig8d_REM.png'])
    savefig(summaryFigure,[dirpath 'PupilGamma_Relationship']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'PupilGamma_Relationship'])
end

end
