function [AnalysisResults] = Fig6_TBD(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Generate figure panel 8 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% Pupil-HbT relationship
resultsStruct = 'Results_PupilHbTRelationship';
load(resultsStruct);
animalIDs = fieldnames(Results_PupilHbTRelationship);
behavFields = {'Awake','NREM','REM'};
% take data from each animal corresponding to the CBV-gamma relationship
 data.HbTRel.catHbT = [];  data.HbTRel.catPupil = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        if isfield(data.HbTRel.catHbT,behavField) == false
             data.HbTRel.catHbT.(behavField) = []; 
             data.HbTRel.catPupil.(behavField) = [];
        end
        data.HbTRel.catHbT.(behavField) = cat(1,data.HbTRel.catHbT.(behavField),Results_PupilHbTRelationship.(animalID).(behavField).HbT);
        data.HbTRel.catPupil.(behavField) = cat(1,data.HbTRel.catPupil.(behavField),Results_PupilHbTRelationship.(animalID).(behavField).Pupil);
    end
end
%% Pupil-Gamma relationship
resultsStruct = 'Results_PupilGammaRelationship';
load(resultsStruct);
animalIDs = fieldnames(Results_PupilGammaRelationship);
behavFields = {'Awake','NREM','REM'};
% take data from each animal corresponding to the CBV-gamma relationship
data.GammaRel.catGamma = [];  data.GammaRel.catPupil = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        if isfield(data.GammaRel.catGamma,behavField) == false
            data.GammaRel.catGamma.(behavField) = [];
            data.GammaRel.catPupil.(behavField) = [];
        end
        data.GammaRel.catGamma.(behavField) = cat(1,data.GammaRel.catGamma.(behavField),Results_PupilGammaRelationship.(animalID).(behavField).Gamma*100);
        data.GammaRel.catPupil.(behavField) = cat(1,data.GammaRel.catPupil.(behavField),Results_PupilGammaRelationship.(animalID).(behavField).Pupil);
    end
end
%% Sleep probability based on pupil mm diameter
resultsStruct = 'Results_SleepProbability';
load(resultsStruct);
diameterAllCatMeans = Results_SleepProbability.diameterCatMeans;
awakeProbPerc = Results_SleepProbability.awakeProbPerc;
nremProbPerc = Results_SleepProbability.nremProbPerc;
remProbPerc = Results_SleepProbability.remProbPerc;
asleepProbPerc = Results_SleepProbability.asleepProbPerc;
%% Sleep model accuracy based on pupil zDiameter alone
resultsStructB = 'Results_PupilSleepModel';
load(resultsStructB);
animalIDs = fieldnames(Results_PupilSleepModel);
data.pupil.holdXlabels = []; data.pupil.holdYlabels = [];
for dd = 1:length(animalIDs)
    animalID = animalIDs{dd,1};
    data.pupil.loss(dd,1) = Results_PupilSleepModel.(animalID).SVM.loss;
    data.pupil.holdXlabels = cat(1,data.pupil.holdXlabels,Results_PupilSleepModel.(animalID).SVM.testXlabels);
    data.pupil.holdYlabels = cat(1,data.pupil.holdYlabels,Results_PupilSleepModel.(animalID).SVM.testYlabels);
end
data.pupil.meanLoss = mean(data.pupil.loss,1);
data.pupil.stdLoss = std(data.pupil.loss,0,1);
%% Sleep model accuracy based on physiology 
resultsStructB = 'Results_PhysioSleepModel';
load(resultsStructB);
animalIDs = fieldnames(Results_PhysioSleepModel);
data.physio.holdXlabels = []; data.physio.holdYlabels = [];
for dd = 1:length(animalIDs)
    animalID = animalIDs{dd,1};
    data.physio.loss(dd,1) = Results_PhysioSleepModel.(animalID).SVM.loss;
    data.physio.holdXlabels = cat(1,data.physio.holdXlabels,Results_PhysioSleepModel.(animalID).SVM.testXlabels);
    data.physio.holdYlabels = cat(1,data.physio.holdYlabels,Results_PhysioSleepModel.(animalID).SVM.testYlabels);
end
data.physio.meanLoss = mean(data.physio.loss,1);
data.physio.stdLoss = std(data.physio.loss,0,1);
%% pupil model coherence
resultsStruct = 'Results_PupilModelCoherence.mat';
load(resultsStruct);
animalIDs = fieldnames(Results_PupilModelCoherence);
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    data.Coherr.pupilf(aa,:) = Results_PupilModelCoherence.(animalID).Pupil.f;
    data.Coherr.pupilC(aa,:) = Results_PupilModelCoherence.(animalID).Pupil.C;
    data.Coherr.physiof(aa,:) = Results_PupilModelCoherence.(animalID).Physio.f;
    data.Coherr.physioC(aa,:) = Results_PupilModelCoherence.(animalID).Physio.C;
end
data.Coherr.meanPupilf = mean(data.Coherr.pupilf,1);
data.Coherr.meanPupilC = mean(data.Coherr.pupilC,1);
data.Coherr.stdPupilC = std(data.Coherr.pupilC,0,1)/sqrt(size(data.Coherr.pupilC,1));
data.Coherr.meanPhysiof = mean(data.Coherr.physiof,1);
data.Coherr.meanPhysioC = mean(data.Coherr.physioC,1);
data.Coherr.stdPhysioC = std(data.Coherr.physioC,0,1)/sqrt(size(data.Coherr.physioC,1));
%% Figure
HbTawakeHist = figure;
h1 = histogram2(data.HbTRel.catPupil.Awake,data.HbTRel.catHbT.Awake,'DisplayStyle','tile','ShowEmptyBins','on','XBinedges',-5:0.025:3,'YBinedges',-25:2.5:125,'Normalization','probability');
h1Vals = h1.Values;
% RGB image for Awake
HbTawakeRGB = figure;
s = pcolor(-4.975:0.025:3,-22.5:2.5:125,h1Vals');
s.FaceColor = 'interp';
set(s,'EdgeColor','none');
n = 50;         
R = linspace(0,1,n);
G = linspace(0,1,n); 
B = linspace(0,1,n);
colormap(flipud([R(:),G(:),B(:)]));
cax = caxis;
caxis([cax(1),cax(2)/1.5])
axis off
h1Frame = getframe(gcf);
h1Img = frame2im(h1Frame);
close(HbTawakeHist)
close(HbTawakeRGB)
% histogram for NREM
HbTnremHist = figure;
h2 = histogram2(data.HbTRel.catPupil.NREM,data.HbTRel.catHbT.NREM,'DisplayStyle','tile','ShowEmptyBins','on','XBinedges',-5:0.025:3,'YBinedges',-25:2.5:125,'Normalization','probability');
h2Vals = h2.Values;
% RGB image for NREM
HbTnremRGB = figure;
s = pcolor(-4.975:0.025:3,-22.5:2.5:125,h2Vals');
s.FaceColor = 'interp';
set(s,'EdgeColor','none');
n = 50;         
R = linspace(0,1,n);
G = linspace(0.4,1,n); 
B = linspace(0,1,n);
colormap(flipud([R(:),G(:),B(:)]));
cax = caxis;
caxis([cax(1),cax(2)/1.5])
axis off
h2Frame = getframe(gcf);
h2Img = frame2im(h2Frame);
close(HbTnremHist)
close(HbTnremRGB)
% histogram for REM
HbTremHist = figure;
h3 = histogram2(data.HbTRel.catPupil.REM,data.HbTRel.catHbT.REM,'DisplayStyle','tile','ShowEmptyBins','on','XBinedges',-5:0.025:3,'YBinedges',-25:2.5:125,'Normalization','probability');
h3Vals = h3.Values;
% RGB image for REM
HbTRemRGB = figure;
s = pcolor(-4.975:0.025:3,-22.5:2.5:125,h3Vals');
s.FaceColor = 'interp';
set(s,'EdgeColor','none');
n = 50;         
R = linspace(1,1,n);
G = linspace(0,1,n); 
B = linspace(1,1,n);
colormap(flipud([R(:),G(:),B(:)]));
cax = caxis;
caxis([cax(1),cax(2)/1.5])
axis off
h3Frame = getframe(gcf);
h3Img = frame2im(h3Frame);
close(HbTremHist)
close(HbTRemRGB)
GammaAwakeHist = figure;
h4 = histogram2(data.HbTRel.catPupil.Awake,data.GammaRel.catGamma.Awake,'DisplayStyle','tile','ShowEmptyBins','on','XBinedges',-5:0.025:3,'YBinedges',-25:2.5:100,'Normalization','probability');
h4Vals = h4.Values;
% RGB image for Awake
GammaAwakeRGB = figure;
s = pcolor(-4.975:0.025:3,-22.5:2.5:100,h4Vals');
s.FaceColor = 'interp';
set(s,'EdgeColor','none');
n = 50;         
R = linspace(0,1,n);
G = linspace(0,1,n); 
B = linspace(0,1,n);
colormap(flipud([R(:),G(:),B(:)]));
cax = caxis;
caxis([cax(1),cax(2)/1.5])
axis off
h4Frame = getframe(gcf);
h4Img = frame2im(h4Frame);
close(GammaAwakeHist)
close(GammaAwakeRGB)
% histogram for NREM
GammaNremHist = figure;
h5 = histogram2(data.GammaRel.catPupil.NREM,data.GammaRel.catGamma.NREM,'DisplayStyle','tile','ShowEmptyBins','on','XBinedges',-5:0.025:3,'YBinedges',-25:2.5:100,'Normalization','probability');
h5Vals = h5.Values;
% RGB image for NREM
GammaNremRGB = figure;
s = pcolor(-4.975:0.025:3,-22.5:2.5:100,h5Vals');
s.FaceColor = 'interp';
set(s,'EdgeColor','none');
n = 50;         
R = linspace(0,1,n);
G = linspace(0.4,1,n); 
B = linspace(0,1,n);
colormap(flipud([R(:),G(:),B(:)]));
cax = caxis;
caxis([cax(1),cax(2)/1.5])
axis off
h5Frame = getframe(gcf);
h5Img = frame2im(h5Frame);
close(GammaNremHist)
close(GammaNremRGB)
% histogram for REM
GammaRemHist = figure;
h6 = histogram2(data.GammaRel.catPupil.REM,data.GammaRel.catGamma.REM,'DisplayStyle','tile','ShowEmptyBins','on','XBinedges',-5:0.025:3,'YBinedges',-25:2.5:100,'Normalization','probability');
h6Vals = h6.Values;
% RGB image for REM
GammaRemRGB = figure;
s = pcolor(-4.975:0.025:3,-22.5:2.5:100,h6Vals');
s.FaceColor = 'interp';
set(s,'EdgeColor','none');
n = 50;         
R = linspace(1,1,n);
G = linspace(0,1,n); 
B = linspace(1,1,n);
colormap(flipud([R(:),G(:),B(:)]));
cax = caxis;
caxis([cax(1),cax(2)/1.5])
axis off
h6Frame = getframe(gcf);
h6Img = frame2im(h6Frame);
close(GammaRemHist)
close(GammaRemRGB)
%% axis for composite images
Fig6A = figure('Name','Figure Panel 6 - Turner et al. 2022','Units','Normalized','OuterPosition',[0,0,1,1]);
subplot(1,2,1)
img = imagesc(-4.975:0.025:3,-22.5:2.5:100,h4Vals');
xlabel('Diameter (z-units)')
ylabel('\DeltaP/P (%)')
title('Pupil-Gamma axis template')
set(gca,'box','off')
axis square
axis xy
delete(img)
subplot(1,2,2)
img = imagesc(-4.975:0.025:3,-22.5:2.5:125,h1Vals');
xlabel('Diameter (z-units)')
ylabel('\Delta[HbT] (\muM)')
title('Pupil-HbT axis template')
set(gca,'box','off')
axis square
axis xy
delete(img)
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    set(Fig6A,'PaperPositionMode','auto');
    print('-vector','-dpdf','-fillpage',[dirpath 'Fig6A_JNeurosci2022'])
    close(Fig6A)
    imwrite(h1Img,[dirpath 'Fig6_HbTAwake.png'])
    imwrite(h2Img,[dirpath 'Fig6_HbTNREM.png'])
    imwrite(h3Img,[dirpath 'Fig6_HbTREM.png'])
    imwrite(h4Img,[dirpath 'Fig6_GammaAwake.png'])
    imwrite(h5Img,[dirpath 'Fig6_GammaNREM.png'])
    imwrite(h6Img,[dirpath 'Fig6_GammaREM.png'])
end
%% Figure
Fig6B = figure('Name','Figure Panel 6 - Turner et al. 2022','Units','Normalized','OuterPosition',[0,0,1,1]);
ax1 = subplot(2,3,1);
edges = -8:0.1:6.5;
yyaxis right
h1 = histogram(diameterAllCatMeans,edges,'Normalization','probability','EdgeColor','k','FaceColor',colors('dark candy apple red'));
ylabel('Probability','rotation',-90,'VerticalAlignment','bottom')
yyaxis left
p1 = plot(edges,sgolayfilt(medfilt1(awakeProbPerc,10,'truncate'),3,17),'-','color',colors('black'),'LineWidth',2);
hold on
p2 = plot(edges,sgolayfilt(medfilt1(nremProbPerc,10,'truncate'),3,17),'-','color',[0,0.4,0],'LineWidth',2);
p3 = plot(edges,sgolayfilt(medfilt1(remProbPerc,10,'truncate'),3,17),'-','color','m','LineWidth',2);
p4 = plot(edges,sgolayfilt(medfilt1(asleepProbPerc,10,'truncate'),3,17),'-','color',colors('royal purple'),'LineWidth',2);
ylabel({'Arousal-state probability (%)'})
xlim([-8,6.5])
ylim([0,100])
legend([p1,p2,p3,p4,h1],'Awake','NREM','REM','Asleep','\DeltaArea','Location','NorthEast')
title('Diameter vs. arousal state')
xlabel('Diameter (z-units)')
axis square
set(gca,'box','off')
set(gca,'TickLength',[0.03,0.03]);
set(h1,'facealpha',0.2);
ax1.TickLength = [0.03,0.03];
ax1.YAxis(1).Color = 'k';
ax1.YAxis(2).Color = colors_eLife2020('dark candy apple red');
%% Gamma
subplot(2,3,2)
gammaPupilImg = imread('GammaPupilStack.png'); % needs made by combining images in ImageJ (Z project min)
imshow(gammaPupilImg)
axis off
title('Pupil-Gamma')
xlabel('Diameter (z-units)') 
ylabel('\DeltaP/P (%)') 
%% HbT
subplot(2,3,3)
HbTPupilImg = imread('HbTPupilStack.png'); % needs made by combining images in ImageJ (Z project min)
imshow(HbTPupilImg)
axis off
title('Pupil-HbT')
xlabel('Diameter (z-units)') 
ylabel('\Delta[HbT] (\muM)') 
%% sleep model confusion matrix
subplot(2,4,5)
cm = confusionchart(data.physio.holdYlabels,data.physio.holdXlabels);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,4])/totalScores))*100,1);
cm.Title = {'Physio SVM',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
%% sleep model confusion matrix
subplot(2,4,6)
cm = confusionchart(data.pupil.holdYlabels,data.pupil.holdXlabels);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,4])/totalScores))*100,1);
cm.Title = {'Pupil SVM',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
%% sleep model 10-fold loss
ax6 = subplot(2,4,7);
s1 = scatter(ones(1,length(data.physio.loss))*1,data.physio.loss,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('black'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.physio.meanLoss,data.physio.stdLoss,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(ones(1,length(data.pupil.loss))*2,data.pupil.loss,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','on','jitterAmount',0.25);
hold on
e2 = errorbar(2,data.pupil.meanLoss,data.pupil.stdLoss,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
title('10-fold cross validation')
ylabel('Loss (mean squared error)')
legend([s1,s2],'Physio mdl','Pupil mdl','Location','NorthWest')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,3])
ylim([0,0.2])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% model coherence
subplot(2,4,8);
s1 = semilogx(data.Coherr.meanPhysiof,data.Coherr.meanPhysioC,'color',colors('black'),'LineWidth',2);
hold on
semilogx(data.Coherr.meanPhysiof,data.Coherr.meanPhysioC + data.Coherr.stdPhysioC,'color',colors('black'),'LineWidth',0.5);
semilogx(data.Coherr.meanPhysiof,data.Coherr.meanPhysioC - data.Coherr.stdPhysioC,'color',colors('black'),'LineWidth',0.5);
s2 = semilogx(data.Coherr.meanPupilf,data.Coherr.meanPupilC,'color',colors('sapphire'),'LineWidth',2);
semilogx(data.Coherr.meanPupilf,data.Coherr.meanPupilC + data.Coherr.stdPupilC,'color',colors('sapphire'),'LineWidth',0.5);
semilogx(data.Coherr.meanPupilf,data.Coherr.meanPupilC - data.Coherr.stdPupilC,'color',colors('sapphire'),'LineWidth',0.5);
x1 = xline(1/30,'color',[0,0.4,0]);
x2 = xline(1/60,'color','m');
title('Model accuracy coherence')
ylabel('Coherence')
xlabel('Freq (Hz)')
legend([s1,s2,x1,x2],'Physio mdl','Pupil mdl','NREM req','REM req')
axis square
% xlim([0.003,0])
ylim([0,1])
set(gca,'box','off')
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    set(Fig6B,'PaperPositionMode','auto');
    print('-vector','-dpdf','-fillpage',[dirpath 'Fig6B_JNeurosci2022'])
end

end
