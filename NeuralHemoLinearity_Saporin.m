function [AnalysisResults] = NeuralHemoLinearity_Saporin(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Generate figure panel 8 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% set-up and process data
%% set-up and process data
animalIDs = {'T141','T155','T156','T157','T142','T144','T159','T172','T150','T165','T166','T177','T179','T186','T187','T188','T189'};
C57BL6J_IDs = {'T141','T155','T156','T157','T186','T187','T188','T189'};
SSP_SAP_IDs = {'T142','T144','T159','T172'};
Blank_SAP_IDs = {'T150','T165','T166','T177','T179'};
treatments = {'C57BL6J','SSP_SAP','Blank_SAP'};
behavFields = {'Awake','NREM','REM'};
%% average coherence during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    % recognize treatment based on animal group
    if ismember(animalID,C57BL6J_IDs) == true
        treatment = 'C57BL6J';
    elseif ismember(animalIDs{1,aa},SSP_SAP_IDs) == true
        treatment = 'SSP_SAP';
    elseif ismember(animalIDs{1,aa},Blank_SAP_IDs) == true
        treatment = 'Blank_SAP';
    end
    % pre-allocate necessary variable fields
    data.(treatment).(behavField).dummCheck = 1;
    if isfield(data.(treatment).(behavField),'LH_gamma') == false
        data.(treatment).(behavField).LH_gamma = {};
        data.(treatment).(behavField).LH_HbT = {};
        data.(treatment).(behavField).RH_gamma = {};
        data.(treatment).(behavField).RH_HbT = {};
    end
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        data.(treatment).(behavField).LH_gamma = cat(2,data.(treatment).(behavField).LH.C,AnalysisResults.(animalID).WhiskHemoCoherence.(behavField).adjLH.C);
        data.(treatment).(behavField).LH_HbT = cat(1,data.(treatment).(behavField).LH.f,AnalysisResults.(animalID).WhiskHemoCoherence.(behavField).adjLH.f);
        data.(treatment).(behavField).RH_gamma = cat(2,data.(treatment).(behavField).RH.C,AnalysisResults.(animalID).WhiskHemoCoherence.(behavField).adjRH.C);
        data.(treatment).(behavField).RH_HbT = cat(1,data.(treatment).(behavField).RH.f,AnalysisResults.(animalID).WhiskHemoCoherence.(behavField).adjRH.f);
    end
end
%% take data from each animal corresponding to the CBV-gamma relationship
catHbT = []; catGam = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        if isfield(catHbT,behavField) == false
            catHbT.(behavField) = []; 
            catGam.(behavField) = [];
        end
        catHbT.(behavField) = cat(1,catHbT.(behavField),AnalysisResults.(animalID).HbTvsGamma.(behavField).HbT);
        catGam.(behavField) = cat(1,catGam.(behavField),AnalysisResults.(animalID).HbTvsGamma.(behavField).Gamma);
    end
end
%% average coherence during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.NeuralHemoCoherence = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        % create the behavior folder for the first iteration of the loop
        if isfield(data.NeuralHemoCoherence,behavField) == false
            data.NeuralHemoCoherence.(behavField) = [];
        end
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.C) == false
                % create the data type folder for the first iteration of the loop
                if isfield(data.NeuralHemoCoherence.(behavField),dataType) == false
                    data.NeuralHemoCoherence.(behavField).(dataType).C = [];
                    data.NeuralHemoCoherence.(behavField).(dataType).f = [];
                    data.NeuralHemoCoherence.(behavField).(dataType).confC = [];
                end
                % concatenate C/f for existing data - exclude any empty sets
                data.NeuralHemoCoherence.(behavField).(dataType).C = cat(2,data.NeuralHemoCoherence.(behavField).(dataType).C,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.C,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjRH.C);
                data.NeuralHemoCoherence.(behavField).(dataType).f = cat(1,data.NeuralHemoCoherence.(behavField).(dataType).f,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.f,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjRH.f);
                data.NeuralHemoCoherence.(behavField).(dataType).confC = cat(1,data.NeuralHemoCoherence.(behavField).(dataType).confC,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.confC,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjRH.confC);
            end
        end
    end
end
% take mean/StD of C/f and determine confC line
for ee = 1:length(behavFields)
    behavField = behavFields{1,ee};
    for ff = 1:length(dataTypes)
        dataType = dataTypes{1,ff};
        data.NeuralHemoCoherence.(behavField).(dataType).meanC = mean(data.NeuralHemoCoherence.(behavField).(dataType).C,2);
        data.NeuralHemoCoherence.(behavField).(dataType).stdC = std(data.NeuralHemoCoherence.(behavField).(dataType).C,0,2);
        data.NeuralHemoCoherence.(behavField).(dataType).meanf = mean(data.NeuralHemoCoherence.(behavField).(dataType).f,1);
        data.NeuralHemoCoherence.(behavField).(dataType).maxConfC = geomean(data.NeuralHemoCoherence.(behavField).(dataType).confC);
        data.NeuralHemoCoherence.(behavField).(dataType).maxConfC_Y = ones(length(data.NeuralHemoCoherence.(behavField).(dataType).meanf),1)*data.NeuralHemoCoherence.(behavField).(dataType).maxConfC;
    end
end
%% Fig. 8
summaryFigure = figure('Name','Fig8 (a-d)');
sgtitle('Figure 8 - Turner et al. 2020')
%% [8a] HbT vs. arousal-state probability
ax1 = subplot(2,3,1);
edges = -35:1:115;
yyaxis right
h1 = histogram(HbTallCatMeans,edges,'Normalization','probability','EdgeColor','k','FaceColor',colors_eLife2020('dark candy apple red'));
ylabel({'5-sec Mean \DeltaHbT','Probability distribution'},'rotation',-90,'VerticalAlignment','bottom')
yyaxis left
p1 = plot(edges,sgolayfilt(medfilt1(awakeProbPerc,10,'truncate'),3,17),'-','color',colors_eLife2020('rich black'),'LineWidth',2);
hold on
p2 = plot(edges,sgolayfilt(medfilt1(nremProbPerc,10,'truncate'),3,17),'-','color',colorNREM,'LineWidth',2);
p3 = plot(edges,sgolayfilt(medfilt1(remProbPerc,10,'truncate'),3,17),'-','color',colorREM,'LineWidth',2);
ylabel({'Arousal-state probability (%)'})
xlim([-35,115])
ylim([0,85])
legend([p1,p2,p3,h1],'Awake','NREM','REM','\DeltaHbT','Location','NorthEast')
title({'[8a] 5-sec mean \Delta[HbT] (\muM)','vs. arousal state probability'})
xlabel({'\Delta[HbT] (\muM)','1 \muM bins'})
axis square
set(gca,'box','off')
set(gca,'TickLength',[0.03,0.03]);
ylim([0,90])
xlim([-35,115])
set(h1,'facealpha',0.2);
ax1.TickLength = [0.03,0.03];
ax1.YAxis(1).Color = 'k';
ax1.YAxis(2).Color = colors_eLife2020('dark candy apple red');
%% [8b] D/D vs. arousal-state probability
ax2 = subplot(2,3,2);
edges = -20:1:50;
yyaxis right
h2 = histogram(TwoPallCatMeans,edges,'Normalization','probability','EdgeColor','k','FaceColor',colors_eLife2020('dark candy apple red'));
ylabel({'5-sec Mean \DeltaD/D (%)','Probability distribution'},'rotation',-90,'VerticalAlignment','bottom')
yyaxis left
p1 = plot(edges,sgolayfilt(medfilt1(TwoPawakeProbPerc,10,'truncate'),3,17),'-','color',colors_eLife2020('rich black'),'LineWidth',2);
hold on
p2 = plot(edges,sgolayfilt(medfilt1(TwoPnremProbPerc,10,'truncate'),3,17),'-','color',colorNREM,'LineWidth',2);
p3 = plot(edges,sgolayfilt(medfilt1(TwoPremProbPerc,10,'truncate'),3,17),'-','color',colorREM,'LineWidth',2);
ylabel({'Arousal-state probability (%)'})
xlim([-20,50])
ylim([0,90])
legend([p1,p2,p3,h2],'Awake','NREM','REM','\DeltaD/D','Location','NorthEast')
title({'[S20a] 5-sec mean \DeltaD/D (%)','vs. arousal state probability'})
xlabel({'\DeltaD/D (%)','1 (%) bins'})
axis square
set(gca,'box','off')
set(gca,'TickLength',[0.03,0.03]);
ylim([0,90])
xlim([-20,50])
set(h2,'facealpha',0.2);
ax2.TickLength = [0.03,0.03];
ax2.YAxis(1).Color = 'k';
ax2.YAxis(2).Color = colors_eLife2020('dark candy apple red');
%% [8c] coherence^2 between HbT and gamma-band power during different arousal-states
ax3 = subplot(2,3,3);
s1 = semilogx(data.NeuralHemoCoherence.Rest.gammaBandPower.meanf,data.NeuralHemoCoherence.Rest.gammaBandPower.meanC.^2,'color',colorRest,'LineWidth',2);
hold on
rectangle('Position',[0.005,0.05,0.1 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
s2 = semilogx(data.NeuralHemoCoherence.NREM.gammaBandPower.meanf,data.NeuralHemoCoherence.NREM.gammaBandPower.meanC.^2,'color',colorNREM,'LineWidth',2);
rectangle('Position',[0.005,0.05,1/30 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
s3 = semilogx(data.NeuralHemoCoherence.REM.gammaBandPower.meanf,data.NeuralHemoCoherence.REM.gammaBandPower.meanC.^2,'color',colorREM,'LineWidth',2);
rectangle('Position',[0.005,0.05,1/60 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
s4 = semilogx(data.NeuralHemoCoherence.Awake.gammaBandPower.meanf,data.NeuralHemoCoherence.Awake.gammaBandPower.meanC.^2,'color',colorAlert,'LineWidth',2);
s5 = semilogx(data.NeuralHemoCoherence.Sleep.gammaBandPower.meanf,data.NeuralHemoCoherence.Sleep.gammaBandPower.meanC.^2,'color',colorAsleep,'LineWidth',2);
s6 = semilogx(data.NeuralHemoCoherence.All.gammaBandPower.meanf,data.NeuralHemoCoherence.All.gammaBandPower.meanC.^2,'color',colorAll,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[8c] Neural-hemo coherence^2','Gamma-band power and \Delta[HbT] \muM (%)'})
legend([s1,s2,s3,s4,s5,s6],'Rest','NREM','REM','Alert','Asleep','All','Location','SouthEast')
axis square
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [8d] histogram images for HbT vs. gamma-band power during different arousal-states
% histogram for Awake
awakeHist = figure;
h1 = histogram2(catGam.Awake,catHbT.Awake,'DisplayStyle','tile','ShowEmptyBins','on','XBinedges',-0.25:0.025:1,'YBinedges',-25:2.5:125,'Normalization','probability');
h1Vals = h1.Values;
% RGB image for Awake
awakeRGB = figure;
s = pcolor(-0.225:0.025:1,-22.5:2.5:125,h1Vals');
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
subplot(2,3,4)
imshow(h1Img)
axis off
title('[8d] rfc-Awake Gamma-[HbT]')
% histogram for NREM
nremHist = figure;
h2 = histogram2(catGam.NREM,catHbT.NREM,'DisplayStyle','tile','ShowEmptyBins','on','XBinedges',-0.25:0.025:1,'YBinedges',-25:2.5:125,'Normalization','probability');
h2Vals = h2.Values;
% RGB image for NREM
nremRGB = figure;
s = pcolor(-0.225:0.025:1,-22.5:2.5:125,h2Vals');
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
subplot(2,3,5)
imshow(h2Img)
axis off
title('[8d] rfc-NREM Gamma-[HbT]')
% histogram for REM
remHist = figure;
h3 = histogram2(catGam.REM,catHbT.REM,'DisplayStyle','tile','ShowEmptyBins','on','XBinedges',-0.25:0.025:1,'YBinedges',-25:2.5:125,'Normalization','probability');
h3Vals = h3.Values;
% RGB image for REM
remRGB = figure;
s = pcolor(-0.225:0.025:1,-22.5:2.5:125,h3Vals');
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
subplot(2,3,6)
imshow(h3Img)
axis off
title('[8d] rfc-REM Gamma-[HbT]')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    imwrite(h1Img,[dirpath 'Fig8d_Awake.png'])
    imwrite(h2Img,[dirpath 'Fig8d_NREM.png'])
    imwrite(h3Img,[dirpath 'Fig8d_REM.png'])
    savefig(summaryFigure,[dirpath 'Fig8']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig8'])
end

end
