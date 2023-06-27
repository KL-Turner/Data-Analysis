function [] = Fig6_nNOS(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
path = [rootFolder delim 'Results_Turner'];
cd(path)
% load file and gather information
exampleProcFile = 'T165_210222_12_45_38_ProcData.mat';
load(exampleProcFile)
exampleSpecFile = 'T165_210222_12_45_38_SpecDataA.mat';
load(exampleSpecFile)
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z,p,k] = butter(4,1/(ProcData.notes.CBVCamSamplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
binWhiskers = ProcData.data.binWhiskerAngle;
% data
HbT = filtfilt(sos,g,ProcData.data.HbT.RH);
% RH gamma baseline for T151 Jan 30
gammaBaseline = 1.016060094087967e-10;
gamma = filtfilt(sos,g,((ProcData.data.cortical_RH.gammaBandPower - gammaBaseline)./gammaBaseline)*100);
% cortical and hippocampal spectrograms
cortNormS = SpecData.cortical_RH.normS.*100;
T = SpecData.cortical_RH.T;
F = SpecData.cortical_RH.F;
% Yvals for behavior Indices
whisking_Yvals = 155*ones(size(binWhiskers));
whiskInds = binWhiskers.*whisking_Yvals;
% set whisk indeces
for x = 1:length(whiskInds)
    if whiskInds(1,x) == 0
        whiskInds(1,x) = NaN;
    end
end
% figure
summaryFigure = figure;
% force sensor and EMG
ax1 = subplot(4,1,1:3);
% hemodynamic and behavioral indeces
s1 = scatter((1:length(binWhiskers))/ProcData.notes.dsFs,whiskInds,'.','MarkerEdgeColor',colors('black'));
hold on;
p1 = plot((1:length(HbT))/ProcData.notes.CBVCamSamplingRate,HbT,'color',colors('black'),'LineWidth',1);
ylabel('\Delta[Hb] (\muM)')
ylim([-50,160])
yyaxis right
p2 = plot((1:length(gamma))/ProcData.notes.CBVCamSamplingRate,gamma,'color',colors('magenta'),'LineWidth',1);
ylim([-300,1300])
legend([p1,p2,s1],'HbT','gamma power','whisking')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([310,370,430,490])
title('RH SIBF')
xlim([310,490])
ax1.YAxis(1).Color = colors('black');
ax1.YAxis(2).Color = colors('magenta');
% cortical electrode spectrogram
ax2 = subplot(4,1,4);
Semilog_ImageSC(T,F,cortNormS,'y')
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
clim([-100,100])
title('Cortical LFP')
xlabel('Time (min)')
ylabel('Frequency (Hz)')
xlim([310,490])
set(gca,'box','off')
xticks([310,370,430,490])
xticklabels({'0','1','2','3'});
% axis properties
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax2Pos(3) = ax1Pos(3);
set(ax2,'position',ax2Pos);
% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures' delim 'Example Files' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'AwakeExample_Ephys']);
    cla(ax2);
    set(ax2,'YLim',[1,99]);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-vector','-dpdf','-fillpage',[dirpath 'AwakeExample_Ephys'])
    close(summaryFigure)
    % spectrogram image
    subplotImgs = figure;
    Semilog_ImageSC(T,F,cortNormS,'y')
    clim([-100,100])
    set(gca,'box','off')
    axis xy
    axis tight
    axis off
    xlim([310,490])
    print('-vector','-dtiffn',[dirpath 'AwakeExample_Ephys_Spectrogram'])
    close(subplotImgs)
end