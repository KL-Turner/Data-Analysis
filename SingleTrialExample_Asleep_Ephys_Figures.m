function [] = SingleTrialExample_Asleep_Ephys_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
path = [rootFolder delim 'Results_Turner'];
cd(path)
% load file and gather information
exampleProcFile = 'T151_210118_15_17_58_ProcData.mat';
% exampleProcFile = 'T195_210518_13_47_08_ProcData.mat';
load(exampleProcFile)
exampleSpecFile = 'T151_210118_15_17_58_SpecDataA.mat';
% exampleSpecFile = 'T195_210518_13_47_08_SpecDataA.mat';
load(exampleSpecFile)
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z,p,k] = butter(4,1/(ProcData.notes.CBVCamSamplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
binWhiskers = ProcData.data.binWhiskerAngle;
% data
HbT = filtfilt(sos,g,ProcData.data.HbT.RH);
% RH gamma baseline for T151 Jan 18
gammaBaseline = 3.513854434241187e-10;
gamma = filtfilt(sos,g,((ProcData.data.cortical_RH.gammaBandPower - gammaBaseline)./gammaBaseline)*100);
% cortical and hippocampal spectrograms
hipNormS = SpecData.hippocampus.normS.*100;
T = SpecData.cortical_LH.T;
F = SpecData.cortical_LH.F;
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
ylim([-500,1000])
legend([p1,p2,s1],'HbT','gamma power','whisking')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([310,370,430,490,550,610,670])
title('RH SIBF')
xlim([310,670])
ax1.YAxis(1).Color = colors('black');
ax1.YAxis(2).Color = colors('magenta');
% hippocampal electrode spectrogram
ax2 = subplot(4,1,4);
Semilog_ImageSC(T,F,hipNormS,'y')
c2 = colorbar;
ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
title('Hippocampal LFP')
xlabel('Time (min)')
ylabel('Frequency (Hz)')
xlim([310,670])
set(gca,'box','off')
xticks([310,370,430,490,550,610,670])
xticklabels({'0','1','2','3','4','5','6'})
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
    savefig(summaryFigure,[dirpath 'AsleepExample_Ephys']);
    cla(ax2);
    set(ax2,'YLim',[1,99]);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-vector','-dpdf','-fillpage',[dirpath 'AsleepExample_Ephys'])
    close(summaryFigure)
    % spectrogram image
    subplotImgs = figure;
    Semilog_ImageSC(T,F,hipNormS,'y')
    caxis([-100,100])
    set(gca,'box','off')
    axis xy
    axis tight
    axis off
    xlim([310,670])
    print('-vector','-dtiffn',[dirpath 'AsleepExample_Ephys_Spectrogram'])
    close(subplotImgs)
end