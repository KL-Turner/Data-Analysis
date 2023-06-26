function [] = SingleTrialExample_Asleep_Ephys_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
path = [rootFolder delim 'Results_Turner'];
cd(path)
% load file and gather information
exampleProcFile_AsleepEphys = 'T165_210223_12_37_03_ProcData.mat';
load(exampleProcFile_AsleepEphys)
exampleSpecFile_AsleepEphys = 'T165_210223_12_37_03_SpecDataA.mat';
load(exampleSpecFile_AsleepEphys)
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z_AsleepEphys,p_AsleepEphys,k_AsleepEphys] = butter(4,1/(ProcData.notes.CBVCamSamplingRate/2),'low');
[sos_AsleepEphys,g_AsleepEphys] = zp2sos(z_AsleepEphys,p_AsleepEphys,k_AsleepEphys);
binWhiskers_AsleepEphys = ProcData.data.binWhiskerAngle;
% data
HbT_AsleepEphys = filtfilt(sos_AsleepEphys,g_AsleepEphys,ProcData.data.HbT.RH);
% RH gamma baseline for T151 Jan 18
gammaBaseline_AsleepEphys = 3.513854434241187e-10;
gamma_AsleepEphys = filtfilt(sos_AsleepEphys,g_AsleepEphys,((ProcData.data.cortical_RH.gammaBandPower - gammaBaseline_AsleepEphys)./gammaBaseline_AsleepEphys)*100);
% cortical and hippocampal spectrograms
hipNormS_AsleepEphys = SpecData.hippocampus.normS.*100;
T_AsleepEphys = SpecData.cortical_LH.T;
F_AsleepEphys = SpecData.cortical_LH.F;
% Yvals for behavior Indices
whisking_Yvals_AsleepEphys = 155*ones(size(binWhiskers_AsleepEphys));
whiskInds_AsleepEphys = binWhiskers_AsleepEphys.*whisking_Yvals_AsleepEphys;
% set whisk indeces
for x = 1:length(whiskInds_AsleepEphys)
    if whiskInds_AsleepEphys(1,x) == 0
        whiskInds_AsleepEphys(1,x) = NaN;
    end
end
% figure
summaryFigure = figure;
% force sensor and EMG
ax1 = subplot(4,1,1:3);
% hemodynamic and behavioral indeces
s1 = scatter((1:length(binWhiskers_AsleepEphys))/ProcData.notes.dsFs,whiskInds_AsleepEphys,'.','MarkerEdgeColor',colors('black'));
hold on;
p1 = plot((1:length(HbT_AsleepEphys))/ProcData.notes.CBVCamSamplingRate,HbT_AsleepEphys,'color',colors('black'),'LineWidth',1);
ylabel('\Delta[Hb] (\muM)')
ylim([-50,160])
yyaxis right
p2 = plot((1:length(gamma_AsleepEphys))/ProcData.notes.CBVCamSamplingRate,gamma_AsleepEphys,'color',colors('magenta'),'LineWidth',1);
ylim([-300,1300])
legend([p1,p2,s1],'HbT','gamma power','whisking')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([450,510,570,630,690,750,810])
title('RH SIBF')
xlim([450,810])
ax1.YAxis(1).Color = colors('black');
ax1.YAxis(2).Color = colors('magenta');
% hippocampal electrode spectrogram
ax2 = subplot(4,1,4);
Semilog_ImageSC(T_AsleepEphys,F_AsleepEphys,hipNormS_AsleepEphys,'y')
c2 = colorbar;
ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
clim([-100,100])
title('Hippocampal LFP')
xlabel('Time (min)')
ylabel('Frequency (Hz)')
xlim([450,810])
set(gca,'box','off')
xticks([450,510,570,630,690,750,810])
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
    Semilog_ImageSC(T_AsleepEphys,F_AsleepEphys,hipNormS_AsleepEphys,'y')
    clim([-100,100])
    set(gca,'box','off')
    axis xy
    axis tight
    axis off
    xlim([450,810])
    print('-vector','-dtiffn',[dirpath 'AsleepExample_Ephys_Spectrogram'])
    close(subplotImgs)
end