function [] = SingleTrialExample_Awake_GCaMP_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
path = [rootFolder delim 'Results_Turner'];
cd(path)
% load file and gather information
exampleProcFile_AwakeGCaMP = 'T233_220206_13_07_08_ProcData.mat';
load(exampleProcFile_AwakeGCaMP)
exampleSpecFile_AwakeGCaMP = 'T233_220206_13_07_08_SpecDataA.mat';
load(exampleSpecFile_AwakeGCaMP)
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z_AwakeGCaMP,p_AwakeGCaMP,k_AwakeGCaMP] = butter(4,1/(ProcData.notes.CBVCamSamplingRate/2),'low');
[sos_AwakeGCaMP,g_AwakeGCaMP] = zp2sos(z_AwakeGCaMP,p_AwakeGCaMP,k_AwakeGCaMP);
binWhiskers_AwakeGCaMP = ProcData.data.binWhiskerAngle;
% data
HbT_AwakeGCaMP = filtfilt(sos_AwakeGCaMP,g_AwakeGCaMP,ProcData.data.HbT.RH);
HbO_AwakeGCaMP = filtfilt(sos_AwakeGCaMP,g_AwakeGCaMP,ProcData.data.HbO.RH);
HbR_AwakeGCaMP = filtfilt(sos_AwakeGCaMP,g_AwakeGCaMP,ProcData.data.HbR.RH);
GCaMP_AwakeGCaMP = filtfilt(sos_AwakeGCaMP,g_AwakeGCaMP,ProcData.data.GCaMP.RH);
% cortical and hippocampal spectrograms
cortNormS_AwakeGCaMP = SpecData.cortical_RH.normS.*100;
T_AwakeGCaMP = SpecData.cortical_LH.T;
F_AwakeGCaMP = SpecData.cortical_LH.F;
% Yvals for behavior Indices
whisking_Yvals_AwakeGCaMP = 155*ones(size(binWhiskers_AwakeGCaMP));
whiskInds_AwakeGCaMP = binWhiskers_AwakeGCaMP.*whisking_Yvals_AwakeGCaMP;
% set whisk indeces
for x = 1:length(whiskInds_AwakeGCaMP)
    if whiskInds_AwakeGCaMP(1,x) == 0
        whiskInds_AwakeGCaMP(1,x) = NaN;
    end
end
% figure
summaryFigure = figure;
% force sensor and EMG
ax1 = subplot(4,1,1:3);
% hemodynamic and behavioral indeces
s1 = scatter((1:length(binWhiskers_AwakeGCaMP))/ProcData.notes.dsFs,whiskInds_AwakeGCaMP,'.','MarkerEdgeColor',colors('black'));
hold on;
p1 = plot((1:length(HbO_AwakeGCaMP))/ProcData.notes.CBVCamSamplingRate,HbO_AwakeGCaMP,'color',colors('blue grotto'),'LineWidth',1);
p2 = plot((1:length(HbR_AwakeGCaMP))/ProcData.notes.CBVCamSamplingRate,HbR_AwakeGCaMP,'color',colors('gold'),'LineWidth',1);
p3 = plot((1:length(HbT_AwakeGCaMP))/ProcData.notes.CBVCamSamplingRate,HbT_AwakeGCaMP,'color',colors('black'),'LineWidth',1);
ylabel('\Delta[Hb] (\muM)')
ylim([-50,160])
yyaxis right
p4 = plot((1:length(GCaMP_AwakeGCaMP))/ProcData.notes.CBVCamSamplingRate,(GCaMP_AwakeGCaMP - 1)*100,'color',colors('fluorescent green'),'LineWidth',1);
ylabel('\DeltaF/F (%)','rotation',-90,'VerticalAlignment','bottom')
ylim([-7,20])
legend([p1,p2,p3,p4,s1],'HbO','HbR','HbT','GCaMP','whisking')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([260,320,380,440])
title('RH SIBF')
xlim([260,440])
ax1.YAxis(1).Color = colors('black');
ax1.YAxis(2).Color = colors('fluorescent green');
% cortical electrode spectrogram
ax2 = subplot(4,1,4);
Semilog_ImageSC(T_AwakeGCaMP,F_AwakeGCaMP,cortNormS_AwakeGCaMP,'y')
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
clim([-100,100])
title('Cortical LFP')
xlabel('Time (min)')
ylabel('Frequency (Hz)')
xlim([260,440])
set(gca,'box','off')
xticks([260,320,380,440])
xticklabels({'0','1','2','3'})
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
    savefig(summaryFigure,[dirpath 'AwakeExample_GCaMP']);
    cla(ax2);
    set(ax2,'YLim',[1,99]);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-vector','-dpdf','-fillpage',[dirpath 'AwakeExample_GCaMP'])
    close(summaryFigure)
    % spectrogram image
    subplotImgs = figure;
    Semilog_ImageSC(T_AwakeGCaMP,F_AwakeGCaMP,cortNormS_AwakeGCaMP,'y')
    clim([-100,100])
    set(gca,'box','off')
    axis xy
    axis tight
    axis off
    xlim([260,440])
    print('-vector','-dtiffn',[dirpath 'AwakeExample_GCaMP_Spectrogram'])
    close(subplotImgs)
end