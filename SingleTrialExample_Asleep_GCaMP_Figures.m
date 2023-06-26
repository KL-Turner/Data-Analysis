function [] = SingleTrialExample_Asleep_GCaMP_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
path = [rootFolder delim 'Results_Turner'];
cd(path)
% load file and gather information
exampleProcFile_AsleepGCaMP = 'T233_220206_14_10_44_ProcData.mat';
load(exampleProcFile_AsleepGCaMP)
exampleSpecFile_AsleepGCaMP = 'T233_220206_14_10_44_SpecDataA.mat';
load(exampleSpecFile_AsleepGCaMP)
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z_AsleepGCaMP,p_AsleepGCaMP,k_AsleepGCaMP] = butter(4,1/(ProcData.notes.CBVCamSamplingRate/2),'low');
[sos_AsleepGCaMP,g_AsleepGCaMP] = zp2sos(z_AsleepGCaMP,p_AsleepGCaMP,k_AsleepGCaMP);
binWhiskers_AsleepGCaMP = ProcData.data.binWhiskerAngle;
% data
HbT_AsleepGCaMP = filtfilt(sos_AsleepGCaMP,g_AsleepGCaMP,ProcData.data.HbT.RH);
HbO_AsleepGCaMP = filtfilt(sos_AsleepGCaMP,g_AsleepGCaMP,ProcData.data.HbO.RH);
HbR_AsleepGCaMP = filtfilt(sos_AsleepGCaMP,g_AsleepGCaMP,ProcData.data.HbR.RH);
GCaMP_AsleepGCaMP = filtfilt(sos_AsleepGCaMP,g_AsleepGCaMP,ProcData.data.GCaMP.RH);
% cortical and hippocampal spectrograms
hipNormS_AsleepGCaMP = SpecData.hippocampus.normS.*100;
T_AsleepGCaMP = SpecData.cortical_LH.T;
F_AsleepGCaMP = SpecData.cortical_LH.F;
% Yvals for behavior Indices
whisking_Yvals_AsleepGCaMP = 155*ones(size(binWhiskers_AsleepGCaMP));
whiskInds_AsleepGCaMP = binWhiskers_AsleepGCaMP.*whisking_Yvals_AsleepGCaMP;
% set whisk indeces
for x = 1:length(whiskInds_AsleepGCaMP)
    if whiskInds_AsleepGCaMP(1,x) == 0
        whiskInds_AsleepGCaMP(1,x) = NaN;
    end
end
% figure
summaryFigure = figure;
% force sensor and EMG
ax1 = subplot(4,1,1:3);
% hemodynamic and behavioral indeces
s1 = scatter((1:length(binWhiskers_AsleepGCaMP))/ProcData.notes.dsFs,whiskInds_AsleepGCaMP,'.','MarkerEdgeColor',colors('black'));
hold on;
p1 = plot((1:length(HbO_AsleepGCaMP))/ProcData.notes.CBVCamSamplingRate,HbO_AsleepGCaMP,'color',colors('blue grotto'),'LineWidth',1);
p2 = plot((1:length(HbR_AsleepGCaMP))/ProcData.notes.CBVCamSamplingRate,HbR_AsleepGCaMP,'color',colors('gold'),'LineWidth',1);
p3 = plot((1:length(HbT_AsleepGCaMP))/ProcData.notes.CBVCamSamplingRate,HbT_AsleepGCaMP,'color',colors('black'),'LineWidth',1);
ylabel('\Delta[Hb] (\muM)')
ylim([-50,160])
yyaxis right
p4 = plot((1:length(GCaMP_AsleepGCaMP))/ProcData.notes.CBVCamSamplingRate,(GCaMP_AsleepGCaMP - 1)*100,'color',colors('fluorescent green'),'LineWidth',1);
ylabel('\DeltaF/F (%)','rotation',-90,'VerticalAlignment','bottom')
ylim([-7,20])
legend([p1,p2,p3,p4,s1],'HbO','HbR','HbT','GCaMP','whisking')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([215,275,335,395,455,515,575])
title('RH SIBF')
xlim([215,575])
ax1.YAxis(1).Color = colors('black');
ax1.YAxis(2).Color = colors('fluorescent green');
% hippocampal electrode spectrogram
ax2 = subplot(4,1,4);
Semilog_ImageSC(T_AsleepGCaMP,F_AsleepGCaMP,hipNormS_AsleepGCaMP,'y')
c2 = colorbar;
ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
clim([-100,100])
title('Hippocampal LFP')
xlabel('Time (min)')
ylabel('Frequency (Hz)')
xlim([215,575])
set(gca,'box','off')
xticks([215,275,335,395,455,515,575])
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
    savefig(summaryFigure,[dirpath 'AsleepExample_GCaMP']);
    cla(ax2);
    set(ax2,'YLim',[1,99]);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-vector','-dpdf','-fillpage',[dirpath 'AsleepExample_GCaMP'])
    close(summaryFigure)
    % spectrogram image
    subplotImgs = figure;
    Semilog_ImageSC(T_AsleepGCaMP,F_AsleepGCaMP,hipNormS_AsleepGCaMP,'y')
    clim([-100,100])
    set(gca,'box','off')
    axis xy
    axis tight
    axis off
    xlim([215,575])
    print('-vector','-dtiffn',[dirpath 'AsleepExample_GCaMP_Spectrogram'])
    close(subplotImgs)
end