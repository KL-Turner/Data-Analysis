function [] = SingleTrialExample_Asleep_GCaMP_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
path = [rootFolder delim 'Results_Turner'];
cd(path)
% load file and gather information
exampleProcFile = 'T233_220206_14_10_44_ProcData.mat';
load(exampleProcFile)
exampleSpecFile = 'T233_220206_14_10_44_SpecDataA.mat';
load(exampleSpecFile)
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z,p,k] = butter(4,1/(ProcData.notes.CBVCamSamplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
binWhiskers = ProcData.data.binWhiskerAngle;
% data
HbT = filtfilt(sos,g,ProcData.data.HbT.RH);
HbO = filtfilt(sos,g,ProcData.data.HbO.RH);
HbR = filtfilt(sos,g,ProcData.data.HbR.RH);
GCaMP = filtfilt(sos,g,ProcData.data.GCaMP.RH);
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
p1 = plot((1:length(HbO))/ProcData.notes.CBVCamSamplingRate,HbO,'color',colors('blue grotto'),'LineWidth',1);
p2 = plot((1:length(HbR))/ProcData.notes.CBVCamSamplingRate,HbR,'color',colors('gold'),'LineWidth',1);
p3 = plot((1:length(HbT))/ProcData.notes.CBVCamSamplingRate,HbT,'color',colors('magenta'),'LineWidth',1);
ylabel('\Delta[Hb] (\muM)')
ylim([-30,160])
yyaxis right
p4 = plot((1:length(GCaMP))/ProcData.notes.CBVCamSamplingRate,(GCaMP - 1)*100,'color',colors('fluorescent green'),'LineWidth',1);
ylabel('\DeltaF/F (%)','rotation',-90,'VerticalAlignment','bottom')
ylim([-7,17])
legend([p1,p2,p3,p4,s1],'HbO','HbR','HbT','GCaMP','whisking')
xlabel('Time (s)')
set(gca,'box','off')
xlim([140,540])
ax1.YAxis(1).Color = colors('black');
ax1.YAxis(2).Color = colors('fluorescent green');
% Hippocampal electrode spectrogram
ax2 = subplot(4,1,4);
Semilog_ImageSC(T,F,hipNormS,'y')
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
xlim([140,540])
set(gca,'box','off')
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
end