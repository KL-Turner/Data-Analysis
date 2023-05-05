function [] = SingleTrialExample_GCaMP_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
path = [rootFolder delim 'Results_Turner'];
cd(path)
% load file and gather information
exampleProcFile = 'T233_220206_11_32_05_ProcData.mat';
load(exampleProcFile)
exampleSpecFile = 'T233_220206_11_32_05_SpecDataA.mat';
load(exampleSpecFile)
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z,p,k] = butter(4,1/(ProcData.notes.CBVCamSamplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
binWhiskers = ProcData.data.binWhiskerAngle;
% stimulations
LPadSol = ProcData.data.stimulations.LPadSol;
RPadSol = ProcData.data.stimulations.RPadSol;
AudSol = ProcData.data.stimulations.AudSol;
% data
HbT = filtfilt(sos,g,ProcData.data.HbT.RH);
HbO = filtfilt(sos,g,ProcData.data.HbO.RH);
HbR = filtfilt(sos,g,ProcData.data.HbR.RH);
GCaMP = filtfilt(sos,g,ProcData.data.GCaMP.RH);
% cortical and hippocampal spectrograms
cortNormS = SpecData.cortical_RH.normS.*100;
hipNormS = SpecData.hippocampus.normS.*100;
T = SpecData.cortical_LH.T;
F = SpecData.cortical_LH.F;
% Yvals for behavior Indices
whisking_Yvals = 1.10*max(HbO)*ones(size(binWhiskers));
LPad_Yvals = 1.50*max(HbO)*ones(size(LPadSol));
RPad_Yvals = 1.50*max(HbO)*ones(size(RPadSol));
Aud_Yvals = 1.50*max(HbO)*ones(size(AudSol));
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
ax1 = subplot(2,1,1);
% hemodynamic and behavioral indeces
s1 = scatter((1:length(binWhiskers))/ProcData.notes.dsFs,whiskInds,'.','MarkerEdgeColor',colors('black'));
hold on
s2 = scatter(LPadSol,LPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','c');
s3 = scatter(RPadSol,RPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','m');
s4 = scatter(AudSol,Aud_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','g');
p1 = plot((1:length(HbO))/ProcData.notes.CBVCamSamplingRate,HbO,'color',colors('blue grotto'),'LineWidth',1);
p2 = plot((1:length(HbR))/ProcData.notes.CBVCamSamplingRate,HbR,'color',colors('gold'),'LineWidth',1);
p3 = plot((1:length(HbT))/ProcData.notes.CBVCamSamplingRate,HbT,'color',colors('magenta'),'LineWidth',1);
ylabel('\Delta[Hb] (\muM)')
ylim([-25,80])
yyaxis right
p4 = plot((1:length(GCaMP))/ProcData.notes.CBVCamSamplingRate,(GCaMP - 1)*100,'color',colors('fluorescent green'),'LineWidth',1);
ylabel('\DeltaF/F (%)','Rotation',90)
ylim([-2,15])
legend([p1,p2,p3,p4,s1,s2,s3,s4],'HbO','HbR','HbT','GCaMP','whisking',',LPad sol','RPad sol','Aud sol')
xlabel('Time (s)')
set(gca,'box','off')
axis tight
xlim([145,325])
% Hippocampal electrode spectrogram
ax2 = subplot(2,1,2);
Semilog_ImageSC(T,F,cortNormS,'y')
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','Rotation',90)
caxis([-50,100])
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
xlim([145,325])
set(gca,'box','off')
% axis properties
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax2Pos(3:4) = ax1Pos(3:4);
set(ax2,'position',ax2Pos);
% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures' delim 'Awake Example GCaMP' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'AwakeExample_GCaMP']);
end