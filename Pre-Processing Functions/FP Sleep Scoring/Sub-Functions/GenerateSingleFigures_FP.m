function [figHandle,ax1,ax2,ax3,ax4,ax5,ax6] = GenerateSingleFigures_FP(procDataFileID,saveFigs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Create a summary figure for a single n minute IOS trial
%________________________________________________________________________________________________________________________

% load file and gather information
load(procDataFileID)
[animalID,fileDate,fileID] = GetFileInfo_FP(procDataFileID);
strDay = ConvertDate_FP(fileDate);
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1,k1] = butter(4,10/(ProcData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(4,1/(ProcData.notes.dsFs/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
% whisker angle
filteredWhiskerAngle = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
binWhiskers = ProcData.data.binWhiskerAngle;
% force sensor
filtForceSensor = filtfilt(sos1,g1,ProcData.data.forceSensor);
binForce = ProcData.data.binForceSensor;
% emg
EMG = ProcData.data.EMG.emg;
% stimulations
LPadSol = ProcData.data.stimulations.LPadSol;
RPadSol = ProcData.data.stimulations.RPadSol;
AudSol = ProcData.data.stimulations.AudSol;
% HbT and GCaMP data
LH_HbT = ProcData.data.HbT.LH;
filtLH_HbT = filtfilt(sos2,g2,LH_HbT);
RH_HbT = ProcData.data.HbT.RH;
filtRH_HbT = filtfilt(sos2,g2,RH_HbT);
LH_GCaMP = ProcData.data.GCaMP7s.LH;
filtLH_GCaMP = filtfilt(sos2,g2,LH_GCaMP);
RH_GCaMP = ProcData.data.GCaMP7s.RH;
filtRH_GCaMP = filtfilt(sos2,g2,RH_GCaMP);
% cortical and hippocampal spectrograms
specDataFile = [animalID '_' fileID '_SpecDataA.mat'];
load(specDataFile,'-mat');
cortical_LHnormS = SpecData.cortical_LH.normS.*100;
cortical_RHnormS = SpecData.cortical_RH.normS.*100;
hippocampusNormS = SpecData.hippocampus.normS.*100;
T = SpecData.cortical_LH.T;
F = SpecData.cortical_LH.F;
% Yvals for behavior Indices
indecesMax = max([filtLH_HbT;filtRH_HbT;filtLH_GCaMP;filtRH_GCaMP]);
whisking_Yvals = 1.10*max(indecesMax)*ones(size(binWhiskers));
force_Yvals = 1.20*max(indecesMax)*ones(size(binForce));
LPad_Yvals = 1.30*max(indecesMax)*ones(size(LPadSol));
RPad_Yvals = 1.30*max(indecesMax)*ones(size(RPadSol));
Aud_Yvals = 1.30*max(indecesMax)*ones(size(AudSol));
forceInds = binForce.*force_Yvals;
whiskInds = binWhiskers.*whisking_Yvals;
% set force indeces
for x = 1:length(forceInds)
    if forceInds(1,x) == 0
        forceInds(1,x) = NaN;
    end
end
% set whisk indeces
for x = 1:length(whiskInds)
    if whiskInds(1,x) == 0
        whiskInds(1,x) = NaN;
    end
end
%% Figure
figHandle = figure;
% force sensor and EMG
ax1 = subplot(7,1,1);
fileID2 = strrep(fileID,'_',' ');
p1 = plot((1:length(filtForceSensor))/ProcData.notes.dsFs,filtForceSensor,'color',colors('north texas green'),'LineWidth',1);
title([animalID ' IOS behavioral characterization and CBV dynamics for ' fileID2])
ylabel('Force Sensor (Volts)')
xlim([0,ProcData.notes.trialDuration_sec])
yyaxis right
p2 = plot((1:length(EMG))/ProcData.notes.dsFs,EMG,'color',colors('deep carrot orange'),'LineWidth',1);
ylabel('EMG (Volts^2)')
xlim([0,ProcData.notes.trialDuration_sec])
legend([p1,p2],'force sensor','EMG')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% Whisker angle and heart rate
ax2 = subplot(7,1,2);
plot((1:length(filteredWhiskerAngle))/ProcData.notes.dsFs,-filteredWhiskerAngle,'color',colors('dark pink'),'LineWidth',1);
ylabel('Angle (deg)')
xlim([0,ProcData.notes.trialDuration_sec])
ylim([-20,60])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% CBV and behavioral indeces
ax3 = subplot(7,1,3);
s1 = scatter((1:length(binForce))/ProcData.notes.dsFs,forceInds,'.','MarkerEdgeColor',colors('north texas green'));
hold on
s2 = scatter((1:length(binWhiskers))/ProcData.notes.dsFs,whiskInds,'.','MarkerEdgeColor',colors('dark pink'));
p5 = plot((1:length(filtLH_HbT))/ProcData.notes.dsFs,filtLH_HbT,'color',colors('dark candy apple red'),'LineWidth',1);
p6 = plot((1:length(filtRH_HbT))/ProcData.notes.dsFs,filtRH_HbT,'color',colors('sapphire'),'LineWidth',1);
ylabel('CBV zScore')
legend([p5,p6,s1,s2],'LH CBV','RH CBV','movement','whisking')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% GCaMP and behavioral indeces
ax4 = subplot(7,1,4);
p7 = plot((1:length(filtLH_GCaMP))/ProcData.notes.dsFs,filtLH_GCaMP,'color',colors('electric purple'),'LineWidth',1);
hold on
p8 = plot((1:length(filtRH_GCaMP))/ProcData.notes.dsFs,filtRH_GCaMP,'color',colors('vegas gold'),'LineWidth',1);
ylabel('GCaMP7s zScore')
legend([p7,p8],'LH GCaMP7s','RH GCaMP7s')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% Left cortical electrode spectrogram
ax5 = subplot(7,1,5);
Semilog_ImageSC(T,F,cortical_LHnormS,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)')
caxis([-100,100])
ylabel('Frequency (Hz)')
set(gca,'Yticklabel','10^1')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
yyaxis right
ylabel('Left cortical LFP')
set(gca,'Yticklabel', [])
% Right cortical electrode spectrogram
ax6 = subplot(7,1,6);
Semilog_ImageSC(T,F,cortical_RHnormS,'y')
axis xy
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)')
caxis([-100,100])
ylabel('Frequency (Hz)')
set(gca,'Yticklabel','10^1')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
yyaxis right
ylabel('Right cortical LFP')
set(gca,'Yticklabel',[])
% Hippocampal electrode spectrogram
ax7 = subplot(7,1,7);
Semilog_ImageSC(T,F,hippocampusNormS,'y')
c7 = colorbar;
ylabel(c7,'\DeltaP/P (%)')
caxis([-100,100])
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'box','off')
yyaxis right
ylabel('Hippocampal LFP')
set(gca,'Yticklabel',[])
% Axes properties
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'x')
ax1Pos = get(ax1,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
ax5Pos(3:4) = ax1Pos(3:4);
ax6Pos(3:4) = ax1Pos(3:4);
ax7Pos(3:4) = ax1Pos(3:4);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
set(ax7,'position',ax7Pos);
%% save the file to directory.
if strcmp(saveFigs,'y') == true
    [pathstr,~,~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Single Trial Figures/'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(figHandle,[dirpath animalID '_' fileID '_SingleTrialFig']);
end

end
