function [figHandle] = GenerateSingleFigures_SVM(procDataFileID, RestingBaselines, baselineType)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: 
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%
%   Last Revised: July 27th, 2019
%________________________________________________________________________________________________________________________

load(procDataFileID)

[animalID, fileDate, fileID] = GetFileInfo_IOS(procDataFileID);
strDay = ConvertDate_IOS(fileDate);

%% BLOCK PURPOSE: Behavior
% Setup butterworth filter coefficients for a 10 Hz lowpass based on the sampling rate (30 Hz).
[B, A] = butter(4, 10/(ProcData.notes.dsFs/2), 'low');

% Whiskers
filteredWhiskerAngle = filtfilt(B, A, ProcData.data.whiskerAngle);
binWhiskers = ProcData.data.binWhiskerAngle;

% Force Sensor
filtForceSensor = filtfilt(B, A, ProcData.data.forceSensor);
binForce = ProcData.data.binForceSensor;

% EMG
EMG = ProcData.data.EMG.emg;

% Heart Rate
heartRate = ProcData.data.heartRate;

% Solenoids
LPadSol = ProcData.data.solenoids.LPadSol;
RPadSol = ProcData.data.solenoids.RPadSol;
AudSol = ProcData.data.solenoids.AudSol;

%% CBV data - normalize and then lowpass filer
% Setup butterworth filter coefficients for a 1 Hz lowpass based on the sampling rate (20 Hz).
[D, C] = butter(4, 1/(ProcData.notes.CBVCamSamplingRate/2), 'low');
LH_CBV = ProcData.data.CBV.adjLH;
RH_CBV = ProcData.data.CBV.adjRH;
normLH_CBV = (LH_CBV - RestingBaselines.(baselineType).CBV.adjLH.(strDay))./(RestingBaselines.(baselineType).CBV.adjLH.(strDay));
filtLH_CBV = (filtfilt(D, C, normLH_CBV))*100;
normRH_CBV = (RH_CBV - RestingBaselines.(baselineType).CBV.adjRH.(strDay))./(RestingBaselines.(baselineType).CBV.adjRH.(strDay));
filtRH_CBV = (filtfilt(D, C, normRH_CBV))*100;

%% Neural data
deltaPower_LH = ProcData.data.cortical_LH.deltaBandPower;
normDeltaPower_LH = (deltaPower_LH - RestingBaselines.(baselineType).cortical_LH.deltaBandPower.(strDay))./(RestingBaselines.(baselineType).cortical_LH.deltaBandPower.(strDay));
filtDeltaPower_LH = filtfilt(D, C, normDeltaPower_LH);

deltaPower_RH = ProcData.data.cortical_RH.deltaBandPower;
normDeltaPower_RH = (deltaPower_RH - RestingBaselines.(baselineType).cortical_RH.deltaBandPower.(strDay))./(RestingBaselines.(baselineType).cortical_RH.deltaBandPower.(strDay));
filtDeltaPower_RH = filtfilt(D, C, normDeltaPower_RH);

gammaPower_LH = ProcData.data.cortical_LH.gammaBandPower;
normGammaPower_LH = (gammaPower_LH - RestingBaselines.(baselineType).cortical_LH.gammaBandPower.(strDay))./(RestingBaselines.(baselineType).cortical_LH.gammaBandPower.(strDay));
filtGammaPower_LH = filtfilt(D, C, normGammaPower_LH);

gammaPower_RH = ProcData.data.cortical_RH.gammaBandPower;
normGammaPower_RH = (gammaPower_RH - RestingBaselines.(baselineType).cortical_RH.gammaBandPower.(strDay))./(RestingBaselines.(baselineType).cortical_RH.gammaBandPower.(strDay));
filtGammaPower_RH = filtfilt(D, C, normGammaPower_RH);

thetaPower_Hipp = ProcData.data.hippocampus.thetaBandPower;
normThetaPower_Hipp = (thetaPower_Hipp - RestingBaselines.(baselineType).hippocampus.thetaBandPower.(strDay))./(RestingBaselines.(baselineType).hippocampus.thetaBandPower.(strDay));
filtThetaPower_Hipp = filtfilt(D, C, normThetaPower_Hipp);

%% Normalized neural spectrogram
specDataFile = [animalID '_' fileID '_SpecData.mat'];
load(specDataFile, '-mat');
cortical_LHnormS = SpecData.cortical_LH.fiveSec.normS;
cortical_RHnormS = SpecData.cortical_RH.fiveSec.normS;
hippocampusNormS = SpecData.hippocampus.fiveSec.normS;
T = SpecData.cortical_LH.fiveSec.T;
F = SpecData.cortical_LH.fiveSec.F;

%% Yvals for behavior Indices
if max(filtLH_CBV) >= max(filtRH_CBV)
    whisking_Yvals = 1.10*max(filtLH_CBV)*ones(size(binWhiskers));
    force_Yvals = 1.20*max(filtLH_CBV)*ones(size(binForce));
    LPad_Yvals = 1.30*max(filtLH_CBV)*ones(size(LPadSol));
    RPad_Yvals = 1.30*max(filtLH_CBV)*ones(size(RPadSol));
    Aud_Yvals = 1.30*max(filtLH_CBV)*ones(size(AudSol));
else
    whisking_Yvals = 1.10*max(filtRH_CBV)*ones(size(binWhiskers));
    force_Yvals = 1.20*max(filtRH_CBV)*ones(size(binForce));
    LPad_Yvals = 1.30*max(filtRH_CBV)*ones(size(LPadSol));
    RPad_Yvals = 1.30*max(filtRH_CBV)*ones(size(RPadSol));
    Aud_Yvals = 1.30*max(filtRH_CBV)*ones(size(AudSol));
end

whiskInds = binWhiskers.*whisking_Yvals;
forceInds = binForce.*force_Yvals;
for x = 1:length(whiskInds)
    if whiskInds(1, x) == 0
        whiskInds(1, x) = NaN;
    end
    
    if forceInds(1, x) == 0
        forceInds(1, x) = NaN;
    end
end

%% Figure
figHandle = figure;

% Force sensor and EMG
ax1 = subplot(6,1,1);
fileID2 = strrep(fileID, '_', ' ');
plot((1:length(filtForceSensor))/ProcData.notes.dsFs, filtForceSensor, 'color', colors_IOS('sapphire'))
title([animalID ' IOS behavioral characterization and CBV dynamics for ' fileID2])
ylabel('Force Sensor (Volts)')
xlim([0 ProcData.notes.trialDuration_sec])
yyaxis right
plot((1:length(EMG))/ProcData.notes.dsFs, EMG, 'color', colors_IOS('deep carrot orange'))
ylabel('EMG (Volts^2)')
xlim([0 ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0, 0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight

% Whisker angle and heart rate
ax2 = subplot(6,1,2);
plot((1:length(filteredWhiskerAngle))/ProcData.notes.dsFs, -filteredWhiskerAngle, 'color', colors_IOS('electric purple'))
ylabel('Angle (deg)')
xlim([0 ProcData.notes.trialDuration_sec])
ylim([-20 60])
yyaxis right
plot((1:length(heartRate)), heartRate, 'color', colors_IOS('dark sea green'), 'LineWidth', 2)
ylabel('Heart Rate (Hz)')
ylim([6 15])
set(gca,'TickLength',[0, 0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight

% CBV and behavioral indeces
ax3 = subplot(6,1,3);
plot((1:length(filtLH_CBV))/ProcData.notes.CBVCamSamplingRate, filtLH_CBV, 'color', colors_IOS('dark candy apple red'))
hold on;
plot((1:length(filtRH_CBV))/ProcData.notes.CBVCamSamplingRate, filtRH_CBV, 'color', colors_IOS('rich black'))

scatter((1:length(binForce))/ProcData.notes.dsFs, forceInds, '.', 'MarkerEdgeColor', colors_IOS('sapphire'));
scatter((1:length(binWhiskers))/ProcData.notes.dsFs, whiskInds, '.', 'MarkerEdgeColor', colors_IOS('electric purple'));
scatter(LPadSol, LPad_Yvals, 'v', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c');
scatter(RPadSol, RPad_Yvals, 'v', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'm');
scatter(AudSol, Aud_Yvals, 'v', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g');

ylabel('% change (\DeltaR/R)')
xlim([0 ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0, 0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight

% yyaxis right
% plot((1:length(filtHbtLH_CBV))/ProcData.notes.CBVCamSamplingRate, filtHbtLH_CBV, 'color', colors_IOS('coral red'))
% hold on;
% plot((1:length(filtHbtRH_CBV))/ProcData.notes.CBVCamSamplingRate, filtHbtRH_CBV, 'color', colors_IOS('battleship grey'))
% ylabel('\DeltaHbT')

xlim([0 ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0, 0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight

% % Neural bands
% ax4 = subplot(6,1,4);
% plot((1:length(filtDeltaPower_LH))/ProcData.notes.CBVCamSamplingRate, filtDeltaPower_LH, 'color', colors_IOS('dark candy apple red'))
% hold on
% plot((1:length(filtDeltaPower_RH))/ProcData.notes.CBVCamSamplingRate, filtDeltaPower_RH, 'color', colors_IOS('coral red'))
% plot((1:length(filtGammaPower_LH))/ProcData.notes.CBVCamSamplingRate, filtGammaPower_LH, 'color', colors_IOS('electric blue'))
% plot((1:length(filtGammaPower_RH))/ProcData.notes.CBVCamSamplingRate, filtGammaPower_RH, 'color', colors_IOS('sapphire'))
% plot((1:length(filtThetaPower_Hipp))/ProcData.notes.CBVCamSamplingRate, filtThetaPower_Hipp, 'color', colors_IOS('vegas gold'))
% ylabel('Normalized Power)')
% xlim([0 ProcData.notes.trialDuration_sec])
% set(gca,'TickLength',[0, 0])
% set(gca,'Xticklabel',[])
% set(gca,'box','off')
% axis tight

% Left cortical electrode spectrogram
ax4 = subplot(6,1,4);
semilog_imagesc_IOS(T, F, cortical_LHnormS, 'y')
axis xy
caxis([-1 2])
ylabel('Frequency (Hz)')
set(gca,'Yticklabel', '10^1')
xlim([0 ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0, 0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
yyaxis right
ylabel('Left cortical LFP')
set(gca,'Yticklabel', [])

% Right cortical electrode spectrogram
ax5 = subplot(6,1,5);
semilog_imagesc_IOS(T, F, cortical_RHnormS, 'y')
axis xy
caxis([-1 2])
ylabel('Frequency (Hz)')
set(gca,'Yticklabel', '10^1')
xlim([0 ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0, 0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
yyaxis right
ylabel('Right cortical LFP')
set(gca,'Yticklabel', [])

% Hippocampal electrode spectrogram
ax6 = subplot(6,1,6);
semilog_imagesc_IOS(T, F, hippocampusNormS, 'y')
caxis([-1 2])
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
xlim([0 ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0, 0])
set(gca,'box','off')
yyaxis right
ylabel('Hippocampal LFP')
set(gca,'Yticklabel', [])

linkaxes([ax1 ax2 ax3 ax4 ax5 ax6], 'x')

end
