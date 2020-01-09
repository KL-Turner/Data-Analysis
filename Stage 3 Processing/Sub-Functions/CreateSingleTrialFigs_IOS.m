
function CreateSingleTrialFigs(procDataFile, RestingBaselines, SpectrogramData)
%___________________________________________________________________________________________________
% Written by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%___________________________________________________________________________________________________
%
%   Purpose: Create a single figure for easy analysis of the whisker, CBV, and neural data for each
%            five minute trial.
%___________________________________________________________________________________________________
%
%   Inputs: Single procDataFile string that will be used to load the RawData and ProcData structures.
%          
%
%   Outputs: Output figure will be a 4x1 subplot that can be saved to directory under the animal ID.
%
%            (4,1,1): Whisker angle with an overlain scatter of the stimuli. Filtered
%                     signal is lowpass filtered at 10 Hz. Added HR.
%
%            (4,1,2:3): Biltaral LH and RH CBV signals, lowpass filtered at 2 Hz. Normalized by the
%            entire day's resting average.
%
%            (4,1,4): LH LFP normalized spectrogram (Bandpass filtered between 0.1 and 100 Hz).
%
%            (4,1,5): RH LFP normalized spectrogram (Bandpass filtered between 0.1 and 100 Hz).
%___________________________________________________________________________________________________

%% This function can be run through a loop of filenames, or called independently for convenience.
% If there is no filename (procDataFile) loaded in, prompt the user to load in a ProcData.mat file
% from the current directory.
if nargin == 0  
    procDataFile = uigetfile('*_ProcData.mat');
end

load(procDataFile); 

% Obtain the animal and file ID information from the raw/procDataFile. It can be either or.
[animal, ~, fileDate, fileID] = GetFileInfo(procDataFile);
strDay = ConvertDate(fileDate);

%% BLOCK PURPOSE: Filter the whisker angle and identify the solenoid timing and location.
% Setup butterworth filter coefficients for a 10 Hz lowpass based on the sampling rate (150 Hz).
[B, A] = butter(4, 10 / (150/2), 'low');
filteredWhiskerAngle = filtfilt(B, A, ProcData.Data.Behavior.whiskers);
whiskers = ProcData.Data.Behavior.whiskers;

% Link closely occuring whisker events for a scatterplot overlain on the LH/RH CBV plot.
whiskingThreshold = 0.1;
dCrit = [75, 0];
whiskerMoves = gt(abs(diff(filteredWhiskerAngle, 2)), whiskingThreshold);
[linkedWhiskerMovements] = LinkBinaryEvents(whiskerMoves, dCrit);
whiskingInds = find(linkedWhiskerMovements);

HeartRate = ProcData.Data.HeartRate.HR;
tR = ProcData.Data.HeartRate.tr;

% Identify the solenoid times from the ProcData file.
solenoidLeftPad = floor(ProcData.Data.Sol.solenoidLeftPad);
solenoidRightPad = floor(ProcData.Data.Sol.solenoidRightPad);
solenoidAuditory = floor(ProcData.Data.Sol.solenoidAuditory);

%% CBV data - normalize and then lowpass filer
LH_CBV = ProcData.Data.CBV.LH;
RH_CBV = ProcData.Data.CBV.RH;

normLH_CBV = (LH_CBV - RestingBaselines.CBV.LH.(strDay)) ./ (RestingBaselines.CBV.LH.(strDay));
normRH_CBV = (RH_CBV - RestingBaselines.CBV.RH.(strDay)) ./ (RestingBaselines.CBV.RH.(strDay));

[D, C] = butter(4, 2 / (30 / 2), 'low');
filtLH_CBV = filtfilt(D, C, normLH_CBV);
filtRH_CBV = filtfilt(D, C, normRH_CBV);

%% Neural spectrograms
for f = 1:length(SpectrogramData.LH.FileIDs)
    specFileID = SpectrogramData.LH.FileIDs{f};
    if strcmp(fileID, specFileID)
        LH_S_norm = SpectrogramData.LH.FiveSec.S_Norm{f};
        RH_S_norm = SpectrogramData.RH.FiveSec.S_Norm{f};
        F = SpectrogramData.LH.FiveSec.F{f};
        T = SpectrogramData.LH.FiveSec.T{f};
    end
end

%% Yvals for behavior Indices
whisking_YVals = 1.10*max(-whiskers)*ones(size(whiskingInds));
LH_solenoid_YVals = 1.20*max(-whiskers)*ones(size(solenoidLeftPad));
RH_solenoid_YVals = 1.20*max(-whiskers)*ones(size(solenoidRightPad));
Aud_solenoid_YVals = 1.20*max(-whiskers)*ones(size(solenoidAuditory));

%% Figure
singleTrialFig = figure;
ax1 = subplot(5,1,1);
plot((1:length(ProcData.Data.Behavior.whiskers)) / ProcData.Notes.downsampledWhiskerSamplingRate,...
    -ProcData.Data.Behavior.whiskers, 'k');
ylabel('Degrees');
hold on;

scatter((whiskingInds / ProcData.Notes.downsampledWhiskerSamplingRate), whisking_YVals, '.k');
scatter(solenoidLeftPad, LH_solenoid_YVals, 'vk', 'MarkerFaceColor', 'b');
scatter(solenoidRightPad, RH_solenoid_YVals, 'vk', 'MarkerFaceColor', 'g');
scatter(solenoidAuditory, Aud_solenoid_YVals, 'vk', 'MarkerFaceColor', 'r');

yyaxis right
plot(tR, HeartRate, 'm');
ylabel('Heart Rate (Hz)');
ylim([6 15]);

title({[animal ' ' fileID ' Single Trial'], 'Behavioral State'});
set(gca, 'Ticklength', [0 0])
legend('Whisker Angle', 'Whisking Events', 'Left Pad Sol', 'Right Pad Sol', 'Auditory Sol', 'Heart Rate', 'Location', 'NorthEast')

ax2 = subplot(5,1,2:3);
timeVector = (1:9000) ./ 30;
plot(timeVector, filtLH_CBV*100, 'k');
hold on;
plot(timeVector, filtRH_CBV*100, 'b');
title('Normalized & Filtered CBV Signals');
ylabel('Reflectance (%)')
axis tight
legend('LH CBV', 'RH CBV', 'Location', 'NorthEast')
set(gca, 'Ticklength', [0 0])

ax4 = subplot(5,1,4);
imagesc(T, F, LH_S_norm);
colormap parula
colorbar
caxis([-1 2])
xlim([0 300])
title('LH Spectrogram, 5 second sliding window with 5,9 tapers')
ylabel('Frequency (Hz)')
axis xy
set(gca, 'Ticklength', [0 0])

ax5 = subplot(5,1,5);
imagesc(T, F, RH_S_norm);
colormap parula
colorbar
caxis([-1 2])
xlim([0 300])
title('RH Spectrogram, 5 second sliding window with 5,9 tapers')
ylabel('Frequency (Hz)')
xlabel('Time (sec)')
axis xy
set(gca, 'Ticklength', [0 0])

linkaxes([ax1 ax2 ax4 ax5], 'x')

%% Save the file to directory.
[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Single Trial Figures/'];

if ~exist(dirpath, 'dir') 
    mkdir(dirpath); 
end

savefig(singleTrialFig, [dirpath animal '_' fileID '_SingleTrialFig']);

end
