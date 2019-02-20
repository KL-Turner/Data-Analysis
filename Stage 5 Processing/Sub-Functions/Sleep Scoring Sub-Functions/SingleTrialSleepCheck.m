function SingleTrialSleepCheck(animal, fileID, windowCamFile, ProcData, SleepEventData)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: 
%________________________________________________________________________________________________________________________
%
%   Inputs:
%          
%
%   Outputs: 
%________________________________________________________________________________________________________________________

%% BLOCK PURPOSE: Find the sleeping times from the SleepEventData for this particular trial
% This loop creates a logical that matches the inputed FileID with all potential sleeping trials
for fileNumber = 1:length(SleepEventData.FileIDs)                   % Loop through each sleeping event
    if strcmp(fileID, SleepEventData.FileIDs{fileNumber, 1})        % If the fileID matches the sleep event, create a 1
        sleepTimeLogical{fileNumber, 1} = 1;
    else
        sleepTimeLogical{fileNumber, 1} = 0;                        % Else create a 0
    end
end

% Now that we have a logical showing the sleeping events that match this specific trial, we want to pull out the sleep times
% for that specific trial
x = 1;                                              % This is the first index point for our new "Sleeping times" cell
for fileNumber = 1:length(sleepTimeLogical)         % Loop through the sleep logical
    if sleepTimeLogical{fileNumber, 1} == 1         % If the logical denotes a 1, aka the fileID matches the sleep logical
        sleepTimes{x, 1} = SleepEventData.BinTimes{fileNumber, 1};  % Pull out the associated bin times
        x = x + 1;                                                  % This adds additional times within the same trial down the new struct
    end
end

strDay = ConvertDate(fileID(1:6));  % Find the fileID's date

%% BLOCK PURPOSE: Filter the whisker angle and identify the solenoid timing and location.
% Setup butterworth filter coefficients for a 10 Hz lowpass based on the sampling rate (150 Hz).
[B, A] = butter(4, 10 / (150/2), 'low');
filteredWhiskerAngle = filtfilt(B, A, ProcData.Data.Behavior.whiskers);

% Link closely occuring whisker events for a scatterplot overlain on the LH/RH CBV plot.
whiskingThreshold = 0.1;
dCrit = [75, 0];
whiskerMoves = gt(abs(diff(filteredWhiskerAngle, 2)), whiskingThreshold);
[linkedWhiskerMovements] = LinkBinaryEvents(whiskerMoves, dCrit);
whiskingInds = find(linkedWhiskerMovements);

% Identify the solenoid times from the ProcData file.
solenoidLeftPad = floor(ProcData.Data.Sol.solenoidLeftPad);
solenoidRightPad = floor(ProcData.Data.Sol.solenoidRightPad);

%% CBV data loaded in the GetBilateralCBV function
baselineType = 'Sleep';
[LH_CBV, RH_CBV] = GetBilateralCBV(windowCamFile, baselineType);
timeVector = (1:length(LH_CBV)) / 30;

[~, LH_tr, ~, LH_HR] = FindHeartRate(ProcData.Data.CBV.LH, ProcData.Notes.CBVCamSamplingRate);
[~, ~, ~, RH_HR] = FindHeartRate(ProcData.Data.CBV.RH, ProcData.Notes.CBVCamSamplingRate);
HR = (LH_HR + RH_HR)/ 2;
[D, C] = butter(4, 2 / (ProcData.Notes.CBVCamSamplingRate / 2), 'low');
HeartRate = filtfilt(D, C, HR);

%% Neural Data
% Delta Power
[H, G] = butter(4, 1 / (ProcData.Notes.deltaBandSamplingRate / 2), 'low');
LH_Delta = ProcData.Data.Neuro.LH.deltaBand_Power;
RH_Delta = ProcData.Data.Neuro.RH.deltaBand_Power;

% Smooth with filter
LH_DeltaNeuro = filtfilt(H, G, LH_Delta); 
RH_DeltaNeuro = filtfilt(H, G, RH_Delta);

% Theta Power
LH_Theta = ProcData.Data.Neuro.LH.thetaBand_Power;
RH_Theta = ProcData.Data.Neuro.RH.thetaBand_Power;

% Smooth with filter
LH_ThetaNeuro = filtfilt(H, G, LH_Theta); 
RH_ThetaNeuro = filtfilt(H, G, RH_Theta);

%% Neural Normalization
DeltaBandBaselines = ls('*_DeltaBandBaselines.mat');
load(DeltaBandBaselines);

LH_baselineDelta = mean(DeltaBandBaselines.LH.whisk.(strDay).Means);
RH_baselineDelta = mean(DeltaBandBaselines.RH.whisk.(strDay).Means);

LH_NormDeltaNeuro = (LH_DeltaNeuro - LH_baselineDelta) / LH_baselineDelta; 
RH_NormDeltaNeuro = (RH_DeltaNeuro - RH_baselineDelta) / RH_baselineDelta; 

ThetaBandBaselines = ls('*_ThetaBandBaselines.mat');
load(ThetaBandBaselines);

LH_baselineTheta = mean(ThetaBandBaselines.LH.whisk.(strDay).Means);
RH_baselineTheta = mean(ThetaBandBaselines.RH.whisk.(strDay).Means);

LH_NormThetaNeuro = (LH_ThetaNeuro - LH_baselineTheta) / LH_baselineTheta; 
RH_NormThetaNeuro = (RH_ThetaNeuro - RH_baselineTheta) / RH_baselineTheta; 

%% Yvals for behavior Indices
if max(LH_CBV) >= max(RH_CBV)
    whisking_YVals = 1.20*max(LH_CBV)*ones(size(whiskingInds));
    sleeping_YVal = 1.20*max(LH_CBV);
else
    whisking_YVals = 1.20*max(RH_CBV)*ones(size(whiskingInds));
    sleeping_YVal = 1.20*max(RH_CBV);
end

solenoidLeftPad_YVals = 1.20*max(ProcData.Data.Behavior.whiskers)*ones(size(solenoidLeftPad));
solenoidRightPad_YVals = 1.20*max(ProcData.Data.Behavior.whiskers)*ones(size(solenoidRightPad));

%% Figure
singleSleepTrial = figure;
ax1 = subplot(5,1,1);
plot((1:length(ProcData.Data.Behavior.whiskers)) / ProcData.Notes.downsampledWhiskerSamplingRate,...
    -ProcData.Data.Behavior.whiskers, 'color', colors('ash grey'));
hold on;
plot((1:length(filteredWhiskerAngle))/ProcData.Notes.downsampledWhiskerSamplingRate,...
    -filteredWhiskerAngle, 'k');
% scatter(solenoidLeftPad, solenoidLeftPad_YVals, 'vb', 'MarkerFaceColor', 'k');
% scatter(solenoidRightPad, solenoidRightPad_YVals, 'vk', 'MarkerFaceColor', 'k');
axis tight
ylabel('Degrees');
yyaxis right
ylim([6 15]);
plot(LH_tr, HeartRate);
ylabel('Heart Rate (Hz)');
title([animal ' ' fileID ' Behavior']);
set(gca, 'Ticklength', [0 0])
legend('Unfiltered Whisker Angle', 'Filtered Whisker Angle', 'HeartRate', 'Location', 'NorthEast')

ax2 = subplot(5,1,2:3);
plot(timeVector, LH_CBV, 'k');
hold on;
plot(timeVector, RH_CBV, 'b');
for sleepT = 1:length(sleepTimes)
    scatter(sleepTimes{sleepT, 1}, (ones(1, length(sleepTimes{sleepT, 1})))*sleeping_YVal, 'm')
end
title('Normalized, Filtered, and Mean-Subtracted CBV Signals');
ylabel('Reflectance (%)')
axis tight
legend('LH CBV', 'RH CBV', 'Sleep', 'Location', 'NorthEast')
set(gca, 'Ticklength', [0 0])

ax4 = subplot(5,1,4);
plot((1:length(LH_Delta)) / ProcData.Notes.deltaBandSamplingRate, LH_NormDeltaNeuro, 'k')
hold on;
plot((1:length(RH_Delta)) / ProcData.Notes.deltaBandSamplingRate, RH_NormDeltaNeuro, 'c')
title('Delta Band Power')
ylabel('Power')
ylim([-5 30])
legend('LH Delta Power', 'RH Delta Power', 'Location', 'NorthEast')
set(gca, 'Ticklength', [0 0])

ax5 = subplot(5,1,5);
plot((1:length(LH_Theta)) / ProcData.Notes.thetaBandSamplingRate, LH_NormThetaNeuro, 'k')
hold on;
plot((1:length(RH_Theta)) / ProcData.Notes.thetaBandSamplingRate, RH_NormThetaNeuro, 'r')
title('Theta Band Power')
ylabel('Power')
xlabel('Time (sec)')
ylim([-5 30])
legend('LH Theta Power', 'RH Theta Power', 'Location', 'NorthEast')
set(gca, 'Ticklength', [0 0])

linkaxes([ax1 ax2 ax4 ax5], 'x')

%% Save the file to directory.
[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Single Trial Sleep Checks/'];

if ~exist(dirpath, 'dir') 
    mkdir(dirpath); 
end

savefig(singleSleepTrial, [dirpath animal '_' fileID '_SingleTrialSleepCheck']);

end
