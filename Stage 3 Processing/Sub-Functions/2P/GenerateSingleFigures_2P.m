function [] = GenerateSingleFigures_2P(mergedDataFiles, RestingBaselines)

for a = 1:size(mergedDataFiles, 1)

    mergedDataFile = mergedDataFiles(a,:);
    load(mergedDataFile)
    
    [animalID, fileDate, fileID, vesselID, imageID] = GetFileInfo2_2P(mergedDataFile);
    strDay = ConvertDate_2P(fileDate);
    
    %% BLOCK PURPOSE: Filter the whisker angle and identify the solenoid timing and location.
    % Setup butterworth filter coefficients for a 10 Hz lowpass based on the sampling rate (30 Hz).
    [B, A] = butter(4, 10/(MergedData.notes.dsFs/2), 'low');
    filteredWhiskerAngle = filtfilt(B, A, MergedData.data.whiskerAngle);
    filtForceSensor = filtfilt(B, A, MergedData.data.forceSensorM);
    binWhiskers = MergedData.data.binWhiskerAngle;
    binForce = MergedData.data.binForceSensorM;
    
    %% CBV data - normalize and then lowpass filer
    % Setup butterworth filter coefficients for a 1 Hz lowpass based on the sampling rate (20 Hz).
    [D, C] = butter(4, 1/(MergedData.notes.p2Fs/2), 'low');
    vesselDiameter = MergedData.data.vesselDiameter;
    normVesselDiameter = (vesselDiameter - RestingBaselines.(vesselID).(strDay).vesselDiameter.baseLine)./(RestingBaselines.(vesselID).(strDay).vesselDiameter.baseLine);
    filtVesselDiameter = (filtfilt(D, C, normVesselDiameter))*100;
    
    %% Normalized neural spectrogram
    specDataFile = [animalID '_' vesselID '_' fileID '_' imageID '_SpecData.mat'];
    load(specDataFile, '-mat');
    normS = SpecData.fiveSec.normS;
    T = SpecData.fiveSec.T;
    F = SpecData.fiveSec.F;
    
    %% Yvals for behavior Indices
    whisking_YVals = 1.10*max(detrend(filtVesselDiameter, 'constant'))*ones(size(binWhiskers));
    force_YVals = 1.20*max(detrend(filtVesselDiameter, 'constant'))*ones(size(binForce));
    
    %% Figure
    singleTrialFig = figure;
    ax1 = subplot(4,1,1);
    plot((1:length(filtForceSensor))/MergedData.notes.dsFs, filtForceSensor, 'color', colors('sapphire'))
    title({[animalID ' Two-photon behavioral characterization and vessel ' vesselID ' diameter changes for ' fileID], 'Force sensor and whisker angle'})
    xlabel('Time (sec)')
    ylabel('Force Sensor (Volts)')
    xlim([0 MergedData.notes.trialDuration_Sec])  
    yyaxis right
    plot((1:length(filteredWhiskerAngle))/MergedData.notes.dsFs, -filteredWhiskerAngle, 'color', colors('ash grey'))
    ylabel('Angle (deg)')
    legend('Force sensor', 'Whisker angle')
    xlim([0 MergedData.notes.trialDuration_Sec])

    ax2 = subplot(4,1,2:3);
    plot((1:length(filtVesselDiameter))/MergedData.notes.p2Fs, detrend(filtVesselDiameter, 'constant'), 'color', colors('dark candy apple red'))
    hold on;
    whiskInds = binWhiskers.*whisking_YVals;
    forceInds = binForce.*force_YVals;
    for x = 1:length(whiskInds)
        if whiskInds(1, x) == 0
            whiskInds(1, x) = NaN;
        end
        
        if forceInds(1, x) == 0
            forceInds(1, x) = NaN;
        end
    end
    scatter((1:length(binForce))/MergedData.notes.dsFs, forceInds, '.', 'MarkerEdgeColor', colors('rich black'));
    scatter((1:length(binWhiskers))/MergedData.notes.dsFs, whiskInds, '.', 'MarkerEdgeColor', colors('sapphire'));
    title('Vessel diameter in response to behaviorial events')
    xlabel('Time (sec)')
    ylabel('% change (diameter)')
    legend('Vessel diameter', 'Binarized movement events', 'binarized whisking events')
    xlim([0 MergedData.notes.trialDuration_Sec])
    
    ax3 = subplot(4,1,4);
    imagesc(T,F,normS)
    axis xy
    colorbar
    caxis([-0.5 0.75])
    linkaxes([ax1 ax2 ax3], 'x')
    title('Hippocampal (LFP) spectrogram')
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    xlim([0 MergedData.notes.trialDuration_Sec])
    pause(1)
    
    %% Save the file to directory.
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Single Trial Figures/'];

    if ~exist(dirpath, 'dir') 
        mkdir(dirpath); 
    end

    savefig(singleTrialFig, [dirpath animalID '_' vesselID '_' fileID '_SingleTrialFig']);
    close(singleTrialFig)
end

end