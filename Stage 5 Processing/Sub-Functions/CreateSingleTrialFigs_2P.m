function CreateSingleTrialFigs_2P(mergedDataFiles, RestingBaselines, SpectrogramData)
%___________________________________________________________________________________________________
% Written by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%___________________________________________________________________________________________________
%
%   Purpose:
%___________________________________________________________________________________________________
%
%   Inputs: 
%          
%
%   Outputs: 
%___________________________________________________________________________________________________

for f = 1:size(mergedDataFiles, 1)
    mergedDataFile = mergedDataFiles(f, :);
    disp(['Analyzing single trial figure ' num2str(f) ' of ' num2str(size(mergedDataFiles, 1)) '...']); disp(' ');
    [animalID, fileDate, fileID, vesselID] = GetFileInfo_2P(mergedDataFile);
    load(mergedDataFile)
    strDay = ConvertDate(fileDate);

    %% BLOCK PURPOSE: Filter the whisker angle and identify the solenoid timing and location.
    % Setup butterworth filter coefficients for a 10 Hz lowpass based on the sampling rate (150 Hz).
    [B, A] = butter(4, 10 / (30/2), 'low');
    filteredWhiskerAngle = filtfilt(B, A, MergedData.Data.Whisker_Angle);
    filteredForceSensor = filtfilt(B, A, MergedData.Data.Force_Sensor_M);
    filteredEMG = filtfilt(B, A, MergedData.Data.EMG);
    binWhiskers = MergedData.Data.binWhisker_Angle;
    binForce = MergedData.Data.binForce_Sensor_M;

    %% CBV data - normalize and then lowpass filer
    Vessel_Diameter = MergedData.Data.Vessel_Diameter;
    normVessel_Diameter = (Vessel_Diameter - RestingBaselines.(vesselID).(strDay).Vessel_Diameter.baseLine) ./ (RestingBaselines.(vesselID).(strDay).Vessel_Diameter.baseLine);
    [D, C] = butter(4, 1/(20/2), 'low');
    filteredVessel_Diameter = (filtfilt(D, C, normVessel_Diameter))*100;

    %% Neural spectrograms
    S = SpectrogramData.FiveSec.S{f, 1};
    S_Norm = SpectrogramData.FiveSec.S_Norm{f, 1};
    T = SpectrogramData.FiveSec.T{f, 1};
    F = SpectrogramData.FiveSec.F{f, 1};
    
    %% Yvals for behavior Indices
    whisking_YVals = 1.10*max(filteredVessel_Diameter)*ones(size(binWhiskers));
    force_YVals = 1.20*max(filteredVessel_Diameter)*ones(size(binForce));

    %% Figure
    singleTrialFig = figure;
    ax1 = subplot(4,1,1);
    plot((1:length(filteredForceSensor))/30, filteredForceSensor, 'color', colors('sapphire'))
    xlim([0 280])
    title({[animalID ' Two-photon behavioral characterization and vessel ' vesselID ' diameter changes for ' fileID], 'Force sensor and EMG'})
    xlabel('Time (sec)')
    ylabel('Force Sensor (Volts)')
    
    yyaxis right
    plot((1:length(filteredEMG))/30, filteredEMG, 'color', colors('harvest gold'));
    ylabel('EMG (Volts)')
    legend('Force Sensor', 'EMG')
    
    ax2 = subplot(4,1,2:3);
    yyaxis right
    plot((1:length(filteredWhiskerAngle))/30, -filteredWhiskerAngle, 'color', colors('ash grey'))
    xlim([0 280])
    ylabel('Angle (deg)')
    ylim([-10 60])
    
    yyaxis left
    plot((1:length(filteredVessel_Diameter))/20, filteredVessel_Diameter, 'color', colors('dark candy apple red'))
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
    scatter((1:length(binForce))/30, forceInds, '.', 'MarkerEdgeColor', colors('rich black'));
    scatter((1:length(binWhiskers))/30, whiskInds, '.', 'MarkerEdgeColor', colors('sapphire'));
    xlim([0 280])
    ylim([(min(filteredVessel_Diameter))-0.1 (max(filteredVessel_Diameter))*1.3])
    title('Vessel diameter in response to behavior events')
    xlabel('Time (sec)')
    ylabel('% change (diameter)')
    legend('Whisker angle', 'Vessel diameter', 'Binarized movement events', 'binarized whisking events')
    
    ax3 = subplot(4,1,4);
    imagesc(T,F,S_Norm)
    axis xy
    caxis([-1 2])
    linkaxes([ax1 ax2 ax3], 'x')
    title('Hippocampal (LFP) spectrogram, caxis([-1 2])')
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')

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
