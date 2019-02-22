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
    
    binWhiskers = MergedData.Data.binWhisker_Angle;
    binForce = MergedData.Data.binForce_Sensor_M;

    %% CBV data - normalize and then lowpass filer
    Vessel_Diameter = MergedData.Data.Vessel_Diameter;
    normVessel_Diameter = (Vessel_Diameter - RestingBaselines.(vesselID).(strDay).Vessel_Diameter.baseLine) ./ (RestingBaselines.(vesselID).(strDay).Vessel_Diameter.baseLine);
    [D, C] = butter(4, 2 / (6 / 2), 'low');
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
    ax1 = subplot(3,1,1);
    plot((1:length(filteredForceSensor))/30, filteredForceSensor, 'color', colors('dark candy apple red'))
    xlim([0 290])
    title({[animalID ' Two-photon behavioral characterization and vessel ' vesselID ' diameter changes for ' fileID], 'Force Sensor'})
    xlabel('Time (sec)')
    ylabel('Force Sensor (Volts)')
    
    ax2 = subplot(3,1,2);
    plot((1:length(filteredWhiskerAngle))/30, filteredWhiskerAngle, 'color', colors('brandeis blue'))
    xlim([0 290])
    title('Whisker Angle')
    xlabel('Time (sec)')
    ylabel('Angle (deg)')
  
    ax3 = subplot(3,1,3);
    plot((1:length(filteredVessel_Diameter))/6.06, filteredVessel_Diameter, 'k')
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
    scatter((1:length(binForce))/30, forceInds, 'filled', 'MarkerFaceColor', colors('dark candy apple red'));
    scatter((1:length(binWhiskers))/30, whiskInds, 'filled', 'MarkerFaceColor', colors('brandeis blue'));
    xlim([0 290])
    ylim([(min(filteredVessel_Diameter))-0.1 (max(filteredVessel_Diameter))*1.3])
    title('Vessel Diameter')
    xlabel('Time (sec)')
    ylabel('Percentage Change (Diameter)')
    legend('Vessel Diameter', 'Movement events', 'Whisking events')

    ax4 = subplot(4,1,4);
    imagesc(T,F,S_Norm)
    axis xy
    colorbar
    caxis([-0.5 1])
    linkaxes([ax1 ax2 ax3 ax4], 'x')
    title('Hippocampal LFP Spectrogram')
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')

    %% Save the file to directory.
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Single Trial Figures/'];

    if ~exist(dirpath, 'dir') 
        mkdir(dirpath); 
    end

    savefig(singleTrialFig, [dirpath animalID '_' vesselID '_' fileID '_SingleTrialFig']);
    close all
end

end
