function [] = GenerateSingleFigures_IOS(procDataFileIDs, RestingBaselines)

for a = 2:size(procDataFileIDs, 1)

    procDataFile = procDataFileIDs(a,:);
    load(procDataFile)
    
    [animalID, fileDate, fileID] = GetFileInfo_IOS(procDataFile);
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
    EMG = ProcData.data.EMG;
    
    % Heart Rate
    heartRate = ProcData.data.heartRate;
    
    % Solenoids
    LPadSol = ProcData.data.solenoids.LPadSol;
    RPadSol = ProcData.data.solenoids.RPadSol;
    AudSol = ProcData.data.solenoids.AudSol;

    
    %% CBV data - normalize and then lowpass filer
    % Setup butterworth filter coefficients for a 1 Hz lowpass based on the sampling rate (20 Hz).
    [D, C] = butter(4, 1/(ProcData.notes.CBVCamSamplingRate/2), 'low');
    LH_CBV = ProcData.data.CBV.LH;
    normLH_CBV = (LH_CBV - RestingBaselines.CBV.LH.(strDay))./(RestingBaselines.CBV.LH.(strDay));
    filtLH_CBV = (filtfilt(D, C, normLH_CBV))*100;
    
    RH_CBV = ProcData.data.CBV.RH;
    normRH_CBV = (RH_CBV - RestingBaselines.CBV.RH.(strDay))./(RestingBaselines.CBV.RH.(strDay));
    filtRH_CBV = (filtfilt(D, C, normRH_CBV))*100;
    
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
    singleTrialFig = figure;
    % Force sensor and EMG
    ax1 = subplot(7,1,1);
    plot((1:length(filtForceSensor))/ProcData.notes.dsFs, filtForceSensor, 'color', colors_IOS('sapphire'))
    title({[animalID ' IOS behavioral characterization and CBV dynamics for ' fileID], 'Force sensor and EMG'})
    ylabel('Force Sensor (Volts)')
    xlim([0 ProcData.notes.trialDuration_sec])
    yyaxis right
    plot((1:length(EMG))/ProcData.notes.dsFs, log10(EMG), 'color', colors_IOS('deep carrot orange'))
    ylabel('EMG (Volts^2)')
    legend('Force sensor', 'EMG')
    xlim([0 ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0, 0])
    set(gca,'box','off')
    axis tight
    
    % Whisker angle and heart rate
    ax2 = subplot(7,1,2);
    plot((1:length(filteredWhiskerAngle))/ProcData.notes.dsFs, -filteredWhiskerAngle, 'color', colors_IOS('electric purple'))
    title('Whisker angle and heart rate')
    ylabel('Angle (deg)')
    xlim([0 ProcData.notes.trialDuration_sec])
    yyaxis right
    plot((1:length(heartRate)), heartRate, 'color', colors_IOS('dark sea green'), 'LineWidth', 2)
    ylabel('Heart Rate (Hz)')
    legend('Whisker angle', 'Heart rate')
    set(gca,'TickLength',[0, 0])
    set(gca,'box','off')
    axis tight
    
    % CBV and behavioral indeces
    ax3 = subplot(7,1,3:4);
    plot((1:length(filtLH_CBV))/ProcData.notes.CBVCamSamplingRate, filtLH_CBV, 'color', colors_IOS('dark candy apple red'))
    hold on;
    plot((1:length(filtRH_CBV))/ProcData.notes.CBVCamSamplingRate, filtRH_CBV, 'color', colors_IOS('rich black'))
    
    scatter((1:length(binForce))/ProcData.notes.dsFs, forceInds, '.', 'MarkerEdgeColor', colors_IOS('sapphire'));
    scatter((1:length(binWhiskers))/ProcData.notes.dsFs, whiskInds, '.', 'MarkerEdgeColor', colors_IOS('electric purple'));
    scatter(LPadSol, LPad_Yvals, 'v', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c');
    scatter(RPadSol, RPad_Yvals, 'v', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'm');
    scatter(AudSol, Aud_Yvals, 'v', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g');
    
    title('Cerebral blood volume in response to behaviorial events')
    ylabel('% change (\DeltaR/R)')
    legend('LH CBV', 'RH CBV', 'Movement', 'Whisking', 'LPadSol', 'RPadSol', 'AudSol')
    xlim([0 ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0, 0])
    set(gca,'box','off')
    axis tight
    
    % Left cortical electrode spectrogram
    ax4 = subplot(7,1,5);
    imagesc(T, F, cortical_LHnormS)
    axis xy
    caxis([-1 2])
    title('LH cortical (LFP) spectrogram')
    ylabel('Frequency (Hz)')
    xlim([0 ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0, 0])
    set(gca,'box','off')
    
    % Right cortical electrode spectrogram
    ax5 = subplot(7,1,6);
    imagesc(T, F, cortical_RHnormS)
    axis xy
    caxis([-1 2])
    title('RH cortical (LFP) spectrogram')
    ylabel('Frequency (Hz)')
    xlim([0 ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0, 0])
    set(gca,'box','off')
    
    % Hippocampal electrode spectrogram
    ax6 = subplot(7,1,7);
    imagesc(T, F, hippocampusNormS)
    axis xy
    caxis([-0.5 0.75])
    title('Hippocampal (LFP) spectrogram')
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    xlim([0 ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0, 0])
    set(gca,'box','off')
    pause(1)
    
    linkaxes([ax1 ax2 ax3 ax4 ax5 ax6], 'x');
    
    %% Save the file to directory.
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Single Trial Figures/'];

    if ~exist(dirpath, 'dir') 
        mkdir(dirpath); 
    end

    savefig(singleTrialFig, [dirpath animalID '_' fileID '_SingleTrialFig']);
    close(singleTrialFig)
end

end