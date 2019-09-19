function [singleTrialFig] = GenerateSingleFigures_IOS(procDataFileIDs, RestingBaselines, baselineType, saveFigs)
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
%   Last Revised: June 30th, 2019    
%________________________________________________________________________________________________________________________

for a = 1:size(procDataFileIDs, 1)

    procDataFile = procDataFileIDs(a,:);
    load(procDataFile)
    disp(['Creating single trial summary figure ' num2str(a) ' of ' num2str(size(procDataFileIDs,1)) '...']); disp(' ')
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
    LH_CBV = ProcData.data.CBV.LH;
    normLH_CBV = (LH_CBV - RestingBaselines.(baselineType).CBV.LH.(strDay))./(RestingBaselines.(baselineType).CBV.LH.(strDay));
    filtLH_CBV = filtfilt(D, C, normLH_CBV)*100;
    
    RH_CBV = ProcData.data.CBV.RH;
    normRH_CBV = (RH_CBV - RestingBaselines.(baselineType).CBV.RH.(strDay))./(RestingBaselines.(baselineType).CBV.RH.(strDay));
    filtRH_CBV = filtfilt(D, C, normRH_CBV)*100;
    
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
    figTemplet = tight_subplot_IOS(6,1,[.005 .005],[.045 .04],[.025 .025]);
    
    % Force sensor and EMG
    axes(figTemplet(1));
    fileID2 = strrep(fileID, '_', ' ');
    plot((1:length(filtForceSensor))/ProcData.notes.dsFs, filtForceSensor, 'color', colors_IOS('sapphire'))
    title([animalID ' IOS behavioral characterization and CBV dynamics for ' fileID2])
    ylabel('Force Sensor (Volts)')
    xlim([0 ProcData.notes.trialDuration_sec])
    yyaxis right
    plot((1:length(EMG))/ProcData.notes.dsFs, EMG, 'color', colors_IOS('deep carrot orange'))
    ylabel('EMG (Volts^2)')
    legend('Force sensor', 'EMG')
    xlim([0 ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0, 0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    
    % Whisker angle and heart rate
    axes(figTemplet(2));
    plot((1:length(filteredWhiskerAngle))/ProcData.notes.dsFs, -filteredWhiskerAngle, 'color', colors_IOS('electric purple'))
    ylabel('Angle (deg)')
    xlim([0 ProcData.notes.trialDuration_sec])
    yyaxis right
    plot((1:length(heartRate)), heartRate, 'color', colors_IOS('dark sea green'), 'LineWidth', 2)
    ylabel('Heart Rate (Hz)')
    legend('Whisker angle', 'Heart rate')
    set(gca,'TickLength',[0, 0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    
    % CBV and behavioral indeces
    axes(figTemplet(3));
    plot((1:length(filtLH_CBV))/ProcData.notes.CBVCamSamplingRate, filtLH_CBV, 'color', colors_IOS('dark candy apple red'))
    hold on;
    plot((1:length(filtRH_CBV))/ProcData.notes.CBVCamSamplingRate, filtRH_CBV, 'color', colors_IOS('rich black'))
    
    scatter((1:length(binForce))/ProcData.notes.dsFs, forceInds, '.', 'MarkerEdgeColor', colors_IOS('sapphire'));
    scatter((1:length(binWhiskers))/ProcData.notes.dsFs, whiskInds, '.', 'MarkerEdgeColor', colors_IOS('electric purple'));
    scatter(LPadSol, LPad_Yvals, 'v', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c');
    scatter(RPadSol, RPad_Yvals, 'v', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'm');
    scatter(AudSol, Aud_Yvals, 'v', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g');
    
    ylabel('% change (\DeltaR/R)')
    legend('LH CBV', 'RH CBV', 'Movement', 'Whisking', 'LPadSol', 'RPadSol', 'AudSol')
    xlim([0 ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0, 0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    
    % Left cortical electrode spectrogram
    axes(figTemplet(4));
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
    axes(figTemplet(5));
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
    axes(figTemplet(6));
    semilog_imagesc_IOS(T, F, hippocampusNormS, 'y')
    caxis([-0.5 0.75])
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    set(gca,'Yticklabel', '10^1')
    set(gca,'Xticklabel', [0 100 200 300 400 500 600 700 800 900])
    xlim([0 ProcData.notes.trialDuration_sec])
    set(gca,'TickLength',[0, 0])
    set(gca,'box','off')
    yyaxis right
    ylabel('Hippocampal LFP')
    set(gca,'Yticklabel', [])
    
    pause(1)
    
    if strcmp(saveFigs, 'y') == true
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

end