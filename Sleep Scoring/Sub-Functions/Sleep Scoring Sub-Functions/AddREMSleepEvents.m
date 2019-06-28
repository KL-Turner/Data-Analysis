function [SleepData] = AddREMSleepEvents(animal, hem, RestingBaselines, SleepData, SpectrogramData)
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
%   Outputs:
%________________________________________________________________________________________________________________________

sleepFileList = unique(SleepData.NREM.FileIDs);
SleepData.REM = [];

for f = 1:length(sleepFileList)
    clear sleepTimes
    clear LH_deltaPower RH_deltaPower LH_thetaPower RH_thetaPower LH_gammaPower RH_gammaPower WhiskerAcceleration HeartRate CBV_LH CBV_RH ElectrodeCBV_LH ElectrodeCBV_RH BinTimes
    clear cellLH_DeltaPower cellRH_DeltaPower cellLH_ThetaPower cellRH_ThetaPower cellLH_GammaPower cellRH_GammaPower cellWhiskerAcceleration cellHeartRate cellLH_CBV cellRH_CBV cellLH_ElectrodeCBV cellRH_ElectrodeCBV cellBinTimes
    clear mat2CellLH_DeltaPower mat2CellRH_DeltaPower mat2CellLH_ThetaPower mat2CellRH_TheatPower mat2CellLH_GammaPower mat2CellRH_GammaPower mat2CellWhiskerAcceleration mat2CellHeartRate mat2CellLH_CBV mat2CellRH_CBV mat2cellLH_ElectrodeCBV mat2cellRH_ElectrodeCBV mat2CellBinTimes
    clear matLH_DeltaPower matRH_DeltaPower matLH_ThetaPower matRH_ThetaPower matLH_GammaPower matRH_GammaPower matWhiskerAcceleration matHeartRate matLH_CBV matRH_CBV matLH_ElectrodeCBV matRH_ElectrodeCBV matBinTimes
    
    sleepFile = sleepFileList{f};
    sleepFileID = ([animal '_' hem '_' sleepFile '_ProcData.mat']);
    load(sleepFileID);
    [~, ~, ~, fileDate] = GetFileInfo(sleepFileID);
    [strDay] = ConvertDate(fileDate);
    
    % This loop creates a logical that matches the inputed FileID with all potential sleeping trials
    for fileNumber = 1:length(SleepData.NREM.FileIDs)   % Loop through each sleeping event
        if strcmp(sleepFile, SleepData.NREM.FileIDs{fileNumber, 1})   % If the fileID matches the sleep event, create a 1
            sleepTimeLogical{fileNumber, 1} = 1;
        else
            sleepTimeLogical{fileNumber, 1} = 0;   % Else create a 0
        end
    end
    
    % Now that we have a logical showing the sleeping events that match this specific trial, we want to pull out the sleep times
    % for that specific trial
    x = 1;   % This is the first index point for our new "Sleeping times" cell
    for fileNumber = 1:length(sleepTimeLogical)   % Loop through the sleep logical
        if sleepTimeLogical{fileNumber, 1} == 1   % If the logical denotes a 1, aka the fileID matches the sleep logical
            sleepTimes{x, 1} = SleepData.NREM.BinTimes{fileNumber, 1};   % Pull out the associated bin times
            x = x + 1;   % This adds additional times within the same trial down the new struct
        end
    end
    
    %% Filter the whisker angle and identify the solenoid timing and location.
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
    
    HR = ProcData.Data.HeartRate.HR;
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
    if max(filtLH_CBV) >= max(filtRH_CBV)
        sleeping_YVal = 1.20*max(filtLH_CBV*100);
        scaleMaxVal = 1.3*max(filtLH_CBV*100);
    else
        sleeping_YVal = 1.20*max(filtRH_CBV*100);
        scaleMaxVal = 1.3*max(filtRH_CBV*100);
    end
    
    whisking_YVals = 1.10*max(-whiskers)*ones(size(whiskingInds));
    LH_solenoid_YVals = 1.20*max(-whiskers)*ones(size(solenoidLeftPad));
    RH_solenoid_YVals = 1.20*max(-whiskers)*ones(size(solenoidRightPad));
    Aud_solenoid_YVals = 1.20*max(-whiskers)*ones(size(solenoidAuditory));
    
    %% Figure
    singleTrialFig = figure;
    ax1 = subplot(5,1,1);
    plot((1:length(ProcData.Data.Behavior.whiskers)) / ProcData.Notes.downsampledWhiskerSamplingRate,...
        -ProcData.Data.Behavior.whiskers, 'k');
    hold on;
    
    scatter((whiskingInds / ProcData.Notes.downsampledWhiskerSamplingRate), whisking_YVals, '.k');
    scatter(solenoidLeftPad, LH_solenoid_YVals, 'vk', 'MarkerFaceColor', 'b');
    scatter(solenoidRightPad, RH_solenoid_YVals, 'vk', 'MarkerFaceColor', 'g');
    scatter(solenoidAuditory, Aud_solenoid_YVals, 'vk', 'MarkerFaceColor', 'r');
    ylabel('Degrees')
    ylim([min(-whiskers), 1.30*max(-whiskers)])
    
    yyaxis right
    plot(tR, HR, 'm');
    ylabel('Heart Rate (Hz)');
    ylim([6 15]);
    
    title({[animal ' ' fileID ' Single Trial'], 'Behavioral State'});
    set(gca, 'Ticklength', [0 0])
    legend('Whisker Angle', 'Whisking Events', 'Left Pad Sol', 'Right Pad Sol', 'Auditory Sol', 'Heart Rate', 'Location', 'NorthEast')
    
    ax2 = subplot(5,1,2:3);
    timeVector = (1:9000) ./ 30;
    plot(timeVector, filtLH_CBV*100, 'k');
    hold on;
    plot(timeVector, filtRH_CBV*100, 'Color', colors('electric purple'));
    for sleepT = 1:length(sleepTimes)
        scatter(sleepTimes{sleepT, 1}, (ones(1, length(sleepTimes{sleepT, 1})))*sleeping_YVal, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c')
    end
    title('Normalized & Filtered CBV Signals');
    ylabel('Reflectance (%)')
    legend('LH CBV', 'RH CBV', 'NREM Sleep', 'Location', 'NorthEast')
    set(gca, 'Ticklength', [0 0])
    axis tight;
    
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
    
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Figures/NREM Single Trial Figures/'];

    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end

    savefig(singleTrialFig, [dirpath animal '_' sleepFile '_NREMSingleTrialFig']);
    
    %% Input the starting and ending points of the REM sleep trial
    while(1)
        REMcase = input(['Does file ID: ' sleepFileID ' contain a period of REM sleep? (y/n): '], 's');
        disp(' ');
        if strcmp(REMcase, 'y')
            startingIndex = input('Input the lagging edge time for the starting bin (from 5:300 at 5 second intervals: ');
            disp(' ')
            binStart = startingIndex / 5;
            
            endingIndex = input('Input the lagging edge time for the ending bin (from 5:30 at 5 second intervals): ');
            disp(' ')
            binEnd = endingIndex / 5;
            
            close(singleTrialFig)
            binTimes = 5*(binStart:binEnd);
            
            indexCount = 1;
            for index = binStart:binEnd   % Loop through the length of Sleep Index, and pull out associated data
                LH_deltaPower{indexCount, 1} = ProcData.Sleep.Parameters.DeltaBand_Power.LH{index, 1};
                RH_deltaPower{indexCount, 1} = ProcData.Sleep.Parameters.DeltaBand_Power.RH{index, 1};
                LH_thetaPower{indexCount, 1} = ProcData.Sleep.Parameters.ThetaBand_Power.LH{index, 1};
                RH_thetaPower{indexCount, 1} = ProcData.Sleep.Parameters.ThetaBand_Power.RH{index, 1};
                LH_gammaPower{indexCount, 1} = ProcData.Sleep.Parameters.GammaBand_Power.LH{index, 1};
                RH_gammaPower{indexCount, 1} = ProcData.Sleep.Parameters.GammaBand_Power.RH{index, 1};
                WhiskerAcceleration{indexCount, 1} = ProcData.Sleep.Parameters.WhiskerAcceleration{index, 1};
                HeartRate{indexCount, 1} = ProcData.Sleep.Parameters.HeartRate{index, 1};
                CBV_LH{indexCount, 1} = ProcData.Sleep.Parameters.CBV.LH{index, 1};
                CBV_RH{indexCount, 1} = ProcData.Sleep.Parameters.CBV.RH{index, 1};
                ElectrodeCBV_LH{indexCount, 1} = ProcData.Sleep.Parameters.CBV.LH_Electrode{index, 1};
                ElectrodeCBV_RH{indexCount, 1} = ProcData.Sleep.Parameters.CBV.RH_Electrode{index, 1};
                indexCount = indexCount + 1;
            end
            
            matLH_DeltaPower = cell2mat(LH_deltaPower);
            arrayLH_DeltaPower = reshape(matLH_DeltaPower', [1, size(matLH_DeltaPower, 2)*size(matLH_DeltaPower, 1)]);
            cellLH_DeltaPower = {arrayLH_DeltaPower};
            
            matRH_DeltaPower = cell2mat(RH_deltaPower);
            arrayRH_DeltaPower = reshape(matRH_DeltaPower', [1, size(matRH_DeltaPower, 2)*size(matRH_DeltaPower, 1)]);
            cellRH_DeltaPower = {arrayRH_DeltaPower};
            
            matLH_ThetaPower = cell2mat(LH_thetaPower);
            arrayLH_ThetaPower = reshape(matLH_ThetaPower', [1, size(matLH_ThetaPower, 2)*size(matLH_ThetaPower, 1)]);
            cellLH_ThetaPower = {arrayLH_ThetaPower};
            
            matRH_ThetaPower = cell2mat(RH_thetaPower);
            arrayRH_ThetaPower = reshape(matRH_ThetaPower', [1, size(matRH_ThetaPower, 2)*size(matRH_ThetaPower, 1)]);
            cellRH_ThetaPower = {arrayRH_ThetaPower};
            
            matLH_GammaPower = cell2mat(LH_gammaPower);
            arrayLH_GammaPower = reshape(matLH_GammaPower', [1, size(matLH_GammaPower, 2)*size(matLH_GammaPower, 1)]);
            cellLH_GammaPower = {arrayLH_GammaPower};
            
            matRH_GammaPower = cell2mat(RH_gammaPower);
            arrayRH_GammaPower = reshape(matRH_GammaPower', [1, size(matRH_GammaPower, 2)*size(matRH_GammaPower, 1)]);
            cellRH_GammaPower = {arrayRH_GammaPower};
            
            for x = 1:length(WhiskerAcceleration)
                targetPoints = size(WhiskerAcceleration{1, 1}, 2);
                if size(WhiskerAcceleration{x, 1}, 2) ~= targetPoints
                    maxLength = size(WhiskerAcceleration{x, 1}, 2);
                    difference = targetPoints - size(WhiskerAcceleration{x, 1}, 2);
                    for y = 1:difference
                        WhiskerAcceleration{x, 1}(maxLength + y) = 0;
                    end
                end
            end
            
            matWhiskerAcceleration = cell2mat(WhiskerAcceleration);
            arrayWhiskerAcceleration = reshape(matWhiskerAcceleration', [1, size(matWhiskerAcceleration, 2)*size(matWhiskerAcceleration, 1)]);
            cellWhiskerAcceleration = {arrayWhiskerAcceleration};
            
            for x = 1:length(HeartRate)
                targetPoints = size(HeartRate{1, 1}, 2);
                if size(HeartRate{x, 1}, 2) ~= targetPoints
                    maxLength = size(HeartRate{x, 1}, 2);
                    difference = targetPoints - size(HeartRate{x, 1}, 2);
                    for y = 1:difference
                        HeartRate{x, 1}(maxLength + y) = mean(HeartRate{x, 1});
                    end
                end
            end
            
            matHeartRate = cell2mat(HeartRate);
            arrayHeartRate = reshape(matHeartRate', [1, size(matHeartRate, 2)*size(matHeartRate, 1)]);
            cellHeartRate = {arrayHeartRate};
            
            matLH_CBV = cell2mat(CBV_LH);
            arrayLH_CBV = reshape(matLH_CBV', [1, size(matLH_CBV, 2)*size(matLH_CBV, 1)]);
            cellLH_CBV = {arrayLH_CBV};
            
            matRH_CBV = cell2mat(CBV_RH);
            arrayRH_CBV = reshape(matRH_CBV', [1, size(matRH_CBV, 2)*size(matRH_CBV, 1)]);
            cellRH_CBV = {arrayRH_CBV};
            
            matLH_ElectrodeCBV = cell2mat(ElectrodeCBV_LH);
            arrayLH_ElectrodeCBV = reshape(matLH_ElectrodeCBV', [1, size(matLH_ElectrodeCBV, 2)*size(matLH_ElectrodeCBV, 1)]);
            cellLH_ElectrodeCBV = {arrayLH_ElectrodeCBV};
            
            matRH_ElectrodeCBV = cell2mat(ElectrodeCBV_RH);
            arrayRH_ElectrodeCBV = reshape(matRH_ElectrodeCBV', [1, size(matRH_ElectrodeCBV, 2)*size(matRH_ElectrodeCBV, 1)]);
            cellRH_ElectrodeCBV = {arrayRH_ElectrodeCBV};
            
            cellBinTimes = {binTimes};
            
            %% Check results
            singleTrialFig = figure;
            ax1 = subplot(5,1,1);
            plot((1:length(ProcData.Data.Behavior.whiskers)) / ProcData.Notes.downsampledWhiskerSamplingRate,...
                -ProcData.Data.Behavior.whiskers, 'k');
            hold on;
            
            scatter((whiskingInds / ProcData.Notes.downsampledWhiskerSamplingRate), whisking_YVals, '.k');
            scatter(solenoidLeftPad, LH_solenoid_YVals, 'vk', 'MarkerFaceColor', 'b');
            scatter(solenoidRightPad, RH_solenoid_YVals, 'vk', 'MarkerFaceColor', 'g');
            scatter(solenoidAuditory, Aud_solenoid_YVals, 'vk', 'MarkerFaceColor', 'r');
            ylabel('Degrees')
            ylim([min(-whiskers), 1.30*max(-whiskers)])
            
            yyaxis right
            plot(tR, HR, 'm');
            ylabel('Heart Rate (Hz)');
            ylim([6 15]);
            
            title({[animal ' ' fileID ' Single Trial'], 'Behavioral State'});
            set(gca, 'Ticklength', [0 0])
            legend('Whisker Angle', 'Whisking Events', 'Left Pad Sol', 'Right Pad Sol', 'Auditory Sol', 'Heart Rate', 'Location', 'NorthEast')
            
            ax2 = subplot(5,1,2:3);
            timeVector = (1:9000) ./ 30;
            plot(timeVector, filtLH_CBV*100, 'k');
            hold on;
            plot(timeVector, filtRH_CBV*100, 'Color', colors('electric purple'));
            scatter(binTimes, (ones(1, length(binTimes))*sleeping_YVal), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors('carrot orange'))
            for sleepT = 1:length(sleepTimes)
                scatter(sleepTimes{sleepT, 1}, (ones(1, length(sleepTimes{sleepT, 1})))*sleeping_YVal, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c')
            end
            
            title('Normalized & Filtered CBV Signals');
            ylabel('Reflectance (%)')
            legend('LH CBV', 'RH CBV', 'REM Sleep', 'NREM Sleep', 'Location', 'NorthEast')
            set(gca, 'Ticklength', [0 0])
            axis tight
            
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
            
            isOK = input(['Are the displayed REM sleep times from bin edges ' num2str(binStart*5) ':' num2str(binEnd*5) ' seconds accurate? (y/n): '], 's');
            disp(' ')
            
            if strcmp(isOK, 'y')
                 [pathstr, ~, ~] = fileparts(cd);
                 dirpath = [pathstr '/Figures/REM Single Trial Figures/'];
                 
                 if ~exist(dirpath, 'dir')
                     mkdir(dirpath);
                 end
                 
                 savefig(singleTrialFig, [dirpath animal '_' sleepFile '_REMSingleTrialFig']);
                
                close(singleTrialFig)
            else
                continue
            end
            
            %% BLOCK PURPOSE: Save the data in the SleepEventData struct
            if isempty(SleepData.REM)  % If the structure is empty, we need a special case to format the struct properly
                SleepData.REM.Data.DeltaBand_Power.LH{1, 1} = cellLH_DeltaPower{1, 1};
                SleepData.REM.Data.DeltaBand_Power.RH{1, 1} = cellRH_DeltaPower{1, 1};
                SleepData.REM.Data.ThetaBand_Power.LH{1, 1} = cellLH_ThetaPower{1, 1};
                SleepData.REM.Data.ThetaBand_Power.RH{1, 1} = cellRH_ThetaPower{1, 1};
                SleepData.REM.Data.GammaBand_Power.LH{1, 1} = cellLH_GammaPower{1, 1};
                SleepData.REM.Data.GammaBand_Power.RH{1, 1} = cellRH_GammaPower{1, 1};
                SleepData.REM.Data.CBV.LH{1, 1} = cellLH_CBV{1, 1};
                SleepData.REM.Data.CBV.RH{1, 1} = cellRH_CBV{1, 1};
                SleepData.REM.Data.CBV.LH_Electrode{1, 1} = cellLH_ElectrodeCBV{1, 1};
                SleepData.REM.Data.CBV.RH_Electrode{1, 1} = cellRH_ElectrodeCBV{1, 1};
                SleepData.REM.Data.WhiskerAcceleration{1, 1} = cellWhiskerAcceleration{1, 1};
                SleepData.REM.Data.HeartRate{1, 1} = cellHeartRate{1, 1};
                SleepData.REM.FileIDs{1, 1} = sleepFile;
                SleepData.REM.BinTimes{1, 1} = cellBinTimes{1, 1};
            else    % If the struct is not empty, add each new iteration after previous data
                SleepData.REM.Data.DeltaBand_Power.LH{size(SleepData.REM.Data.DeltaBand_Power.LH, 1) + 1, 1} = cellLH_DeltaPower{1, 1};
                SleepData.REM.Data.DeltaBand_Power.RH{size(SleepData.REM.Data.DeltaBand_Power.RH, 1) + 1, 1} = cellRH_DeltaPower{1, 1};
                SleepData.REM.Data.ThetaBand_Power.LH{size(SleepData.REM.Data.ThetaBand_Power.LH, 1) + 1, 1} = cellLH_ThetaPower{1, 1};
                SleepData.REM.Data.ThetaBand_Power.RH{size(SleepData.REM.Data.ThetaBand_Power.RH, 1) + 1, 1} = cellRH_ThetaPower{1, 1};
                SleepData.REM.Data.GammaBand_Power.LH{size(SleepData.REM.Data.GammaBand_Power.LH, 1) + 1, 1} = cellLH_GammaPower{1, 1};
                SleepData.REM.Data.GammaBand_Power.RH{size(SleepData.REM.Data.GammaBand_Power.RH, 1) + 1, 1} = cellRH_GammaPower{1, 1};
                SleepData.REM.Data.CBV.LH{size(SleepData.REM.Data.CBV.LH, 1) + 1, 1} = cellLH_CBV{1, 1};
                SleepData.REM.Data.CBV.RH{size(SleepData.REM.Data.CBV.RH, 1) + 1, 1} = cellRH_CBV{1, 1};
                SleepData.REM.Data.CBV.LH_Electrode{size(SleepData.REM.Data.CBV.LH_Electrode, 1) + 1, 1} = cellLH_ElectrodeCBV{1, 1};
                SleepData.REM.Data.CBV.RH_Electrode{size(SleepData.REM.Data.CBV.RH_Electrode, 1) + 1, 1} = cellRH_ElectrodeCBV{1, 1};
                SleepData.REM.Data.WhiskerAcceleration{size(SleepData.REM.Data.WhiskerAcceleration, 1) + 1, 1} = cellWhiskerAcceleration{1, 1};
                SleepData.REM.Data.HeartRate{size(SleepData.REM.Data.HeartRate, 1) + 1, 1} = cellHeartRate{1, 1};
                SleepData.REM.FileIDs{size(SleepData.REM.FileIDs, 1) + 1, 1} = sleepFile;
                SleepData.REM.BinTimes{size(SleepData.REM.BinTimes, 1) + 1, 1} = cellBinTimes{1, 1};
            end
            break
        elseif strcmp(REMcase, 'n')
            disp('Loading next file...'); disp(' ')
            close(singleTrialFig)
            break
        end
    end
end

save([animal '_SleepData.mat'], 'SleepData');

end
