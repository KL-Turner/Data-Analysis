function CorrectLabVIEWOffset(labviewDataFiles, mscanDataFiles)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner 
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
%
%   Last Revised: February 20th, 2019    
%________________________________________________________________________________________________________________________

for f = 3:size(mscanDataFiles, 1)
    %% Find offset between the two force sensor signals using the cross correlation
    disp(['Correcting offset in file number ' num2str(f) ' of ' num2str(size(mscanDataFiles, 1)) '...']); disp(' '); 
    mscanDataFile = mscanDataFiles(f, :);
    load(mscanDataFile);
    labviewDataFile = labviewDataFiles(f, :);
    load(labviewDataFile)
    
    [animalID, hem, fileDate, fileID] = GetFileInfo(labviewDataFile);
    imageID = MScanData.Notes.imageID;

    analogSamplingRate = LabVIEWData.Notes.analogSamplingRate;
    whiskerCamSamplingRate = LabVIEWData.Notes.whiskerCamSamplingRate;
    dsSamplingRate = LabVIEWData.Notes.downsampledForceSensorSamplingRate;
    vesselSamplingRate = floor(MScanData.Notes.frameRate);
    trialDuration = LabVIEWData.Notes.trialDuration_Seconds;
    
    labviewForce = detrend(LabVIEWData.Data.dsForce_Sensor_L, 'constant');
    mscanForce = detrend(MScanData.Data.dsForce_Sensor_M, 'constant');
    analog_labviewForce = detrend(LabVIEWData.Data.Force_Sensor, 'constant');
    analog_mscanForce = detrend(MScanData.Data.MScan_Force_Sensor, 'constant');
    
    maxLag = 30*dsSamplingRate;
    analog_MaxLag = 30*analogSamplingRate;
    [analog_r, analog_lags] = xcorr(analog_labviewForce, analog_mscanForce, analog_MaxLag);
    [r, lags] = xcorr(labviewForce, mscanForce, maxLag);
    [~, analog_index] = max(analog_r);
    [~, index] = max(r);
    offset = lags(index);
    analog_offset = analog_lags(analog_index);
    analog_forceOffset = round(abs(analog_offset)/analogSamplingRate);
    analog_whiskerOffset = round(abs(analog_offset)/whiskerCamSamplingRate);
    dsOffset = round(dsSamplingRate*(abs(offset)/dsSamplingRate));
    
    if offset > 0
        analog_forceShift = analog_labviewForce(analog_forceOffset:end);
        analog_whiskerShift = LabVIEWData.Data.WhiskerAngle(analog_whiskerOffset:end);
        dsForceShift = labviewForce(offset:end);
        dsWhiskShift = LabVIEWData.Data.dsWhisker_Angle(offset:end);
        binForceShift = LabVIEWData.Data.binForce_Sensor_L(offset:end);
        binWhiskShift = LabVIEWData.Data.binWhisker_Angle(offset:end);
    elseif offset <= 0
        analog_fpad = zeros(1, abs(analog_forceOffset));
        analog_wpad = zeros(1, abs(analog_whiskerOffset));
        pad = zeros(1, abs(dsOffset));
        analog_forceShift = horzcat(analog_fpad, analog_labviewForce);
        analog_whiskerShift = horzcat(analog_wpad, LabVIEWData.Data.WhiskerAngle);
        dsForceShift = horzcat(pad, labviewForce);
        dsWhiskShift = horzcat(pad, LabVIEWData.Data.dsWhisker_Angle);
        binForceShift = horzcat(pad, LabVIEWData.Data.binForce_Sensor_L);
        binWhiskShift = horzcat(pad, LabVIEWData.Data.binWhisker_Angle);
    end
    
    %% Original and Corrected figure
    corrOffset = figure;
    ax1 = subplot(3,1,1);
    plot((1:length(mscanForce))/dsSamplingRate, mscanForce, 'k')
    hold on;
    plot((1:length(labviewForce))/dsSamplingRate, labviewForce, 'r')
    title({[animalID ' ' fileID ' ' imageID ' force sensor data'], 'Offset correction between MScan and LabVIEW DAQ'})
    legend('Original MScan', 'Original LabVIEW')
    ylabel('A.U.')
    xlabel('Time (sec)')
    set(gca, 'Ticklength', [0 0])
    axis tight
    
    ax2 = subplot(3,1,2); %#ok<NASGU>
    plot(analog_lags/dsSamplingRate, analog_r, 'k')
    title('Cross Correlation between the two signals')
    ylabel('Correlation (A.U.)')
    xlabel('Lag (sec)')
    set(gca, 'Ticklength', [0 0])
    axis tight
    
    ax3 = subplot(3,1,3);
    plot((1:length(mscanForce))/dsSamplingRate, mscanForce, 'k')
    hold on;
    plot((1:length(dsForceShift))/dsSamplingRate, dsForceShift, 'b')
    title({'Shifted correction between MScan and LabVIEW DAQ', ['Offset value: ' num2str(offset) ' samples or ~' num2str(offset/dsSamplingRate) ' seconds']})
    legend('Original MScan', 'Shifted LabVIEW')
    ylabel('A.U.')
    xlabel('Time (sec)')
    set(gca, 'Ticklength', [0 0])
    axis tight
    linkaxes([ax1 ax3], 'x')
    
    %% Save the file to directory.
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Offset Correction/'];
    
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    
    saveas(corrOffset, [dirpath animalID '_' fileID '_' imageID '_CorrectedOffset'], 'tiff');
    close all
    
    %% Apply correction to the data, and trim excess time
    frontCut = 5;
    endCut = 5;
    
    mscanAnalogSampleDiff = analogSamplingRate*trialDuration - length(MScanData.Data.MScan_Force_Sensor);
    mscanAnalogCut = endCut*analogSamplingRate - mscanAnalogSampleDiff;
    
    mscan_dsAnalogSampleDiff = dsSamplingRate*trialDuration - length(MScanData.Data.dsForce_Sensor_M);
    mscan_dsAnalogCut = endCut*dsSamplingRate - mscan_dsAnalogSampleDiff;
    
    mscan_binForceSampleDiff = dsSamplingRate*trialDuration - length(MScanData.Data.binForce_Sensor_M);
    mscan_binForceCut = endCut*dsSamplingRate - mscan_binForceSampleDiff;
    
    labview_AnalogSampleDiff = analogSamplingRate*trialDuration - length(analog_forceShift);
    labview_AnalogCut = endCut*analogSamplingRate - labview_AnalogSampleDiff;
    
    labview_WhiskerSampleDiff = whiskerCamSamplingRate*trialDuration - length(analog_whiskerShift);
    labview_WhiskerCut = endCut*whiskerCamSamplingRate - labview_WhiskerSampleDiff;
    
    labview_dsWhiskSamplingDiff = dsSamplingRate*trialDuration - length(dsWhiskShift);
    labview_dsWhiskCut = endCut*dsSamplingRate - labview_dsWhiskSamplingDiff;
    
    labview_dsForceSamplingDiff = dsSamplingRate*trialDuration - length(dsForceShift);
    labview_dsForceCut = endCut*dsSamplingRate - labview_dsForceSamplingDiff;
    
    labview_binForceSampleDiff = dsSamplingRate*trialDuration - length(binForceShift);
    labview_binForceCut = endCut*dsSamplingRate - labview_binForceSampleDiff;
    
    labview_binWhiskSamplingDiff = dsSamplingRate*trialDuration - length(binWhiskShift);
    labview_binWhiskCut = endCut*dsSamplingRate - labview_binWhiskSamplingDiff;
    
    MScanData.Data.MScan_Force_Sensor = MScanData.Data.MScan_Force_Sensor(frontCut*analogSamplingRate:end - (mscanAnalogCut + 1))';
    MScanData.Data.MScan_Neural_Data = MScanData.Data.MScan_Neural_Data(frontCut*analogSamplingRate:end - (mscanAnalogCut + 1))';
    MScanData.Data.Vessel_Diameter = MScanData.Data.Vessel_Diameter(frontCut*vesselSamplingRate:end - (endCut*vesselSamplingRate + 1));
    MScanData.Data.Vessel_RawDiameter = MScanData.Data.Vessel_RawDiameter(frontCut*vesselSamplingRate:end - (endCut*vesselSamplingRate + 1));
    MScanData.Data.Vessel_TempDiameter = MScanData.Data.Vessel_TempDiameter(frontCut*vesselSamplingRate:end - (endCut*vesselSamplingRate + 1));
    MScanData.Data.MUA_Power = MScanData.Data.MUA_Power(frontCut*dsSamplingRate:end - (mscan_dsAnalogCut + 1))';
    MScanData.Data.GammaBand_Power = MScanData.Data.GammaBand_Power(frontCut*dsSamplingRate:end - (mscan_dsAnalogCut + 1))';
    MScanData.Data.BetaBand_Power = MScanData.Data.BetaBand_Power(frontCut*dsSamplingRate:end - (mscan_dsAnalogCut + 1))';
    MScanData.Data.AlphaBand_Power = MScanData.Data.AlphaBand_Power(frontCut*dsSamplingRate:end - (mscan_dsAnalogCut + 1))';
    MScanData.Data.ThetaBand_Power = MScanData.Data.ThetaBand_Power(frontCut*dsSamplingRate:end - (mscan_dsAnalogCut + 1))';
    MScanData.Data.DeltaBand_Power = MScanData.Data.DeltaBand_Power(frontCut*dsSamplingRate:end - (mscan_dsAnalogCut + 1))';
    MScanData.Data.dsForce_Sensor_M = MScanData.Data.dsForce_Sensor_M(frontCut*dsSamplingRate:end - (mscan_dsAnalogCut + 1))';
    MScanData.Data.binForce_Sensor_M = MScanData.Data.binForce_Sensor_M(frontCut*dsSamplingRate:end - (mscan_binForceCut + 1))';

    LabVIEWData.Data.Force_Sensor = analog_forceShift(frontCut*analogSamplingRate:end - (labview_AnalogCut + 1));
    LabVIEWData.Data.WhiskerAngle = analog_whiskerShift(frontCut*whiskerCamSamplingRate:end - (labview_WhiskerCut + 1));
    LabVIEWData.Data.dsWhisker_Angle = dsWhiskShift(frontCut*dsSamplingRate:end - (labview_dsWhiskCut + 1));
    LabVIEWData.Data.binWhisker_Angle = binWhiskShift(frontCut*dsSamplingRate:end - (labview_binWhiskCut + 1));
    LabVIEWData.Data.dsForce_Sensor_L = dsForceShift(frontCut*dsSamplingRate:end - (labview_dsForceCut + 1));
    LabVIEWData.Data.binForce_Sensor_L = binForceShift(frontCut*dsSamplingRate:end - (labview_binForceCut + 1));

    disp('Updating MScanData and LabVIEW Files...'); disp(' ')
    save([animalID '_' fileDate '_' imageID '_MScanData'], 'MScanData') 
    save([animalID '_' hem '_' fileID '_LabVIEWData'], 'LabVIEWData')
 
end
