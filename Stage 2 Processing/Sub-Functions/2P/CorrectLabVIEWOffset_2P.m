function CorrectLabVIEWOffset_2P(labviewDataFiles, mscanDataFiles, trimTime)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: MScan triggers the LabVIEW acquisition program to start recording, but there is a slight (~ 1 second) lag
%            associated with the MScan data. The force sensor is duplicated, and this function serves to correct that
%            offset by finding the peak in the cross correlation, and shifting the LabVIEW signals based on the number of
%            lags. The beginning/end of all signals are then snipped appropriately after shifting.
%________________________________________________________________________________________________________________________
%
%   Inputs: List of MScan and LabVIEW data files, along with the number of seconds to be trimmed off after shifting.
%
%   Outputs: Saves updated fields to the MScan and LabVIEW data structures.
%
%   Last Revised: February 29th, 2019
%________________________________________________________________________________________________________________________

for a = 1:size(mscanDataFiles,1)
    %% Find offset between the two force sensor signals using the cross correlation
    mscanDataFile = mscanDataFiles(a, :);
    load(mscanDataFile);
    labviewDataFile = MScanData.notes.labviewFileID;
    load(labviewDataFile)
    if MScanData.notes.checklist.offsetCorrect == false
        disp(['Correcting offset in file number ' num2str(a) ' of ' num2str(size(mscanDataFiles, 1)) '...']); disp(' ');
        [animalID, hem, fileDate, fileID] = GetFileInfo_2P(labviewDataFile);
        imageID = MScanData.notes.imageID;
        vesselID = MScanData.notes.vesselID;

        analogSamplingRate = LabVIEWData.notes.analogSamplingRate_Hz;
        whiskerCamSamplingRate = LabVIEWData.notes.whiskerCamSamplingRate_Hz;
        dsFs = MScanData.notes.downSampledFs;
        vesselSamplingRate = floor(MScanData.notes.frameRate);
        trialDuration = LabVIEWData.notes.trialDuration_Seconds;
        
        labviewForce = detrend(LabVIEWData.data.dsForceSensorL, 'constant');
        mscanForce = detrend(MScanData.data.dsForceSensorM, 'constant');
        analog_labviewForce = detrend(LabVIEWData.data.forceSensor, 'constant');
        analog_mscanForce = detrend(MScanData.data.forceSensor, 'constant');
        
        maxLag = 30*dsFs;
        analog_MaxLag = 30*analogSamplingRate;
        [analog_r, analog_lags] = xcorr(analog_labviewForce, analog_mscanForce, analog_MaxLag);
        [r, lags] = xcorr(labviewForce, mscanForce, maxLag);
        [~, analog_index] = max(analog_r);
        [~, index] = max(r);
        offset = lags(index);
        analog_offset = analog_lags(analog_index);
        analog_forceOffset = round(abs(analog_offset)/analogSamplingRate);
        analog_whiskerOffset = round(abs(analog_offset)/whiskerCamSamplingRate);
        dsOffset = round(dsFs*(abs(offset)/dsFs));
        disp(['LabVIEW trailed MScan by ' num2str(-offset/dsFs) ' seconds.']); disp(' ')
        
        if offset > 0
            analog_forceShift = analog_labviewForce(analog_forceOffset:end);
            analog_whiskerShift = LabVIEWData.data.whiskerAngle(analog_whiskerOffset:end);
            dsForceShift = labviewForce(offset:end);
            dsWhiskShift = LabVIEWData.data.dsWhiskerAngle(offset:end);
            binForceShift = LabVIEWData.data.binForceSensorL(offset:end);
            binWhiskShift = LabVIEWData.data.binWhiskerAngle(offset:end);
        elseif offset <= 0
            analog_fpad = zeros(1, abs(analog_forceOffset));
            analog_wpad = zeros(1, abs(analog_whiskerOffset));
            pad = zeros(1, abs(dsOffset));
            analog_forceShift = horzcat(analog_fpad, analog_labviewForce);
            analog_whiskerShift = horzcat(analog_wpad, LabVIEWData.data.whiskerAngle);
            dsForceShift = horzcat(pad, labviewForce);
            dsWhiskShift = horzcat(pad, LabVIEWData.data.dsWhiskerAngle);
            binForceShift = horzcat(pad, LabVIEWData.data.binForceSensorL);
            binWhiskShift = horzcat(pad, LabVIEWData.data.binWhiskerAngle);
        end
        
        corrOffset = figure;
        ax1 = subplot(3,1,1);
        plot((1:length(mscanForce))/dsFs, mscanForce, 'k')
        hold on;
        plot((1:length(labviewForce))/dsFs, labviewForce, 'r')
        title({[animalID ' ' fileID ' ' imageID ' force sensor data'], 'Offset correction between MScan and LabVIEW DAQ'})
        legend('Original MScan', 'Original LabVIEW')
        ylabel('A.U.')
        xlabel('Time (sec)')
        set(gca, 'Ticklength', [0 0])
        axis tight
        
        ax2 = subplot(3,1,2); %#ok<NASGU>
        plot(analog_lags/dsFs, analog_r, 'k')
        title('Cross Correlation between the two signals')
        ylabel('Correlation (A.U.)')
        xlabel('Lag (sec)')
        set(gca, 'Ticklength', [0 0])
        axis tight
        
        ax3 = subplot(3,1,3);
        plot((1:length(mscanForce))/dsFs, mscanForce, 'k')
        hold on;
        plot((1:length(dsForceShift))/dsFs, dsForceShift, 'b')
        title({'Shifted correction between MScan and LabVIEW DAQ', ['Offset value: ' num2str(offset) ' samples or ~' num2str(offset/dsFs) ' seconds']})
        legend('Original MScan', 'Shifted LabVIEW')
        ylabel('A.U.')
        xlabel('Time (sec)')
        set(gca, 'Ticklength', [0 0])
        axis tight
        linkaxes([ax1 ax3], 'x')
        
        %% Apply correction to the data, and trim excess time
        frontCut = trimTime;
        endCut = trimTime;
        
        mscanAnalogSampleDiff = analogSamplingRate*trialDuration - length(MScanData.data.forceSensor);
        mscanAnalogCut = endCut*analogSamplingRate - mscanAnalogSampleDiff;
        
        mscan_dsAnalogSampleDiff = dsFs*trialDuration - length(MScanData.data.dsForceSensorM);
        mscan_dsAnalogCut = endCut*dsFs - mscan_dsAnalogSampleDiff;
        
        mscan_binForceSampleDiff = dsFs*trialDuration - length(MScanData.data.binForceSensorM);
        mscan_binForceCut = endCut*dsFs - mscan_binForceSampleDiff;
        
        labview_AnalogSampleDiff = analogSamplingRate*trialDuration - length(analog_forceShift);
        labview_AnalogCut = endCut*analogSamplingRate - labview_AnalogSampleDiff;
        
        labview_WhiskerSampleDiff = whiskerCamSamplingRate*trialDuration - length(analog_whiskerShift);
        labview_WhiskerCut = endCut*whiskerCamSamplingRate - labview_WhiskerSampleDiff;
        
        labview_dsWhiskSamplingDiff = dsFs*trialDuration - length(dsWhiskShift);
        labview_dsWhiskCut = endCut*dsFs - labview_dsWhiskSamplingDiff;
        
        labview_dsForceSamplingDiff = dsFs*trialDuration - length(dsForceShift);
        labview_dsForceCut = endCut*dsFs - labview_dsForceSamplingDiff;
        
        labview_binForceSampleDiff = dsFs*trialDuration - length(binForceShift);
        labview_binForceCut = endCut*dsFs - labview_binForceSampleDiff;
        
        labview_binWhiskSamplingDiff = dsFs*trialDuration - length(binWhiskShift);
        labview_binWhiskCut = endCut*dsFs - labview_binWhiskSamplingDiff;
        
        MScanData.data.forceSensor_trim = MScanData.data.forceSensor(frontCut*analogSamplingRate:end - (mscanAnalogCut + 1))';
        MScanData.data.rawNeuralData_trim = MScanData.data.rawNeuralData(frontCut*analogSamplingRate:end - (mscanAnalogCut + 1))';
        MScanData.data.EMG_trim = MScanData.data.EMG(frontCut*analogSamplingRate:end - (mscanAnalogCut + 1))';
        MScanData.data.vesselDiameter_trim = MScanData.data.vesselDiameter(frontCut*vesselSamplingRate:end - (endCut*vesselSamplingRate + 1));
        MScanData.data.rawVesselDiameter_trim = MScanData.data.rawVesselDiameter(frontCut*vesselSamplingRate:end - (endCut*vesselSamplingRate + 1));
        MScanData.data.tempVesselDiameter_trim = MScanData.data.tempVesselDiameter(frontCut*vesselSamplingRate:end - (endCut*vesselSamplingRate + 1));
        MScanData.data.muaPower_trim = MScanData.data.muaPower(frontCut*dsFs:end - (mscan_dsAnalogCut + 1))';
        MScanData.data.gammaPower_trim = MScanData.data.gammaPower(frontCut*dsFs:end - (mscan_dsAnalogCut + 1))';
        MScanData.data.betaPower_trim = MScanData.data.betaPower(frontCut*dsFs:end - (mscan_dsAnalogCut + 1))';
        MScanData.data.alphaPower_trim = MScanData.data.alphaPower(frontCut*dsFs:end - (mscan_dsAnalogCut + 1))';
        MScanData.data.thetaPower_trim = MScanData.data.thetaPower(frontCut*dsFs:end - (mscan_dsAnalogCut + 1))';
        MScanData.data.deltaPower_trim = MScanData.data.deltaPower(frontCut*dsFs:end - (mscan_dsAnalogCut + 1))';
        MScanData.data.dsForceSensorM_trim = MScanData.data.dsForceSensorM(frontCut*dsFs:end - (mscan_dsAnalogCut + 1))';
        MScanData.data.binForceSensorM_trim = MScanData.data.binForceSensorM(frontCut*dsFs:end - (mscan_binForceCut + 1))';
        MScanData.data.dsEMG_trim = MScanData.data.dsEMG(frontCut*dsFs:end - (mscan_dsAnalogCut + 1))';
        MScanData.notes.checklist.offsetCorrect = true;
        
        LabVIEWData.data.forceSensor_trim = analog_forceShift(frontCut*analogSamplingRate:end - (labview_AnalogCut + 1));
        LabVIEWData.data.whiskerAngle_trim = analog_whiskerShift(frontCut*whiskerCamSamplingRate:end - (labview_WhiskerCut + 1));
        LabVIEWData.data.dsWhiskerAngle_trim = dsWhiskShift(frontCut*dsFs:end - (labview_dsWhiskCut + 1));
        LabVIEWData.data.binWhiskerAngle_trim = binWhiskShift(frontCut*dsFs:end - (labview_binWhiskCut + 1));
        LabVIEWData.data.dsForceSensorL_trim = dsForceShift(frontCut*dsFs:end - (labview_dsForceCut + 1));
        LabVIEWData.data.binForceSensorL_trim = binForceShift(frontCut*dsFs:end - (labview_binForceCut + 1));
        LabVIEWData.notes.checklist.offsetCorrect = true;
        LabVIEWData.notes.trimTime = trimTime;
        LabVIEWData.notes.trialDuration_Seconds_trim = LabVIEWData.notes.trialDuration_Seconds - 2*trimTime;
        
        disp('Updating MScanData and LabVIEW Files...'); disp(' ')
        save([animalID '_' fileDate '_' vesselID '_' imageID '_MScanData'], 'MScanData')
        save([animalID '_' hem '_' fileID '_LabVIEWData'], 'LabVIEWData')
    else
        disp(['Offset in ' mscanDataFile ' and ' labviewDataFile ' has already been corrected. Continuing...']); disp(' ');
    end
end

end

