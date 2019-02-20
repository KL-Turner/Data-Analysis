function Correct_LabVIEW_Offset(combDataFiles)
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

for f = 1:size(combDataFiles, 1)
    %% Find offset between the two force sensor signals using the cross correlation
    disp(['Correcting offset in file number ' num2str(f) ' of ' num2str(size(combDataFiles, 1)) '...']); disp(' '); 
    combDataFile = combDataFiles(f, :);
    load(combDataFile);
    [animalID, ~, fileID, imageID] = GetFileInfo2(combDataFile);

    analogSamplingRate = CombData.Notes.LabVIEW.analogSamplingRate;
    whiskerSamplingRate = CombData.Notes.LabVIEW.whiskerCamSamplingRate;
    vesselSamplingRate = floor(CombData.Notes.MScan.frameRate);
    
    labviewForce = detrend(CombData.Data.Force_Sensor_L, 'constant');
    mscanForce = detrend(CombData.Data.Force_Sensor_M, 'constant');
    
    maxLag = 30*analogSamplingRate;
    [r, lags] = xcorr(labviewForce, mscanForce, maxLag);
    [~, index] = max(r);
    offset = lags(index);
    whiskerOffset = round(whiskerSamplingRate*(abs(offset)/analogSamplingRate));
    
    if offset > 0
        labviewShift = labviewForce(offset:end);
        whiskerShift = CombData.Data.Whisker_Angle(whiskerOffset:end);
    elseif offset <= 0
        fpad = zeros(1, abs(offset));
        wpad = zeros(1, abs(whiskerOffset));
        labviewShift = horzcat(fpad, labviewForce);
        whiskerShift = horzcat(wpad, CombData.Data.Whisker_Angle);
    end
    
    %% Original and Corrected figure
    corrOffset = figure;
    ax1 = subplot(3,1,1);
    plot((1:length(mscanForce))/analogSamplingRate, mscanForce, 'k')
    hold on;
    plot((1:length(labviewForce))/analogSamplingRate, labviewForce, 'r')
    title({[animalID ' ' fileID ' ' imageID ' force sensor data'], 'Offset correction between MScan and LabVIEW DAQ'})
    legend('Original MScan', 'Original LabVIEW')
    ylabel('A.U.')
    xlabel('Time (sec)')
    set(gca, 'Ticklength', [0 0])
    axis tight
    
    ax2 = subplot(3,1,2); %#ok<NASGU>
    plot(lags/analogSamplingRate, r, 'k')
    title('Cross Correlation between the two signals')
    ylabel('Correlation (A.U.)')
    xlabel('Lag (sec)')
    set(gca, 'Ticklength', [0 0])
    axis tight
    
    ax3 = subplot(3,1,3);
    plot((1:length(mscanForce))/analogSamplingRate, mscanForce, 'k')
    hold on;
    plot((1:length(labviewShift))/analogSamplingRate, labviewShift, 'b')
    title({'Shifted correction between MScan and LabVIEW DAQ', ['Offset value: ' num2str(offset) ' samples or ~' num2str(offset/analogSamplingRate) ' seconds']})
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
    
    saveas(corrOffset, [dirpath animalID '_' imageID '_' fileID '_CorrectedOffset'], 'tiff');
    close all
    
    %% Apply correction to the data, and trim excess time
    frontCut = 5;
    endCut = 5;
    
    mscanAnalogSampleDiff = analogSamplingRate*CombData.Notes.LabVIEW.trialDuration_Seconds - length(CombData.Data.Neural_Data);
    mscanAnalogCut = endCut*analogSamplingRate - mscanAnalogSampleDiff;
    
    labviewAnalogSampleDiff = analogSamplingRate*CombData.Notes.LabVIEW.trialDuration_Seconds - length(labviewShift);
    labviewAnalogCut = endCut*analogSamplingRate - labviewAnalogSampleDiff;
    
    labviewWhiskerSampleDiff = whiskerSamplingRate*CombData.Notes.LabVIEW.trialDuration_Seconds - length(whiskerShift);
    labviewWhiskerCut = endCut*whiskerSamplingRate - labviewWhiskerSampleDiff;
    
    CombData.Data.Vessel_Diameter = CombData.Data.Vessel_Diameter(frontCut*vesselSamplingRate:end - (endCut*vesselSamplingRate + 1));
    CombData.Data.Force_Sensor_M = CombData.Data.Force_Sensor_M(frontCut*analogSamplingRate:end - (mscanAnalogCut + 1));
    CombData.Data.Neural_Data = CombData.Data.Neural_Data(frontCut*analogSamplingRate:end - (mscanAnalogCut + 1));
    
    CombData.Data.Whisker_Angle = whiskerShift(frontCut*whiskerSamplingRate:end - (labviewWhiskerCut + 1));
    CombData.Data.Force_Sensor_L = labviewShift(frontCut*analogSamplingRate:end - (labviewAnalogCut + 1));
    
     disp('Updating CombData File...'); disp(' ')
     save([animalID '_' fileID '_' imageID '_CombData'], 'CombData')
     
end
