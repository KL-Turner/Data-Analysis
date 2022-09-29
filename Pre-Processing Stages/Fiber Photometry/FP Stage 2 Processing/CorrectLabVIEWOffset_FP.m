function [] = CorrectLabVIEWOffset_FP(trialDataFileIDs,rawDataFileIDs,trimTime)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Doric triggers the LabVIEW acquisition program to start recording, but there is a slight (~ 1 second) lag
%            associated with the Doric data. The force sensor is duplicated, and this function serves to correct that
%            offset by finding the peak in the cross correlation, and shifting the LabVIEW signals based on the number of
%            lags. The beginning/end of all signals are then snipped appropriately after shifting.
%________________________________________________________________________________________________________________________

for aa = 1:size(trialDataFileIDs,1)
    %% find offset between the two force sensor signals using the cross correlation
    trialDataFileID = trialDataFileIDs(aa,:);
    load(trialDataFileID);
    rawDataFileID = rawDataFileIDs(aa,:);
    load(rawDataFileID)
    disp(['Correcting offset in file number ' num2str(aa) ' of ' num2str(size(trialDataFileIDs, 1)) '...']); disp(' ');
    [animalID,~,fileID] = GetFileInfo_FP(rawDataFileID);
    analogSamplingRate = RawData.notes.analogSamplingRate;
    whiskCamSamplingRate = RawData.notes.whiskCamSamplingRate;
    doricSamplingRate = trialData.notes.samplingRate;
    whiskerCamSamplingRate = RawData.notes.whiskCamSamplingRate;
    trialDuration = RawData.notes.trialDuration_long;
    doricTrialDuration = trialDuration - 1;
    dsFs = 100; % Hz
    labviewSyncSignal = detrend(resample(RawData.data.sync_long,dsFs,analogSamplingRate),'constant');
    labviewSyncSignal = labviewSyncSignal(dsFs*trialDuration/4:dsFs*trialDuration/2);
    fiberSyncSignal = detrend(resample(trialData.syncData,dsFs,doricSamplingRate),'constant');
    fiberSyncSignal = fiberSyncSignal(round(dsFs*doricTrialDuration/4):round(dsFs*doricTrialDuration/2));
    sampleDiff = length(labviewSyncSignal) - length(fiberSyncSignal);
    fiberSyncSignal = vertcat(fiberSyncSignal,zeros(sampleDiff,1));
    % dsFs xcorr
    maxLag = 10*dsFs;
    [r,lags] = xcorr(labviewSyncSignal,fiberSyncSignal,maxLag);
    [~,index] = max(r(1:(length(r) - 1)/2) - dsFs/10);
    offset = lags(index);
    offsetTime = abs(offset/dsFs);
    disp(['LabVIEW trailed Doric by ' num2str(offsetTime) ' seconds.']); disp(' ')
    % offset
    dsOffset = round(dsFs*(abs(offset)/dsFs));
    dsFs_pad = zeros(1,abs(dsOffset));
    labviewSyncShift = horzcat(dsFs_pad,labviewSyncSignal);
    corrOffset = figure;
    sgtitle([rawDataFileID ' cross ' trialDataFileID])
    ax1 = subplot(3,1,1);
    plot((1:length(fiberSyncSignal))/dsFs,fiberSyncSignal,'b')
    hold on;
    plot((1:length(labviewSyncSignal))/dsFs,labviewSyncSignal,'r')
    title({[animalID ' ' fileID ' sync channel data'],'offset correction between Doric and LabVIEW DAQ'})
    legend('Original Doric','Original LabVIEW')
    ylabel('A.U.')
    xlabel('Time (sec)')
    set(gca,'Ticklength',[0,0])
    axis tight
    ax2 = subplot(3,1,2); %#ok<NASGU>
    plot(lags/dsFs,r,'k')
    title('Cross Correlation between the two signals')
    ylabel('Correlation (A.U.)')
    xlabel('Lag (sec)')
    set(gca,'Ticklength',[0,0])
    axis tight
    ax3 = subplot(3,1,3);
    plot((1:length(fiberSyncSignal))/dsFs,fiberSyncSignal,'k')
    hold on;
    plot((1:length(labviewSyncShift))/dsFs,labviewSyncShift,'b')
    title({'Shifted correction between Doric and LabVIEW DAQ',['Offset value: ' num2str(offset) ' samples or ~' num2str(offset/dsFs) ' seconds']})
    legend('Original Doric','Shifted LabVIEW')
    ylabel('A.U.')
    xlabel('Time (sec)')
    set(gca,'Ticklength',[0,0])
    axis tight
    linkaxes([ax1,ax3],'x')
    %% apply doric correction to the data and trim excess time
    doricSampleDiff = doricSamplingRate*(trialDuration) - length(trialData.syncData);
    doricCut = floor(trimTime*doricSamplingRate - doricSampleDiff);
    doricFields = {'RH_405','RH_465','RH_560','LH_405','LH_465','LH_560'};
    for cc = 1:length(doricFields)
        subfields = fieldnames(trialData.Opticaldata.(doricFields{1,cc}));
        for dd = 1:length(subfields)
            RawData.data.(doricFields{1,cc}).(char(subfields(dd,1))) = trialData.Opticaldata.(doricFields{1,cc}).(char(subfields(dd,1)))(floor(trimTime*doricSamplingRate):end - (doricCut + 1));
        end
    end
    RawData.data.fiberSync = trialData.syncData(floor(trimTime*doricSamplingRate):end - (doricCut + 1));
    %% apply labview correction to the data and trim excess time
    labviewFields = {'cortical_LH','cortical_RH','hippocampus','EMG','forceSensor','stimulations','sync'};
    for cc = 1:length(labviewFields)
        labviewAnalogPad = zeros(1,round(offsetTime*analogSamplingRate));
        labviewAnalogShift = horzcat(labviewAnalogPad,RawData.data.([labviewFields{1,cc} '_long']));
        labviewAnalogSampleDiff = analogSamplingRate*trialDuration - length(labviewAnalogShift);
        labviewAnalogCut = trimTime*analogSamplingRate - labviewAnalogSampleDiff;
        RawData.data.(labviewFields{1,cc}) = labviewAnalogShift(floor(trimTime*analogSamplingRate):end - (labviewAnalogCut + 1));
    end
    labviewWhiskerPad = zeros(1,round(offsetTime*whiskCamSamplingRate));
    labviewWhiskerShift = horzcat(labviewWhiskerPad,RawData.data.whiskerAngle_long);
    labviewWhiskerSampleDiff = whiskCamSamplingRate*trialDuration - length(labviewWhiskerShift);
    labviewWhiskerCut = trimTime*whiskerCamSamplingRate - labviewWhiskerSampleDiff;
    RawData.data.whiskerAngle = labviewWhiskerShift(floor(trimTime*whiskCamSamplingRate):end - (labviewWhiskerCut + 1));
    RawData.notes.offsetCorrect = true;
    RawData.notes.trimTime = trimTime;
    RawData.notes.trialDuration_sec = trialDuration - 2*trimTime;
    %% check shift
    labviewSyncSignal_2 = resample(detrend(RawData.data.fiberSync,'constant'),dsFs,doricSamplingRate);
    fiberSyncSignal_2 = resample(detrend(RawData.data.sync,'constant'),dsFs,analogSamplingRate);
    checkShift = figure;
    p1  =  plot((1:length(fiberSyncSignal_2))/dsFs,fiberSyncSignal_2 ,'b');
    hold on;
    p2 = plot((1:length(labviewSyncSignal_2))/dsFs,labviewSyncSignal_2 ,'r');
    legend([p1,p2],'Doric','LabVIEW')
    title([rawDataFileID ' cross ' trialDataFileID])
    axis tight
    %% save files
    RawData.notes.doric = trialData.notes;
    save(rawDataFileID,'RawData','-v7.3')
    %% Save the file to directory.
    [pathstr,~,~] = fileparts(cd);
    dirpath = [pathstr '/Figures/XCorr Shift/'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(corrOffset,[dirpath animalID '_' fileID '_XCorrShift']);
    close(corrOffset)
    dirpath = [pathstr '/Figures/Shift Check/'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(checkShift,[dirpath animalID '_' fileID '_ShiftCheck']);
    close(checkShift)
end

end
