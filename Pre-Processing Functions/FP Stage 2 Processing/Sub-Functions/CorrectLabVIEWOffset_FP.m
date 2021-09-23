function [] = CorrectLabVIEWOffset_FP(fiberDataFileIDs,rawDataFileIDs,trimTime)
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

for aa = 1:size(fiberDataFileIDs,1)
    %% find offset between the two force sensor signals using the cross correlation
    fiberDataFileID = fiberDataFileIDs(aa,:);
    load(fiberDataFileID);
    rawDataFileID = rawDataFileIDs(aa,:);
    load(rawDataFileID)
    disp(['Correcting offset in file number ' num2str(aa) ' of ' num2str(size(fiberDataFileIDs, 1)) '...']); disp(' ');
    [animalID,~,fileID] = GetFileInfo_FP(rawDataFileID);
    analogSamplingRate = RawData.notes.analogSamplingRate;
    whiskCamSamplingRate = RawData.notes.whiskCamSamplingRate;
    doricSamplingRate = FiberData.notes.samplingRate;
    whiskerCamSamplingRate = RawData.notes.whiskCamSamplingRate;
    trialDuration = RawData.notes.trialDuration_sec;
    dsFs = 100; % Hz
    labviewSyncSignal = resample(detrend(RawData.data.sync,'constant'),dsFs,analogSamplingRate);
    fiberSyncSignal = resample(detrend(FiberData.data.sync,'constant'),dsFs,doricSamplingRate);
    sampleDiff = length(labviewSyncSignal) - length(fiberSyncSignal);
    fiberSyncSignal = vertcat(fiberSyncSignal,zeros(sampleDiff,1));
    % dsFs xcorr
    maxLag = 10*dsFs;
    [r,lags] = xcorr(labviewSyncSignal,fiberSyncSignal,maxLag);
    [~,index] = max(r(1:(length(r) - 1)/2));
    offset = lags(index);
    offsetTime = abs(offset/dsFs);
    disp(['LabVIEW trailed Doric by ' num2str(offsetTime) ' seconds.']); disp(' ')
    % offset
    dsOffset = round(dsFs*(abs(offset)/dsFs));
    dsFs_pad = zeros(1,abs(dsOffset));
    labviewSyncShift = horzcat(dsFs_pad,labviewSyncSignal);
    corrOffset = figure;
    ax1 = subplot(3,1,1);
    plot((1:length(fiberSyncSignal))/dsFs,fiberSyncSignal,'k')
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
    doricSampleDiff = doricSamplingRate*trialDuration - length(FiberData.data.sync);
    doricCut = floor(trimTime*doricSamplingRate - doricSampleDiff);
    doricFields = {'RH_405','RH_465','RH_560','LH_405','LH_465','LH_560'};
    for cc = 1:length(doricFields)
        subfields = fieldnames(FiberData.data.(doricFields{1,cc}));
        for dd = 1:length(subfields)
            RawData.data.(doricFields{1,cc}).([char(subfields(dd,1)) '_trim']) = FiberData.data.(doricFields{1,cc}).(char(subfields(dd,1)))(floor(trimTime*doricSamplingRate):end - (doricCut + 1));
        end
    end
    RawData.data.sync_trim = FiberData.data.sync(floor(trimTime*doricSamplingRate):end - (doricCut + 1));
    %% apply labview correction to the data and trim excess time
    labviewAnalogPad = zeros(1,round(offsetTime*analogSamplingRate));
    labviewAnalogShift = horzcat(labviewAnalogPad,RawData.data.sync);
    labviewAnalogSampleDiff = analogSamplingRate*trialDuration - length(labviewAnalogShift);
    labviewAnalogCut = trimTime*analogSamplingRate - labviewAnalogSampleDiff;
    labviewFields = {'cortical_LH','cortical_RH','hippocampus','EMG','forceSensor','stimulations','LH_560'};
    for cc = 1:length(labviewFields)
        RawData.data.([labviewFields{1,cc} '_trim']) = labviewAnalogShift(floor(trimTime*analogSamplingRate):end - (labviewAnalogCut + 1));
    end
    labviewWhiskerPad = zeros(1,round(offsetTime*whiskCamSamplingRate));
    labviewWhiskerShift = horzcat(labviewWhiskerPad,RawData.data.whiskerAngle);
    labviewWhiskerSampleDiff = whiskCamSamplingRate*trialDuration - length(labviewWhiskerShift);
    labviewWhiskerCut = trimTime*whiskerCamSamplingRate - labviewWhiskerSampleDiff;
    RawData.data.whiskerAngle_trim = labviewWhiskerShift(floor(trimTime*whiskCamSamplingRate):end - (labviewWhiskerCut + 1));
    RawData.notes.offsetCorrect = true;
    RawData.notes.trimTime = trimTime;
    RawData.notes.trialDuration_sec_trim = trialDuration - 2*trimTime;
    %% check shift
    labviewSyncSignal_2 = resample(detrend(RawData.data.LH_560_trim,'constant'),dsFs,analogSamplingRate);
    fiberSyncSignal_2 = resample(detrend(RawData.data.sync_trim,'constant'),dsFs,doricSamplingRate);
    checkShift = figure;
    p1  =  plot((1:length(fiberSyncSignal_2))/dsFs,fiberSyncSignal_2 ,'k');
    hold on;
    p2 = plot((1:length(labviewSyncSignal_2))/dsFs,labviewSyncSignal_2 ,'b');
    legend([p1,p2],'Doric','LabVIEW')
    axis tight
    %% save files
    RawData.notes.doric = FiberData.notes;
    save(rawDataFileID,'RawData','-v7.3')
    %% Save the file to directory.
    [pathstr,~,~] = fileparts(cd);
    dirpath = [pathstr '/Figures/XCorr Shift/'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(corrOffset,[dirpath animalID '_' fileID '_XCorrShift']);
    close(corrOffset)
    savefig(checkShift,[dirpath animalID '_' fileID '_ShiftCheck']);
    close(checkShift)
end

end
