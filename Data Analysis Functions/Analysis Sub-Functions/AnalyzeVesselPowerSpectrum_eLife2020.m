function [AnalysisResults] = AnalyzeVesselPowerSpectrum_eLife2020(animalID,saveFigs,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the spectral power of arteriole diameter D/D (2PLSM)
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs = {'T115','T116','T117','T118','T125','T126'};
modelType = 'Manual';
params.minTime.Rest = 10;
params.minTime.NREM = 30;
params.minTime.REM = 60;
%% only run analysis for valid animal IDs
if any(strcmp(animalIDs,animalID))
    dataLocation = [rootFolder '/' animalID '/2P Data/'];
    cd(dataLocation)
    % character list of all MergedData files
    mergedDirectory = dir('*_MergedData.mat');
    mergedDataFiles = {mergedDirectory.name}';
    mergedDataFileIDs = char(mergedDataFiles);
    % character list of all TrainingData files
    trainingDirectory = dir('*_TrainingData.mat');
    trainingDataFiles = {trainingDirectory.name}';
    trainingDataFileIDs = char(trainingDataFiles);
    % find and load EventData.mat struct
    eventDataFileStruct = dir('*_EventData.mat');
    eventDataFile = {eventDataFileStruct.name}';
    eventDataFileID = char(eventDataFile);
    load(eventDataFileID)
    % find and load RestData.mat struct
    restDataFileStruct = dir('*_RestData.mat');
    restDataFile = {restDataFileStruct.name}';
    restDataFileID = char(restDataFile);
    load(restDataFileID)
    % find and load manual baseline event information
    manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
    manualBaselineFile = {manualBaselineFileStruct.name}';
    manualBaselineFileID = char(manualBaselineFile);
    load(manualBaselineFileID)
    % find and load RestingBaselines.mat strut
    baselineDataFileStruct = dir('*_RestingBaselines.mat');
    baselineDataFile = {baselineDataFileStruct.name}';
    baselineDataFileID = char(baselineDataFile);
    load(baselineDataFileID)
    % find and load SleepData.mat strut
    sleepDataFileStruct = dir('*_SleepData.mat');
    sleepDataFile = {sleepDataFileStruct.name}';
    sleepDataFileID = char(sleepDataFile);
    load(sleepDataFileID)
    % lowpass filter
    samplingRate = RestData.vesselDiameter.data.samplingRate;
    [z,p,k] = butter(4,1/(samplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    % criteria for resting
    RestCriteria.Fieldname = {'durations'};
    RestCriteria.Comparison = {'gt'};
    RestCriteria.Value = {params.minTime.Rest};
    %% analyze power spectra during periods of rest
    % pull data from RestData.mat structure
    [restLogical] = FilterEvents_2P_eLife2020(RestData.vesselDiameter.data,RestCriteria);
    restLogical = logical(restLogical);
    restingData = RestData.vesselDiameter.data.data(restLogical,:);
    restFileIDs = RestData.vesselDiameter.data.fileIDs(restLogical,:);
    restVesselIDs = RestData.vesselDiameter.data.vesselIDs(restLogical,:);
    restEventTimes = RestData.vesselDiameter.data.eventTimes(restLogical,:);
    restDurations = RestData.vesselDiameter.data.durations(restLogical,:);
    % keep only the data that occurs within the manually-approved awake regions
    [finalRestData,finalRestFileIDs,finalRestVesselIDs,~,~] = RemoveInvalidData_2P_eLife2020(restingData,restFileIDs,restVesselIDs,restDurations,restEventTimes,ManualDecisions);
    % filter, detrend, and truncate data to minimum length to match events
    for aa = 1:length(finalRestData)
        restStrDay = ConvertDate_2P_eLife2020(finalRestFileIDs{aa,1}(1:6));
        if length(finalRestData{aa,1}) < params.minTime.Rest*samplingRate
            restChunkSampleDiff = params.minTime.Rest*samplingRate - length(finalRestData{aa,1});
            restPad = (ones(1,restChunkSampleDiff))*finalRestData{aa,1}(end);
            arrayRestData = horzcat(finalRestData{aa,1},restPad); %#ok<*AGROW>
            normRestData = (arrayRestData - RestingBaselines.vesselDiameter.data.(finalRestVesselIDs{aa,1}).(restStrDay))./RestingBaselines.vesselDiameter.data.(finalRestVesselIDs{aa,1}).(restStrDay);
            procRestData{aa,1} = filtfilt(sos,g,detrend(normRestData,'constant'));
        else
            normRestData = (finalRestData{aa,1}(1:(params.minTime.Rest*samplingRate)) - RestingBaselines.manualSelection.vesselDiameter.data.(finalRestVesselIDs{aa,1}).(restStrDay))./RestingBaselines.manualSelection.vesselDiameter.data.(finalRestVesselIDs{aa,1}).(restStrDay);
            procRestData{aa,1} = filtfilt(sos,g,detrend(normRestData,'constant'));
        end
    end
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
    reshapedRestData = zeros(length(procRestData{1,1}),length(procRestData));
    for bb = 1:length(procRestData)
        reshapedRestData(:,bb) = procRestData{bb,1};
    end
    % remove veins from the artery list
    uniqueRestVesselIDs = unique(finalRestVesselIDs);
    dd = 1;
    for cc = 1:length(uniqueRestVesselIDs)
        if strcmp(uniqueRestVesselIDs{cc,1}(1),'V') == false
            restArterioleIDs{dd,1} = uniqueRestVesselIDs{cc,1};
            dd = dd + 1;
        end
    end
    % split the data based on different arteries
    for ee = 1:length(restArterioleIDs)
        restArterioleID = restArterioleIDs{ee,1};
        gg = 1;
        for ff = 1:length(finalRestVesselIDs)
            if strcmp(restArterioleID,finalRestVesselIDs{ff,1}) == true
                restArterioleData.(restArterioleID)(:,gg) = reshapedRestData(:,ff);
                gg = gg + 1;
            end
        end
    end
    % parameters for mtspectrumc - information available in function
    params.tapers = [1,1];   % Tapers [n, 2n - 1]
    params.pad = 1;
    params.Fs = samplingRate;
    params.fpass = [0,0.5];   % Pass band [0, nyquist]
    params.trialave = 1;
    params.err = [2,0.05];
    % calculate the power spectra of the desired signals
    for hh = 1:length(restArterioleIDs)
        restVID = restArterioleIDs{hh,1};
        [rest_S{hh,1},rest_f{hh,1},rest_sErr{hh,1}] = mtspectrumc_eLife2020(restArterioleData.(restVID),params);
        % save results
        AnalysisResults.(animalID).PowerSpectra.Rest.(restVID).S = rest_S{hh,1};
        AnalysisResults.(animalID).PowerSpectra.Rest.(restVID).f = rest_f{hh,1};
        AnalysisResults.(animalID).PowerSpectra.Rest.(restVID).sErr = rest_sErr{hh,1};
        % save figures if desired
        if strcmp(saveFigs,'y') == true
            RestPower = figure;
            loglog(rest_f{hh,1},rest_S{hh,1},'k')
            hold on;
            loglog(rest_f{hh,1},rest_sErr{hh,1},'color',colors_eLife2020('battleship grey'))
            xlabel('Freq (Hz)');
            ylabel('Power');
            title([animalID  ' ' restVID ' vessel diameter power during awake rest']);
            set(gca,'Ticklength',[0,0]);
            xlim([0,1])
            axis square
            set(gca,'box','off')
            [pathstr,~,~] = fileparts(cd);
            dirpath = [pathstr '/Figures/Vessel Power Spectrum/'];
            if ~exist(dirpath,'dir')
                mkdir(dirpath);
            end
            savefig(RestPower,[dirpath animalID '_' restVID '_RestPowerSpectra']);
            close(RestPower)
        end
    end
    %% analyze power spectra during periods of alert
    awakeData = [];
    for aa = 1:size(mergedDataFileIDs,1)
        mergedDataFileID = mergedDataFileIDs(aa,:);
        [~,~,fileDate,~,~,vesselID] = GetFileInfo2_2P_eLife2020(mergedDataFileID);
        if strcmp(vesselID(1),'V') == false
            load(mergedDataFileID,'-mat')
            trainingDataFileID = trainingDataFileIDs(aa,:);
            load(trainingDataFileID,'-mat')
            strDay = ConvertDate_2P_eLife2020(fileDate);
            if sum(strcmp(trainingTable.behavState,'Not Sleep')) > 144
                if isfield(awakeData,vesselID) == false
                    awakeData.(vesselID) = [];
                    binWhisking.(vesselID) = [];
                end
                % filter and detrend data
                vesselDiam = MergedData.data.vesselDiameter.data;
                normVesselDiam = filtfilt(sos,g,detrend((vesselDiam - RestingBaselines.manualSelection.vesselDiameter.data.(vesselID).(strDay))./ RestingBaselines.manualSelection.vesselDiameter.data.(vesselID).(strDay),'constant'));
                awakeData.(vesselID) = horzcat(awakeData.(vesselID),normVesselDiam');
                binWhisk = MergedData.data.binWhiskerAngle;
                [linkedBinarizedWhiskers] = LinkBinaryEvents_2P_eLife2020(gt(binWhisk,0),[round(30/3),0]);
                binWhiskingPercent = sum(linkedBinarizedWhiskers)/length(linkedBinarizedWhiskers)*100 ;
                binWhisking.(vesselID) = horzcat(binWhisking.(vesselID),binWhiskingPercent);
            end
        end
    end
    % parameters for mtspectrumc - information available in function
    params.tapers = [10,19];   % Tapers [n, 2n - 1]
    params.pad = 1;
    params.Fs = samplingRate;
    params.fpass = [0,0.5];   % Pass band [0, nyquist]
    params.trialave = 1;
    params.err = [2,0.05];
    if isempty(awakeData) == false
        awakeDataVesselIDs = fieldnames(awakeData);
        for bb = 1:length(awakeDataVesselIDs)
            awakeDataVID = awakeDataVesselIDs{bb,1};
            [awakeData_S{bb,1},awakeData_f{bb,1},awakeData_sErr{bb,1}] = mtspectrumc_eLife2020(awakeData.(awakeDataVID),params);
            % save results
            AnalysisResults.(animalID).PowerSpectra.Awake.(awakeDataVID).whiskingPerc = mean(binWhisking.(awakeDataVID));
            AnalysisResults.(animalID).PowerSpectra.Awake.(awakeDataVID).S = awakeData_S{bb,1};
            AnalysisResults.(animalID).PowerSpectra.Awake.(awakeDataVID).f = awakeData_f{bb,1};
            AnalysisResults.(animalID).PowerSpectra.Awake.(awakeDataVID).sErr = awakeData_sErr{bb,1};
            % save figures if desired
            if strcmp(saveFigs,'y') == true
                awakeDataPower = figure;
                loglog(awakeData_f{bb,1},awakeData_S{bb,1},'k')
                hold on;
                loglog(awakeData_f{bb,1},awakeData_sErr{bb,1},'color',colors_eLife2020('battleship grey'))
                xlabel('Freq (Hz)');
                ylabel('Power');
                title([animalID  ' ' awakeDataVID ' vessel diameter power during all awake data']);
                set(gca,'Ticklength',[0,0]);
                xlim([0,1])
                axis square
                set(gca,'box','off')
                [pathstr,~,~] = fileparts(cd);
                dirpath = [pathstr '/Figures/Vessel Power Spectrum/'];
                if ~exist(dirpath,'dir')
                    mkdir(dirpath);
                end
                savefig(awakeDataPower,[dirpath animalID '_' awakeDataVID '_AllDataPowerSpectra']);
                close(awakeDataPower)
            end
        end
    else
        % save results
        AnalysisResults.(animalID).PowerSpectra.Awake = [];
    end
    %% analyze power spectra during periods of all data
    allData = [];
    for aa = 1:size(mergedDataFileIDs,1)
        mergedDataFileID = mergedDataFileIDs(aa,:);
        [~,~,fileDate,~,~,vesselID] = GetFileInfo2_2P_eLife2020(mergedDataFileID);
        if strcmp(vesselID(1),'V') == false
            load(mergedDataFileID,'-mat')
            strDay = ConvertDate_2P_eLife2020(fileDate);
            if isfield(allData,vesselID) == false
                allData.(vesselID) = [];
                binWhisking.(vesselID) = [];
            end
            % filter and detrend data
            vesselDiam = MergedData.data.vesselDiameter.data;
            normVesselDiam = filtfilt(sos,g,detrend((vesselDiam - RestingBaselines.manualSelection.vesselDiameter.data.(vesselID).(strDay))./ RestingBaselines.manualSelection.vesselDiameter.data.(vesselID).(strDay),'constant'));
            allData.(vesselID) = horzcat(allData.(vesselID),normVesselDiam');
            binWhisk = MergedData.data.binWhiskerAngle;
            [linkedBinarizedWhiskers] = LinkBinaryEvents_2P_eLife2020(gt(binWhisk,0),[round(30/3),0]);
            binWhiskingPercent = sum(linkedBinarizedWhiskers)/length(linkedBinarizedWhiskers)*100 ;
            binWhisking.(vesselID) = horzcat(binWhisking.(vesselID),binWhiskingPercent);
        end
    end
    % parameters for mtspectrumc - information available in function
    params.tapers = [10,19];   % Tapers [n, 2n - 1]
    params.pad = 1;
    params.Fs = samplingRate;
    params.fpass = [0,0.5];   % Pass band [0, nyquist]
    params.trialave = 1;
    params.err = [2,0.05];
    if isempty(allData) == false
        allDataVesselIDs = fieldnames(allData);
        for bb = 1:length(allDataVesselIDs)
            allDataVID = allDataVesselIDs{bb,1};
            [allData_S{bb,1},allData_f{bb,1},allData_sErr{bb,1}] = mtspectrumc_eLife2020(allData.(allDataVID),params);
            % save results
            AnalysisResults.(animalID).PowerSpectra.All.(allDataVID).whiskingPerc = mean(binWhisking.(allDataVID));
            AnalysisResults.(animalID).PowerSpectra.All.(allDataVID).S = allData_S{bb,1};
            AnalysisResults.(animalID).PowerSpectra.All.(allDataVID).f = allData_f{bb,1};
            AnalysisResults.(animalID).PowerSpectra.All.(allDataVID).sErr = allData_sErr{bb,1};
            % save figures if desired
            if strcmp(saveFigs,'y') == true
                allDataPower = figure;
                loglog(allData_f{bb,1},allData_S{bb,1},'k')
                hold on;
                loglog(allData_f{bb,1},allData_sErr{bb,1},'color',colors_eLife2020('battleship grey'))
                xlabel('Freq (Hz)');
                ylabel('Power');
                title([animalID  ' ' allDataVID ' vessel diameter power during all awake data']);
                set(gca,'Ticklength',[0,0]);
                xlim([0,1])
                axis square
                set(gca,'box','off')
                [pathstr,~,~] = fileparts(cd);
                dirpath = [pathstr '/Figures/Vessel Power Spectrum/'];
                if ~exist(dirpath,'dir')
                    mkdir(dirpath);
                end
                savefig(allDataPower,[dirpath animalID '_' allDataVID '_AllDataPowerSpectra']);
                close(allDataPower)
            end
        end
    end
    %% analyze power spectra during periods of NREM
    if isfield(SleepData.(modelType),'NREM') == true
        % pull data from SleepData.mat structure
        nremData = SleepData.(modelType).NREM.data.vesselDiameter.data;
        nremVesselIDs = SleepData.(modelType).NREM.VesselIDs;
        % filter, detrend, and truncate data to minimum length to match events
        for aa = 1:length(nremData)
            nremData{aa,1} = filtfilt(sos,g,detrend(nremData{aa,1}(1:(params.minTime.NREM*samplingRate)),'constant'));
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        nremReshape = zeros(length(nremData{1,1}),length(nremData));
        for bb = 1:length(nremData)
            nremReshape(:,bb) = nremData{bb,1};
        end
        % remove veins from the artery list
        uniqueNREMVesselIDs = unique(nremVesselIDs);
        dd = 1;
        for cc = 1:length(uniqueNREMVesselIDs)
            if strcmp(uniqueNREMVesselIDs{cc,1}(1),'V') == false
                nremArterioleIDs{dd,1} = uniqueNREMVesselIDs{cc,1};
                dd = dd + 1;
            end
        end
        % split the data based on different arteries
        for ee = 1:length(nremArterioleIDs)
            nremArterioleID = nremArterioleIDs{ee,1};
            gg = 1;
            for ff = 1:length(nremVesselIDs)
                if strcmp(nremArterioleID,nremVesselIDs{ff,1}) == true
                    nremArterioleData.(nremArterioleID)(:,gg) = nremReshape(:,ff);
                    gg = gg + 1;
                end
            end
        end
        % parameters for mtspectrumc - information available in function
        params.tapers = [3,5];   % Tapers [n, 2n - 1]
        params.pad = 1;
        params.Fs = samplingRate;
        params.fpass = [0,0.5];   % Pass band [0, nyquist]
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the power spectra of the desired signals
        for bb = 1:length(nremArterioleIDs)
            nremVID = nremArterioleIDs{bb,1};
            [nrem_S{bb,1},nrem_f{bb,1},nrem_sErr{bb,1}] = mtspectrumc_eLife2020(nremArterioleData.(nremVID),params);
            % save results
            AnalysisResults.(animalID).PowerSpectra.NREM.(nremVID).S = nrem_S{bb,1};
            AnalysisResults.(animalID).PowerSpectra.NREM.(nremVID).f = nrem_f{bb,1};
            AnalysisResults.(animalID).PowerSpectra.NREM.(nremVID).sErr = nrem_sErr{bb,1};
            % save figures if desired
            if strcmp(saveFigs,'y') == true
                nremPower = figure;
                loglog(nrem_f{bb,1},nrem_S{bb,1},'k')
                hold on;
                loglog(nrem_f{bb,1},nrem_sErr{bb,1},'color',colors_eLife2020('battleship grey'))
                xlabel('Freq (Hz)');
                ylabel('Power');
                title([animalID  ' ' nremVID ' vessel diameter power during NREM']);
                set(gca,'Ticklength',[0,0]);
                xlim([0,1])
                axis square
                set(gca,'box','off')
                [pathstr,~,~] = fileparts(cd);
                dirpath = [pathstr '/Figures/Vessel Power Spectrum/'];
                if ~exist(dirpath,'dir')
                    mkdir(dirpath);
                end
                savefig(nremPower,[dirpath animalID '_' nremVID '_NREMPowerSpectra']);
                close(nremPower)
            end
        end
    end
    %% analyze power spectra during periods of REM
    if isfield(SleepData.(modelType),'REM') == true
        % pull data from SleepData.mat structure
        remData = SleepData.(modelType).REM.data.vesselDiameter.data;
        remVesselIDs = SleepData.(modelType).REM.VesselIDs;
        % filter, detrend, and truncate data to minimum length to match events
        for aa = 1:length(remData)
            remDataA{aa,1} = filtfilt(sos,g,detrend(remData{aa,1}(1:(params.minTime.REM*samplingRate)),'constant'));
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        remReshape = zeros(length(remDataA{1,1}),length(remDataA));
        for bb = 1:length(remDataA)
            remReshape(:,bb) = remDataA{bb,1};
        end
        % remove veins from the artery list
        uniqueREMVesselIDs = unique(remVesselIDs);
        dd = 1;
        for cc = 1:length(uniqueREMVesselIDs)
            if strcmp(uniqueREMVesselIDs{cc,1}(1),'V') == false
                remArterioleIDs{dd,1} = uniqueREMVesselIDs{cc,1};
                dd = dd + 1;
            end
        end
        % split the data based on different arteries
        for ee = 1:length(remArterioleIDs)
            remArterioleID = remArterioleIDs{ee,1};
            gg = 1;
            for ff = 1:length(remVesselIDs)
                if strcmp(remArterioleID,remVesselIDs{ff,1}) == true
                    remArterioleData.(remArterioleID)(:,gg) = remReshape(:,ff);
                    gg = gg + 1;
                end
            end
        end
        % parameters for mtspectrumc - information available in function
        params.tapers = [5,9];   % Tapers [n, 2n - 1]
        params.pad = 1;
        params.Fs = samplingRate;
        params.fpass = [0,0.5];   % Pass band [0, nyquist]
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the power spectra of the desired signals
        for bb = 1:length(remArterioleIDs)
            remVID = remArterioleIDs{bb,1};
            [rem_S{bb,1},rem_f{bb,1},rem_sErr{bb,1}] = mtspectrumc_eLife2020(remArterioleData.(remVID),params);
            % save results
            AnalysisResults.(animalID).PowerSpectra.REM.(remVID).S = rem_S{bb,1};
            AnalysisResults.(animalID).PowerSpectra.REM.(remVID).f = rem_f{bb,1};
            AnalysisResults.(animalID).PowerSpectra.REM.(remVID).sErr = rem_sErr{bb,1};
            % save figures if desired
            if strcmp(saveFigs,'y') == true
                remPower = figure;
                loglog(rem_f{bb,1},rem_S{bb,1},'k')
                hold on;
                loglog(rem_f{bb,1},rem_sErr{bb,1},'color',colors_eLife2020('battleship grey'))
                xlabel('Freq (Hz)');
                ylabel('Power');
                title([animalID  ' ' remVID ' vessel diameter power during REM']);
                set(gca,'Ticklength',[0,0]);
                xlim([0,1])
                axis square
                set(gca,'box','off')
                [pathstr,~,~] = fileparts(cd);
                dirpath = [pathstr '/Figures/Vessel Power Spectrum/'];
                if ~exist(dirpath,'dir')
                    mkdir(dirpath);
                end
                savefig(remPower,[dirpath animalID '_' remVID '_REMPowerSpectra']);
                close(remPower)
            end
        end
    end
end
% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults')

end
