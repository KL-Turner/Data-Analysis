function [AnalysisResults] = AnalyzePowerSpectrum_IOS(dataType, params, AnalysisResults)
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

% list of all Procdata.mat files
procdataFileStruct = dir('*_Procdata.mat');
procdataFiles = {procdataFileStruct.name}';
procdataFileIDs = char(procdataFiles);

% find and load RestData.mat struct
restdataFileStruct = dir('*_Restdata.mat');
restdataFile = {restdataFileStruct.name}';
restdataFileID = char(restdataFile);
load(restdataFileID)

% find and load SleepData.mat strut
sleepdataFileStruct = dir('*_Sleepdata.mat');
sleepdataFile = {sleepdataFileStruct.name}';
sleepdataFileID = char(sleepdataFile);
load(sleepdataFileID)

% Parameters for coherencyc_IOS - information available in function
samplingRate = RestData.CBV.LH.CBVCamSamplingRate;
params.tapers = [3 5];   % Tapers [n, 2n - 1]
params.pad = 1;
params.Fs = samplingRate;   % Sampling Rate
params.fpass = [0 1];   % Pass band [0, nyquist] 
params.trialave = 1;
params.err = [2 0.05];

fileBreaks = strfind(restdataFileID, '_');
animalID = restdataFileID(1:fileBreaks(1)-1);
trialDurationMin = 15;   % min

%% Load in relevant data from the RestData struct:
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};

PuffCriteria.Fieldname = {'puffDistances'};
PuffCriteria.Comparison = {'gt'};
PuffCriteria.Value = {5};

% Use the RestCriteria we specified earlier to find all resting events that are greater than the criteria
if strcmp(dataType, 'CBV')
    [restLogical] = FilterEvents_IOS(RestData.(dataType).LH, RestCriteria);
    [puffLogical] = FilterEvents_IOS(RestData.(dataType).LH, PuffCriteria);
    combRestLogical = logical(restLogical.*puffLogical);
    allRestFiles = RestData.(dataType).LH.fileIDs(combRestLogical, :);
    LH_allRestingData = RestData.(dataType).LH.NormData(combRestLogical, :);
    RH_allRestingData = RestData.(dataType).RH.NormData(combRestLogical, :);
else
    [restLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType), RestCriteria);
    [puffLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType), PuffCriteria);
    combRestLogical = logical(restLogical.*puffLogical);
    allRestFiles = RestData.cortical_LH.(dataType).fileIDs(combRestLogical, :);
    LH_allRestingData =RestData.cortical_LH.(dataType).NormData(combRestLogical, :);
    RH_allRestingData = RestData.cortical_RH.(dataType).NormData(combRestLogical, :);
    Hip_allRestingData = RestData.hippocampus.(dataType).NormData(combRestLogical, :);
end
    
uniqueDays = GetUniqueDays_IOS(allRestFiles);   % Find the unique days of imaging
uniqueFiles = unique(allRestFiles);   % Find the unique files from the filelist. This removes duplicates
% since most files have more than one resting event
numberOfFiles = length(unique(allRestFiles));   % Find the number of unique files
fileTarget = params.targetMinutes/trialDurationMin;

% Loop through each unique day in order to create a logical to filter the file list so that it only includes the first
% x number of files that fall within the targetMinutes requirement
for uD = 1:length(uniqueDays)
    day = uniqueDays(uD);
    x = 1;
    for nOF = 1:numberOfFiles
        file = uniqueFiles(nOF);
        fileID = file{1}(1:6);
        if strcmp(day, fileID) && x <= fileTarget
            filtLogical{uD, 1}(nOF, 1) = 1;
            x = x + 1;
        else
            filtLogical{uD, 1}(nOF, 1) = 0;
        end
    end
end

% Combine the 3 logicals so that it reflects the first "x" number of files from each day
finalLogical = any(sum(cell2mat(filtLogical'), 2), 2);

% Now that the appropriate files from each day are identified, loop through each file name with respect to the original
% list of all relevant files, only keeping the ones that fall within the first targetMinutes of each day.
filtRestFiles = uniqueFiles(finalLogical, :);
for rF = 1:length(allRestFiles)
    logic = strcmp(allRestFiles{rF}, filtRestFiles);
    logicSum = sum(logic);
    if logicSum == 1
        fileFilter(rF, 1) = 1;
    else
        fileFilter(rF, 1) = 0;
    end
end

finalFileFilter = logical(fileFilter);

LH_finalRestData = LH_allRestingData(finalFileFilter,:);
RH_finalRestData = RH_allRestingData(finalFileFilter,:);
if strcmp(dataType, 'CBV') == false
    Hip_finalRestData = Hip_allRestingData(finalFileFilter,:);
end

% Take only the first n number of samples (for example, 10 seconds worth) based on the minimum resting length. Filter that data below 2 Hz (CBV)
% or 0.5 Hz (Gamma Band Power) and then detrend it.
[B, A] = butter(4,1/(samplingRate/2), 'low');
for ii = 1:length(LH_finalRestData)
    LH_finalRestData{ii, 1} = detrend(filtfilt(B, A, LH_finalRestData{ii, 1}(1:(params.minTime.Rest*samplingRate))), 'constant');
    RH_finalRestData{ii, 1} = detrend(filtfilt(B, A, RH_finalRestData{ii, 1}(1:(params.minTime.Rest*samplingRate))), 'constant');
    if strcmp(dataType, 'CBV') == false
        Hip_finalRestData{ii, 1} = detrend(filtfilt(B, A, Hip_finalRestData{ii, 1}(1:(params.minTime.Rest*samplingRate))), 'constant');
    end
end

% Input data as time(1st dimension, vertical) by trials (2nd dimension, horizontally) for coherencyc_IOS and calc. coherence using chronux toolbox.
% Pre-allocate sizes
LH_restData = zeros(length(LH_finalRestData{1, 1}), length(LH_finalRestData));
RH_restData = zeros(length(RH_finalRestData{1, 1}), length(RH_finalRestData));
if strcmp(dataType, 'CBV') == false
    Hip_restData = zeros(length(Hip_finalRestData{1, 1}), length(Hip_finalRestData));
end

for n = 1:length(LH_finalRestData)
    LH_restData(:, n) = LH_finalRestData{n, 1};
    RH_restData(:, n) = RH_finalRestData{n, 1};
    if strcmp(dataType, 'CBV') == false
        Hip_restData(:, n) = Hip_finalRestData{n, 1};
    end
end


% Calculate the power spectrum of the desired signals and save those signals in a the comparison data structure.
disp(['Analyzing the power spectrum of the LH RestData ' dataType ' signal power...']); disp(' ')
[LH_rest_S, LH_rest_f, LH_rest_sErr] = mtspectrumc_IOS(LH_restData, params);
disp(['Analyzing the power spectrum of the RH RestData ' dataType ' signal power...']); disp(' ')
[RH_rest_S, RH_rest_f, RH_rest_sErr] = mtspectrumc_IOS(RH_restData, params);
if strcmp(dataType, 'CBV') == false
    disp(['Analyzing the power spectrum of the Hippocampal RestData ' dataType ' signal power...']); disp(' ')
    [Hip_rest_S, Hip_rest_f, Hip_rest_sErr] = mtspectrumc_IOS(Hip_restData, params);
end

% nboot = 1000;
% LH_restCI = bootci(nboot, @mtspectrumc2, LH_restData', params);
% RH_restCI = bootci(nboot, @mtspectrumc2, RH_restData', params);
% if strcmp(dataType, 'CBV') == false
%     Hip_restCI = bootci(nboot, @mtspectrumc2, Hip_restData', params);
% end

LH_RestPower = figure;
loglog(LH_rest_f, LH_rest_S)
hold on;
loglog(LH_rest_f, LH_rest_sErr)
% hold on;
% loglog(LH_rest_f, restLH_CI)
xlabel('Frequency (Hz)');
ylabel('Power');
title([animalID  ' LH ' dataType ' Power during awake rest']);
set(gca, 'Ticklength', [0 0]);
legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Location', 'Southeast');
set(legend, 'FontSize', 6);
xlim([0 1])
axis square

RH_RestPower = figure;
loglog(RH_rest_f, RH_rest_S)
hold on;
loglog(RH_rest_f, RH_rest_sErr)
% hold on;
% loglog(RH_rest_f, restRH_CI)
xlabel('Frequency (Hz)');
ylabel('Power');
title([animalID  ' RH ' dataType ' Power during awake rest']);
set(gca, 'Ticklength', [0 0]);
legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Location', 'Southeast');
set(legend, 'FontSize', 6);
xlim([0 1])
axis square

if strcmp(dataType, 'CBV') == false    
    Hip_RestPower = figure;
    loglog(Hip_rest_f, Hip_rest_S)
    hold on;
    loglog(Hip_rest_f, Hip_rest_sErr)
    % hold on;
    % loglog(Hip_rest_f, restHip_CI)
    xlabel('Frequency (Hz)');
    ylabel('Power');
    title([animalID  ' Hippocampal ' dataType ' Power during awake rest']);
    set(gca, 'Ticklength', [0 0]);
    legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Location', 'Southeast');
    set(legend, 'FontSize', 6);
    xlim([0 1])
    axis square
end

AnalysisResults.PowerSpectra.Rest.(dataType).LH.S = LH_rest_S;
AnalysisResults.PowerSpectra.Rest.(dataType).LH.f = LH_rest_f;
AnalysisResults.PowerSpectra.Rest.(dataType).LH.sErr = LH_rest_sErr;

AnalysisResults.PowerSpectra.Rest.(dataType).RH.S = RH_rest_S;
AnalysisResults.PowerSpectra.Rest.(dataType).RH.f = RH_rest_f;
AnalysisResults.PowerSpectra.Rest.(dataType).RH.sErr = RH_rest_sErr;

if strcmp(dataType, 'CBV') == false
    AnalysisResults.PowerSpectra.Rest.(dataType).Hip.S = Hip_rest_S;
    AnalysisResults.PowerSpectra.Rest.(dataType).Hip.f = Hip_rest_f;
    AnalysisResults.PowerSpectra.Rest.(dataType).Hip.sErr = Hip_rest_sErr;
end

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Analysis Power Spectra/'];

if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end

savefig(LH_RestPower, [dirpath animalID '_Rest_LH_' dataType '_PowerSpectra']);
savefig(RH_RestPower, [dirpath animalID '_Rest_RH_' dataType '_PowerSpectra']);
if strcmp(dataType, 'CBV') == false
    savefig(Hip_RestPower, [dirpath animalID '_Rest_Hippocampal_' dataType '_PowerSpectra']);
end

%% Load in relevant data from SleepData
if strcmp(dataType, 'CBV')
    LH_nremData = SleepData.NREM.data.(dataType).LH;
    RH_nremData = SleepData.NREM.data.(dataType).RH;
else
    LH_nremData = SleepData.NREM.data.cortical_LH.(dataType);
    RH_nremData = SleepData.NREM.data.cortical_RH.(dataType);
    Hip_nremData = SleepData.NREM.data.hippocampus.(dataType);
end

for ii = 1:length(LH_nremData)
    LH_nremData{ii, 1} = detrend(LH_nremData{ii, 1}(1:(params.minTime.NREM*samplingRate)), 'constant');
    RH_nremData{ii, 1} = detrend(RH_nremData{ii, 1}(1:(params.minTime.NREM*samplingRate)), 'constant');
    if strcmp(dataType, 'CBV') == false
        Hip_nremData{ii, 1} = detrend(Hip_nremData{ii, 1}(1:(params.minTime.NREM*samplingRate)), 'constant');
    end
end

LH_nrem = zeros(length(LH_nremData{1, 1}), length(LH_nremData));
RH_nrem = zeros(length(RH_nremData{1, 1}), length(RH_nremData));
if strcmp(dataType, 'CBV') == false
    Hip_nrem = zeros(length(Hip_nremData{1, 1}), length(Hip_nremData));
end

for n = 1:length(LH_nremData)
    LH_nrem(:, n) = LH_nremData{n, 1};
    RH_nrem(:, n) = RH_nremData{n, 1};
    if strcmp(dataType, 'CBV') == false
        Hip_nrem(:, n) = Hip_nremData{n, 1};
    end
end

% Calculate the power spectrum of the desired signals and save those signals in a the comparison data structure.
disp(['Analyzing the power spectrum of the LH NREM ' dataType ' signal power...']); disp(' ')
[LH_nrem_S, LH_nrem_f, LH_nrem_sErr] = mtspectrumc_IOS(LH_nrem, params);
disp(['Analyzing the power spectrum of the RH NREM ' dataType ' signal power...']); disp(' ')
[RH_nrem_S, RH_nrem_f, RH_nrem_sErr] = mtspectrumc_IOS(RH_nrem, params);
if strcmp(dataType, 'CBV') == false
    disp(['Analyzing the power spectrum of the Hippocampal NREM ' dataType ' signal power...']); disp(' ')
    [Hip_nrem_S, Hip_nrem_f, Hip_nrem_sErr] = mtspectrumc_IOS(Hip_nrem, params);
end

% nboot = 1000;
% LH_nremCI = bootci(nboot, @mtspectrumc2, LH_nrem', params);
% RH_nremCI = bootci(nboot, @mtspectrumc2, RH_nrem', params);
% if strcmp(dataType, 'CBV') == false
%     Hip_nremCI = bootci(nboot, @mtspectrumc2, Hip_nrem', params);
% end

LH_nremPower = figure;
loglog(LH_nrem_f, LH_nrem_S)
hold on;
loglog(LH_nrem_f, LH_nrem_sErr)
% hold on;
% loglog(LH_nrem_f, nremLH_CI)
xlabel('Frequency (Hz)');
ylabel('Power');
title([animalID  ' LH ' dataType ' Power during NREM sleep']);
set(gca, 'Ticklength', [0 0]);
legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Location', 'Southeast');
set(legend, 'FontSize', 6);
xlim([0 1])
axis square

RH_nremPower = figure;
loglog(RH_nrem_f, RH_nrem_S)
hold on;
loglog(RH_nrem_f, RH_nrem_sErr)
% hold on;
% loglog(RH_nrem_f, nremRH_CI)
xlabel('Frequency (Hz)');
ylabel('Power');
title([animalID  ' RH ' dataType ' Power during NREM sleep']);
set(gca, 'Ticklength', [0 0]);
legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Location', 'Southeast');
set(legend, 'FontSize', 6);
xlim([0 1])
axis square

if strcmp(dataType, 'CBV') == false    
    Hip_nremPower = figure;
    loglog(Hip_nrem_f, Hip_nrem_S)
    hold on;
    loglog(Hip_nrem_f, Hip_nrem_sErr)
    % hold on;
    % loglog(Hip_nrem_f, nremHip_CI)
    xlabel('Frequency (Hz)');
    ylabel('Power');
    title([animalID  ' Hippocampal ' dataType ' Power during NREM sleep']);
    set(gca, 'Ticklength', [0 0]);
    legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Location', 'Southeast');
    set(legend, 'FontSize', 6);
    xlim([0 1])
    axis square
end

AnalysisResults.PowerSpectra.NREM.(dataType).LH.S = LH_nrem_S;
AnalysisResults.PowerSpectra.NREM.(dataType).LH.f = LH_nrem_f;
AnalysisResults.PowerSpectra.NREM.(dataType).LH.sErr = LH_nrem_sErr;

AnalysisResults.PowerSpectra.NREM.(dataType).RH.S = RH_nrem_S;
AnalysisResults.PowerSpectra.NREM.(dataType).RH.f = RH_nrem_f;
AnalysisResults.PowerSpectra.NREM.(dataType).RH.sErr = RH_nrem_sErr;

if strcmp(dataType, 'CBV') == false
    AnalysisResults.PowerSpectra.NREM.(dataType).Hip.S = Hip_nrem_S;
    AnalysisResults.PowerSpectra.NREM.(dataType).Hip.f = Hip_nrem_f;
    AnalysisResults.PowerSpectra.NREM.(dataType).Hip.sErr = Hip_nrem_sErr;
end

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Analysis Power Spectra/'];

if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end

savefig(LH_nremPower, [dirpath animalID '_NREM_LH_' dataType '_PowerSpectra']);
savefig(RH_nremPower, [dirpath animalID '_NREM_RH_' dataType '_PowerSpectra']);
if strcmp(dataType, 'CBV') == false
    savefig(Hip_nremPower, [dirpath animalID '_NREM_Hippocampal_' dataType '_PowerSpectra']);
end

%% REM
if strcmp(dataType, 'CBV')
    LH_remData = SleepData.REM.data.(dataType).LH;
    RH_remData = SleepData.REM.data.(dataType).RH;
else
    LH_remData = SleepData.REM.data.cortical_LH.(dataType);
    RH_remData = SleepData.REM.data.cortical_RH.(dataType);
    Hip_remData = SleepData.REM.data.hippocampus.(dataType);
end

for ii = 1:length(LH_remData)
    LH_remData{ii, 1} = detrend(LH_remData{ii, 1}(1:(params.minTime.REM*samplingRate)), 'constant');
    RH_remData{ii, 1} = detrend(RH_remData{ii, 1}(1:(params.minTime.REM*samplingRate)), 'constant');
    if strcmp(dataType, 'CBV') == false
        Hip_remData{ii, 1} = detrend(Hip_remData{ii, 1}(1:(params.minTime.REM*samplingRate)), 'constant');
    end
end

LH_rem = zeros(length(LH_remData{1, 1}), length(LH_remData));
RH_rem = zeros(length(RH_remData{1, 1}), length(RH_remData));
if strcmp(dataType, 'CBV') == false
    Hip_rem = zeros(length(Hip_remData{1, 1}), length(Hip_remData));
end

for n = 1:length(LH_remData)
    LH_rem(:, n) = LH_remData{n, 1};
    RH_rem(:, n) = RH_remData{n, 1};
    if strcmp(dataType, 'CBV') == false
        Hip_rem(:, n) = Hip_remData{n, 1};
    end
end

% Calculate the power spectrum of the desired signals and save those signals in a the comparison data structure.
disp(['Analyzing the power spectrum of the LH REM ' dataType ' signal power...']); disp(' ')
[LH_rem_S, LH_rem_f, LH_rem_sErr] = mtspectrumc_IOS(LH_rem, params);
disp(['Analyzing the power spectrum of the RH REM ' dataType ' signal power...']); disp(' ')
[RH_rem_S, RH_rem_f, RH_rem_sErr] = mtspectrumc_IOS(RH_rem, params);
if strcmp(dataType, 'CBV') == false
    disp(['Analyzing the power spectrum of the Hippocampal REM ' dataType ' signal power...']); disp(' ')
    [Hip_rem_S, Hip_rem_f, Hip_rem_sErr] = mtspectrumc_IOS(Hip_rem, params);
end

% nboot = 1000;
% LH_remCI = bootci(nboot, @mtspectrumc2, LH_rem', params);
% RH_remCI = bootci(nboot, @mtspectrumc2, RH_rem', params);
% if strcmp(dataType, 'CBV') == false
%     Hip_remCI = bootci(nboot, @mtspectrumc2, Hip_rem', params);
% end

LH_remPower = figure;
loglog(LH_rem_f, LH_rem_S)
hold on;
loglog(LH_rem_f, LH_rem_sErr)
% hold on;
% loglog(LH_rem_f, remLH_CI)
xlabel('Frequency (Hz)');
ylabel('Power');
title([animalID  ' LH ' dataType ' Power during REM sleep']);
set(gca, 'Ticklength', [0 0]);
legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Location', 'Southeast');
set(legend, 'FontSize', 6);
xlim([0 1])
axis square

RH_remPower = figure;
loglog(RH_rem_f, RH_rem_S)
hold on;
loglog(RH_rem_f, RH_rem_sErr)
% hold on;
% loglog(RH_rem_f, remRH_CI)
xlabel('Frequency (Hz)');
ylabel('Power');
title([animalID  ' RH ' dataType ' Power during REM sleep']);
set(gca, 'Ticklength', [0 0]);
legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Location', 'Southeast');
set(legend, 'FontSize', 6);
xlim([0 1])
axis square

if strcmp(dataType, 'CBV') == false    
    Hip_remPower = figure;
    loglog(Hip_rem_f, Hip_rem_S)
    hold on;
    loglog(Hip_rem_f, Hip_rem_sErr)
    % hold on;
    % loglog(Hip_rem_f, remHip_CI)
    xlabel('Frequency (Hz)');
    ylabel('Power');
    title([animalID  ' Hippocampal ' dataType ' Power during REM sleep']);
    set(gca, 'Ticklength', [0 0]);
    legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Location', 'Southeast');
    set(legend, 'FontSize', 6);
    xlim([0 1])
    axis square
end

AnalysisResults.PowerSpectra.REM.(dataType).LH.S = LH_rem_S;
AnalysisResults.PowerSpectra.REM.(dataType).LH.f = LH_rem_f;
AnalysisResults.PowerSpectra.REM.(dataType).LH.sErr = LH_rem_sErr;

AnalysisResults.PowerSpectra.REM.(dataType).RH.S = RH_rem_S;
AnalysisResults.PowerSpectra.REM.(dataType).RH.f = RH_rem_f;
AnalysisResults.PowerSpectra.REM.(dataType).RH.sErr = RH_rem_sErr;

if strcmp(dataType, 'CBV') == false
    AnalysisResults.PowerSpectra.REM.(dataType).Hip.S = Hip_rem_S;
    AnalysisResults.PowerSpectra.REM.(dataType).Hip.f = Hip_rem_f;
    AnalysisResults.PowerSpectra.REM.(dataType).Hip.sErr = Hip_rem_sErr;
end

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Analysis Power Spectra/'];

if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end

savefig(LH_remPower, [dirpath animalID '_REM_LH_' dataType '_PowerSpectra']);
savefig(RH_remPower, [dirpath animalID '_REM_RH_' dataType '_PowerSpectra']);
if strcmp(dataType, 'CBV') == false
    savefig(Hip_remPower, [dirpath animalID '_REM_Hippocampal_' dataType '_PowerSpectra']);
end

%% Load in relevant data from all ProcdataFiles
% for pDF = 1:size(procdataFiles, 1)
%     procdataFile = procdataFiles(pDF, :);
%     load(procdataFile);
%     [~, ~, fileDate, fileID] = GetFileInfo_IOS(procdataFile);
%     strDay = ConvertDate(fileDate);
%     
%     disp(['Adding data from Procdata file ' num2str(pDF) ' of ' num2str(size(procdataFiles, 1))]); disp(' ')
%     LH_CBV{pDF, 1} = (Procdata.data.CBV.LH - RestingBaselines.(baselineType).CBV.LH.(strDay)) / RestingBaselines.CBV.LH.(strDay);
%     RH_CBV{pDF, 1} = (Procdata.data.CBV.RH - RestingBaselines.(baselineType).CBV.RH.(strDay)) / RestingBaselines.CBV.RH.(strDay);
%     LH_GAM{pDF, 1} = (Procdata.data.GammaBand_Power.LH - RestingBaselines.GammaBand_Power.LH.(strDay)) / RestingBaselines.GammaBand_Power.LH.(strDay);
%     RH_GAM{pDF, 1} = (Procdata.data.GammaBand_Power.RH - RestingBaselines.GammaBand_Power.RH.(strDay)) / RestingBaselines.GammaBand_Power.RH.(strDay);
% end
% 
% for ii = 1:length(LH_CBV)
%     LH_CBV{ii, 1} = detrend(filtfilt(B, A, LH_CBV{ii, 1}), 'constant');
%     RH_CBV{ii, 1} = detrend(filtfilt(B, A, RH_CBV{ii, 1}), 'constant');
%     LH_GAM{ii, 1} = detrend(filtfilt(D, C, LH_GAM{ii, 1}), 'constant');
%     RH_GAM{ii, 1} = detrend(filtfilt(D, C, RH_GAM{ii, 1}), 'constant');
% end
% 
% % Input data as time(1st dimension, vertical) by trials (2nd dimension, horizontally) for coherencyc_IOS and calc. coherence using chronux toolbox.
% % Pre-allocate sizes
% LH_allCBV = zeros(length(LH_CBV{1, 1}), length(LH_CBV));
% RH_allCBV = zeros(length(RH_CBV{1, 1}), length(RH_CBV));
% LH_allGAM = zeros(length(LH_GAM{1, 1}), length(LH_GAM));
% RH_allGAM = zeros(length(RH_GAM{1, 1}), length(RH_GAM));
% 
% for n = 1:length(LH_CBV)
%     LH_allCBV(:, n) = LH_CBV{n, 1};
%     RH_allCBV(:, n) = RH_CBV{n, 1};
%     LH_allGAM(:, n) = LH_GAM{n, 1};
%     RH_allGAM(:, n) = RH_GAM{n, 1};
% end
% 
% % Calculate the Coherency between desired signals and save those signals in a the comparison data structure.
% % disp('Analyzing the Coherence between all L/R CBV signals...'); disp(' ')
% % [C_allCBV, phi_allCBV, S12_allCBV, S1_allCBV, S2_allCBV, f_allCBV, confC_allCBV, phiSTD_allCBV, cErr_allCBV] = coherencyc_IOS(LH_allCBV, RH_allCBV, params);
% % disp('Analyzing the Coherence between all L/R gamma band power signals...'); disp(' ')
% % [C_allGAM, phi_allGAM, S12_allGAM, S1_allGAM, S2_allGAM, f_allGAM, confC_allGAM, phiSTD_allGAM, cErr_allGAM] = coherencyc_IOS(LH_allGAM, RH_allGAM, params);
% % 
% % allCBV_CI = bootci(nboot, @coherencyc2, LH_allCBV', RH_allCBV', params);
% % allGAM_CI = bootci(nboot, @coherencyc2, LH_allGAM', RH_allGAM', params);
% % 
% % allCBVCoherence = figure;
% % loglog(f_allCBV, C_allCBV)
% % hold on;
% % loglog(f_allCBV, cErr_allCBV)
% % hold on;
% % loglog(f_allCBV, allCBV_CI)
% % xlabel('Frequency (Hz)');
% % ylabel('Magnitude of Coherence');
% % title([animalID  ' L/R CBV Coherence for All data']);
% % set(gca, 'Ticklength', [0 0]);
% % legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
% % set(legend, 'FontSize', 6);
% % ylim([0 1])
% % xlim([0 1])
% % axis square
% % 
% % allGAMCoherence = figure;
% % loglog(f_allGAM, C_allGAM)
% % hold on;
% % loglog(f_allGAM, cErr_allGAM)
% % hold on;
% % loglog(f_allGAM, allGAM_CI)
% % xlabel('Frequency (Hz)');
% % ylabel('Magnitude of Coherence');
% % title([animalID  ' L/R Gamma Coherence for All data']);
% % set(gca, 'Ticklength', [0 0]);
% % legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
% % set(legend, 'FontSize', 6);
% % ylim([0 1])
% % xlim([0 1])
% % axis square
% % 
% % Comparisondata.Coherence.Alldata.CBV.C = C_allCBV;
% % Comparisondata.Coherence.Alldata.CBV.f = f_allCBV;
% % Comparisondata.Coherence.Alldata.CBV.cErr = cErr_allCBV;
% % Comparisondata.Coherence.Alldata.GAM.C = C_allGAM;
% % Comparisondata.Coherence.Alldata.GAM.f = f_allGAM;
% % Comparisondata.Coherence.Alldata.GAM.cErr = cErr_allGAM;
% 
% [pathstr, ~, ~] = fileparts(cd);
% dirpath = [pathstr '/Figures/Coherence/'];
% 
% if ~exist(dirpath, 'dir')
%     mkdir(dirpath);
% end
% 
% savefig(allCBVCoherence, [dirpath animalID '_Alldata_' 'CBVCoherece']);
% savefig(allGAMCoherence, [dirpath animalID '_Alldata_' 'GAMCoherece']);
% 
% 

%% save results
save([animalID '_AnalysisResults.mat'], 'AnalysisResults');

end