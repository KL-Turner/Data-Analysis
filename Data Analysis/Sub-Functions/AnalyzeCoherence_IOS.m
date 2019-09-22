function [AnalysisResults] = AnalyzeCoherence_IOS(dataType, params, AnalysisResults)
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

% Take only the first n number of samples (for example, 10 seconds worth) based on the minimum resting length. Filter that data below 2 Hz (CBV)
% or 0.5 Hz (Gamma Band Power) and then detrend it.
[B, A] = butter(4,1/(samplingRate/2), 'low');
for ii = 1:length(LH_finalRestData)
    LH_finalRestData{ii, 1} = detrend(filtfilt(B, A, LH_finalRestData{ii, 1}(1:(params.minTime.Rest*samplingRate))), 'constant');
    RH_finalRestData{ii, 1} = detrend(filtfilt(B, A, RH_finalRestData{ii, 1}(1:(params.minTime.Rest*samplingRate))), 'constant');
end

% Input data as time(1st dimension, vertical) by trials (2nd dimension, horizontally) for coherencyc_IOS and calc. coherence using chronux toolbox.
% Pre-allocate sizes
LH_restData = zeros(length(LH_finalRestData{1, 1}), length(LH_finalRestData));
RH_restData = zeros(length(RH_finalRestData{1, 1}), length(RH_finalRestData));

for n = 1:length(LH_finalRestData)
    LH_restData(:, n) = LH_finalRestData{n, 1};
    RH_restData(:, n) = RH_finalRestData{n, 1};
end

% Calculate the Coherency between desired signals and save those signals in a the comparison data structure.
disp(['Analyzing the resting coherence between L/R ' dataType ' signals...']); disp(' ')
[C_RestData, ~, ~, ~, ~, f_RestData, ~, ~, cErr_RestData] = coherencyc_IOS(LH_restData, RH_restData, params);

% nboot = 1000;
% rest_CI = bootci(nboot, @coherencyc2, LH_restData', RH_restData', params);

RestCoherence = figure;
plot(f_RestData, C_RestData)
hold on;
plot(f_RestData, cErr_RestData)
% hold on;
% plot(f_RestCBV, restCBV_CI)
xlabel('Frequency (Hz)');
ylabel('Magnitude of Coherence');
title([animalID  ' L/R ' dataType ' Coherence for Resting data']);
set(gca, 'Ticklength', [0 0]);
legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
set(legend, 'FontSize', 6);
ylim([0 1])
xlim([0 1])
axis square

AnalysisResults.Coherence.Rest.(dataType).C = C_RestData;
AnalysisResults.Coherence.Rest.(dataType).f = f_RestData;
AnalysisResults.Coherence.Rest.(dataType).cErr = cErr_RestData;

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Analysis Coherence/'];

if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end

savefig(RestCoherence, [dirpath animalID '_Rest_' dataType '_Coherence']);

%% Load in relevant data from SleepData
if strcmp(dataType, 'CBV')
    LH_nremData = SleepData.NREM.data.(dataType).LH;
    RH_nremData = SleepData.NREM.data.(dataType).RH;
else
    LH_nremData = SleepData.NREM.data.cortical_LH.(dataType);
    RH_nremData = SleepData.NREM.data.cortical_RH.(dataType);
end

for ii = 1:length(LH_nremData)
    LH_nremData{ii, 1} = detrend(LH_nremData{ii, 1}(1:(params.minTime.NREM*samplingRate)), 'constant');
    RH_nremData{ii, 1} = detrend(RH_nremData{ii, 1}(1:(params.minTime.NREM*samplingRate)), 'constant');
end

LH_nrem = zeros(length(LH_nremData{1, 1}), length(LH_nremData));
RH_nrem = zeros(length(RH_nremData{1, 1}), length(RH_nremData));

for n = 1:length(LH_nremData)
    LH_nrem(:, n) = LH_nremData{n, 1};
    RH_nrem(:, n) = RH_nremData{n, 1};
end

disp(['Analyzing the NREM coherence between L/R ' dataType ' signals...']); disp(' ')
[C_nrem, ~, ~, ~, ~, f_nrem, ~, ~, cErr_nrem] = coherencyc_IOS(LH_nrem, RH_nrem, params);

% nboot = 1000;
% nrem_CI = bootci(nboot, @coherencyc2, LH_nrem', RH_nrem', params);

nremCoherence = figure;
plot(f_nrem, C_nrem)
hold on;
plot(f_nrem, cErr_nrem)
% hold on;
% plot(f_nremCBV, nremCBV_CI)
xlabel('Frequency (Hz)');
ylabel('Magnitude of Coherence');
title([animalID  ' L/R ' dataType ' Coherence for NREM data']);
set(gca, 'Ticklength', [0 0]);
legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
set(legend, 'FontSize', 6);
ylim([0 1])
xlim([0 1])
axis square

AnalysisResults.Coherence.NREM.(dataType).C = C_nrem;
AnalysisResults.Coherence.NREM.(dataType).f = f_nrem;
AnalysisResults.Coherence.NREM.(dataType).cErr = cErr_nrem;

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Analysis Coherence/'];

if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end

savefig(nremCoherence, [dirpath animalID '_NREM_' dataType '_Coherence']);

%% REM
if strcmp(dataType, 'CBV')
    LH_remData = SleepData.REM.data.(dataType).LH;
    RH_remData = SleepData.REM.data.(dataType).RH;
else
    LH_remData = SleepData.REM.data.cortical_LH.(dataType);
    RH_remData = SleepData.REM.data.cortical_RH.(dataType);
end

for ii = 1:length(LH_remData)
    LH_remData{ii, 1} = detrend(LH_remData{ii, 1}(1:(params.minTime.REM*samplingRate)), 'constant');
    RH_remData{ii, 1} = detrend(RH_remData{ii, 1}(1:(params.minTime.REM*samplingRate)), 'constant');
end

LH_rem = zeros(length(LH_remData{1, 1}), length(LH_remData));
RH_rem = zeros(length(RH_remData{1, 1}), length(RH_remData));

for n = 1:length(LH_remData)
    LH_rem(:, n) = LH_remData{n, 1};
    RH_rem(:, n) = RH_remData{n, 1};
end

disp(['Analyzing the NREM coherence between L/R ' dataType ' signals...']); disp(' ')
[C_rem, ~, ~, ~, ~, f_rem, ~, ~, cErr_rem] = coherencyc_IOS(LH_rem, RH_rem, params);

% nboot = 1000;
% rem_CI = bootci(nboot, @coherencyc2, LH_rem', RH_rem', params);

remCoherence = figure;
plot(f_rem, C_rem)
hold on;
plot(f_rem, cErr_rem)
% hold on;
% plot(f_rem, rem_CI)
xlabel('Frequency (Hz)');
ylabel('Magnitude of Coherence');
title([animalID  ' L/R ' dataType ' Coherence for REM data']);
set(gca, 'Ticklength', [0 0]);
legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
set(legend, 'FontSize', 6);
ylim([0 1])
xlim([0 1])
axis square

AnalysisResults.Coherence.REM.(dataType).C = C_rem;
AnalysisResults.Coherence.REM.(dataType).f = f_rem;
AnalysisResults.Coherence.REM.(dataType).cErr = cErr_rem;

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Analysis Coherence/'];

if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end

savefig(remCoherence, [dirpath animalID '_REM_' dataType '_Coherence']);

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
% % plot(f_allCBV, C_allCBV)
% % hold on;
% % plot(f_allCBV, cErr_allCBV)
% % hold on;
% % plot(f_allCBV, allCBV_CI)
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
% % plot(f_allGAM, C_allGAM)
% % hold on;
% % plot(f_allGAM, cErr_allGAM)
% % hold on;
% % plot(f_allGAM, allGAM_CI)
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