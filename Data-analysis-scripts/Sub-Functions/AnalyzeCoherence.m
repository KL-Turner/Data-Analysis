function [ComparisonData] = AnalyzeCoherence(animal, params, procDataFiles, RestingBaselines, RestData, SleepData, ComparisonData)
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

% Parameters for coherencyc - information available in function
samplingRate = RestData.CBV.LH.CBVCamSamplingRate;
params.tapers = [3 5];   % Tapers [n, 2n - 1]
params.pad = 1;
params.Fs = samplingRate;   % Sampling Rate
params.fpass = [0 1];   % Pass band [0, nyquist] 
params.trialave = 1;
params.err = [2 0.05];

%% Load in relevant data from the RestData struct:
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};

PuffCriteria.Fieldname = {'puffDistances'};
PuffCriteria.Comparison = {'gt'};
PuffCriteria.Value = {5};

% Use the RestCriteria we specified earlier to find all resting events that are greater than the criteria
[restLogical] = FilterEvents(RestData.CBV.LH, RestCriteria);   % Output is a logical
[puffLogical] = FilterEvents(RestData.CBV.LH, PuffCriteria);   % Output is a logical
combRestLogical = logical(restLogical.*puffLogical);
allRestFiles = RestData.CBV.LH.fileIDs(combRestLogical, :);   % Overall logical for all resting file names that meet criteria

LH_allRestingCBVData = RestData.CBV.LH.NormData(combRestLogical, :);   % Pull out data from all those resting files that meet criteria
RH_allRestingCBVData = RestData.CBV.RH.NormData(combRestLogical, :);   % Pull out data from all those resting files that meet criteria
LH_allRestingGAMData = RestData.GammaBand_Power.LH.NormData(combRestLogical, :);   % Pull out data from all those resting files that meet criteria
RH_allRestingGAMData = RestData.GammaBand_Power.RH.NormData(combRestLogical, :);   % Pull out data from all those resting files that meet criteria

uniqueDays = GetUniqueDays(RestData.CBV.LH.fileIDs);   % Find the unique days of imaging
uniqueFiles = unique(RestData.CBV.LH.fileIDs);   % Find the unique files from the filelist. This removes duplicates
% since most files have more than one resting event
numberOfFiles = length(unique(RestData.CBV.LH.fileIDs));   % Find the number of unique files
fileTarget = params.targetMinutes / 5;   % Divide that number of unique files by 5 (minutes) to get the number of files that
% corresponds to the desired targetMinutes

% Loop through each unique day in order to create a logical to filter the file list so that it only includes the first
% x number of files that fall within the targetMinutes requirement
for uD = 1:length(uniqueDays)
    day = uniqueDays(uD);
    x = 1;
    for nOF = 1:numberOfFiles
        file = uniqueFiles(nOF);
        fileID = file{1}(1:6);
        if strcmp(params.Infusion, 'y')
            if strcmp(day, fileID) && x > fileTarget
                filtLogical{uD, 1}(nOF, 1) = 1;
            else
                filtLogical{uD, 1}(nOF, 1) = 0;
                x = x + 1;
            end
        else
            if strcmp(day, fileID) && x <= fileTarget
                filtLogical{uD, 1}(nOF, 1) = 1;
                x = x + 1;
            else
                filtLogical{uD, 1}(nOF, 1) = 0;
            end
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

LH_finalRestCBVData = LH_allRestingCBVData(finalFileFilter, :);
RH_finalRestCBVData = RH_allRestingCBVData(finalFileFilter, :);
LH_finalRestGAMData = LH_allRestingGAMData(finalFileFilter, :);
RH_finalRestGAMData = RH_allRestingGAMData(finalFileFilter, :);

% Take only the first n number of samples (for example, 10 seconds worth) based on the minimum resting length. Filter that data below 2 Hz (CBV)
% or 0.5 Hz (Gamma Band Power) and then detrend it.
[B, A] = butter(4, 2 / (30 / 2), 'low');
[D, C] = butter(4, 0.5 / (30 / 2), 'low');
for ii = 1:length(LH_finalRestCBVData)
    LH_finalRestCBVData{ii, 1} = detrend(filtfilt(B, A, LH_finalRestCBVData{ii, 1}(1:(params.minTime.Rest*samplingRate))), 'constant');
    RH_finalRestCBVData{ii, 1} = detrend(filtfilt(B, A, RH_finalRestCBVData{ii, 1}(1:(params.minTime.Rest*samplingRate))), 'constant');
    LH_finalRestGAMData{ii, 1} = detrend(filtfilt(D, C, LH_finalRestGAMData{ii, 1}(1:(params.minTime.Rest*samplingRate))), 'constant');
    RH_finalRestGAMData{ii, 1} = detrend(filtfilt(D, C, RH_finalRestGAMData{ii, 1}(1:(params.minTime.Rest*samplingRate))), 'constant');
end

% Input data as time(1st dimension, vertical) by trials (2nd dimension, horizontally) for coherencyc and calc. coherence using chronux toolbox.
% Pre-allocate sizes
LH_restCBV = zeros(length(LH_finalRestCBVData{1, 1}), length(LH_finalRestCBVData));
RH_restCBV = zeros(length(RH_finalRestCBVData{1, 1}), length(RH_finalRestCBVData));
LH_restGAM = zeros(length(LH_finalRestGAMData{1, 1}), length(LH_finalRestGAMData));
RH_restGAM = zeros(length(RH_finalRestGAMData{1, 1}), length(RH_finalRestGAMData));

for n = 1:length(LH_finalRestCBVData)
    LH_restCBV(:, n) = LH_finalRestCBVData{n, 1};
    RH_restCBV(:, n) = RH_finalRestCBVData{n, 1};
    LH_restGAM(:, n) = LH_finalRestGAMData{n, 1};
    RH_restGAM(:, n) = RH_finalRestGAMData{n, 1};
end

% Calculate the Coherency between desired signals and save those signals in a the comparison data structure.
disp('Analyzing the Coherence between the RestData L/R CBV signals...'); disp(' ')
[C_RestCBV, phi_RestCBV, S12_RestCBV, S1_RestCBV, S2_RestCBV, f_RestCBV, confC_RestCBV, phiSTD_RestCBV, cErr_RestCBV] = coherencyc(LH_restCBV, RH_restCBV, params);
disp('Analyzing the Coherence between the RestData L/R gamma band power signals...'); disp(' ')
[C_RestGAM, phi_RestGAM, S12_RestGAM, S1_RestGAM, S2_RestGAM, f_RestGAM, confC_RestGAM, phiSTD_RestGAM, cErr_RestGAM] = coherencyc(LH_restGAM, RH_restGAM, params);

nboot = 1000;
restCBV_CI = bootci(nboot, @coherencyc2, LH_restCBV', RH_restCBV', params);
restGAM_CI = bootci(nboot, @coherencyc2, LH_restGAM', RH_restGAM', params);

RestCBVCoherence = figure;
plot(f_RestCBV, C_RestCBV)
hold on;
plot(f_RestCBV, cErr_RestCBV)
hold on;
plot(f_RestCBV, restCBV_CI)
xlabel('Frequency (Hz)');
ylabel('Magnitude of Coherence');
title([animal  ' L/R CBV Coherence for Resting Data']);
set(gca, 'Ticklength', [0 0]);
legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
set(legend, 'FontSize', 6);
ylim([0 1])
xlim([0 1])
axis square

RestGAMCoherence = figure;
plot(f_RestGAM, C_RestGAM)
hold on;
plot(f_RestGAM, cErr_RestGAM)
hold on;
plot(f_RestGAM, restGAM_CI)
xlabel('Frequency (Hz)');
ylabel('Magnitude of Coherence');
title([animal  ' L/R Gamma Coherence for Resting Data']);
set(gca, 'Ticklength', [0 0]);
legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
set(legend, 'FontSize', 6);
ylim([0 1])
xlim([0 1])
axis square

ComparisonData.Coherence.Rest.CBV.C = C_RestCBV;
ComparisonData.Coherence.Rest.CBV.f = f_RestCBV;
ComparisonData.Coherence.Rest.CBV.cErr = cErr_RestCBV;
ComparisonData.Coherence.Rest.GAM.C = C_RestGAM;
ComparisonData.Coherence.Rest.GAM.f = f_RestGAM;
ComparisonData.Coherence.Rest.GAM.cErr = cErr_RestGAM;

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Coherence/'];

if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end

savefig(RestCBVCoherence, [dirpath animal '_Rest_' 'CBVCoherece']);
savefig(RestGAMCoherence, [dirpath animal '_Rest_' 'GAMCoherece']);

%% Load in relevant data from SleepData
if ~isempty(SleepData)
    LH_nremCBV = SleepData.NREM.Data.CBV.LH;
    RH_nremCBV = SleepData.NREM.Data.CBV.RH;
    LH_nremGAM = SleepData.NREM.Data.GammaBand_Power.LH;
    RH_nremGAM = SleepData.NREM.Data.GammaBand_Power.RH;
    
    for ii = 1:length(LH_nremCBV)
        LH_nremCBV{ii, 1} = detrend(LH_nremCBV{ii, 1}(1:(params.minTime.NREM*samplingRate)), 'constant');
        RH_nremCBV{ii, 1} = detrend(RH_nremCBV{ii, 1}(1:(params.minTime.NREM*samplingRate)), 'constant');
        LH_nremGAM{ii, 1} = detrend(LH_nremGAM{ii, 1}(1:(params.minTime.NREM*samplingRate)), 'constant');
        RH_nremGAM{ii, 1} = detrend(RH_nremGAM{ii, 1}(1:(params.minTime.NREM*samplingRate)), 'constant');
    end
    
    LH_nrem_CBV = zeros(length(LH_nremCBV{1, 1}), length(LH_nremCBV));
    RH_nrem_CBV = zeros(length(RH_nremCBV{1, 1}), length(RH_nremCBV));
    LH_nrem_GAM = zeros(length(LH_nremGAM{1, 1}), length(LH_nremGAM));
    RH_nrem_GAM = zeros(length(RH_nremGAM{1, 1}), length(RH_nremGAM));
    
    for n = 1:length(LH_nremCBV)
        LH_nrem_CBV(:, n) = LH_nremCBV{n, 1};
        RH_nrem_CBV(:, n) = RH_nremCBV{n, 1};
        LH_nrem_GAM(:, n) = LH_nremGAM{n, 1};
        RH_nrem_GAM(:, n) = RH_nremGAM{n, 1};
    end
    
    disp('Analyzing the Coherence between the NREM SleepData L/R CBV signals...'); disp(' ')
    [C_nremCBV, phi_nremCBV, S12_nremCBV, S1_nremCBV, S2_nremCBV, f_nremCBV, confC_nremCBV, phiSTD_nremCBV, cErr_nremCBV] = coherencyc(LH_nrem_CBV, RH_nrem_CBV, params);
    disp('Analyzing the Coherence between the NREM SleepData L/R gamma band power signals...'); disp(' ')
    [C_nremGAM, phi_nremGAM, S12_nremGAM, S1_nremGAM, S2_nremGAM, f_nremGAM, confC_nremGAM, phiSTD_nremGAM, cErr_nremGAM] = coherencyc(LH_nrem_GAM, RH_nrem_GAM, params);
    
    nremCBV_CI = bootci(nboot, @coherencyc2, LH_nrem_CBV', RH_nrem_CBV', params);
    nremGAM_CI = bootci(nboot, @coherencyc2, LH_nrem_GAM', RH_nrem_GAM', params);
    
    nremCBVCoherence = figure;
    plot(f_nremCBV, C_nremCBV)
    hold on;
    plot(f_nremCBV, cErr_nremCBV)
    hold on;
    plot(f_nremCBV, nremCBV_CI)
    xlabel('Frequency (Hz)');
    ylabel('Magnitude of Coherence');
    title([animal  ' L/R CBV Coherence for NREM Data']);
    set(gca, 'Ticklength', [0 0]);
    legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
    set(legend, 'FontSize', 6);
    ylim([0 1])
    xlim([0 1])
    axis square
    
    nremGAMCoherence = figure;
    plot(f_nremGAM, C_nremGAM)
    hold on;
    plot(f_nremGAM, cErr_nremGAM)
    hold on;
    plot(f_nremCBV, nremGAM_CI)
    xlabel('Frequency (Hz)');
    ylabel('Magnitude of Coherence');
    title([animal  ' L/R Gamma Coherence for NREM Data']);
    set(gca, 'Ticklength', [0 0]);
    legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
    set(legend, 'FontSize', 6);
    ylim([0 1])
    xlim([0 1])
    axis square
    
    ComparisonData.Coherence.NREM.CBV.C = C_nremCBV;
    ComparisonData.Coherence.NREM.CBV.f = f_nremCBV;
    ComparisonData.Coherence.NREM.CBV.cErr = cErr_nremCBV;
    ComparisonData.Coherence.NREM.GAM.C = C_nremGAM;
    ComparisonData.Coherence.NREM.GAM.f = f_nremGAM;
    ComparisonData.Coherence.NREM.GAM.cErr = cErr_nremGAM;
    
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Coherence/'];
    
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    
    savefig(nremCBVCoherence, [dirpath animal '_NREM_' 'CBVCoherece']);
    savefig(nremGAMCoherence, [dirpath animal '_NREM_' 'GAMCoherece']);

    %% REM Sleep
    if ~isempty(SleepData.REM)
        LH_remCBV = SleepData.REM.Data.CBV.LH;
        RH_remCBV = SleepData.REM.Data.CBV.RH;
        LH_remGAM = SleepData.REM.Data.GammaBand_Power.LH;
        RH_remGAM = SleepData.REM.Data.GammaBand_Power.RH;
        
        for ii = 1:length(LH_remCBV)
            LH_remCBV{ii, 1} = detrend(LH_remCBV{ii, 1}(1:(params.minTime.REM*samplingRate)), 'constant');
            RH_remCBV{ii, 1} = detrend(RH_remCBV{ii, 1}(1:(params.minTime.REM*samplingRate)), 'constant');
            LH_remGAM{ii, 1} = detrend(LH_remGAM{ii, 1}(1:(params.minTime.REM*samplingRate)), 'constant');
            RH_remGAM{ii, 1} = detrend(RH_remGAM{ii, 1}(1:(params.minTime.REM*samplingRate)), 'constant');
        end
        
        LH_rem_CBV = zeros(length(LH_remCBV{1, 1}), length(LH_remCBV));
        RH_rem_CBV = zeros(length(RH_remCBV{1, 1}), length(RH_remCBV));
        LH_rem_GAM = zeros(length(LH_remGAM{1, 1}), length(LH_remGAM));
        RH_rem_GAM = zeros(length(RH_remGAM{1, 1}), length(RH_remGAM));
        
        for n = 1:length(LH_remCBV)
            LH_rem_CBV(:, n) = LH_remCBV{n, 1};
            RH_rem_CBV(:, n) = RH_remCBV{n, 1};
            LH_rem_GAM(:, n) = LH_remGAM{n, 1};
            RH_rem_GAM(:, n) = RH_remGAM{n, 1};
        end
        
        
        disp('Analyzing the Coherence between the REM SleepData L/R CBV signals...'); disp(' ')
        [C_remCBV, phi_remCBV, S12_remCBV, S1_remCBV, S2_remCBV, f_remCBV, confC_remCBV, phiSTD_remCBV, cErr_remCBV] = coherencyc(LH_rem_CBV, RH_rem_CBV, params);
        disp('Analyzing the Coherence between the REM SleepData L/R gamma band power signals...'); disp(' ')
        [C_remGAM, phi_remGAM, S12_remGAM, S1_remGAM, S2_remGAM, f_remGAM, confC_remGAM, phiSTD_remGAM, cErr_remGAM] = coherencyc(LH_rem_GAM, RH_rem_GAM, params);
        
        remCBV_CI = bootci(nboot, @coherencyc2, LH_rem_CBV', RH_rem_CBV', params);
        remGAM_CI = bootci(nboot, @coherencyc2, LH_rem_GAM', RH_rem_GAM', params);
        
        remCBVCoherence = figure;
        plot(f_remCBV, C_remCBV)
        hold on;
        plot(f_remCBV, cErr_remCBV)
        hold on;
        plot(f_remCBV, remCBV_CI)
        xlabel('Frequency (Hz)');
        ylabel('Magnitude of Coherence');
        title([animal  ' L/R CBV Coherence for REM Data']);
        set(gca, 'Ticklength', [0 0]);
        legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
        set(legend, 'FontSize', 6);
        ylim([0 1])
        xlim([0 1])
        axis square
        
        remGAMCoherence = figure;
        plot(f_remGAM, C_remGAM)
        hold on;
        plot(f_remGAM, cErr_remGAM)
        hold on;
        plot(f_remGAM, remGAM_CI)
        xlabel('Frequency (Hz)');
        ylabel('Magnitude of Coherence');
        title([animal  ' L/R Gamma Coherence for REM Data']);
        set(gca, 'Ticklength', [0 0]);
        legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
        set(legend, 'FontSize', 6);
        ylim([0 1])
        xlim([0 1])
        axis square
        
        ComparisonData.Coherence.REM.CBV.C = C_remCBV;
        ComparisonData.Coherence.REM.CBV.f = f_remCBV;
        ComparisonData.Coherence.REM.CBV.cErr = cErr_remCBV;
        ComparisonData.Coherence.REM.GAM.C = C_remGAM;
        ComparisonData.Coherence.REM.GAM.f = f_remGAM;
        ComparisonData.Coherence.REM.GAM.cErr = cErr_remGAM;
        
        [pathstr, ~, ~] = fileparts(cd);
        dirpath = [pathstr '/Figures/Coherence/'];
        
        if ~exist(dirpath, 'dir')
            mkdir(dirpath);
        end
        
        savefig(remCBVCoherence, [dirpath animal '_REM_' 'CBVCoherece']);
        savefig(remGAMCoherence, [dirpath animal '_REM_' 'GAMCoherece']);
        
    end
end

%% Load in relevant data from all ProcDataFiles
for pDF = 1:size(procDataFiles, 1)
    procDataFile = procDataFiles(pDF, :);
    load(procDataFile);
    [~, ~, fileDate, fileID] = GetFileInfo_IOS(procDataFile);
    strDay = ConvertDate(fileDate);
    
    disp(['Adding data from ProcData file ' num2str(pDF) ' of ' num2str(size(procDataFiles, 1))]); disp(' ')
    LH_CBV{pDF, 1} = (ProcData.Data.CBV.LH - RestingBaselines.CBV.LH.(strDay)) / RestingBaselines.CBV.LH.(strDay);
    RH_CBV{pDF, 1} = (ProcData.Data.CBV.RH - RestingBaselines.CBV.RH.(strDay)) / RestingBaselines.CBV.RH.(strDay);
    LH_GAM{pDF, 1} = (ProcData.Data.GammaBand_Power.LH - RestingBaselines.GammaBand_Power.LH.(strDay)) / RestingBaselines.GammaBand_Power.LH.(strDay);
    RH_GAM{pDF, 1} = (ProcData.Data.GammaBand_Power.RH - RestingBaselines.GammaBand_Power.RH.(strDay)) / RestingBaselines.GammaBand_Power.RH.(strDay);
end

for ii = 1:length(LH_CBV)
    LH_CBV{ii, 1} = detrend(filtfilt(B, A, LH_CBV{ii, 1}), 'constant');
    RH_CBV{ii, 1} = detrend(filtfilt(B, A, RH_CBV{ii, 1}), 'constant');
    LH_GAM{ii, 1} = detrend(filtfilt(D, C, LH_GAM{ii, 1}), 'constant');
    RH_GAM{ii, 1} = detrend(filtfilt(D, C, RH_GAM{ii, 1}), 'constant');
end

% Input data as time(1st dimension, vertical) by trials (2nd dimension, horizontally) for coherencyc and calc. coherence using chronux toolbox.
% Pre-allocate sizes
LH_allCBV = zeros(length(LH_CBV{1, 1}), length(LH_CBV));
RH_allCBV = zeros(length(RH_CBV{1, 1}), length(RH_CBV));
LH_allGAM = zeros(length(LH_GAM{1, 1}), length(LH_GAM));
RH_allGAM = zeros(length(RH_GAM{1, 1}), length(RH_GAM));

for n = 1:length(LH_CBV)
    LH_allCBV(:, n) = LH_CBV{n, 1};
    RH_allCBV(:, n) = RH_CBV{n, 1};
    LH_allGAM(:, n) = LH_GAM{n, 1};
    RH_allGAM(:, n) = RH_GAM{n, 1};
end

% Calculate the Coherency between desired signals and save those signals in a the comparison data structure.
% disp('Analyzing the Coherence between all L/R CBV signals...'); disp(' ')
% [C_allCBV, phi_allCBV, S12_allCBV, S1_allCBV, S2_allCBV, f_allCBV, confC_allCBV, phiSTD_allCBV, cErr_allCBV] = coherencyc(LH_allCBV, RH_allCBV, params);
% disp('Analyzing the Coherence between all L/R gamma band power signals...'); disp(' ')
% [C_allGAM, phi_allGAM, S12_allGAM, S1_allGAM, S2_allGAM, f_allGAM, confC_allGAM, phiSTD_allGAM, cErr_allGAM] = coherencyc(LH_allGAM, RH_allGAM, params);
% 
% allCBV_CI = bootci(nboot, @coherencyc2, LH_allCBV', RH_allCBV', params);
% allGAM_CI = bootci(nboot, @coherencyc2, LH_allGAM', RH_allGAM', params);
% 
% allCBVCoherence = figure;
% plot(f_allCBV, C_allCBV)
% hold on;
% plot(f_allCBV, cErr_allCBV)
% hold on;
% plot(f_allCBV, allCBV_CI)
% xlabel('Frequency (Hz)');
% ylabel('Magnitude of Coherence');
% title([animal  ' L/R CBV Coherence for All Data']);
% set(gca, 'Ticklength', [0 0]);
% legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
% set(legend, 'FontSize', 6);
% ylim([0 1])
% xlim([0 1])
% axis square
% 
% allGAMCoherence = figure;
% plot(f_allGAM, C_allGAM)
% hold on;
% plot(f_allGAM, cErr_allGAM)
% hold on;
% plot(f_allGAM, allGAM_CI)
% xlabel('Frequency (Hz)');
% ylabel('Magnitude of Coherence');
% title([animal  ' L/R Gamma Coherence for All Data']);
% set(gca, 'Ticklength', [0 0]);
% legend('Coherence', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
% set(legend, 'FontSize', 6);
% ylim([0 1])
% xlim([0 1])
% axis square
% 
% ComparisonData.Coherence.AllData.CBV.C = C_allCBV;
% ComparisonData.Coherence.AllData.CBV.f = f_allCBV;
% ComparisonData.Coherence.AllData.CBV.cErr = cErr_allCBV;
% ComparisonData.Coherence.AllData.GAM.C = C_allGAM;
% ComparisonData.Coherence.AllData.GAM.f = f_allGAM;
% ComparisonData.Coherence.AllData.GAM.cErr = cErr_allGAM;

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Coherence/'];

if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end

savefig(allCBVCoherence, [dirpath animal '_AllData_' 'CBVCoherece']);
savefig(allGAMCoherence, [dirpath animal '_AllData_' 'GAMCoherece']);


save([animal '_ComparisonData.mat'], 'ComparisonData');

end