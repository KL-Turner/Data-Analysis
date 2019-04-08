function [ComparisonData] = AnalyzePowerSpectrum(animal, params, procDataFiles, RestingBaselines, RestData, SleepData, ComparisonData)
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

% Parameters for mtspectrumc - information available in function
samplingRate = RestData.CBV.LH.CBVCamSamplingRate;
params.tapers = [1 1];   % Tapers [n, 2n - 1]
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

% Calculate the power spectrum of the desired signals and save those signals in a the comparison data structure.
disp('Analyzing the power spectrum of the LH RestData CBV signal...'); disp(' ')
[LH_restCBV_S, LH_restCBV_f, LH_restCBV_sErr] = mtspectrumc(LH_restCBV, params);
disp('Analyzing the power spectrum of the RH RestData CBV signal...'); disp(' ')
[RH_restCBV_S, RH_restCBV_f, RH_restCBV_sErr] = mtspectrumc(RH_restCBV, params);

disp('Analyzing the power spectrum of the LH RestData Gamma Band signal...'); disp(' ')
[LH_restGAM_S, LH_restGAM_f, LH_restGAM_sErr] = mtspectrumc(LH_restGAM, params);
disp('Analyzing the power spectrum of the RH RestData Gamma Band signal...'); disp(' ')
[RH_restGAM_S, RH_restGAM_f, RH_restGAM_sErr] = mtspectrumc(RH_restGAM, params);

nboot = 1000;
restLHCBV_CI = bootci(nboot, @mtspectrumc2, LH_restCBV', params);
restRHCBV_CI = bootci(nboot, @mtspectrumc2, RH_restCBV', params);
restLHGAM_CI = bootci(nboot, @mtspectrumc2, LH_restGAM', params);
restRHGAM_CI = bootci(nboot, @mtspectrumc2, RH_restGAM', params);

% Rest CBV fig
RestCBVPowerSpec = figure;
ax1 = subplot(1,2,1);
loglog(LH_restCBV_f, LH_restCBV_S)
hold on; 
loglog(LH_restCBV_f, LH_restCBV_sErr)
hold on;
loglog(LH_restCBV_f, restLHCBV_CI)
xlabel('Frequency (Hz)')
ylabel('Power')
title([animal  ' LH Resting CBV Power Spectrum'])
legend('Power Spectra', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
set(gca, 'Ticklength', [0 0])
xlim([0 1])
ylim([1e-7 1e-2])
axis square

ax2 = subplot(1,2,2);
loglog(RH_restCBV_f, RH_restCBV_S)
hold on; 
loglog(RH_restCBV_f, RH_restCBV_sErr)
hold on;
loglog(RH_restCBV_f, restRHCBV_CI)
xlabel('Frequency (Hz)')
ylabel('Power')
title([animal  ' RH Resting CBV Power Spectrum'])
legend('Power Spectra', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
set(gca, 'Ticklength', [0 0])
xlim([0 1])
axis square
linkaxes([ax1 ax2], 'y')

% Rest Gam fig
RestGammaPowerSpec = figure;
ax3 = subplot(1,2,1);
loglog(LH_restGAM_f, LH_restGAM_S)
hold on; 
loglog(LH_restGAM_f, LH_restGAM_sErr)
loglog(LH_restGAM_f, restLHGAM_CI)
xlabel('Frequency (Hz)')
ylabel('Power')
title([animal  ' LH Resting Gamma Power Spectrum'])
legend('Power Spectra', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
set(gca, 'Ticklength', [0 0])
xlim([0 1])
ylim([1e-4 1e0])
axis square

ax4 = subplot(1,2,2);
loglog(RH_restGAM_f, RH_restGAM_S)
hold on;
loglog(RH_restGAM_f, RH_restGAM_sErr)
loglog(RH_restGAM_f, restRHGAM_CI)
xlabel('Frequency (Hz)')
ylabel('Power')
title([animal  ' RH Resting Gamma Power Spectrum'])
legend('Power Spectra', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
set(gca, 'Ticklength', [0 0])
xlim([0 1])
axis square
linkaxes([ax3 ax4], 'y')

ComparisonData.PowerSpectrum.Rest.CBV.LH.S = LH_restCBV_S;
ComparisonData.PowerSpectrum.Rest.CBV.LH.f = LH_restCBV_f;
ComparisonData.PowerSpectrum.Rest.CBV.LH.sErr = LH_restCBV_sErr;
ComparisonData.PowerSpectrum.Rest.CBV.RH.S = RH_restCBV_S;
ComparisonData.PowerSpectrum.Rest.CBV.RH.f = RH_restCBV_f;
ComparisonData.PowerSpectrum.Rest.CBV.RH.sErr = RH_restCBV_sErr;
ComparisonData.PowerSpectrum.Rest.GAM.LH.S = LH_restGAM_S;
ComparisonData.PowerSpectrum.Rest.GAM.LH.f = LH_restGAM_f;
ComparisonData.PowerSpectrum.Rest.GAM.LH.sErr = LH_restGAM_sErr;
ComparisonData.PowerSpectrum.Rest.GAM.RH.S = RH_restGAM_S;
ComparisonData.PowerSpectrum.Rest.GAM.RH.f = RH_restGAM_f;
ComparisonData.PowerSpectrum.Rest.GAM.RH.sErr = RH_restGAM_sErr;
save([animal '_ComparisonData.mat'], 'ComparisonData');

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Power Spectrum/'];

if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end

savefig(RestCBVPowerSpec, [dirpath animal '_Rest_CBVPowerSpec']);
savefig(RestGammaPowerSpec, [dirpath animal '_Rest_GammaPowerSpec']);

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
    
    % Calculate the power spectrum of the desired signals and save those signals in a the comparison data structure.
    disp('Analyzing the power spectrum of the LH NREM SleepData CBV signal...'); disp(' ')
    [LH_nremCBV_S, LH_nremCBV_f, LH_nremCBV_sErr] = mtspectrumc(LH_nrem_CBV, params);
    disp('Analyzing the power spectrum of the RH NREM SleepData CBV signal...'); disp(' ')
    [RH_nremCBV_S, RH_nremCBV_f, RH_nremCBV_sErr] = mtspectrumc(RH_nrem_CBV, params);
    
    disp('Analyzing the power spectrum of the LH NREM SleepData Gamma Band signal...'); disp(' ')
    [LH_nremGAM_S, LH_nremGAM_f, LH_nremGAM_sErr] = mtspectrumc(LH_nrem_GAM, params);
    disp('Analyzing the power spectrum of the RH NREM SleepData Gamma Band signal...'); disp(' ')
    [RH_nremGAM_S, RH_nremGAM_f, RH_nremGAM_sErr] = mtspectrumc(RH_nrem_GAM, params);
    
    nremLHCBV_CI = bootci(nboot, @mtspectrumc2, LH_nrem_CBV', params);
    nremRHCBV_CI = bootci(nboot, @mtspectrumc2, RH_nrem_CBV', params);
    nremLHGAM_CI = bootci(nboot, @mtspectrumc2, LH_nrem_GAM', params);
    nremRHGAM_CI = bootci(nboot, @mtspectrumc2, RH_nrem_GAM', params);
    
    % Rest CBV fig
    nremCBVPowerSpec = figure;
    ax5 = subplot(1,2,1);
    loglog(LH_nremCBV_f, LH_nremCBV_S)
    hold on;
    loglog(LH_nremCBV_f, LH_nremCBV_sErr)
    loglog(LH_nremCBV_f, nremLHCBV_CI)
    xlabel('Frequency (Hz)')
    ylabel('Power')
    title([animal  ' LH NREM CBV Power Spectrum'])
    legend('Power Spectra', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
    set(gca, 'Ticklength', [0 0])
    xlim([0 1])
    ylim([1e-7 1e-2])
    axis square
    
    ax6 = subplot(1,2,2);
    loglog(RH_nremCBV_f, RH_nremCBV_S)
    hold on; 
    loglog(RH_nremCBV_f, RH_nremCBV_sErr)
    loglog(RH_nremCBV_f, nremRHCBV_CI)
    xlabel('Frequency (Hz)')
    ylabel('Power')
    title([animal  ' RH NREM CBV Power Spectrum'])
    legend('Power Spectra', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
    set(gca, 'Ticklength', [0 0])
    xlim([0 1])
    axis square
    linkaxes([ax5 ax6], 'y')
    
    % Rest Gam fig
    nremGammaPowerSpec = figure;
    ax7 = subplot(1,2,1);
    loglog(LH_nremGAM_f, LH_nremGAM_S)
    hold on; 
    loglog(RH_nremGAM_f, RH_nremGAM_sErr)
    loglog(RH_nremGAM_f, nremRHGAM_CI)
    xlabel('Frequency (Hz)')
    ylabel('Power')
    title([animal  ' LH NREM Gamma Power Spectrum'])
    legend('Power Spectra', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
    set(gca, 'Ticklength', [0 0])
    xlim([0 1])
    ylim([1e-4 1e0])
    axis square
    
    ax8 = subplot(1,2,2);
    loglog(RH_nremGAM_f, RH_nremGAM_S)
    hold on; 
    loglog(RH_nremGAM_f, RH_nremGAM_sErr)
    loglog(RH_nremGAM_f, nremRHGAM_CI)
    xlabel('Frequency (Hz)')
    ylabel('Power')
    title([animal  ' RH NREM Gamma Power Spectrum'])
    legend('Power Spectra', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
    set(gca, 'Ticklength', [0 0])
    xlim([0 1])
    axis square
    linkaxes([ax7 ax8], 'y')
    
    ComparisonData.PowerSpectrum.NREM.CBV.LH.S = LH_nremCBV_S;
    ComparisonData.PowerSpectrum.NREM.CBV.LH.f = LH_nremCBV_f;
    ComparisonData.PowerSpectrum.NREM.CBV.LH.sErr = LH_nremCBV_sErr;
    ComparisonData.PowerSpectrum.NREM.CBV.RH.S = RH_nremCBV_S;
    ComparisonData.PowerSpectrum.NREM.CBV.RH.f = RH_nremCBV_f;
    ComparisonData.PowerSpectrum.NREM.CBV.RH.sErr = RH_nremCBV_sErr;
    ComparisonData.PowerSpectrum.NREM.GAM.LH.S = LH_nremGAM_S;
    ComparisonData.PowerSpectrum.NREM.GAM.LH.f = LH_nremGAM_f;
    ComparisonData.PowerSpectrum.NREM.GAM.LH.sErr = LH_nremGAM_sErr;
    ComparisonData.PowerSpectrum.NREM.GAM.RH.S = RH_nremGAM_S;
    ComparisonData.PowerSpectrum.NREM.GAM.RH.f = RH_nremGAM_f;
    ComparisonData.PowerSpectrum.NREM.GAM.RH.sErr = RH_nremGAM_sErr;
    save([animal '_ComparisonData.mat'], 'ComparisonData');
    
    [pathstr, ~, ~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Power Spectrum/'];
    
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    
    savefig(nremCBVPowerSpec, [dirpath animal '_NREM_CBVPowerSpec']);
    savefig(nremGammaPowerSpec, [dirpath animal '_NREM_GammaPowerSpec']);
    
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
        
        
        % Calculate the power spectrum of the desired signals and save those signals in a the comparison data structure.
        disp('Analyzing the power spectrum of the LH REM SleepData CBV signal...'); disp(' ')
        [LH_remCBV_S, LH_remCBV_f, LH_remCBV_sErr] = mtspectrumc(LH_rem_CBV, params);
        disp('Analyzing the power spectrum of the RH REM SleepData CBV signal...'); disp(' ')
        [RH_remCBV_S, RH_remCBV_f, RH_remCBV_sErr] = mtspectrumc(RH_rem_CBV, params);
        
        disp('Analyzing the power spectrum of the LH REM SleepData Gamma Band signal...'); disp(' ')
        [LH_remGAM_S, LH_remGAM_f, LH_remGAM_sErr] = mtspectrumc(LH_rem_GAM, params);
        disp('Analyzing the power spectrum of the RH REM SleepData Gamma Band signal...'); disp(' ')
        [RH_remGAM_S, RH_remGAM_f, RH_remGAM_sErr] = mtspectrumc(RH_rem_GAM, params);
        
        % Rest CBV fig
        remCBVPowerSpec = figure;
        ax9 = subplot(1,2,1);
        loglog(LH_remCBV_f, LH_remCBV_S)
        xlabel('Frequency (Hz)')
        ylabel('Power')
        title([animal  ' LH REM CBV Power Spectrum'])
        legend('Power Spectra', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
        set(gca, 'Ticklength', [0 0])
        xlim([0 1])
        axis square
        
        ax10 = subplot(1,2,2);
        loglog(RH_remCBV_f, RH_remCBV_S)
        xlabel('Frequency (Hz)')
        ylabel('Power')
        title([animal  ' RH REM CBV Power Spectrum'])
        legend('Power Spectra', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
        set(gca, 'Ticklength', [0 0])
        xlim([0 1])
        axis square
        linkaxes([ax9 ax10], 'y')
        
        % Rest Gam fig
        remGammaPowerSpec = figure;
        ax11 = subplot(1,2,1);
        loglog(LH_remGAM_f, LH_remGAM_S)
        xlabel('Frequency (Hz)')
        ylabel('Power')
        title([animal  ' LH NREM Gamma Power Spectrum'])
        set(gca, 'Ticklength', [0 0])
        xlim([0 1])
        axis square
        
        ax12 = subplot(1,2,2);
        loglog(RH_remGAM_f, RH_remGAM_S)
        xlabel('Frequency (Hz)')
        ylabel('Power')
        title([animal  ' RH REM Gamma Power Spectrum'])
        set(gca, 'Ticklength', [0 0])
        xlim([0 1])
        axis square
        linkaxes([ax11 ax12], 'y')
        
        ComparisonData.PowerSpectrum.REM.CBV.LH.S = LH_remCBV_S;
        ComparisonData.PowerSpectrum.REM.CBV.LH.f = LH_remCBV_f;
        ComparisonData.PowerSpectrum.REM.CBV.LH.sErr = LH_remCBV_sErr;
        ComparisonData.PowerSpectrum.REM.CBV.RH.S = RH_remCBV_S;
        ComparisonData.PowerSpectrum.REM.CBV.RH.f = RH_remCBV_f;
        ComparisonData.PowerSpectrum.REM.CBV.RH.sErr = RH_remCBV_sErr;
        ComparisonData.PowerSpectrum.REM.GAM.LH.S = LH_remGAM_S;
        ComparisonData.PowerSpectrum.REM.GAM.LH.f = LH_remGAM_f;
        ComparisonData.PowerSpectrum.REM.GAM.LH.sErr = LH_remGAM_sErr;
        ComparisonData.PowerSpectrum.REM.GAM.RH.S = RH_remGAM_S;
        ComparisonData.PowerSpectrum.REM.GAM.RH.f = RH_remGAM_f;
        ComparisonData.PowerSpectrum.REM.GAM.RH.sErr = RH_remGAM_sErr;
        save([animal '_ComparisonData.mat'], 'ComparisonData');
        
        [pathstr, ~, ~] = fileparts(cd);
        dirpath = [pathstr '/Figures/Power Spectrum/'];
        
        if ~exist(dirpath, 'dir')
            mkdir(dirpath);
        end
        
        savefig(remCBVPowerSpec, [dirpath animal '_REM_CBVPowerSpec']);
        savefig(remGammaPowerSpec, [dirpath animal '_REM_GammaPowerSpec']);
        
    end
end

% %% Load in relevant data from all ProcDataFiles
% for pDF = 1:size(procDataFiles, 1)
%     procDataFile = procDataFiles(pDF, :);
%     load(procDataFile);
%     [~, ~, fileDate, fileID] = GetFileInfo_IOS(procDataFile);
%     strDay = ConvertDate(fileDate);
%     
%     disp(['Adding data from ProcData file ' num2str(pDF) ' of ' num2str(size(procDataFiles, 1))]); disp(' ')
%     LH_CBV{pDF, 1} = (ProcData.Data.CBV.LH - RestingBaselines.CBV.LH.(strDay)) / RestingBaselines.CBV.LH.(strDay);
%     RH_CBV{pDF, 1} = (ProcData.Data.CBV.RH - RestingBaselines.CBV.RH.(strDay)) / RestingBaselines.CBV.RH.(strDay);
%     LH_GAM{pDF, 1} = (ProcData.Data.GammaBand_Power.LH - RestingBaselines.GammaBand_Power.LH.(strDay)) / RestingBaselines.GammaBand_Power.LH.(strDay);
%     RH_GAM{pDF, 1} = (ProcData.Data.GammaBand_Power.RH - RestingBaselines.GammaBand_Power.RH.(strDay)) / RestingBaselines.GammaBand_Power.RH.(strDay);
% end
% 
% for ii = 1:length(LH_CBV)
%     LH_CBV{ii, 1} = detrend(filtfilt(B, A, LH_CBV{ii, 1}), 'constant');
%     RH_CBV{ii, 1} = detrend(filtfilt(B, A, RH_CBV{ii, 1}), 'constant');
%     LH_GAM{ii, 1} = detrend(filtfilt(D, C, LH_GAM{ii, 1}), 'constant');
%     RH_GAM{ii, 1} = detrend(filtfilt(D, C, RH_GAM{ii, 1}), 'constant');
% end
% 
% % Input data as time(1st dimension, vertical) by trials (2nd dimension, horizontally) for coherencyc and calc. coherence using chronux toolbox.
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
% % Calculate the power spectrum of the desired signals and save those signals in a the comparison data structure.
% disp('Analyzing the power spectrum of the LH all data CBV signal...'); disp(' ')
% [LH_allCBV_S, LH_allCBV_f, LH_allCBV_sErr] = mtspectrumc(LH_allCBV, params);
% disp('Analyzing the power spectrum of the RH all data CBV signal...'); disp(' ')
% [RH_allCBV_S, RH_allCBV_f, RH_allCBV_sErr] = mtspectrumc(RH_allCBV, params);
% 
% disp('Analyzing the power spectrum of the LH all data Gamma Band signal...'); disp(' ')
% [LH_allGAM_S, LH_allGAM_f, LH_allGAM_sErr] = mtspectrumc(LH_allGAM, params);
% disp('Analyzing the power spectrum of the RH all data Gamma Band signal...'); disp(' ')
% [RH_allGAM_S, RH_allGAM_f, RH_allGAM_sErr] = mtspectrumc(RH_allGAM, params);
% 
% allLHCBV_CI = bootci(nboot, @mtspectrumc2, LH_allCBV', params);
% allRHCBV_CI = bootci(nboot, @mtspectrumc2, RH_allCBV', params);
% 
% % Rest CBV fig
% allCBVPowerSpec = figure;
% ax13 = subplot(1,2,1);
% loglog(LH_allCBV_f, LH_allCBV_S)
% hold on;
% loglog(LH_allCBV_f, LH_allCBV_sErr)
% loglog(LH_allCBV_f, allLHCBV_CI)
% xlabel('Frequency (Hz)')
% ylabel('Power')
% title([animal  ' LH all CBV Power Spectrum'])
% legend('Power Spectra', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
% set(gca, 'Ticklength', [0 0])
% xlim([0 1])
% axis square
% 
% ax14 = subplot(1,2,2);
% loglog(RH_allCBV_f, RH_allCBV_S)
% hold on;
% loglog(RH_allCBV_f, RH_allCBV_sErr)
% loglog(RH_allCBV_f, allRHCBV_CI)
% xlabel('Frequency (Hz)')
% ylabel('Power')
% title([animal  ' RH all CBV Power Spectrum'])
% legend('Power Spectra', 'Jackknife Lower', 'JackknifeUpper', 'Bootstrap Lower', 'Bootstrap Upper', 'Location', 'Southeast');
% set(gca, 'Ticklength', [0 0])
% xlim([0 1])
% axis square
% linkaxes([ax13 ax14], 'y')
% 
% % Rest Gam fig
% allGammaPowerSpec = figure;
% ax15 = subplot(1,2,1);
% loglog(LH_allGAM_f, LH_allGAM_S)
% xlabel('Frequency (Hz)')
% ylabel('Power')
% title([animal  ' LH all Gamma Power Spectrum'])
% set(gca, 'Ticklength', [0 0])
% xlim([0 1])
% axis square
% 
% ax16 = subplot(1,2,2);
% loglog(RH_allGAM_f, RH_allGAM_S)
% xlabel('Frequency (Hz)')
% ylabel('Power')
% title([animal  ' RH all Gamma Power Spectrum'])
% set(gca, 'Ticklength', [0 0])
% xlim([0 1])
% axis square
% linkaxes([ax15 ax16], 'y')
% 
% ComparisonData.PowerSpectrum.AllData.CBV.LH.S = LH_allCBV_S;
% ComparisonData.PowerSpectrum.AllData.CBV.LH.f = LH_allCBV_f;
% ComparisonData.PowerSpectrum.AllData.CBV.LH.sErr = LH_allCBV_sErr;
% ComparisonData.PowerSpectrum.AllData.CBV.RH.S = RH_allCBV_S;
% ComparisonData.PowerSpectrum.AllData.CBV.RH.f = RH_allCBV_f;
% ComparisonData.PowerSpectrum.AllData.CBV.RH.sErr = RH_allCBV_sErr;
% ComparisonData.PowerSpectrum.AllData.GAM.LH.S = LH_allGAM_S;
% ComparisonData.PowerSpectrum.AllData.GAM.LH.f = LH_allGAM_f;
% ComparisonData.PowerSpectrum.AllData.GAM.LH.sErr = LH_allGAM_sErr;
% ComparisonData.PowerSpectrum.AllData.GAM.RH.S = RH_allGAM_S;
% ComparisonData.PowerSpectrum.AllData.GAM.RH.f = RH_allGAM_f;
% ComparisonData.PowerSpectrum.AllData.GAM.RH.sErr = RH_allGAM_sErr;
% save([animal '_ComparisonData.mat'], 'ComparisonData');

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Power Spectrum/'];

if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end

savefig(allCBVPowerSpec, [dirpath animal '_AllData_CBVPowerSpec']);
savefig(allGammaPowerSpec, [dirpath animal '_AllData_GammaPowerSpec']);


save([animal '_ComparisonData.mat'], 'ComparisonData');

end



