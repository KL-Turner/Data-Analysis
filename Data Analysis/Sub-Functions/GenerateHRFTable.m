function [ComparisonData] = GenerateHRFTable(procDataFiles, RestingBaselines, SleepData, ComparisonData)
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

LH_HRF = ComparisonData.HRFs.LH.dataStruct_out.the_best_kernel;
RH_HRF = ComparisonData.HRFs.RH.dataStruct_out.the_best_kernel;
HRF_t = ComparisonData.HRFs.LH.dataStruct_out.kernel_t;
LH_r = NaN(size(procDataFiles, 1), 1);
RH_r = NaN(size(procDataFiles, 1), 1);
LH_r2 = NaN(size(procDataFiles, 1), 1);
RH_r2 = NaN(size(procDataFiles, 1), 1);

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Single Trial HRF Figures/'];

if ~exist(dirpath, 'dir') 
    mkdir(dirpath); 
end

for f = 1:size(procDataFiles)
    disp(['Analyzing file ' num2str(f) ' of ' num2str(size(procDataFiles, 1)) '...']); disp(' ')
    % Handle file number, load each successive file
    procDataFile = procDataFiles(f, :);
    load(procDataFile);
    [animal, ~, fileDate, fileID] = GetFileInfo_IOS(procDataFile);
    fileIDs{f, 1} = fileID;
    strDay = ConvertDate(fileDate);
    
    %% CBV data - normalize by rest, lowpass filer, detrend
    LH_CBV = ProcData.Data.CBV.LH;
    RH_CBV = ProcData.Data.CBV.RH;
    LH_normCBV = (LH_CBV - RestingBaselines.CBV.LH.(strDay)) ./ (RestingBaselines.CBV.LH.(strDay));
    RH_normCBV = (RH_CBV - RestingBaselines.CBV.RH.(strDay)) ./ (RestingBaselines.CBV.RH.(strDay));
    [B, A] = butter(4, 2 / (30 / 2), 'low');   % 2 Hz lowpass
    LH_filtCBV = detrend(filtfilt(B, A, LH_normCBV), 'constant');
    RH_filtCBV = detrend(filtfilt(B, A, RH_normCBV), 'constant');
    
    %% Gamma band power data - normalize by rest, lowpass filer, detrend
    LH_gammaPower = ProcData.Data.GammaBand_Power.LH;
    RH_gammaPower = ProcData.Data.GammaBand_Power.RH;
    LH_normGamma = (LH_gammaPower - RestingBaselines.GammaBand_Power.LH.(strDay)) ./ (RestingBaselines.GammaBand_Power.LH.(strDay));
    RH_normGamma = (RH_gammaPower - RestingBaselines.GammaBand_Power.RH.(strDay)) ./ (RestingBaselines.GammaBand_Power.RH.(strDay));
    LH_filtGamma = detrend(filtfilt(B, A, LH_normGamma), 'constant');
    RH_filtGamma = detrend(filtfilt(B, A, RH_normGamma), 'constant');
    
    %% R % R2 values based on predictions
    LH_CBV_prediction = conv(LH_filtGamma, LH_HRF, 'full');
    LH_CBV_prediction = LH_CBV_prediction(1:length(LH_filtCBV));
    RH_CBV_prediction = conv(RH_filtGamma, RH_HRF, 'full');
    RH_CBV_prediction = RH_CBV_prediction(1:length(RH_filtCBV));
    
    LH_r_tbl = corrcoef(LH_filtCBV, LH_CBV_prediction);
    RH_r_tbl = corrcoef(RH_filtCBV, RH_CBV_prediction);
    LH_r(f, 1) = LH_r_tbl(2, 1);
    RH_r(f, 1) = RH_r_tbl(2, 1);
    
    % Error Variance
    LH_SSE = sum((LH_filtCBV - LH_CBV_prediction).^2);
    RH_SSE = sum((RH_filtCBV - RH_CBV_prediction).^2);
    
    % Total Variance
    LH_SST = sum((LH_filtCBV - (ones(size(LH_filtCBV, 1), 1)*mean(LH_filtCBV))).^2);
    RH_SST = sum((RH_filtCBV - (ones(size(RH_filtCBV, 1), 1)*mean(RH_filtCBV))).^2);
    
    % Check that the sum of the residuals is small compared to SSE + SSR
    LH_r2(f, 1) = ones(size(LH_SSE)) - LH_SSE ./ LH_SST;
    RH_r2(f, 1) = ones(size(RH_SSE)) - RH_SSE ./ RH_SST;
    
    %% Figure
    fig = figure;
    ax1 = subplot(2,1,1);
    plot((1:length(LH_filtCBV)) / 30, LH_filtCBV*100, 'k')
    hold on;
    plot((1:length(LH_CBV_prediction)) / 30, LH_CBV_prediction*100, 'r');
    title({[animal ' ' fileID ' HRF predictions'], ['r = ' num2str(LH_r(f,1)) ', r^2 = ' num2str(LH_r2(f,1))]})
    xlabel('Time (sec)')
    ylabel('Reflectance (%)')
    
%     yyaxis right
%     plot((1:length(LH_filtGamma)) / 30, LH_filtGamma);
%     ylabel('Normalized Gamma Power')
%     legend('LH CBV', 'LH Prediction', 'Gamma Power', 'Location', 'Northeast')
    
    ax2 = subplot(2,1,2);
    plot((1:length(RH_filtCBV)) / 30, RH_filtCBV*100, 'b')
    hold on;
    plot((1:length(RH_CBV_prediction)) / 30, RH_CBV_prediction*100, 'm');
    title(['r = ' num2str(RH_r(f,1)) ', r^2 = ' num2str(RH_r2(f,1))])
    xlabel('Time (sec)')
    ylabel('Reflectance (%)')
    
%     yyaxis right
%     plot((1:length(RH_filtGamma)) / 30, RH_filtGamma);
%     ylabel('Normalized Gamma Power')
%     legend('RH CBV', 'RH Prediction', 'Gamma Power', 'Location', 'Northeast')
%     linkaxes([ax1 ax2], 'xy');
    
    savefig(fig, [dirpath animal '_' fileID '_SingleTrialHRFFig']);
    close all
end

sleepLogical = zeros(size(procDataFiles, 1), 1);
sleepFiles = unique(SleepData.NREM.FileIDs);
for f = 1:size(procDataFiles, 1)
    procDataFile = procDataFiles(f, :);
    [~, ~, ~, procDataFileID] = GetFileInfo_IOS(procDataFile);
    for ff = 1:size(sleepFiles, 1)
        sleepFileID = char(sleepFiles(ff, 1));
        if strcmp(sleepFileID, procDataFileID)
            sleepLogical(f, 1) = 1;
        end
    end
end

ComparisonData.HRFs.Table = table(fileIDs, sleepLogical, LH_r, LH_r2, RH_r, RH_r2);
save([animal '_ComparisonData.mat'], 'ComparisonData');

for a = 1:length(sleepLogical)
    if sleepLogical(a, 1) == 1
        awakeLogical(a, 1) = 0;
    else
        awakeLogical(a, 1) = 1;
    end
end

sleepLogical = logical(sleepLogical);
awakeLogical = logical(awakeLogical);

LH_r_awake = LH_r(awakeLogical, :);
RH_r_awake = RH_r(awakeLogical, :);
LH_r_sleep = LH_r(sleepLogical, :);
RH_r_sleep = RH_r(sleepLogical, :);

LH_r2_awake = LH_r2(awakeLogical, :);
RH_r2_awake = RH_r2(awakeLogical, :);
LH_r2_sleep = LH_r2(sleepLogical, :);
RH_r2_sleep = RH_r2(sleepLogical, :);

hist = figure;
plot(HRF_t, LH_HRF, 'k')
hold on
plot(HRF_t, RH_HRF, 'g')
title('L/R HRFs from first 30 minutes of data')
xlabel('Time, (sec)')
ylabel('A.U.')
legend('LH HRF', 'RH HRF')

edges = -1:0.05:1;
figure;
ax1 = subplot(2,2,1);
h1 = histogram(LH_r_awake, edges, 'Normalization', 'probability');
hold on
h2 = histogram(LH_r_sleep, edges, 'Normalization', 'probability');
title('Left hemisphere R distribution')
xlabel('R val')
ylabel('Probability')
legend('Awake Trial R', 'Sleep Trial R')

ax2 = subplot(2,2,2);
h3 = histogram(RH_r_awake, edges, 'Normalization', 'probability');
hold on
h4 = histogram(RH_r_sleep, edges, 'Normalization', 'probability');
title('Right hemisphere R distribution')
xlabel('R^2 val')
ylabel('Probability')
legend('Awake Trial R', 'Sleep Trial R')
linkaxes([ax1 ax2], 'xy')

ax3 = subplot(2,2,3);
h5 = histogram(LH_r2_awake, edges, 'Normalization', 'probability');
hold on
h6 = histogram(LH_r2_sleep, edges, 'Normalization', 'probability');
title('Left hemisphere R^2 distribution')
xlabel('R^2 val')
ylabel('Probability')
legend('Awake Trial R^2', 'Sleep Trial R^2')

ax4 = subplot(2,2,4);
h7 = histogram(RH_r2_awake, edges, 'Normalization', 'probability');
hold on
h8 = histogram(RH_r2_sleep, edges, 'Normalization', 'probability');
linkaxes([ax1 ax2], 'xy')
title('Right hemisphere R^2 distribution')
xlabel('R^2 val')
ylabel('Probability')
legend('Awake Trial R^2', 'Sleep Trial R^2')
linkaxes([ax3 ax4], 'xy')

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/HRF Summary/'];

if ~exist(dirpath, 'dir') 
    mkdir(dirpath); 
end

savefig(hist, [dirpath animal '_' fileID '_HRFsummary']);

end
