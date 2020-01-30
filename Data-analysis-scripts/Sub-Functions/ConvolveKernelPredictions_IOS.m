function [AnalysisResults] = ConvolveKernelPredictions_IOS(procDataFiles, RestingBaselines, SleepData, AnalysisResults)
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

LH_HRF = AnalysisResults.HRFs.LH.dataStruct_out.the_best_kernel;
RH_HRF = AnalysisResults.HRFs.RH.dataStruct_out.the_best_kernel;
HRF_t = AnalysisResults.HRFs.LH.dataStruct_out.kernel_t;
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

