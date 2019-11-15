function [] = DataAnalysis_IOS()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Step through various different data analysis and save the results in a summary structure.
%________________________________________________________________________________________________________________________

%% BLOCK PURPOSE: [0] Load the script's necessary variables and data structures.
% Clear the workspace variables and command window.
disp('Analyzing Block [0] Preparing the workspace and loading variables.'); disp(' ')

% Load or create the AnalysisResults.mat structure into the Workspace
resultsDataFileStruct = dir('*_AnalysisResults.mat');
resultsDataFile = {resultsDataFileStruct.name}';
resultsDataFileID = char(resultsDataFile);
if exist(resultsDataFileID)
    load(resultsDataFileID)
else
    AnalysisResults = [];
end

% %% BLOCK PURPOSE: [1] Stimulus and whisking evoked averages
% disp('Analyzing Block [1] Analyzing the whisking-evoked and stimulus-evoked hemodynamic and neural responses.'); disp(' ')
% evoked_dataTypes = {'adjLH','adjRH'};
% [AnalysisResults] = AnalyzeEvokedResponses_IOS(evoked_dataTypes,AnalysisResults);
% 
% %% BLOCK PURPOSE: [2] Cross correlation
% disp('Analyzing Block [2] Analzying the cross-correlation between hemodynamics and neural data.'); disp(' ')
% xcorr_dataTypes = {'adjLH','adjRH'};
% params.minTime.Rest = 10;   % seconds
% params.minTime.NREM = 30;   % seconds
% params.minTime.REM = 30;   % seconds
% [AnalysisResults] = AnalyzeXCorr_IOS(xcorr_dataTypes,params,AnalysisResults);
%  
% %% BLOCK PURPOSE: [3] Coherence
% disp('Analyzing Block [3] Analyzing the coherence between L/R hemodynamic and neural data.'); disp(' ')
% coherr_dataTypes = {'CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
% params.minTime.Rest = 10;   % seconds
% params.minTime.NREM = 30;   % seconds
% params.minTime.REM = 30;   % seconds
% [AnalysisResults] = AnalyzeCoherence_IOS(coherr_dataTypes,params, AnalysisResults);
% 
% %% BLOCK PURPOSE: [4] Power Spectra
% disp('Analyzing Block [4] Analyzing the power spectra of hemodynamic and neural data.'); disp(' ')
% powerspec_dataTypes =  {'CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
% params.minTime.Rest = 10;   % seconds
% params.minTime.NREM = 30;   % seconds
% params.minTime.REM = 30;   % seconds
% [AnalysisResults] = AnalyzePowerSpectrum_IOS(powerspec_dataTypes,params,AnalysisResults);
% 
% %% BLOCK PURPOSE: [5] Pearson's correlation coefficient
% disp('Analyzing Block [5] Analyzing the Pearson''s correlation coefficient between bilateral hemodynamic and neural data.'); disp(' ')
% corrCoeff_dataTypes =  {'CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
% params.minTime.Rest = 10;   % seconds
% [AnalysisResults] = AnalyzeCorrCoeffs_IOS(corrCoeff_dataTypes,params,AnalysisResults);

% %% BLOCK PURPOSE: [6] Mean CBV values
% disp('Analyzing Block [6] Analyzing the mean CBV during different behaviors.'); disp(' ')
% params.minTime.Rest = 10;   % seconds
% [AnalysisResults] = AnalyzeMeanCBV_IOS(params,AnalysisResults);
% 
% %% BLOCK PURPOSE: [7] Mean heart rate values
% disp('Analyzing Block [7] Analyzing the mean heart rate during different behaviors.'); disp(' ')
% params.minTime.Rest = 10;   % seconds
% [AnalysisResults] = AnalyzeMeanHeartRate_IOS(params,AnalysisResults);

%% BLOCK PURPOSE: [8] Hemodynamic response functions
% disp('Analyzing Block [8] Analyzing the hemodynamic response function and predictability of awake data.'); disp(' ')
% hemDataTypes = {'adjLH','adjRH'};
% neuralBands =  {'gammaBandPower','muaPower'};
% params.minTime.Rest = 10;   % seconds
% for d = 1:length(hemDataTypes)
%     hemDataType = hemDataTypes{1,d};
%     for e = 1:length(neuralBands)
%         neuralBand = neuralBands{1,e};
%         [AnalysisResults] = AnalyzeAwakeHRF_IOS(hemDataType,neuralBand,AnalysisResults);
% %         [AnalysisResults] = PredictHemodynamicChanges_IOS(params,hemDataType,neuralBand,AnalysisResults);
%     end
% end

%%
%% Figure 2d - Gamma-band HRF
clearvars -except Stats
clc
% animals = {'T108'};
CBVTypes = {'adjLH','adjRH'};
Behaviors = {'Contra','VW','Rest'};
% Behaviors = {'VW'};
HRFLims = [0 5];
% GamHRF = figure;
% set(gcf,'name','Figure 2d','numbertitle','off')
% ColorOrd = ['b','y','c'];
% display('Generating figure 2d...')
display('Calculating HRFs based on Gamma-band power...')
for q = 1:length(CBVTypes)
    CBVType = CBVTypes{1,q};
    for B = 1:length(Behaviors)
        Beh = Behaviors{B};
        display([num2str(B) ' of ' num2str(length(Behaviors)) ' behaviors...'])
        %         for a = 1:length(animals)
        %             display(sprintf(['\t' num2str(a) ' of ' num2str(length(animals)) ' animals...']));
        %         prevdir = cd([animals{a} filesep]);
        CalculateHRFDeconvolution_IOS('gammaBandPower',CBVType,Beh);
        %             TimeLims = HRFs.timevec>=HRFLims(1) & HRFs.timevec<=HRFLims(2);
        %             if a==1
        %                 AllHRFs.Gamma.(Beh).HRFs = NaN*ones(length(animals),sum(TimeLims));
        %                 AllHRFs.Gamma.(Beh).GammaHRFs = NaN*ones(length(animals),sum(TimeLims));
        %             end
        %             AllHRFs.Gamma.(Beh).HRFs(a,:) = HRFs.HRF(TimeLims);
        %             AllHRFs.Gamma.(Beh).GammaHRFs(a,:) = HRFs.GammaHRF;
        %             %         cd(prevdir)
        % %         end
        %         AllHRFs.Gamma.(Beh).Timevec = HRFs.timevec(TimeLims);
        %         AllHRFs.Gamma.(Beh).samplingRate = HRFs.samplingRate;
        %     mHRFs = mean(AllHRFs.Gamma.(Beh).HRFs);
        %         mHRFs = AllHRFs.Gamma.(Beh).HRFs;
        %     figure(GamHRF);
        %     plot(HRFs.timevec(TimeLims),mHRFs,ColorOrd(B));
        %     hold on;
    end
end
% figure;
% plot(AllHRFs.Gamma.Rest.GammaHRFs)
% ylabel('HRF Amplitude (A.U.)');
% xlabel('HRF Time (s)');
% legend({'Sensory Evoked','Whisking','Rest'},'location','southeast')
% title('Gamma-based HRFs')
% save('HRFs.mat','AllHRFs'); % save the result
% pause(0.001);

disp('Data Analysis - Complete.'); disp(' ')

