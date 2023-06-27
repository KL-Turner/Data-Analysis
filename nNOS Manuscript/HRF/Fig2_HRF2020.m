function [AnalysisResults] = Fig2_HRF2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S2 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% arousal-state colors
colorRest = [(0/256),(166/256),(81/256)];
colorWhisk = [(31/256),(120/256),(179/256)];
colorStim = [(255/256),(28/256),(206/256)];
colorNREM = [(191/256),(0/256),(255/256)];
colorREM = [(254/256),(139/256),(0/256)];
colorAll = [(183/256),(115/256),(51/256)];
%% data types
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
solenoidNames = {'LPadSol','RPadSol','AudSol'};
compDataTypes = {'Ipsi','Contra','Auditory'};
whiskDataTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
dataTypes = {'cortLH','cortRH'};
behavFields = {'Rest','NREM','REM'};
neuralBands = {'gammaBandPower','muaPower'};
behaviors = {'Contra','Whisk','Rest','NREM','REM','All'};
samplingRate = 30;
%% set-up and process data
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for dd = 1:length(solenoidNames)
            solenoidName = solenoidNames{1,dd};
            data.Stim.(dataType).(solenoidName).count(:,aa) = AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoidName).count;
            data.Stim.(dataType).(solenoidName).HbT(:,aa) = AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoidName).CBV_HbT.HbT;
            data.Stim.(dataType).(solenoidName).CBV(:,aa) = AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoidName).CBV.CBV;
            data.Stim.(dataType).(solenoidName).cortMUA(:,aa) = AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoidName).MUA.corticalData;
            data.Stim.(dataType).(solenoidName).hipMUA(:,aa) = AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoidName).MUA.hippocampalData;
            data.Stim.(dataType).(solenoidName).cortGam(:,aa) = AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoidName).Gam.corticalData;
            data.Stim.(dataType).(solenoidName).hipGam(:,aa) = AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoidName).Gam.hippocampalData;
            data.Stim.(dataType).(solenoidName).timeVector(:,aa) = AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoidName).timeVector;
            data.Stim.(dataType).(solenoidName).cortS(:,:,aa) = AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoidName).LFP.corticalS;
            data.Stim.(dataType).(solenoidName).cortS_Gam(:,:,aa) = AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoidName).LFP.corticalS(49:end,20:23);
            data.Stim.(dataType).(solenoidName).hipS(:,:,aa) = AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoidName).LFP.hippocampalS;
            data.Stim.(dataType).(solenoidName).hipS_Gam(:,:,aa) = AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoidName).LFP.hippocampalS(49:end,20:23);
            data.Stim.(dataType).(solenoidName).T(:,aa) = AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoidName).LFP.T;
            data.Stim.(dataType).(solenoidName).F(:,aa) = AnalysisResults.Evoked.(animalID).Stim.(dataType).(solenoidName).LFP.F;
        end
    end
end
% concatenate the data from the contra and ipsi data
data.Stim.Contra.count = cat(2,data.Stim.cortLH.RPadSol.count,data.Stim.cortRH.LPadSol.count);
data.Stim.Contra.HbT = cat(2,data.Stim.cortLH.RPadSol.HbT,data.Stim.cortRH.LPadSol.HbT);
data.Stim.Contra.CBV = cat(2,data.Stim.cortLH.RPadSol.CBV,data.Stim.cortRH.LPadSol.CBV);
data.Stim.Contra.cortMUA = cat(2,data.Stim.cortLH.RPadSol.cortMUA,data.Stim.cortRH.LPadSol.cortMUA);
data.Stim.Contra.hipMUA = data.Stim.cortRH.RPadSol.hipMUA;
data.Stim.Contra.cortGam = cat(2,data.Stim.cortLH.RPadSol.cortGam,data.Stim.cortRH.LPadSol.cortGam);
data.Stim.Contra.hipGam = data.Stim.cortRH.RPadSol.hipGam;
data.Stim.Contra.timeVector = cat(2,data.Stim.cortLH.RPadSol.timeVector,data.Stim.cortRH.LPadSol.timeVector);
data.Stim.Contra.cortS = cat(3,data.Stim.cortLH.RPadSol.cortS,data.Stim.cortRH.LPadSol.cortS);
data.Stim.Contra.cortS_Gam = cat(3,data.Stim.cortLH.RPadSol.cortS_Gam,data.Stim.cortRH.LPadSol.cortS_Gam);
data.Stim.Contra.hipS = data.Stim.cortRH.RPadSol.hipS;
data.Stim.Contra.hipS_Gam = data.Stim.cortRH.RPadSol.hipS_Gam;
data.Stim.Contra.T = cat(2,data.Stim.cortLH.RPadSol.T,data.Stim.cortRH.LPadSol.T);
data.Stim.Contra.F = cat(2,data.Stim.cortLH.RPadSol.F,data.Stim.cortRH.LPadSol.F);
data.Stim.Ipsi.count = cat(2,data.Stim.cortLH.LPadSol.count,data.Stim.cortRH.RPadSol.count);
data.Stim.Ipsi.HbT = cat(2,data.Stim.cortLH.LPadSol.HbT,data.Stim.cortRH.RPadSol.HbT);
data.Stim.Ipsi.CBV = cat(2,data.Stim.cortLH.LPadSol.CBV,data.Stim.cortRH.RPadSol.CBV);
data.Stim.Ipsi.cortMUA = cat(2,data.Stim.cortLH.LPadSol.cortMUA,data.Stim.cortRH.RPadSol.cortMUA);
data.Stim.Ipsi.hipMUA = data.Stim.cortRH.LPadSol.hipMUA;
data.Stim.Ipsi.cortGam = cat(2,data.Stim.cortLH.LPadSol.cortGam,data.Stim.cortRH.RPadSol.cortGam);
data.Stim.Ipsi.hipGam = data.Stim.cortRH.LPadSol.hipGam;
data.Stim.Ipsi.timeVector = cat(2,data.Stim.cortLH.LPadSol.timeVector,data.Stim.cortRH.RPadSol.timeVector);
data.Stim.Ipsi.cortS = cat(3,data.Stim.cortLH.LPadSol.cortS,data.Stim.cortRH.RPadSol.cortS);
data.Stim.Ipsi.cortS_Gam = cat(3,data.Stim.cortLH.LPadSol.cortS_Gam,data.Stim.cortRH.RPadSol.cortS_Gam);
data.Stim.Ipsi.hipS = data.Stim.cortRH.LPadSol.hipS;
data.Stim.Ipsi.hipS_Gam = data.Stim.cortRH.LPadSol.hipS_Gam;
data.Stim.Ipsi.T = cat(2,data.Stim.cortLH.LPadSol.T,data.Stim.cortRH.RPadSol.T);
data.Stim.Ipsi.F = cat(2,data.Stim.cortLH.LPadSol.F,data.Stim.cortRH.RPadSol.F);
data.Stim.Auditory.count = cat(2,data.Stim.cortLH.AudSol.count,data.Stim.cortRH.AudSol.count);
data.Stim.Auditory.HbT = cat(2,data.Stim.cortLH.AudSol.HbT,data.Stim.cortRH.AudSol.HbT);
data.Stim.Auditory.CBV = cat(2,data.Stim.cortLH.AudSol.CBV,data.Stim.cortRH.AudSol.CBV);
data.Stim.Auditory.cortMUA = cat(2,data.Stim.cortLH.AudSol.cortMUA,data.Stim.cortRH.AudSol.cortMUA);
data.Stim.Auditory.hipMUA = data.Stim.cortRH.AudSol.hipMUA;
data.Stim.Auditory.cortGam = cat(2,data.Stim.cortLH.AudSol.cortGam,data.Stim.cortRH.AudSol.cortGam);
data.Stim.Auditory.hipGam = data.Stim.cortRH.AudSol.hipGam;
data.Stim.Auditory.timeVector = cat(2,data.Stim.cortLH.AudSol.timeVector,data.Stim.cortRH.AudSol.timeVector);
data.Stim.Auditory.cortS = cat(3,data.Stim.cortLH.AudSol.cortS,data.Stim.cortRH.AudSol.cortS);
data.Stim.Auditory.cortS_Gam = cat(3,data.Stim.cortLH.AudSol.cortS_Gam,data.Stim.cortRH.AudSol.cortS_Gam);
data.Stim.Auditory.hipS = data.Stim.cortRH.AudSol.hipS;
data.Stim.Auditory.hipS_Gam = data.Stim.cortRH.AudSol.hipS_Gam;
data.Stim.Auditory.T = cat(2,data.Stim.cortLH.AudSol.T,data.Stim.cortRH.AudSol.T);
data.Stim.Auditory.F = cat(2,data.Stim.cortLH.AudSol.F,data.Stim.cortRH.AudSol.F);
% take the averages of each field through the proper dimension
for ff = 1:length(compDataTypes)
    compDataType = compDataTypes{1,ff};
    data.Stim.(compDataType).mean_Count = mean(data.Stim.(compDataType).count,2);
    data.Stim.(compDataType).std_Count = std(data.Stim.(compDataType).count,0,2);
    data.Stim.(compDataType).mean_HbT = mean(data.Stim.(compDataType).HbT,2);
    data.Stim.(compDataType).std_HbT = std(data.Stim.(compDataType).HbT,0,2);
    data.Stim.(compDataType).mean_CBV = mean(data.Stim.(compDataType).CBV,2);
    data.Stim.(compDataType).std_CBV = std(data.Stim.(compDataType).CBV,0,2);
    data.Stim.(compDataType).mean_CortMUA = mean(data.Stim.(compDataType).cortMUA,2);
    data.Stim.(compDataType).std_CortMUA = std(data.Stim.(compDataType).cortMUA,0,2);
    data.Stim.(compDataType).mean_HipMUA = mean(data.Stim.(compDataType).hipMUA,2);
    data.Stim.(compDataType).std_HipMUA = std(data.Stim.(compDataType).hipMUA,0,2);
    data.Stim.(compDataType).mean_CortGam = mean(data.Stim.(compDataType).cortGam,2);
    data.Stim.(compDataType).std_CortGam = std(data.Stim.(compDataType).cortGam,0,2);
    data.Stim.(compDataType).mean_HipGam = mean(data.Stim.(compDataType).hipGam,2);
    data.Stim.(compDataType).std_HipGam = std(data.Stim.(compDataType).hipGam,0,2);
    data.Stim.(compDataType).mean_timeVector = mean(data.Stim.(compDataType).timeVector,2);
    data.Stim.(compDataType).mean_CortS = mean(data.Stim.(compDataType).cortS,3).*100;
    data.Stim.(compDataType).mean_CortS_Gam = mean(mean(mean(data.Stim.(compDataType).cortS_Gam.*100,2),1),3);
    data.Stim.(compDataType).std_CortS_Gam = std(mean(mean(data.Stim.(compDataType).cortS_Gam.*100,2),1),0,3);
    data.Stim.(compDataType).mean_HipS = mean(data.Stim.(compDataType).hipS,3).*100;
    data.Stim.(compDataType).mean_HipS_Gam = mean(mean(mean(data.Stim.(compDataType).hipS_Gam.*100,2),1),3);
    data.Stim.(compDataType).std_HipS_Gam = std(mean(mean(data.Stim.(compDataType).hipS_Gam.*100,2),1),0,3);
    data.Stim.(compDataType).mean_T = mean(data.Stim.(compDataType).T,2);
    data.Stim.(compDataType).mean_F = mean(data.Stim.(compDataType).F,2);
end
%% set-up and process data
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for cc = 1:length(whiskDataTypes)
        whiskDataType = whiskDataTypes{1,cc};
        % left cortical
        data.Whisk.(whiskDataType).cortLH.HbT(:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortLH.(whiskDataType).CBV_HbT.HbT;
        data.Whisk.(whiskDataType).cortLH.CBV(:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortLH.(whiskDataType).CBV.CBV;
        data.Whisk.(whiskDataType).cortLH.cortMUA(:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortLH.(whiskDataType).MUA.corticalData;
        data.Whisk.(whiskDataType).cortLH.cortGam(:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortLH.(whiskDataType).Gam.corticalData;
        data.Whisk.(whiskDataType).cortLH.cortS(:,:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortLH.(whiskDataType).LFP.corticalS;
        data.Whisk.(whiskDataType).cortLH.cortS_Gam(:,:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortLH.(whiskDataType).LFP.corticalS(49:end,20:23);
        data.Whisk.(whiskDataType).cortLH.cortT(:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortLH.(whiskDataType).LFP.T;
        data.Whisk.(whiskDataType).cortLH.cortF(:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortLH.(whiskDataType).LFP.F;
        % right cortical
        data.Whisk.(whiskDataType).cortRH.HbT(:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortRH.(whiskDataType).CBV_HbT.HbT;
        data.Whisk.(whiskDataType).cortRH.CBV(:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortRH.(whiskDataType).CBV.CBV;
        data.Whisk.(whiskDataType).cortRH.cortMUA(:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortRH.(whiskDataType).MUA.corticalData;
        data.Whisk.(whiskDataType).cortRH.cortGam(:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortRH.(whiskDataType).Gam.corticalData;
        data.Whisk.(whiskDataType).cortRH.cortS(:,:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortRH.(whiskDataType).LFP.corticalS;
        data.Whisk.(whiskDataType).cortRH.cortS_Gam(:,:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortRH.(whiskDataType).LFP.corticalS(49:end,20:23);
        data.Whisk.(whiskDataType).cortRH.cortT(:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortRH.(whiskDataType).LFP.T;
        data.Whisk.(whiskDataType).cortRH.cortF(:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortRH.(whiskDataType).LFP.F;
        % hippocampal
        data.Whisk.(whiskDataType).Hip.hipMUA(:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortLH.(whiskDataType).MUA.hippocampalData;
        data.Whisk.(whiskDataType).Hip.hipGam(:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortLH.(whiskDataType).Gam.hippocampalData;
        data.Whisk.(whiskDataType).Hip.hipS(:,:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortLH.(whiskDataType).LFP.hippocampalS;
        data.Whisk.(whiskDataType).Hip.hipS_Gam(:,:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortLH.(whiskDataType).LFP.hippocampalS(49:end,20:23);
        data.Whisk.(whiskDataType).Hip.hipT(:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortLH.(whiskDataType).LFP.T;
        data.Whisk.(whiskDataType).Hip.hipF(:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortLH.(whiskDataType).LFP.F;
        % time vector
        data.Whisk.(whiskDataType).timeVector(:,aa) = AnalysisResults.Evoked.(animalID).Whisk.cortLH.(whiskDataType).timeVector;
    end
end
% concatenate the data.Whisk from the contra and ipsi data.Whisk
for ee = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,ee};
    data.Whisk.(whiskDataType).HbT = cat(2,data.Whisk.(whiskDataType).cortLH.HbT,data.Whisk.(whiskDataType).cortRH.HbT);
    data.Whisk.(whiskDataType).CBV = cat(2,data.Whisk.(whiskDataType).cortLH.CBV,data.Whisk.(whiskDataType).cortRH.CBV);
    data.Whisk.(whiskDataType).cortMUA = cat(2,data.Whisk.(whiskDataType).cortLH.cortMUA,data.Whisk.(whiskDataType).cortRH.cortMUA);
    data.Whisk.(whiskDataType).cortGam = cat(2,data.Whisk.(whiskDataType).cortLH.cortGam,data.Whisk.(whiskDataType).cortRH.cortGam);
    data.Whisk.(whiskDataType).cortS = cat(3,data.Whisk.(whiskDataType).cortLH.cortS,data.Whisk.(whiskDataType).cortRH.cortS);
    data.Whisk.(whiskDataType).cortS_Gam = cat(3,data.Whisk.(whiskDataType).cortLH.cortS_Gam,data.Whisk.(whiskDataType).cortRH.cortS_Gam);
    data.Whisk.(whiskDataType).cortT = cat(2,data.Whisk.(whiskDataType).cortLH.cortT,data.Whisk.(whiskDataType).cortRH.cortT);
    data.Whisk.(whiskDataType).cortF = cat(2,data.Whisk.(whiskDataType).cortLH.cortF,data.Whisk.(whiskDataType).cortRH.cortF);
end
% concatenate the data.Whisk from the contra and ipsi data.Whisk
for ee = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,ee};
    data.Whisk.(whiskDataType).meanHbT = mean(data.Whisk.(whiskDataType).HbT,2);
    data.Whisk.(whiskDataType).stdHbT = std(data.Whisk.(whiskDataType).HbT,0,2);
    data.Whisk.(whiskDataType).meanCBV = mean(data.Whisk.(whiskDataType).CBV,2);
    data.Whisk.(whiskDataType).stdCBV = std(data.Whisk.(whiskDataType).CBV,0,2);
    data.Whisk.(whiskDataType).meanCortMUA = mean(data.Whisk.(whiskDataType).cortMUA,2);
    data.Whisk.(whiskDataType).stdCortMUA = std(data.Whisk.(whiskDataType).cortMUA,0,2);
    data.Whisk.(whiskDataType).meanCortGam = mean(data.Whisk.(whiskDataType).cortGam,2);
    data.Whisk.(whiskDataType).stdCortGam = std(data.Whisk.(whiskDataType).cortGam,0,2);
    data.Whisk.(whiskDataType).meanCortS = mean(data.Whisk.(whiskDataType).cortS,3).*100;
    data.Whisk.(whiskDataType).mean_CortS_Gam = mean(mean(mean(data.Whisk.(whiskDataType).cortS_Gam.*100,2),1),3);
    data.Whisk.(whiskDataType).std_CortS_Gam = std(mean(mean(data.Whisk.(whiskDataType).cortS_Gam.*100,2),1),0,3);
    data.Whisk.(whiskDataType).meanCortT = mean(data.Whisk.(whiskDataType).cortT,2);
    data.Whisk.(whiskDataType).meanCortF = mean(data.Whisk.(whiskDataType).cortF,2);
    data.Whisk.(whiskDataType).meanHipMUA = mean(data.Whisk.(whiskDataType).Hip.hipMUA,2);
    data.Whisk.(whiskDataType).stdHipMUA = std(data.Whisk.(whiskDataType).Hip.hipMUA,0,2);
    data.Whisk.(whiskDataType).meanHipGam = mean(data.Whisk.(whiskDataType).Hip.hipGam,2);
    data.Whisk.(whiskDataType).stdHipGam = std(data.Whisk.(whiskDataType).Hip.hipGam,0,2);
    data.Whisk.(whiskDataType).meanHipS = mean(data.Whisk.(whiskDataType).Hip.hipS,3).*100;
    data.Whisk.(whiskDataType).mean_HipS_Gam = mean(mean(mean(data.Whisk.(whiskDataType).Hip.hipS_Gam.*100,2),1),3);
    data.Whisk.(whiskDataType).std_HipS_Gam = std(mean(mean(data.Whisk.(whiskDataType).Hip.hipS_Gam.*100,2),1),0,3);
    data.Whisk.(whiskDataType).meanHipT = mean(data.Whisk.(whiskDataType).Hip.hipT,2);
    data.Whisk.(whiskDataType).meanHipF = mean(data.Whisk.(whiskDataType).Hip.hipF,2);
    data.Whisk.(whiskDataType).meanTimeVector = mean(data.Whisk.(whiskDataType).timeVector(:,aa),2);
end
%% set-up and process data
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        data.XCorr.(behavField).cortLH.HbTvLFPxcVals(:,:,aa) = AnalysisResults.XCorr.(animalID).(behavField).cortLH.HbTvLFPxcVals;
        data.XCorr.(behavField).cortLH.LFP_lags(:,:,aa) = AnalysisResults.XCorr.(animalID).(behavField).cortLH.LFP_lags;
        data.XCorr.(behavField).cortLH.F(:,:,aa) = AnalysisResults.XCorr.(animalID).(behavField).cortLH.F;
        data.XCorr.(behavField).cortRH.HbTvLFPxcVals(:,:,aa) = AnalysisResults.XCorr.(animalID).(behavField).cortRH.HbTvLFPxcVals;
        data.XCorr.(behavField).cortRH.LFP_lags(:,:,aa) = AnalysisResults.XCorr.(animalID).(behavField).cortRH.LFP_lags;
        data.XCorr.(behavField).cortRH.F(:,:,aa) = AnalysisResults.XCorr.(animalID).(behavField).cortRH.F;
        data.XCorr.(behavField).cortLH.HbTvMUAxcVals(:,aa) = AnalysisResults.XCorr.(animalID).(behavField).cortLH.HbTvMUAxcVals;
        data.XCorr.(behavField).cortLH.HbTvMUAxcVals_std(:,aa) = AnalysisResults.XCorr.(animalID).(behavField).cortLH.HbTvMUAxcVals_std;
        data.XCorr.(behavField).cortLH.HbTvGAMxcVals(:,aa) = AnalysisResults.XCorr.(animalID).(behavField).cortLH.HbTvGAMxcVals;
        data.XCorr.(behavField).cortLH.HbTvGAMxcVals_std(:,aa) = AnalysisResults.XCorr.(animalID).(behavField).cortLH.HbTvGAMxcVals_std;
        data.XCorr.(behavField).cortLH.MUA_lags(:,aa) = AnalysisResults.XCorr.(animalID).(behavField).cortLH.MUA_lags;
        data.XCorr.(behavField).cortLH.GAM_lags(:,aa) = AnalysisResults.XCorr.(animalID).(behavField).cortLH.GAM_lags;
        data.XCorr.(behavField).cortRH.HbTvMUAxcVals(:,aa) = AnalysisResults.XCorr.(animalID).(behavField).cortRH.HbTvMUAxcVals;
        data.XCorr.(behavField).cortRH.HbTvMUAxcVals_std(:,aa) = AnalysisResults.XCorr.(animalID).(behavField).cortRH.HbTvMUAxcVals_std;
        data.XCorr.(behavField).cortRH.HbTvGAMxcVals(:,aa) = AnalysisResults.XCorr.(animalID).(behavField).cortRH.HbTvGAMxcVals;
        data.XCorr.(behavField).cortRH.HbTvGAMxcVals_std(:,aa) = AnalysisResults.XCorr.(animalID).(behavField).cortRH.HbTvGAMxcVals_std;
        data.XCorr.(behavField).cortRH.MUA_lags(:,aa) = AnalysisResults.XCorr.(animalID).(behavField).cortRH.MUA_lags;
        data.XCorr.(behavField).cortRH.GAM_lags(:,aa) = AnalysisResults.XCorr.(animalID).(behavField).cortRH.GAM_lags;
        data.XCorr.(behavField).animalID{aa,1} = animalID;
        data.XCorr.(behavField).behavior{aa,1} = behavField;
        data.XCorr.(behavField).LH{aa,1} = 'LH';
        data.XCorr.(behavField).RH{aa,1} = 'RH';
    end
end
% concatenate the data.XCorr from the left and right hemispheres
for dd = 1:length(behavFields)
    behavField = behavFields{1,dd};
    data.XCorr.(behavField).cat_HbTvLFPxcVals = cat(3,data.XCorr.(behavField).cortLH.HbTvLFPxcVals,data.XCorr.(behavField).cortRH.HbTvLFPxcVals);
    data.XCorr.(behavField).cat_LFP_lags = cat(3,data.XCorr.(behavField).cortLH.LFP_lags,data.XCorr.(behavField).cortRH.LFP_lags);
    data.XCorr.(behavField).cat_LFP_F = cat(3,data.XCorr.(behavField).cortLH.F,data.XCorr.(behavField).cortRH.F);
    data.XCorr.(behavField).cat_HbTvMUAxcVals = cat(2,data.XCorr.(behavField).cortLH.HbTvMUAxcVals,data.XCorr.(behavField).cortRH.HbTvMUAxcVals);
    data.XCorr.(behavField).cat_HbTvGAMxcVals = cat(2,data.XCorr.(behavField).cortLH.HbTvGAMxcVals,data.XCorr.(behavField).cortRH.HbTvGAMxcVals);
    data.XCorr.(behavField).cat_MUA_lags = cat(2,data.XCorr.(behavField).cortLH.MUA_lags,data.XCorr.(behavField).cortRH.MUA_lags);
    data.XCorr.(behavField).cat_GAM_lags = cat(2,data.XCorr.(behavField).cortLH.GAM_lags,data.XCorr.(behavField).cortRH.GAM_lags);
end
% take the averages of each field through the proper dimension
for ff = 1:length(behavFields)
    behavField = behavFields{1,ff};
    data.XCorr.(behavField).meanHbTvLFPxcVals = mean(data.XCorr.(behavField).cat_HbTvLFPxcVals,3);
    data.XCorr.(behavField).meanLFP_lags = mean(data.XCorr.(behavField).cat_LFP_lags,3);
    data.XCorr.(behavField).meanLFP_F = mean(data.XCorr.(behavField).cat_LFP_F,3);
    data.XCorr.(behavField).meanHbTvMUAxcVals = mean(data.XCorr.(behavField).cat_HbTvMUAxcVals,2);
    data.XCorr.(behavField).stdHbTvMUAxcVals = std(data.XCorr.(behavField).cat_HbTvMUAxcVals,0,2);
    data.XCorr.(behavField).meanHbTvGAMxcVals = mean(data.XCorr.(behavField).cat_HbTvGAMxcVals,2);
    data.XCorr.(behavField).stdHbTvGAMxcVals = std(data.XCorr.(behavField).cat_HbTvGAMxcVals,0,2);
    data.XCorr.(behavField).meanMUA_lags = mean(data.XCorr.(behavField).cat_MUA_lags,2);
    data.XCorr.(behavField).meanGAM_lags = mean(data.XCorr.(behavField).cat_GAM_lags,2);
end
%% set-up and process data.Kernel
data.Kernel.GammaHbT.gamma = [];
data.Kernel.GammaHbT.HbT = [];
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for bb = 1:length(neuralBands)
        neuralBand = neuralBands{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            % pull HRFs from AnalysisResults.mat structure - dc shift each IR function by offset
            data.Kernel.IR.(neuralBand).(behavior).cortLH(aa,:) = AnalysisResults.HRFs.(animalID).(neuralBand).cortLH.(behavior).IR_function_short - AnalysisResults.HRFs.(animalID).(neuralBand).cortLH.(behavior).IR_function_short(1);
            data.Kernel.IR.(neuralBand).(behavior).cortRH(aa,:) = AnalysisResults.HRFs.(animalID).(neuralBand).cortRH.(behavior).IR_function_short - AnalysisResults.HRFs.(animalID).(neuralBand).cortRH.(behavior).IR_function_short(1);
            data.Kernel.IR_gamma.(neuralBand).(behavior).cortLH(aa,:) = AnalysisResults.HRFs.(animalID).(neuralBand).cortLH.(behavior).IR_gammaFunction;
            data.Kernel.IR_gamma.(neuralBand).(behavior).cortRH(aa,:) = AnalysisResults.HRFs.(animalID).(neuralBand).cortRH.(behavior).IR_gammaFunction;
        end
    end
end
% concatenate the data.Kernel from left and right into a single data.Kernel set
for dd = 1:length(neuralBands)
    neuralBand = neuralBands{1,dd};
    for ee = 1:length(behaviors)
        behavior = behaviors{1,ee};
        % impulse
        data.Kernel.IR.(neuralBand).(behavior).cat = cat(1,data.Kernel.IR.(neuralBand).(behavior).cortLH,data.Kernel.IR.(neuralBand).(behavior).cortRH);
        data.Kernel.IR_gamma.(neuralBand).(behavior).cat = cat(1,data.Kernel.IR_gamma.(neuralBand).(behavior).cortLH,data.Kernel.IR_gamma.(neuralBand).(behavior).cortRH);
    end
end
% mean and std of each arousal-state
for ii = 1:length(neuralBands)
    neuralBand = neuralBands{1,ii};
    for jj = 1:length(behaviors)
        behavior = behaviors{1,jj};
        % IR mean
        data.Kernel.IR.(neuralBand).(behavior).mean = mean(data.Kernel.IR.(neuralBand).(behavior).cat,1);
    end
end
%% gamma HRF based on impulse deconvolution
for kk = 1:length(neuralBands)
    neuralBand = neuralBands{1,kk};
    for ll = 1:length(behaviors)
        behavior = behaviors{1,ll};
        [peak,peakIndex] = max(data.Kernel.IR.(neuralBand).(behavior).mean);
        peakTime = peakIndex/samplingRate;
        threeQuarterMax = max(data.Kernel.IR.(neuralBand).(behavior).mean)/(4/3);
        index1 = find(data.Kernel.IR.(neuralBand).(behavior).mean >= threeQuarterMax,1,'first');
        % find where the data.Kernel last rises above half the max.
        index2 = find(data.Kernel.IR.(neuralBand).(behavior).mean >= threeQuarterMax,1,'last');
        threeQuarterWidth = (index2 - index1 + 1)/samplingRate; % FWHM in indexes.
        initVals = [peak,peakTime,threeQuarterWidth];
        % create gamma function based on impulse values
        t = 0:1/samplingRate:5;
        IR_a = ((initVals(2)/initVals(3))^2*8*log10(2));
        IR_beta = ((initVals(3)^2)/initVals(2)/8/log10(2));
        data.Kernel.IR_gamma.(neuralBand).(behavior).repFunc = initVals(1)*(t/initVals(2)).^IR_a.*exp((t - initVals(2))/(-1*IR_beta));
    end
end
%% Fig. 2
summaryFigure = figure('Name','Fig2 (A-H)');
sgtitle('Figure 2 - Turner et al. 2020')
%% [2A - top] cortical MUA contra stim
ax1 = subplot(3,4,1);
% p1 = plot(data.Stim.Contra.mean_timeVector,data.Stim.Contra.mean_CortMUA,'color',colors_HRF2020('rich black'),'LineWidth',2);
hold on
% plot(data.Stim.Contra.mean_timeVector,data.Stim.Contra.mean_CortMUA + data.Stim.Contra.std_CortMUA,'color',colors_HRF2020('battleship grey'),'LineWidth',0.5)
% plot(data.Stim.Contra.mean_timeVector,data.Stim.Contra.mean_CortMUA - data.Stim.Contra.std_CortMUA,'color',colors_HRF2020('battleship grey'),'LineWidth',0.5)
p2 = plot(data.Stim.Contra.mean_timeVector,data.Stim.Contra.mean_CortGam,'color',colors_HRF2020('dark cyan'),'LineWidth',2);
hold on
p1 = plot(data.Stim.Contra.mean_timeVector,data.Stim.Contra.mean_CortMUA,'color',colors_HRF2020('rich black'),'LineWidth',2);
% plot(data.Stim.Contra.mean_timeVector,data.Stim.Contra.mean_CortGam + data.Stim.Contra.std_CortGam,'color',colors_HRF2020('cyan'),'LineWidth',0.5)
% plot(data.Stim.Contra.mean_timeVector,data.Stim.Contra.mean_CortGam - data.Stim.Contra.std_CortGam,'color',colors_HRF2020('cyan'),'LineWidth',0.5)
title('[2A top] Stim cortical neural')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2],'MUA','Gamma')
axis square
set(gca,'box','off')
xlim([-2,5])
ylim([-100,1200])
ax1.TickLength = [0.03,0.03];
%% [2B - top] HbT contra stim
ax5 = subplot(3,4,5);
plot(data.Stim.Contra.mean_timeVector,data.Stim.Contra.mean_HbT,'color',colors_HRF2020('rich black'),'LineWidth',2)
hold on
plot(data.Stim.Contra.mean_timeVector,data.Stim.Contra.mean_HbT + data.Stim.Contra.std_HbT,'color',colors_HRF2020('battleship grey'),'LineWidth',0.5)
plot(data.Stim.Contra.mean_timeVector,data.Stim.Contra.mean_HbT - data.Stim.Contra.std_HbT,'color',colors_HRF2020('battleship grey'),'LineWidth',0.5)
title('[2A bottom] Stim \Delta[HbT] (\muM)')
ylabel('\Delta[HbT] (\muM)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
xlim([-2,5])
ylim([-10,25])
ax5.TickLength = [0.03,0.03];
%% [2A - bottom] moderate whisks cortical MUA
ax2 = subplot(3,4,2);
plot(data.Whisk.IntermediateWhisks.meanTimeVector,data.Whisk.IntermediateWhisks.meanCortMUA,'color',colors_HRF2020('rich black'),'LineWidth',2);
hold on
% plot(data.Whisk.IntermediateWhisks.meanTimeVector,data.Whisk.IntermediateWhisks.meanCortMUA + data.Whisk.IntermediateWhisks.stdCortMUA,'color',colors_HRF2020('battleship grey'),'LineWidth',0.5)
% plot(data.Whisk.IntermediateWhisks.meanTimeVector,data.Whisk.IntermediateWhisks.meanCortMUA - data.Whisk.IntermediateWhisks.stdCortMUA,'color',colors_HRF2020('battleship grey'),'LineWidth',0.5)
plot(data.Whisk.IntermediateWhisks.meanTimeVector,data.Whisk.IntermediateWhisks.meanCortGam,'color',colors_HRF2020('dark cyan'),'LineWidth',2);
% plot(data.Whisk.IntermediateWhisks.meanTimeVector,data.Whisk.IntermediateWhisks.meanCortGam + data.Whisk.IntermediateWhisks.stdCortGam,'color',colors_HRF2020('cyan'),'LineWidth',0.5)
% plot(data.Whisk.IntermediateWhisks.meanTimeVector,data.Whisk.IntermediateWhisks.meanCortGam - data.Whisk.IntermediateWhisks.stdCortGam,'color',colors_HRF2020('cyan'),'LineWidth',0.5)
title('[2B top] Whisk cortical neural')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
xlim([-2,5])
ylim([-10,120])
ax2.TickLength = [0.03,0.03];
%% [2B - bottom] moderate whisks HbT
ax6 = subplot(3,4,6);
plot(data.Whisk.IntermediateWhisks.meanTimeVector,data.Whisk.IntermediateWhisks.meanHbT,'color',colors_HRF2020('rich black'),'LineWidth',2);
hold on
plot(data.Whisk.IntermediateWhisks.meanTimeVector,data.Whisk.IntermediateWhisks.meanHbT + data.Whisk.IntermediateWhisks.stdHbT,'color',colors_HRF2020('battleship grey'),'LineWidth',0.5)
plot(data.Whisk.IntermediateWhisks.meanTimeVector,data.Whisk.IntermediateWhisks.meanHbT - data.Whisk.IntermediateWhisks.stdHbT,'color',colors_HRF2020('battleship grey'),'LineWidth',0.5)
title('[2B bottom] Whisk \Delta[HbT] (\muM)')
ylabel('\Delta[HbT] (\muM)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
xlim([-2,5])
ylim([-10,25])
ax6.TickLength = [0.03,0.03];
%% [2C] Gamma-HbT XCorr
freq = 30;
lag = 5;
ax1 = subplot(3,4,[3,7]);
plot(data.XCorr.Rest.meanGAM_lags,data.XCorr.Rest.meanHbTvGAMxcVals,'color',colorRest,'LineWidth',2)
hold on
plot(data.XCorr.Rest.meanGAM_lags,data.XCorr.Rest.meanHbTvGAMxcVals + data.XCorr.Rest.stdHbTvGAMxcVals,'color',colorRest,'LineWidth',0.5)
plot(data.XCorr.Rest.meanGAM_lags,data.XCorr.Rest.meanHbTvGAMxcVals - data.XCorr.Rest.stdHbTvGAMxcVals,'color',colorRest,'LineWidth',0.5)
plot(data.XCorr.NREM.meanGAM_lags,data.XCorr.NREM.meanHbTvGAMxcVals,'color',colorNREM,'LineWidth',2)
plot(data.XCorr.NREM.meanGAM_lags,data.XCorr.NREM.meanHbTvGAMxcVals + data.XCorr.NREM.stdHbTvGAMxcVals,'color',colorNREM,'LineWidth',0.5)
plot(data.XCorr.NREM.meanGAM_lags,data.XCorr.NREM.meanHbTvGAMxcVals - data.XCorr.NREM.stdHbTvGAMxcVals,'color',colorNREM,'LineWidth',0.5)
plot(data.XCorr.REM.meanGAM_lags,data.XCorr.REM.meanHbTvGAMxcVals,'color',colorREM,'LineWidth',2)
plot(data.XCorr.REM.meanGAM_lags,data.XCorr.REM.meanHbTvGAMxcVals + data.XCorr.REM.stdHbTvGAMxcVals,'color',colorREM,'LineWidth',0.5)
plot(data.XCorr.REM.meanGAM_lags,data.XCorr.REM.meanHbTvGAMxcVals - data.XCorr.REM.stdHbTvGAMxcVals,'color',colorREM,'LineWidth',0.5)
title('[2C] Gam-[HbT] XCorr')
xticks([-lag*freq,-lag*freq/2,0,lag*freq/2,lag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-lag*freq,lag*freq])
xlabel('Lags (s)')
ylabel({'Corr. coefficient';'Gamma vs. \Delta[HbT] (\muM)'})
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
ylim([-0.1,0.5])
%% [2D] MUA-HbT XCorr
ax1 = subplot(3,4,[4,8]);
plot(data.XCorr.Rest.meanMUA_lags,data.XCorr.Rest.meanHbTvMUAxcVals,'color',colorRest,'LineWidth',2)
hold on
plot(data.XCorr.Rest.meanMUA_lags,data.XCorr.Rest.meanHbTvMUAxcVals + data.XCorr.Rest.stdHbTvMUAxcVals,'color',colorRest,'LineWidth',0.5)
plot(data.XCorr.Rest.meanMUA_lags,data.XCorr.Rest.meanHbTvMUAxcVals - data.XCorr.Rest.stdHbTvMUAxcVals,'color',colorRest,'LineWidth',0.5)
plot(data.XCorr.NREM.meanMUA_lags,data.XCorr.NREM.meanHbTvMUAxcVals,'color',colorNREM,'LineWidth',2)
plot(data.XCorr.NREM.meanMUA_lags,data.XCorr.NREM.meanHbTvMUAxcVals + data.XCorr.NREM.stdHbTvMUAxcVals,'color',colorNREM,'LineWidth',0.5)
plot(data.XCorr.NREM.meanMUA_lags,data.XCorr.NREM.meanHbTvMUAxcVals - data.XCorr.NREM.stdHbTvMUAxcVals,'color',colorNREM,'LineWidth',0.5)
plot(data.XCorr.REM.meanMUA_lags,data.XCorr.REM.meanHbTvMUAxcVals,'color',colorREM,'LineWidth',2)
plot(data.XCorr.REM.meanMUA_lags,data.XCorr.REM.meanHbTvMUAxcVals + data.XCorr.REM.stdHbTvMUAxcVals,'color',colorREM,'LineWidth',0.5)
plot(data.XCorr.REM.meanMUA_lags,data.XCorr.REM.meanHbTvMUAxcVals - data.XCorr.REM.stdHbTvMUAxcVals,'color',colorREM,'LineWidth',0.5)
title('[2D] MUA-[HbT] XCorr')
xticks([-lag*freq,-lag*freq/2,0,lag*freq/2,lag*freq])
xticklabels({'-5','-2.5','0','2.5','5'})
xlim([-lag*freq,lag*freq])
xlabel('Lags (s)')
ylabel({'Corr. coefficient';'MUA vs. \Delta[HbT] (\muM)'})
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
ylim([-0.1,0.5])
%% [2E] gamma-band impulse function
subplot(3,4,9);
p1 = plot(t,data.Kernel.IR.gammaBandPower.Rest.mean,'color',colorRest,'LineWidth',2);
hold on
p2 = plot(t,data.Kernel.IR.gammaBandPower.Whisk.mean,'color',colorWhisk,'LineWidth',2);
p3 = plot(t,data.Kernel.IR.gammaBandPower.Contra.mean,'color',colorStim,'LineWidth',2);
p4 = plot(t,data.Kernel.IR.gammaBandPower.NREM.mean,'color',colorNREM,'LineWidth',2);
p5 = plot(t,data.Kernel.IR.gammaBandPower.REM.mean,'color',colorREM,'LineWidth',2);
p6 = plot(t,data.Kernel.IR.gammaBandPower.All.mean,'color',colorAll,'LineWidth',2);
title('[2E] Gamma [30-100 Hz] IRF')
ylabel('Amplitude (a.u.)')
xlabel('Time (s)')
legend([p1,p2,p3,p4,p5,p6],'Rest','Whisk','Stim','NREM','REM','All','Location','NorthEast')
axis square
xlim([0,5])
ylim([-0.1,0.8])
set(gca,'box','off')
%% [2F] gamma-band gamma function
subplot(3,4,10);
plot(t,data.Kernel.IR_gamma.gammaBandPower.Rest.repFunc,'color',colorRest,'LineWidth',2);
hold on
plot(t,data.Kernel.IR_gamma.gammaBandPower.Whisk.repFunc,'color',colorWhisk,'LineWidth',2);
plot(t,data.Kernel.IR_gamma.gammaBandPower.Contra.repFunc,'color',colorStim,'LineWidth',2);
plot(t,data.Kernel.IR_gamma.gammaBandPower.NREM.repFunc,'color',colorNREM,'LineWidth',2);
plot(t,data.Kernel.IR_gamma.gammaBandPower.REM.repFunc,'color',colorREM,'LineWidth',2);
plot(t,data.Kernel.IR_gamma.gammaBandPower.All.repFunc,'color',colorAll,'LineWidth',2);
title('[2F] Gamma-band [30-100 Hz] GF')
ylabel('Amplitude (a.u.)')
xlabel('Time (s)')
axis square
xlim([0,5])
ylim([-0.1,0.8])
set(gca,'box','off')
%% [2G] MUA impulse function
subplot(3,4,11);
plot(t,data.Kernel.IR.muaPower.Rest.mean,'color',colorRest,'LineWidth',2);
hold on
plot(t,data.Kernel.IR.muaPower.Whisk.mean,'color',colorWhisk,'LineWidth',2);
plot(t,data.Kernel.IR.muaPower.Contra.mean,'color',colorStim,'LineWidth',2);
plot(t,data.Kernel.IR.muaPower.NREM.mean,'color',colorNREM,'LineWidth',2);
plot(t,data.Kernel.IR.muaPower.REM.mean,'color',colorREM,'LineWidth',2);
plot(t,data.Kernel.IR.muaPower.All.mean,'color',colorAll,'LineWidth',2);
title('[2G] MUA [300-3000 Hz] IRF')
ylabel('Amplitude (a.u.)')
xlabel('Time (s)')
axis square
xlim([0,5])
ylim([-0.2,1.5])
set(gca,'box','off')
%% [2H] MUA gamma function
subplot(3,4,12);
plot(t,data.Kernel.IR_gamma.muaPower.Rest.repFunc,'color',colorRest,'LineWidth',2);
hold on
plot(t,data.Kernel.IR_gamma.muaPower.Whisk.repFunc,'color',colorWhisk,'LineWidth',2);
plot(t,data.Kernel.IR_gamma.muaPower.Contra.repFunc,'color',colorStim,'LineWidth',2);
plot(t,data.Kernel.IR_gamma.muaPower.NREM.repFunc,'color',colorNREM,'LineWidth',2);
plot(t,data.Kernel.IR_gamma.muaPower.REM.repFunc,'color',colorREM,'LineWidth',2);
plot(t,data.Kernel.IR_gamma.muaPower.All.repFunc,'color',colorAll,'LineWidth',2);
title('[2H] MUA [300-3000 Hz] GF')
ylabel('Amplitude (a.u.)')
xlabel('Time (s)')
axis square
xlim([0,5])
ylim([-0.2,1.5])
set(gca,'box','off')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig2']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig2'])
    %     %% text diary
    %     diaryFile = [dirpath 'Fig1-S2_Statistics.txt'];
    %     if exist(diaryFile,'file') == 2
    %         delete(diaryFile)
    %     end
    %     diary(diaryFile)
    %     diary on
    %     % text values
    %     disp('======================================================================================================================')
    %     disp('[1-S2] Text values for gamma/HbT/reflectance changes')
    %     disp('======================================================================================================================')
    %     disp('----------------------------------------------------------------------------------------------------------------------')
    %     % cortical MUA/LFP
    %     [~,index] = max(data.Stim.Contra.mean_CortMUA);
    %     disp(['Contra stim Cort gamma MUA P/P (%): ' num2str(round(data.Stim.Contra.mean_CortMUA(index),1)) ' +/- ' num2str(round(data.Stim.Contra.std_CortMUA(index),1))]); disp(' ')
    %     [~,index] = max(data.Stim.Ipsi.mean_CortMUA);
    %     disp(['Ipsil stim Cort gamma MUA P/P (%): ' num2str(round(data.Stim.Ipsi.mean_CortMUA(index),1)) ' +/- ' num2str(round(data.Stim.Ipsi.std_CortMUA(index),1))]); disp(' ')
    %     [~,index] = max(data.Stim.Auditory.mean_CortMUA);
    %     disp(['Audit stim Cort gamma MUA P/P (%): ' num2str(round(data.Stim.Auditory.mean_CortMUA(index),1)) ' +/- ' num2str(round(data.Stim.Auditory.std_CortMUA(index),1))]); disp(' ')
    %     % cortical LFP
    %     disp(['Contra stim Cort gamma LFP P/P (%): ' num2str(round(data.Stim.Contra.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.Stim.Contra.std_CortS_Gam,1))]); disp(' ')
    %     disp(['Ipsil stim Cort gamma LFP P/P (%): ' num2str(round(data.Stim.Ipsi.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.Stim.Ipsi.std_CortS_Gam,1))]); disp(' ')
    %     disp(['Audit stim Cort gamma LFP P/P (%): ' num2str(round(data.Stim.Auditory.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.Stim.Auditory.std_CortS_Gam,1))]); disp(' ')
    %     % hippocampal MUA
    %     [~,index] = max(data.Stim.Contra.mean_HipMUA);
    %     disp(['Contra stim Hip gamma MUA P/P (%): ' num2str(round(data.Stim.Contra.mean_HipMUA(index),1)) ' +/- ' num2str(round(data.Stim.Contra.std_HipMUA(index),1))]); disp(' ')
    %     [~,index] = max(data.Stim.Ipsi.mean_HipMUA);
    %     disp(['Ipsil stim Hip gamma MUA P/P (%): ' num2str(round(data.Stim.Ipsi.mean_HipMUA(index),1)) ' +/- ' num2str(round(data.Stim.Ipsi.std_HipMUA(index),1))]); disp(' ')
    %     [~,index] = max(data.Stim.Auditory.mean_HipMUA);
    %     disp(['Audit stim Hip gamma MUA P/P (%): ' num2str(round(data.Stim.Auditory.mean_HipMUA(index),1)) ' +/- ' num2str(round(data.Stim.Auditory.std_HipMUA(index),1))]); disp(' ')
    %     % hippocampal LFP
    %     disp(['Contra stim Hip gamma LFP P/P (%): ' num2str(round(data.Stim.Contra.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Stim.Contra.std_HipS_Gam,1))]); disp(' ')
    %     disp(['Ipsil stim Hip gamma LFP P/P (%): ' num2str(round(data.Stim.Ipsi.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Stim.Ipsi.std_HipS_Gam,1))]); disp(' ')
    %     disp(['Auditory stim Hip gamma LFP P/P (%): ' num2str(round(data.Stim.Auditory.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Stim.Auditory.std_HipS_Gam,1))]); disp(' ')
    %     % HbT
    %     [~,index] = max(data.Stim.Contra.mean_HbT);
    %     disp(['Contra stim [HbT] (uM): ' num2str(round(data.Stim.Contra.mean_HbT(index),1)) ' +/- ' num2str(round(data.Stim.Contra.std_HbT(index),1))]); disp(' ')
    %     [~,index] = max(data.Stim.Ipsi.mean_HbT);
    %     disp(['Ipsil stim [HbT] (uM): ' num2str(round(data.Stim.Ipsi.mean_HbT(index),1)) ' +/- ' num2str(round(data.Stim.Ipsi.std_HbT(index),1))]); disp(' ')
    %     [~,index] = max(data.Stim.Auditory.mean_HbT);
    %     disp(['Audit stim [HbT] (uM): ' num2str(round(data.Stim.Auditory.mean_HbT(index),1)) ' +/- ' num2str(round(data.Stim.Auditory.std_HbT(index),1))]); disp(' ')
    %     % R/R
    %     [~,index] = min(data.Stim.Contra.mean_CBV);
    %     disp(['Contra stim refl R/R (%): ' num2str(round(data.Stim.Contra.mean_CBV(index),1)) ' +/- ' num2str(round(data.Stim.Contra.std_CBV(index),1))]); disp(' ')
    %     [~,index] = min(data.Stim.Ipsi.mean_CBV);
    %     disp(['Ipsil stim refl R/R (%): ' num2str(round(data.Stim.Ipsi.mean_CBV(index),1)) ' +/- ' num2str(round(data.Stim.Ipsi.std_CBV(index),1))]); disp(' ')
    %     [~,index] = min(data.Stim.Auditory.mean_CBV);
    %     disp(['Audit stim refl R/R (%): ' num2str(round(data.Stim.Auditory.mean_CBV(index),1)) ' +/- ' num2str(round(data.Stim.Auditory.std_CBV(index),1))]); disp(' ')
    %     disp('----------------------------------------------------------------------------------------------------------------------')
    %     diary off
    %
    %
    %
    %
    %
    %     %%
    %      % text values
    %     disp('======================================================================================================================')
    %     disp('[1-S3] Text values for gamma/HbT changes')
    %     disp('======================================================================================================================')
    %     disp('----------------------------------------------------------------------------------------------------------------------')
    %      % cortical MUA/LFP
    %     [~,index] = max(data.Whisk.ShortWhisks.meanCortMUA);
    %     disp(['Brief whisk Cort gamma MUA P/P (%): ' num2str(round(data.Whisk.ShortWhisks.meanCortMUA(index),1)) ' +/- ' num2str(round(data.Whisk.ShortWhisks.stdCortMUA(index),1))]); disp(' ')
    %     [~,index] = max(data.Whisk.IntermediateWhisks.meanCortMUA);
    %     disp(['Moderate whisk Cort gamma MUA P/P (%): ' num2str(round(data.Whisk.IntermediateWhisks.meanCortMUA(index),1)) ' +/- ' num2str(round(data.Whisk.IntermediateWhisks.stdCortMUA(index),1))]); disp(' ')
    %     [~,index] = max(data.Whisk.LongWhisks.meanCortMUA);
    %     disp(['Extended whisk Cort gamma MUA P/P (%): ' num2str(round(data.Whisk.LongWhisks.meanCortMUA(index),1)) ' +/- ' num2str(round(data.Whisk.LongWhisks.stdCortMUA(index),1))]); disp(' ')
    %     % cortical LFP
    %     disp(['Brief whisk Cort gamma LFP P/P (%): ' num2str(round(data.Whisk.ShortWhisks.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.Whisk.ShortWhisks.std_CortS_Gam,1))]); disp(' ')
    %     disp(['Moderate whisk Cort gamma LFP P/P (%): ' num2str(round(data.Whisk.IntermediateWhisks.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.Whisk.IntermediateWhisks.std_CortS_Gam,1))]); disp(' ')
    %     disp(['Extended whisk Cort gamma LFP P/P (%): ' num2str(round(data.Whisk.LongWhisks.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.Whisk.LongWhisks.std_CortS_Gam,1))]); disp(' ')
    %     % hippocampal MUA
    %     [~,index] = max(data.Whisk.ShortWhisks.meanHipMUA);
    %     disp(['Brief whisk Hip gamma MUA P/P (%): ' num2str(round(data.Whisk.ShortWhisks.meanHipMUA(index),1)) ' +/- ' num2str(round(data.Whisk.ShortWhisks.stdHipMUA(index),1))]); disp(' ')
    %     [~,index] = max(data.Whisk.IntermediateWhisks.meanHipMUA);
    %     disp(['Moderate whisk Hip gamma MUA P/P (%): ' num2str(round(data.Whisk.IntermediateWhisks.meanHipMUA(index),1)) ' +/- ' num2str(round(data.Whisk.IntermediateWhisks.stdHipMUA(index),1))]); disp(' ')
    %     [~,index] = max(data.Whisk.LongWhisks.meanHipMUA);
    %     disp(['Extended whisk Hip gamma MUA P/P (%): ' num2str(round(data.Whisk.LongWhisks.meanHipMUA(index),1)) ' +/- ' num2str(round(data.Whisk.LongWhisks.stdHipMUA(index),1))]); disp(' ')
    %     % hippocampal LFP
    %     disp(['Brief whisk Hip gamma LFP P/P (%): ' num2str(round(data.Whisk.ShortWhisks.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Whisk.ShortWhisks.std_HipS_Gam,1))]); disp(' ')
    %     disp(['Moderate whisk Hip gamma LFP P/P (%): ' num2str(round(data.Whisk.IntermediateWhisks.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Whisk.IntermediateWhisks.std_HipS_Gam,1))]); disp(' ')
    %     disp(['Extended whisk Hip gamma LFP P/P (%): ' num2str(round(data.Whisk.LongWhisks.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Whisk.LongWhisks.std_HipS_Gam,1))]); disp(' ')
    %     % HbT
    %     [~,index] = max(data.Whisk.ShortWhisks.meanHbT);
    %     disp(['Brief whisk [HbT] (uM): ' num2str(round(data.Whisk.ShortWhisks.meanHbT(index),1)) ' +/- ' num2str(round(data.Whisk.ShortWhisks.stdHbT(index),1))]); disp(' ')
    %     [~,index] = max(data.Whisk.IntermediateWhisks.meanHbT);
    %     disp(['Moderate whisk [HbT] (uM): ' num2str(round(data.Whisk.IntermediateWhisks.meanHbT(index),1)) ' +/- ' num2str(round(data.Whisk.IntermediateWhisks.stdHbT(index),1))]); disp(' ')
    %     [~,index] = max(data.Whisk.LongWhisks.meanHbT);
    %     disp(['Extended whisk [HbT] (uM): ' num2str(round(data.Whisk.LongWhisks.meanHbT(index),1)) ' +/- ' num2str(round(data.Whisk.LongWhisks.stdHbT(index),1))]); disp(' ')
    %     % R/R
    %     [~,index] = min(data.Whisk.ShortWhisks.meanCBV);
    %     disp(['Brief whisk refl R/R (%): ' num2str(round(data.Whisk.ShortWhisks.meanCBV(index),1)) ' +/- ' num2str(round(data.Whisk.ShortWhisks.stdCBV(index),1))]); disp(' ')
    %     [~,index] = min(data.Whisk.IntermediateWhisks.meanCBV);
    %     disp(['Moderate whisk refl R/R (%): ' num2str(round(data.Whisk.IntermediateWhisks.meanCBV(index),1)) ' +/- ' num2str(round(data.Whisk.IntermediateWhisks.stdCBV(index),1))]); disp(' ')
    %     [~,index] = min(data.Whisk.LongWhisks.meanCBV);
    %     disp(['Extended whisk refl R/R (%): ' num2str(round(data.Whisk.LongWhisks.meanCBV(index),1)) ' +/- ' num2str(round(data.Whisk.LongWhisks.stdCBV(index),1))]); disp(' ')
    %     disp('----------------------------------------------------------------------------------------------------------------------')
    %     diary off
end

end
