function [ComparisonData] = AnalyzeCorrCoeffs(animal, params, procDataFiles, RestingBaselines, RestData, SleepData, ComparisonData)
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

samplingRate = RestData.CBV.LH.CBVCamSamplingRate;

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

for n = 1:length(LH_finalRestCBVData)
    restCBV_CC = corrcoef(LH_finalRestCBVData{n, 1}, RH_finalRestCBVData{n, 1});
    restCBV_R(n, 1) = restCBV_CC(2, 1);
    restGAM_CC = corrcoef(LH_finalRestGAMData{n, 1}, RH_finalRestGAMData{n, 1});
    restGAM_R(n, 1) = restGAM_CC(2, 1);
end

nboot = 1000;
restCBVCC_CI = bootci(nboot, @mean, restCBV_R);
restGAMCC_CI = bootci(nboot, @mean, restGAM_R);

restCBV_CC_mean = mean(restCBV_R);
restCBV_CC_std = std(restCBV_R);
disp([' Resting CBV correlation coefficient is  ' num2str(restCBV_CC_mean) ' with a STD of: ' num2str(restCBV_CC_std)]); disp(' ')

restGAM_CC_mean = mean(restGAM_R);
restGAM_CC_std = std(restGAM_R);
disp([' Resting Gamma correlation coefficient is  ' num2str(restGAM_CC_mean) ' with a STD of: ' num2str(restGAM_CC_std)]); disp(' ')

ComparisonData.CorrCoeff.Rest.CBV.mean = restCBV_CC_mean;
ComparisonData.CorrCoeff.Rest.CBV.std = restCBV_CC_std;
ComparisonData.CorrCoeff.Rest.GAM.mean = restGAM_CC_mean;
ComparisonData.CorrCoeff.Rest.GAM.std = restGAM_CC_std;

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
    
    for n = 1:length(LH_nremCBV)
        nremCBV_CC = corrcoef(LH_nremCBV{n, 1}, RH_nremCBV{n, 1});
        nremCBV_R(n, 1) = nremCBV_CC(2, 1);
        nremGAM_CC = corrcoef(LH_nremGAM{n, 1}, RH_nremGAM{n, 1});
        nremGAM_R(n, 1) = nremGAM_CC(2, 1);
    end
    
    nremCBVCC_CI = bootci(nboot, @mean, nremCBV_R);
    nremGAMCC_CI = bootci(nboot, @mean, nremGAM_R);

    nremCBV_CC_mean = mean(nremCBV_R);
    nremCBV_CC_std = std(nremCBV_R);
    disp([' NREM CBV correlation coefficient is  ' num2str(nremCBV_CC_mean) ' with a STD of: ' num2str(nremCBV_CC_std)]); disp(' ')
    
    nremGAM_CC_mean = mean(nremGAM_R);
    nremGAM_CC_std = std(nremGAM_R);
    disp([' NREM Gamma correlation coefficient is  ' num2str(nremGAM_CC_mean) ' with a STD of: ' num2str(nremGAM_CC_std)]); disp(' ')
    
    ComparisonData.CorrCoeff.NREM.CBV.mean = nremCBV_CC_mean;
    ComparisonData.CorrCoeff.NREM.CBV.std = nremCBV_CC_std;
    ComparisonData.CorrCoeff.NREM.GAM.mean = nremGAM_CC_mean;
    ComparisonData.CorrCoeff.NREM.GAM.std = nremGAM_CC_std;
    
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
        
        for n = 1:length(LH_remCBV)
            remCBV_CC = corrcoef(LH_remCBV{n, 1}, RH_remCBV{n, 1});
            remCBV_R(n, 1) = remCBV_CC(2, 1);
            remGAM_CC = corrcoef(LH_remGAM{n, 1}, RH_remGAM{n, 1});
            remGAM_R(n, 1) = remGAM_CC(2, 1);
        end
        
        remCBV_CC_mean = mean(remCBV_R);
        remCBV_CC_std = std(remCBV_R);
        disp([' REM CBV correlation coefficient is  ' num2str(remCBV_CC_mean) ' with a STD of: ' num2str(remCBV_CC_std)]); disp(' ')
        
        remGAM_CC_mean = mean(remGAM_R);
        remGAM_CC_std = std(remGAM_R);
        disp([' REM Gamma correlation coefficient is  ' num2str(remGAM_CC_mean) ' with a STD of: ' num2str(remGAM_CC_std)]); disp(' ')
        
        ComparisonData.CorrCoeff.REM.CBV.mean = remCBV_CC_mean;
        ComparisonData.CorrCoeff.REM.CBV.std = remCBV_CC_std;
        ComparisonData.CorrCoeff.REM.GAM.mean = remGAM_CC_mean;
        ComparisonData.CorrCoeff.REM.GAM.std = remGAM_CC_std;
    end
end

% %% Load in relevant data from all ProcDataFiles
% for pDF = 1:size(procDataFiles, 1)
%     procDataFile = procDataFiles(pDF, :);
%     load(procDataFile);
%     [~, ~, fileDate, fileID] = GetFileInfo_IOS(procDataFile);
%     strDay = ConvertDate(fileDate);
%     
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
% for n = 1:length(LH_CBV)
%     allCBV_CC = corrcoef(LH_CBV{n, 1}, RH_CBV{n, 1});
%     allCBV_R(n, 1) = allCBV_CC(2, 1);
%     allGAM_CC = corrcoef(LH_GAM{n, 1}, RH_GAM{n, 1});
%     allGAM_R(n, 1) = allGAM_CC(2, 1);
% end
% 
% allCBV_CC_mean = mean(allCBV_R);
% allCBV_CC_std = std(allCBV_R);
% disp([' all data CBV correlation coefficient is  ' num2str(allCBV_CC_mean) ' with a STD of: ' num2str(allCBV_CC_std)]); disp(' ')
% 
% allGAM_CC_mean = mean(allGAM_R);
% allGAM_CC_std = std(allGAM_R);
% disp([' all data Gamma correlation coefficient is  ' num2str(allGAM_CC_mean) ' with a STD of: ' num2str(allGAM_CC_std)]); disp(' ')
% 
% ComparisonData.CorrCoeff.AllData.CBV.mean = allCBV_CC_mean;
% ComparisonData.CorrCoeff.AllData.CBV.std = allCBV_CC_std;
% ComparisonData.CorrCoeff.AllData.GAM.mean = allGAM_CC_mean;
% ComparisonData.CorrCoeff.AllData.GAM.std = allGAM_CC_std;

figure;
errorbar(1, restCBV_CC_mean, abs(restCBV_CC_mean - restCBVCC_CI(1)), abs(restCBV_CC_mean - restCBVCC_CI(2)), 'o')
hold on
errorbar(2, restGAM_CC_mean, abs(restGAM_CC_mean - restGAMCC_CI(1)), abs(restGAM_CC_mean - restGAMCC_CI(2)), 'o')
errorbar(3, nremCBV_CC_mean, abs(nremCBV_CC_mean - nremCBVCC_CI(1)), abs(nremCBV_CC_mean - nremCBVCC_CI(2)), 'o')
errorbar(4, nremGAM_CC_mean, abs(nremGAM_CC_mean - nremGAMCC_CI(1)), abs(nremGAM_CC_mean - nremGAMCC_CI(2)), 'o')
xlim([0 5])
ylim([0 1])
title('Pearsons correlation coefficient with bootstrap confidence intervals')
ylabel('Correlation coefficient')
legend('Rest CBV', 'Rest Gamma', 'NREM CBV', 'NREM Gamma')
set(gca,'xticklabel',{[]})

save([animal '_ComparisonData.mat'], 'ComparisonData');

end