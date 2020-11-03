function [AnalysisResults] = AnalyzeTwoPSleepProbability_eLife2020(TwoP_animalIDs,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the probability of arousal-state classification based on arteriole D/D changes (2PLSM)
%________________________________________________________________________________________________________________________

%% function parameters
allCatLabels = [];
allCatMeans = [];
%% extract data from each animal's sleep scoring results
for aa = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,aa};
    dataLoc = [rootFolder '/' animalID '/2P Data/'];
    cd(dataLoc)
    % character list of all MergedData files
    mergedDirectory = dir('*_MergedData.mat');
    mergedDataFiles = {mergedDirectory.name}';
    mergedDataFileIDs = char(mergedDataFiles);
    % Character list of all TrainingData files
    trainingDirectory = dir('*_TrainingData.mat');
    trainingDataFiles = {trainingDirectory.name}';
    trainingDataFileIDs = char(trainingDataFiles);
    % resting baselines structure
    baselinesFile = [animalID '_RestingBaselines.mat'];
    load(baselinesFile)
    binSize = 5;   % seconds
    numBins = 180;   % 15 minutes with 5 sec bins
    samplingRate = 5;   % Hz
    samplesPerBin = binSize*samplingRate;
    for bb = 1:size(mergedDataFileIDs,1)
        mergedDataFileID = mergedDataFileIDs(bb,:);
        trainingDataFileID = trainingDataFileIDs(bb,:);
        load(mergedDataFileID,'-mat')
        load(trainingDataFileID,'-mat')
        [~,~,fileDate,~,~,vesselID] = GetFileInfo2_2P_eLife2020(mergedDataFileID);
        if strcmp(vesselID(1),'V') == false
            strDay = ConvertDate_2P_eLife2020(fileDate);
            for cc = 1:numBins
                if cc == 1
                    binSamples = MergedData.data.vesselDiameter.data(1:samplesPerBin);
                    normBinSamples = (binSamples - RestingBaselines.manualSelection.vesselDiameter.data.(vesselID).(strDay))./RestingBaselines.manualSelection.vesselDiameter.data.(vesselID).(strDay);
                else
                    binSamples = MergedData.data.vesselDiameter.data((cc - 1)*samplesPerBin + 1:cc*samplesPerBin);
                    normBinSamples = (binSamples - RestingBaselines.manualSelection.vesselDiameter.data.(vesselID).(strDay))./RestingBaselines.manualSelection.vesselDiameter.data.(vesselID).(strDay);
                end
                allCatMeans = cat(1,allCatMeans,mean(normBinSamples)*100);
            end
            allCatLabels = cat(1,allCatLabels,trainingTable.behavState);
        end
    end
end
%% change arousal-state labels to numbers
zAwake = 1; zNrem = 1; zRem = 1;
for zz = 1:length(allCatLabels)
    if strcmp(allCatLabels{zz,1},'Not Sleep') == true
        awakeTwoP(zAwake,1) = allCatMeans(zz,1); %#ok<*NASGU>
        zAwake = zAwake + 1;
    elseif strcmp(allCatLabels{zz,1},'NREM Sleep') == true
        nremTwoP(zNrem,1) = allCatMeans(zz,1);
        zNrem = zNrem + 1;
    elseif strcmp(allCatLabels{zz,1},'REM Sleep') == true
        remTwoP(zRem,1) = allCatMeans(zz,1);
        zRem = zRem + 1;
    end
end
%% put each mean and scoring label into a cell
minTwoP = floor(min(allCatMeans));
maxTwoP = ceil(max(allCatMeans));
awakeBins = minTwoP:1:maxTwoP;
cutDown = abs(minTwoP - (-20));
cutUp = maxTwoP - 50;
probBinLabels = cell(length(minTwoP:1:maxTwoP),1);
probBinMeans = cell(length(minTwoP:1:maxTwoP),1);
discBins = discretize(allCatMeans,awakeBins);
for dd = 1:length(discBins)
    probBinLabels{discBins(dd),1} = cat(1,probBinLabels{discBins(dd),1},{allCatLabels(dd,1)});
    probBinMeans{discBins(dd),1} = cat(1,probBinMeans{discBins(dd),1},{allCatMeans(dd,1)});
end
% condense the left edges of the histogram bins to -35:1:120
cutDownLabels = [];
cutDownMeans = [];
for ee = 1:cutDown
    cutDownLabels = cat(1,cutDownLabels,probBinLabels{ee,1});
    cutDownMeans = cat(1,cutDownMeans,probBinMeans{ee,1});
end
% condense the right edges of the histogram to -35:1:120
cutUpLabels = [];
cutUpMeans = [];
for ff = 1:cutUp + 1
    cutUpLabels = cat(1,cutUpLabels,probBinLabels{end - ff,1});
    cutUpMeans = cat(1,cutUpMeans,probBinMeans{end - ff,1});
end
% reconstruct array of labels based on new edges
finCatLabels = cat(1,{cutDownLabels},probBinLabels(cutDown + 1:end - (cutUp + 2)),{cutUpLabels});
% strcmp the bins and if the bin is asleep (NREM/REM) set to 0, else set 1
for gg = 1:length(finCatLabels)
    for hh = 1:length(finCatLabels{gg,1})
        if strcmp(finCatLabels{gg,1}{hh,1}{1,1},'Not Sleep') == true
            awakeProbEvents{gg,1}(hh,1) = 1; %#ok<*AGROW>
        else
            awakeProbEvents{gg,1}(hh,1) = 0;
        end
    end
end
% strcmp the bins and if the bin is not in NREM (Awake/REM) set to 0, else set 1
for gg = 1:length(finCatLabels)
    for hh = 1:length(finCatLabels{gg,1})
        if strcmp(finCatLabels{gg,1}{hh,1}{1,1},'NREM Sleep') == true
            nremProbEvents{gg,1}(hh,1) = 1;
        else
            nremProbEvents{gg,1}(hh,1) = 0;
        end
    end
end
% strcmp the bins and if the bin is not in REM (Awake/NREM) set to 0, else set 1
for gg = 1:length(finCatLabels)
    for hh = 1:length(finCatLabels{gg,1})
        if strcmp(finCatLabels{gg,1}{hh,1}{1,1},'REM Sleep') == true
            remProbEvents{gg,1}(hh,1) = 1;
        else
            remProbEvents{gg,1}(hh,1) = 0;
        end
    end
end
% take probability of each bin
for ii = 1:length(awakeProbEvents)
    awakeProbPerc(ii,1) = sum(awakeProbEvents{ii,1})/length(awakeProbEvents{ii,1})*100;
    nremProbPerc(ii,1) = sum(nremProbEvents{ii,1})/length(nremProbEvents{ii,1})*100;
    remProbPerc(ii,1) = sum(remProbEvents{ii,1})/length(remProbEvents{ii,1})*100;
end
% save results
AnalysisResults.TwoPSleepProbability.TwoPCatMeans = allCatMeans;
AnalysisResults.TwoPSleepProbability.awakeProbPerc = awakeProbPerc;
AnalysisResults.TwoPSleepProbability.nremProbPerc = nremProbPerc;
AnalysisResults.TwoPSleepProbability.remProbPerc = remProbPerc;
% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults')

end
