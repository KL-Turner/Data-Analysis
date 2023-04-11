function [SleepData] = CreateSleepData_GCaMP_IOS(NREMsleepTime,REMsleepTime,modelName,TrainingFiles,SleepData)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: This function uses the sleep logicals in each ProcData file to find periods where there are 60 seconds of
%          consecutive ones within the sleep logical (12 or more). If a ProcData file's sleep logical contains one or
%          more of these 60 second periods,each of those bins is gathered from the data and put into the SleepEventData.mat
%          struct along with the file's name.
%________________________________________________________________________________________________________________________

% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
if strcmp(modelName,'Manual') == true
    cc = 1;
    % reduce file list to those with the training dates
    for aa = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(aa,:);
        [~,fileDate,~] = GetFileInfo_IOS(procDataFileID);
        if strcmp(fileDate,TrainingFiles.day1) == true || strcmp(fileDate,TrainingFiles.day2) == true
            trainingFileList(cc,:) = procDataFileID;
            cc = cc + 1;
        end
    end
    procDataFileIDs = trainingFileList;
end
% create NREM sleep scored data structure.
% identify sleep epochs and place in SleepEventData.mat structure
sleepBins = NREMsleepTime/5;
for aa = 1:size(procDataFileIDs,1) % loop through the list of ProcData files
    clearvars -except aa procDataFileIDs sleepBins NREMsleepTime REMsleepTime modelName SleepData startingDirectory
    procDataFileID = procDataFileIDs(aa,:); % pull character string associated with the current file
    load(procDataFileID); % load in procDataFile associated with character string
    [~,~,fileID] = GetFileInfo_IOS(procDataFileID); % gather file info
    nremLogical = ProcData.sleep.logicals.(modelName).nremLogical; % logical - ones denote potential sleep epoches (5 second bins)
    targetTime = ones(1,sleepBins); % target time
    sleepIndex = find(conv(nremLogical,targetTime) >= sleepBins) - (sleepBins - 1); % find the periods of time where there are at least 11 more
    % 5 second epochs following. This is not the full list.
    if isempty(sleepIndex) % if sleepIndex is empty,skip this file
        % skip file
    else
        sleepCriteria = (0:(sleepBins - 1)); % this will be used to fix the issue in sleepIndex
        fixedSleepIndex = unique(sleepIndex + sleepCriteria); % sleep Index now has the proper time stamps from sleep logical
        for indexCount = 1:length(fixedSleepIndex) % loop through the length of sleep Index,and pull out associated data
            % cortex
            RH_deltaPower{indexCount,1} = ProcData.sleep.parameters.cortical_RH.deltaBandPower{fixedSleepIndex(indexCount),1};
            RH_thetaPower{indexCount,1} = ProcData.sleep.parameters.cortical_RH.thetaBandPower{fixedSleepIndex(indexCount),1};
            RH_alphaPower{indexCount,1} = ProcData.sleep.parameters.cortical_RH.alphaBandPower{fixedSleepIndex(indexCount),1};
            RH_betaPower{indexCount,1} = ProcData.sleep.parameters.cortical_RH.betaBandPower{fixedSleepIndex(indexCount),1};
            RH_gammaPower{indexCount,1} = ProcData.sleep.parameters.cortical_RH.gammaBandPower{fixedSleepIndex(indexCount),1};
            RH_muaPower{indexCount,1} = ProcData.sleep.parameters.cortical_RH.muaPower{fixedSleepIndex(indexCount),1};
            % hippocampus
            Hip_deltaPower{indexCount,1} = ProcData.sleep.parameters.hippocampus.deltaBandPower{fixedSleepIndex(indexCount),1};
            Hip_thetaPower{indexCount,1} = ProcData.sleep.parameters.hippocampus.thetaBandPower{fixedSleepIndex(indexCount),1};
            Hip_alphaPower{indexCount,1} = ProcData.sleep.parameters.hippocampus.alphaBandPower{fixedSleepIndex(indexCount),1};
            Hip_betaPower{indexCount,1} = ProcData.sleep.parameters.hippocampus.betaBandPower{fixedSleepIndex(indexCount),1};
            Hip_gammaPower{indexCount,1} = ProcData.sleep.parameters.hippocampus.gammaBandPower{fixedSleepIndex(indexCount),1};
            Hip_muaPower{indexCount,1} = ProcData.sleep.parameters.hippocampus.muaPower{fixedSleepIndex(indexCount),1};
            % CBV
            LH_CBV{indexCount,1} = ProcData.sleep.parameters.CBV.LH{fixedSleepIndex(indexCount),1};
            RH_CBV{indexCount,1} = ProcData.sleep.parameters.CBV.RH{fixedSleepIndex(indexCount),1};
            fLH_CBV{indexCount,1} = ProcData.sleep.parameters.CBV.fLH{fixedSleepIndex(indexCount),1};
            fRH_CBV{indexCount,1} = ProcData.sleep.parameters.CBV.fRH{fixedSleepIndex(indexCount),1};
            % HbT
            LH_HbT{indexCount,1} = ProcData.sleep.parameters.HbT.LH{fixedSleepIndex(indexCount),1};
            RH_HbT{indexCount,1} = ProcData.sleep.parameters.HbT.RH{fixedSleepIndex(indexCount),1};
            fLH_HbT{indexCount,1} = ProcData.sleep.parameters.HbT.fLH{fixedSleepIndex(indexCount),1};
            fRH_HbT{indexCount,1} = ProcData.sleep.parameters.HbT.fRH{fixedSleepIndex(indexCount),1};
            % GCaMP
            LH_GCaMP{indexCount,1} = ProcData.sleep.parameters.GCaMP.LH{fixedSleepIndex(indexCount),1};
            RH_GCaMP{indexCount,1} = ProcData.sleep.parameters.GCaMP.RH{fixedSleepIndex(indexCount),1};
            fLH_GCaMP{indexCount,1} = ProcData.sleep.parameters.GCaMP.fLH{fixedSleepIndex(indexCount),1};
            fRH_GCaMP{indexCount,1} = ProcData.sleep.parameters.GCaMP.fRH{fixedSleepIndex(indexCount),1};
            % HbR
            LH_HbR{indexCount,1} = ProcData.sleep.parameters.HbR.LH{fixedSleepIndex(indexCount),1};
            RH_HbR{indexCount,1} = ProcData.sleep.parameters.HbR.RH{fixedSleepIndex(indexCount),1};
            fLH_HbR{indexCount,1} = ProcData.sleep.parameters.HbR.fLH{fixedSleepIndex(indexCount),1};
            fRH_HbR{indexCount,1} = ProcData.sleep.parameters.HbR.fRH{fixedSleepIndex(indexCount),1};
            % HbO
            LH_HbO{indexCount,1} = ProcData.sleep.parameters.HbO.LH{fixedSleepIndex(indexCount),1};
            RH_HbO{indexCount,1} = ProcData.sleep.parameters.HbO.RH{fixedSleepIndex(indexCount),1};
            fLH_HbO{indexCount,1} = ProcData.sleep.parameters.HbO.fLH{fixedSleepIndex(indexCount),1};
            fRH_HbO{indexCount,1} = ProcData.sleep.parameters.HbO.fRH{fixedSleepIndex(indexCount),1};
            % whiskers
            whiskerAcceleration{indexCount,1} = ProcData.sleep.parameters.whiskerAcceleration{fixedSleepIndex(indexCount),1};
            binTimes{indexCount,1} = 5*fixedSleepIndex(indexCount);
        end
        % find if there are numerous sleep periods
        indexBreaks = find(fixedSleepIndex(2:end) - fixedSleepIndex(1:end - 1) > 1);
        % if there is only one period of sleep in this file and not multiple
        if isempty(indexBreaks)
            % RH delta
            matRH_DeltaPower = cell2mat(RH_deltaPower);
            arrayRH_DeltaPower = reshape(matRH_DeltaPower',[1,size(matRH_DeltaPower,2)*size(matRH_DeltaPower,1)]);
            cellRH_DeltaPower = {arrayRH_DeltaPower};
            % RH theta
            matRH_ThetaPower = cell2mat(RH_thetaPower);
            arrayRH_ThetaPower = reshape(matRH_ThetaPower',[1,size(matRH_ThetaPower,2)*size(matRH_ThetaPower,1)]);
            cellRH_ThetaPower = {arrayRH_ThetaPower};
            % RH alpha
            matRH_AlphaPower = cell2mat(RH_alphaPower);
            arrayRH_AlphaPower = reshape(matRH_AlphaPower',[1,size(matRH_AlphaPower,2)*size(matRH_AlphaPower,1)]);
            cellRH_AlphaPower = {arrayRH_AlphaPower};
            % RH beta
            matRH_BetaPower = cell2mat(RH_betaPower);
            arrayRH_BetaPower = reshape(matRH_BetaPower',[1,size(matRH_BetaPower,2)*size(matRH_BetaPower,1)]);
            cellRH_BetaPower = {arrayRH_BetaPower};
            % RH gamma
            matRH_GammaPower = cell2mat(RH_gammaPower);
            arrayRH_GammaPower = reshape(matRH_GammaPower',[1,size(matRH_GammaPower,2)*size(matRH_GammaPower,1)]);
            cellRH_GammaPower = {arrayRH_GammaPower};
            % RH MUA
            matRH_MUAPower = cell2mat(RH_muaPower);
            arrayRH_MUAPower = reshape(matRH_MUAPower',[1,size(matRH_MUAPower,2)*size(matRH_MUAPower,1)]);
            cellRH_MUAPower = {arrayRH_MUAPower};
            % hip delta
            matHip_DeltaPower = cell2mat(Hip_deltaPower);
            arrayHip_DeltaPower = reshape(matHip_DeltaPower',[1,size(matHip_DeltaPower,2)*size(matHip_DeltaPower,1)]);
            cellHip_DeltaPower = {arrayHip_DeltaPower};
            % hip theta
            matHip_ThetaPower = cell2mat(Hip_thetaPower);
            arrayHip_ThetaPower = reshape(matHip_ThetaPower',[1,size(matHip_ThetaPower,2)*size(matHip_ThetaPower,1)]);
            cellHip_ThetaPower = {arrayHip_ThetaPower};
            % hip alpha
            matHip_AlphaPower = cell2mat(Hip_alphaPower);
            arrayHip_AlphaPower = reshape(matHip_AlphaPower',[1,size(matHip_AlphaPower,2)*size(matHip_AlphaPower,1)]);
            cellHip_AlphaPower = {arrayHip_AlphaPower};
            % hip beta
            matHip_BetaPower = cell2mat(Hip_betaPower);
            arrayHip_BetaPower = reshape(matHip_BetaPower',[1,size(matHip_BetaPower,2)*size(matHip_BetaPower,1)]);
            cellHip_BetaPower = {arrayHip_BetaPower};
            % hip gamma
            matHip_GammaPower = cell2mat(Hip_gammaPower);
            arrayHip_GammaPower = reshape(matHip_GammaPower',[1,size(matHip_GammaPower,2)*size(matHip_GammaPower,1)]);
            cellHip_GammaPower = {arrayHip_GammaPower};
            % hip MUA
            matHip_MUAPower = cell2mat(Hip_muaPower);
            arrayHip_MUAPower = reshape(matHip_MUAPower',[1,size(matHip_MUAPower,2)*size(matHip_MUAPower,1)]);
            cellHip_MUAPower = {arrayHip_MUAPower};
            % whisker acceleration
            for x = 1:length(whiskerAcceleration)
                targetPoints = size(whiskerAcceleration{1,1},2);
                if size(whiskerAcceleration{x,1},2) ~= targetPoints
                    maxLength = size(whiskerAcceleration{x,1},2);
                    difference = targetPoints - size(whiskerAcceleration{x,1},2);
                    for y = 1:difference
                        whiskerAcceleration{x,1}(maxLength + y) = 0;
                    end
                end
            end
            matWhiskerAcceleration = cell2mat(whiskerAcceleration);
            arrayWhiskerAcceleration = reshape(matWhiskerAcceleration',[1,size(matWhiskerAcceleration,2)*size(matWhiskerAcceleration,1)]);
            cellWhiskerAcceleration = {arrayWhiskerAcceleration};
            % LH CBV
            matLH_CBV = cell2mat(LH_CBV);
            arrayLH_CBV = reshape(matLH_CBV',[1,size(matLH_CBV,2)*size(matLH_CBV,1)]);
            cellLH_CBV = {arrayLH_CBV};
            % RH CBV
            matRH_CBV = cell2mat(RH_CBV);
            arrayRH_CBV = reshape(matRH_CBV',[1,size(matRH_CBV,2)*size(matRH_CBV,1)]);
            cellRH_CBV = {arrayRH_CBV};
            % f LH CBV
            matfLH_CBV = cell2mat(fLH_CBV);
            arrayfLH_CBV = reshape(matfLH_CBV',[1,size(matfLH_CBV,2)*size(matfLH_CBV,1)]);
            cellfLH_CBV = {arrayfLH_CBV};
            % f RH CBV
            matfRH_CBV = cell2mat(fRH_CBV);
            arrayfRH_CBV = reshape(matfRH_CBV',[1,size(matfRH_CBV,2)*size(matfRH_CBV,1)]);
            cellfRH_CBV = {arrayfRH_CBV};
            % LH HbT
            matLH_HbT = cell2mat(LH_HbT);
            arrayLH_HbT = reshape(matLH_HbT',[1,size(matLH_HbT,2)*size(matLH_HbT,1)]);
            cellLH_HbT = {arrayLH_HbT};
            % RH HbT
            matRH_HbT = cell2mat(RH_HbT);
            arrayRH_HbT = reshape(matRH_HbT',[1,size(matRH_HbT,2)*size(matRH_HbT,1)]);
            cellRH_HbT = {arrayRH_HbT};
            % f LH HbT
            matfLH_HbT = cell2mat(fLH_HbT);
            arrayfLH_HbT = reshape(matfLH_HbT',[1,size(matfLH_HbT,2)*size(matfLH_HbT,1)]);
            cellfLH_HbT = {arrayfLH_HbT};
            % f RH HbT
            matfRH_HbT = cell2mat(fRH_HbT);
            arrayfRH_HbT = reshape(matfRH_HbT',[1,size(matfRH_HbT,2)*size(matfRH_HbT,1)]);
            cellfRH_HbT = {arrayfRH_HbT};
            % LH GCaMP
            matLH_GCaMP = cell2mat(LH_GCaMP);
            arrayLH_GCaMP = reshape(matLH_GCaMP',[1,size(matLH_GCaMP,2)*size(matLH_GCaMP,1)]);
            cellLH_GCaMP = {arrayLH_GCaMP};
            % RH GCaMP
            matRH_GCaMP = cell2mat(RH_GCaMP);
            arrayRH_GCaMP = reshape(matRH_GCaMP',[1,size(matRH_GCaMP,2)*size(matRH_GCaMP,1)]);
            cellRH_GCaMP = {arrayRH_GCaMP};
            % f LH GCaMP
            matfLH_GCaMP = cell2mat(fLH_GCaMP);
            arrayfLH_GCaMP = reshape(matfLH_GCaMP',[1,size(matfLH_GCaMP,2)*size(matfLH_GCaMP,1)]);
            cellfLH_GCaMP = {arrayfLH_GCaMP};
            % f RH GCaMP
            matfRH_GCaMP = cell2mat(fRH_GCaMP);
            arrayfRH_GCaMP = reshape(matfRH_GCaMP',[1,size(matfRH_GCaMP,2)*size(matfRH_GCaMP,1)]);
            cellfRH_GCaMP = {arrayfRH_GCaMP};
            % LH HbR
            matLH_HbR = cell2mat(LH_HbR);
            arrayLH_HbR = reshape(matLH_HbR',[1,size(matLH_HbR,2)*size(matLH_HbR,1)]);
            cellLH_HbR = {arrayLH_HbR};
            % RH HbR
            matRH_HbR = cell2mat(RH_HbR);
            arrayRH_HbR = reshape(matRH_HbR',[1,size(matRH_HbR,2)*size(matRH_HbR,1)]);
            cellRH_HbR = {arrayRH_HbR};
            % f LH HbR
            matfLH_HbR = cell2mat(fLH_HbR);
            arrayfLH_HbR = reshape(matfLH_HbR',[1,size(matfLH_HbR,2)*size(matfLH_HbR,1)]);
            cellfLH_HbR = {arrayfLH_HbR};
            % f RH HbR
            matfRH_HbR = cell2mat(fRH_HbR);
            arrayfRH_HbR = reshape(matfRH_HbR',[1,size(matfRH_HbR,2)*size(matfRH_HbR,1)]);
            cellfRH_HbR = {arrayfRH_HbR};
            % LH HbO
            matLH_HbO = cell2mat(LH_HbO);
            arrayLH_HbO = reshape(matLH_HbO',[1,size(matLH_HbO,2)*size(matLH_HbO,1)]);
            cellLH_HbO = {arrayLH_HbO};
            % RH HbO
            matRH_HbO = cell2mat(RH_HbO);
            arrayRH_HbO = reshape(matRH_HbO',[1,size(matRH_HbO,2)*size(matRH_HbO,1)]);
            cellRH_HbO = {arrayRH_HbO};
            % f LH HbO
            matfLH_HbO = cell2mat(fLH_HbO);
            arrayfLH_HbO = reshape(matfLH_HbO',[1,size(matfLH_HbO,2)*size(matfLH_HbO,1)]);
            cellfLH_HbO = {arrayfLH_HbO};
            % f RH HbO
            matfRH_HbO = cell2mat(fRH_HbO);
            arrayfRH_HbO = reshape(matfRH_HbO',[1,size(matfRH_HbO,2)*size(matfRH_HbO,1)]);
            cellfRH_HbO = {arrayfRH_HbO};
            % bin times
            matBinTimes = cell2mat(binTimes);
            arrayBinTimes = reshape(matBinTimes',[1,size(matBinTimes,2)*size(matBinTimes,1)]);
            cellBinTimes = {arrayBinTimes};
        else
            count = length(fixedSleepIndex);
            holdIndex = zeros(1,(length(indexBreaks) + 1));
            for indexCounter = 1:length(indexBreaks) + 1
                if indexCounter == 1
                    holdIndex(indexCounter) = indexBreaks(indexCounter);
                elseif indexCounter == length(indexBreaks) + 1
                    holdIndex(indexCounter) = count - indexBreaks(indexCounter - 1);
                else
                    holdIndex(indexCounter)= indexBreaks(indexCounter) - indexBreaks(indexCounter - 1);
                end
            end
            % go through each matrix counter
            splitCounter = 1:length(RH_deltaPower);
            convertedMat2Cell = mat2cell(splitCounter',holdIndex);
            for matCounter = 1:length(convertedMat2Cell)
                % cortex
                mat2CellRH_DeltaPower{matCounter,1} = RH_deltaPower(convertedMat2Cell{matCounter,1});
                mat2CellRH_ThetaPower{matCounter,1} = RH_thetaPower(convertedMat2Cell{matCounter,1});
                mat2CellRH_AlphaPower{matCounter,1} = RH_alphaPower(convertedMat2Cell{matCounter,1});
                mat2CellRH_BetaPower{matCounter,1} = RH_betaPower(convertedMat2Cell{matCounter,1});
                mat2CellRH_GammaPower{matCounter,1} = RH_gammaPower(convertedMat2Cell{matCounter,1});
                mat2CellRH_MUAPower{matCounter,1} = RH_muaPower(convertedMat2Cell{matCounter,1});
                % hippocampus
                mat2CellHip_DeltaPower{matCounter,1} = Hip_deltaPower(convertedMat2Cell{matCounter,1});
                mat2CellHip_ThetaPower{matCounter,1} = Hip_thetaPower(convertedMat2Cell{matCounter,1});
                mat2CellHip_AlphaPower{matCounter,1} = Hip_alphaPower(convertedMat2Cell{matCounter,1});
                mat2CellHip_BetaPower{matCounter,1} = Hip_betaPower(convertedMat2Cell{matCounter,1});
                mat2CellHip_GammaPower{matCounter,1} = Hip_gammaPower(convertedMat2Cell{matCounter,1});
                mat2CellHip_MUAPower{matCounter,1} = Hip_muaPower(convertedMat2Cell{matCounter,1});
                % CBV
                mat2CellLH_CBV{matCounter,1} = LH_CBV(convertedMat2Cell{matCounter,1});
                mat2CellRH_CBV{matCounter,1} = RH_CBV(convertedMat2Cell{matCounter,1});
                mat2CellfLH_CBV{matCounter,1} = fLH_CBV(convertedMat2Cell{matCounter,1});
                mat2CellfRH_CBV{matCounter,1} = fRH_CBV(convertedMat2Cell{matCounter,1});
                % HbT
                mat2CellLH_HbT{matCounter,1} = LH_HbT(convertedMat2Cell{matCounter,1});
                mat2CellRH_HbT{matCounter,1} = RH_HbT(convertedMat2Cell{matCounter,1});
                mat2CellfLH_HbT{matCounter,1} = fLH_HbT(convertedMat2Cell{matCounter,1});
                mat2CellfRH_HbT{matCounter,1} = fRH_HbT(convertedMat2Cell{matCounter,1});
                % GCaMP
                mat2CellLH_GCaMP{matCounter,1} = LH_GCaMP(convertedMat2Cell{matCounter,1});
                mat2CellRH_GCaMP{matCounter,1} = RH_GCaMP(convertedMat2Cell{matCounter,1});
                mat2CellfLH_GCaMP{matCounter,1} = fLH_GCaMP(convertedMat2Cell{matCounter,1});
                mat2CellfRH_GCaMP{matCounter,1} = fRH_GCaMP(convertedMat2Cell{matCounter,1});
                % HbR
                mat2CellLH_HbR{matCounter,1} = LH_HbR(convertedMat2Cell{matCounter,1});
                mat2CellRH_HbR{matCounter,1} = RH_HbR(convertedMat2Cell{matCounter,1});
                mat2CellfLH_HbR{matCounter,1} = fLH_HbR(convertedMat2Cell{matCounter,1});
                mat2CellfRH_HbR{matCounter,1} = fRH_HbR(convertedMat2Cell{matCounter,1});
                % HbO
                mat2CellLH_HbO{matCounter,1} = LH_HbO(convertedMat2Cell{matCounter,1});
                mat2CellRH_HbO{matCounter,1} = RH_HbO(convertedMat2Cell{matCounter,1});
                mat2CellfLH_HbO{matCounter,1} = fLH_HbO(convertedMat2Cell{matCounter,1});
                mat2CellfRH_HbO{matCounter,1} = fRH_HbO(convertedMat2Cell{matCounter,1});
                % whiskers
                mat2CellWhiskerAcceleration{matCounter,1} = whiskerAcceleration(convertedMat2Cell{matCounter,1});
                % bin counts
                mat2CellBinTimes{matCounter,1} = binTimes(convertedMat2Cell{matCounter,1});
            end
            % go through each cell counter
            for cellCounter = 1:length(mat2CellRH_DeltaPower)
                % RH delta
                matRH_DeltaPower = cell2mat(mat2CellRH_DeltaPower{cellCounter,1});
                arrayRH_DeltaPower = reshape(matRH_DeltaPower',[1,size(matRH_DeltaPower,2)*size(matRH_DeltaPower,1)]);
                cellRH_DeltaPower{cellCounter,1} = arrayRH_DeltaPower;
                % RH theta
                matRH_ThetaPower = cell2mat(mat2CellRH_ThetaPower{cellCounter,1});
                arrayRH_ThetaPower = reshape(matRH_ThetaPower',[1,size(matRH_ThetaPower,2)*size(matRH_ThetaPower,1)]);
                cellRH_ThetaPower{cellCounter,1} = arrayRH_ThetaPower;
                % RH alpha
                matRH_AlphaPower = cell2mat(mat2CellRH_AlphaPower{cellCounter,1});
                arrayRH_AlphaPower = reshape(matRH_AlphaPower',[1,size(matRH_AlphaPower,2)*size(matRH_AlphaPower,1)]);
                cellRH_AlphaPower{cellCounter,1} = arrayRH_AlphaPower;
                % RH beta
                matRH_BetaPower = cell2mat(mat2CellRH_BetaPower{cellCounter,1});
                arrayRH_BetaPower = reshape(matRH_BetaPower',[1,size(matRH_BetaPower,2)*size(matRH_BetaPower,1)]);
                cellRH_BetaPower{cellCounter,1} = arrayRH_BetaPower;
                % RH gamma
                matRH_GammaPower = cell2mat(mat2CellRH_GammaPower{cellCounter,1});
                arrayRH_GammaPower = reshape(matRH_GammaPower',[1,size(matRH_GammaPower,2)*size(matRH_GammaPower,1)]);
                cellRH_GammaPower{cellCounter,1} = arrayRH_GammaPower;
                % RH MUA
                matRH_MUAPower = cell2mat(mat2CellRH_MUAPower{cellCounter,1});
                arrayRH_MUAPower = reshape(matRH_MUAPower',[1,size(matRH_MUAPower,2)*size(matRH_MUAPower,1)]);
                cellRH_MUAPower{cellCounter,1} = arrayRH_MUAPower;
                % hip delta
                matHip_DeltaPower = cell2mat(mat2CellHip_DeltaPower{cellCounter,1});
                arrayHip_DeltaPower = reshape(matHip_DeltaPower',[1,size(matHip_DeltaPower,2)*size(matHip_DeltaPower,1)]);
                cellHip_DeltaPower{cellCounter,1} = arrayHip_DeltaPower;
                % hip theta
                matHip_ThetaPower = cell2mat(mat2CellHip_ThetaPower{cellCounter,1});
                arrayHip_ThetaPower = reshape(matHip_ThetaPower',[1,size(matHip_ThetaPower,2)*size(matHip_ThetaPower,1)]);
                cellHip_ThetaPower{cellCounter,1} = arrayHip_ThetaPower;
                % hip alpha
                matHip_AlphaPower = cell2mat(mat2CellHip_AlphaPower{cellCounter,1});
                arrayHip_AlphaPower = reshape(matHip_AlphaPower',[1,size(matHip_AlphaPower,2)*size(matHip_AlphaPower,1)]);
                cellHip_AlphaPower{cellCounter,1} = arrayHip_AlphaPower;
                % hip beta
                matHip_BetaPower = cell2mat(mat2CellHip_BetaPower{cellCounter,1});
                arrayHip_BetaPower = reshape(matHip_BetaPower',[1,size(matHip_BetaPower,2)*size(matHip_BetaPower,1)]);
                cellHip_BetaPower{cellCounter,1} = arrayHip_BetaPower;
                % hip gamma
                matHip_GammaPower = cell2mat(mat2CellHip_GammaPower{cellCounter,1});
                arrayHip_GammaPower = reshape(matHip_GammaPower',[1,size(matHip_GammaPower,2)*size(matHip_GammaPower,1)]);
                cellHip_GammaPower{cellCounter,1} = arrayHip_GammaPower;
                % hip MUA
                matHip_MUAPower = cell2mat(mat2CellHip_MUAPower{cellCounter,1});
                arrayHip_MUAPower = reshape(matHip_MUAPower',[1,size(matHip_MUAPower,2)*size(matHip_MUAPower,1)]);
                cellHip_MUAPower{cellCounter,1} = arrayHip_MUAPower;
                % LH CBV
                matLH_CBV = cell2mat(mat2CellLH_CBV{cellCounter,1});
                arrayLH_CBV = reshape(matLH_CBV',[1,size(matLH_CBV,2)*size(matLH_CBV,1)]);
                cellLH_CBV{cellCounter,1} = arrayLH_CBV;
                % RH CBV
                matRH_CBV = cell2mat(mat2CellRH_CBV{cellCounter,1});
                arrayRH_CBV = reshape(matRH_CBV',[1,size(matRH_CBV,2)*size(matRH_CBV,1)]);
                cellRH_CBV{cellCounter,1} = arrayRH_CBV;
                % f LH CBV
                matfLH_CBV = cell2mat(mat2CellfLH_CBV{cellCounter,1});
                arrayfLH_CBV = reshape(matfLH_CBV',[1,size(matfLH_CBV,2)*size(matfLH_CBV,1)]);
                cellfLH_CBV{cellCounter,1} = arrayfLH_CBV;
                % f RH CBV
                matfRH_CBV = cell2mat(mat2CellfRH_CBV{cellCounter,1});
                arrayfRH_CBV = reshape(matfRH_CBV',[1,size(matfRH_CBV,2)*size(matfRH_CBV,1)]);
                cellfRH_CBV{cellCounter,1} = arrayfRH_CBV;
                % LH HbT
                matLH_HbT = cell2mat(mat2CellLH_HbT{cellCounter,1});
                arrayLH_HbT = reshape(matLH_HbT',[1,size(matLH_HbT,2)*size(matLH_HbT,1)]);
                cellLH_HbT{cellCounter,1} = arrayLH_HbT;
                % RH HbT
                matRH_HbT = cell2mat(mat2CellRH_HbT{cellCounter,1});
                arrayRH_HbT = reshape(matRH_HbT',[1,size(matRH_HbT,2)*size(matRH_HbT,1)]);
                cellRH_HbT{cellCounter,1} = arrayRH_HbT;
                % f LH HbT
                matfLH_HbT = cell2mat(mat2CellfLH_HbT{cellCounter,1});
                arrayfLH_HbT = reshape(matfLH_HbT',[1,size(matfLH_HbT,2)*size(matfLH_HbT,1)]);
                cellfLH_HbT{cellCounter,1} = arrayfLH_HbT;
                % f RH HbT
                matfRH_HbT = cell2mat(mat2CellfRH_HbT{cellCounter,1});
                arrayfRH_HbT = reshape(matfRH_HbT',[1,size(matfRH_HbT,2)*size(matfRH_HbT,1)]);
                cellfRH_HbT{cellCounter,1} = arrayfRH_HbT;
                % LH GCaMP
                matLH_GCaMP = cell2mat(mat2CellLH_GCaMP{cellCounter,1});
                arrayLH_GCaMP = reshape(matLH_GCaMP',[1,size(matLH_GCaMP,2)*size(matLH_GCaMP,1)]);
                cellLH_GCaMP{cellCounter,1} = arrayLH_GCaMP;
                % RH GCaMP
                matRH_GCaMP = cell2mat(mat2CellRH_GCaMP{cellCounter,1});
                arrayRH_GCaMP = reshape(matRH_GCaMP',[1,size(matRH_GCaMP,2)*size(matRH_GCaMP,1)]);
                cellRH_GCaMP{cellCounter,1} = arrayRH_GCaMP;
                % f LH GCaMP
                matfLH_GCaMP = cell2mat(mat2CellfLH_GCaMP{cellCounter,1});
                arrayfLH_GCaMP = reshape(matfLH_GCaMP',[1,size(matfLH_GCaMP,2)*size(matfLH_GCaMP,1)]);
                cellfLH_GCaMP{cellCounter,1} = arrayfLH_GCaMP;
                % f RH GCaMP
                matfRH_GCaMP = cell2mat(mat2CellfRH_GCaMP{cellCounter,1});
                arrayfRH_GCaMP = reshape(matfRH_GCaMP',[1,size(matfRH_GCaMP,2)*size(matfRH_GCaMP,1)]);
                cellfRH_GCaMP{cellCounter,1} = arrayfRH_GCaMP;
                % LH HbR
                matLH_HbR = cell2mat(mat2CellLH_HbR{cellCounter,1});
                arrayLH_HbR = reshape(matLH_HbR',[1,size(matLH_HbR,2)*size(matLH_HbR,1)]);
                cellLH_HbR{cellCounter,1} = arrayLH_HbR;
                % RH HbR
                matRH_HbR = cell2mat(mat2CellRH_HbR{cellCounter,1});
                arrayRH_HbR = reshape(matRH_HbR',[1,size(matRH_HbR,2)*size(matRH_HbR,1)]);
                cellRH_HbR{cellCounter,1} = arrayRH_HbR;
                % f LH HbR
                matfLH_HbR = cell2mat(mat2CellfLH_HbR{cellCounter,1});
                arrayfLH_HbR = reshape(matfLH_HbR',[1,size(matfLH_HbR,2)*size(matfLH_HbR,1)]);
                cellfLH_HbR{cellCounter,1} = arrayfLH_HbR;
                % f RH HbR
                matfRH_HbR = cell2mat(mat2CellfRH_HbR{cellCounter,1});
                arrayfRH_HbR = reshape(matfRH_HbR',[1,size(matfRH_HbR,2)*size(matfRH_HbR,1)]);
                cellfRH_HbR{cellCounter,1} = arrayfRH_HbR;
                % LH HbO
                matLH_HbO = cell2mat(mat2CellLH_HbO{cellCounter,1});
                arrayLH_HbO = reshape(matLH_HbO',[1,size(matLH_HbO,2)*size(matLH_HbO,1)]);
                cellLH_HbO{cellCounter,1} = arrayLH_HbO;
                % RH HbO
                matRH_HbO = cell2mat(mat2CellRH_HbO{cellCounter,1});
                arrayRH_HbO = reshape(matRH_HbO',[1,size(matRH_HbO,2)*size(matRH_HbO,1)]);
                cellRH_HbO{cellCounter,1} = arrayRH_HbO;
                % f LH HbO
                matfLH_HbO = cell2mat(mat2CellfLH_HbO{cellCounter,1});
                arrayfLH_HbO = reshape(matfLH_HbO',[1,size(matfLH_HbO,2)*size(matfLH_HbO,1)]);
                cellfLH_HbO{cellCounter,1} = arrayfLH_HbO;
                % f RH HbO
                matfRH_HbO = cell2mat(mat2CellfRH_HbO{cellCounter,1});
                arrayfRH_HbO = reshape(matfRH_HbO',[1,size(matfRH_HbO,2)*size(matfRH_HbO,1)]);
                cellfRH_HbO{cellCounter,1} = arrayfRH_HbO;
                % whisker acceleration
                for x = 1:size(mat2CellWhiskerAcceleration{cellCounter,1},1)
                    targetPoints = size(mat2CellWhiskerAcceleration{cellCounter,1}{1,1},2);
                    if size(mat2CellWhiskerAcceleration{cellCounter,1}{x,1},2) ~= targetPoints
                        maxLength = size(mat2CellWhiskerAcceleration{cellCounter,1}{x,1},2);
                        difference = targetPoints - size(mat2CellWhiskerAcceleration{cellCounter,1}{x,1},2);
                        for y = 1:difference
                            mat2CellWhiskerAcceleration{cellCounter,1}{x,1}(maxLength + y) = 0;
                        end
                    end
                end
                matWhiskerAcceleration = cell2mat(mat2CellWhiskerAcceleration{cellCounter,1});
                arrayWhiskerAcceleration = reshape(matWhiskerAcceleration',[1,size(matWhiskerAcceleration,2)*size(matWhiskerAcceleration,1)]);
                cellWhiskerAcceleration{cellCounter,1} = arrayWhiskerAcceleration;
                % bin times
                matBinTimes = cell2mat(mat2CellBinTimes{cellCounter,1});
                arrayBinTimes = reshape(matBinTimes',[1,size(matBinTimes,2)*size(matBinTimes,1)]);
                cellBinTimes{cellCounter,1} = arrayBinTimes;
            end
        end
        % save the data in the SleepEventData struct
        if isfield(SleepData,(modelName)) == false % if the structure is empty we need a special case to format the struct properly
            for cellLength = 1:size(cellRH_DeltaPower,2) % loop through however many sleep epochs this file has
                % cortex
                SleepData.(modelName).NREM.data.cortical_RH.deltaBandPower{cellLength,1} = cellRH_DeltaPower{1,1};
                SleepData.(modelName).NREM.data.cortical_RH.thetaBandPower{cellLength,1} = cellRH_ThetaPower{1,1};
                SleepData.(modelName).NREM.data.cortical_RH.alphaBandPower{cellLength,1} = cellRH_AlphaPower{1,1};
                SleepData.(modelName).NREM.data.cortical_RH.betaBandPower{cellLength,1} = cellRH_BetaPower{1,1};
                SleepData.(modelName).NREM.data.cortical_RH.gammaBandPower{cellLength,1} = cellRH_GammaPower{1,1};
                SleepData.(modelName).NREM.data.cortical_RH.muaPower{cellLength,1} = cellRH_MUAPower{1,1};
                % hippocampus
                SleepData.(modelName).NREM.data.hippocampus.deltaBandPower{cellLength,1} = cellHip_DeltaPower{1,1};
                SleepData.(modelName).NREM.data.hippocampus.thetaBandPower{cellLength,1} = cellHip_ThetaPower{1,1};
                SleepData.(modelName).NREM.data.hippocampus.alphaBandPower{cellLength,1} = cellHip_AlphaPower{1,1};
                SleepData.(modelName).NREM.data.hippocampus.betaBandPower{cellLength,1} = cellHip_BetaPower{1,1};
                SleepData.(modelName).NREM.data.hippocampus.gammaBandPower{cellLength,1} = cellHip_GammaPower{1,1};
                SleepData.(modelName).NREM.data.hippocampus.muaPower{cellLength,1} = cellHip_MUAPower{1,1};
                % CBV
                SleepData.(modelName).NREM.data.CBV.LH{cellLength,1} = cellLH_CBV{1,1};
                SleepData.(modelName).NREM.data.CBV.RH{cellLength,1} = cellRH_CBV{1,1};
                SleepData.(modelName).NREM.data.CBV.fLH{cellLength,1} = cellfLH_CBV{1,1};
                SleepData.(modelName).NREM.data.CBV.fRH{cellLength,1} = cellfRH_CBV{1,1};
                % HbT
                SleepData.(modelName).NREM.data.HbT.LH{cellLength,1} = cellLH_HbT{1,1};
                SleepData.(modelName).NREM.data.HbT.RH{cellLength,1} = cellRH_HbT{1,1};
                SleepData.(modelName).NREM.data.HbT.fLH{cellLength,1} = cellfLH_HbT{1,1};
                SleepData.(modelName).NREM.data.HbT.fRH{cellLength,1} = cellfRH_HbT{1,1};
                % GCaMP
                SleepData.(modelName).NREM.data.GCaMP.LH{cellLength,1} = cellLH_GCaMP{1,1};
                SleepData.(modelName).NREM.data.GCaMP.RH{cellLength,1} = cellRH_GCaMP{1,1};
                SleepData.(modelName).NREM.data.GCaMP.fLH{cellLength,1} = cellfLH_GCaMP{1,1};
                SleepData.(modelName).NREM.data.GCaMP.fRH{cellLength,1} = cellfRH_GCaMP{1,1};
                % HbR
                SleepData.(modelName).NREM.data.HbR.LH{cellLength,1} = cellLH_HbR{1,1};
                SleepData.(modelName).NREM.data.HbR.RH{cellLength,1} = cellRH_HbR{1,1};
                SleepData.(modelName).NREM.data.HbR.fLH{cellLength,1} = cellfLH_HbR{1,1};
                SleepData.(modelName).NREM.data.HbR.fRH{cellLength,1} = cellfRH_HbR{1,1};
                % HbO
                SleepData.(modelName).NREM.data.HbO.LH{cellLength,1} = cellLH_HbO{1,1};
                SleepData.(modelName).NREM.data.HbO.RH{cellLength,1} = cellRH_HbO{1,1};
                SleepData.(modelName).NREM.data.HbO.fLH{cellLength,1} = cellfLH_HbO{1,1};
                SleepData.(modelName).NREM.data.HbO.fRH{cellLength,1} = cellfRH_HbO{1,1};
                % whiskers
                SleepData.(modelName).NREM.data.whiskerAcceleration{cellLength,1} = cellWhiskerAcceleration{1,1};
                % file IDs & bin times
                SleepData.(modelName).NREM.FileIDs{cellLength,1} = fileID;
                SleepData.(modelName).NREM.BinTimes{cellLength,1} = cellBinTimes{1,1};
            end
        else % if the struct is not empty,add each new iteration after previous data
            for cellLength = 1:size(cellRH_DeltaPower,1) % loop through however many sleep epochs this file has
                % cortex
                SleepData.(modelName).NREM.data.cortical_RH.deltaBandPower{size(SleepData.(modelName).NREM.data.cortical_RH.deltaBandPower,1) + 1,1} = cellRH_DeltaPower{cellLength,1};
                SleepData.(modelName).NREM.data.cortical_RH.thetaBandPower{size(SleepData.(modelName).NREM.data.cortical_RH.thetaBandPower,1) + 1,1} = cellRH_ThetaPower{cellLength,1};
                SleepData.(modelName).NREM.data.cortical_RH.alphaBandPower{size(SleepData.(modelName).NREM.data.cortical_RH.alphaBandPower,1) + 1,1} = cellRH_AlphaPower{cellLength,1};
                SleepData.(modelName).NREM.data.cortical_RH.betaBandPower{size(SleepData.(modelName).NREM.data.cortical_RH.betaBandPower,1) + 1,1} = cellRH_BetaPower{cellLength,1};
                SleepData.(modelName).NREM.data.cortical_RH.gammaBandPower{size(SleepData.(modelName).NREM.data.cortical_RH.gammaBandPower,1) + 1,1} = cellRH_GammaPower{cellLength,1};
                SleepData.(modelName).NREM.data.cortical_RH.muaPower{size(SleepData.(modelName).NREM.data.cortical_RH.muaPower,1) + 1,1} = cellRH_MUAPower{cellLength,1};
                % hippocampus
                SleepData.(modelName).NREM.data.hippocampus.deltaBandPower{size(SleepData.(modelName).NREM.data.hippocampus.deltaBandPower,1) + 1,1} = cellHip_DeltaPower{cellLength,1};
                SleepData.(modelName).NREM.data.hippocampus.thetaBandPower{size(SleepData.(modelName).NREM.data.hippocampus.thetaBandPower,1) + 1,1} = cellHip_ThetaPower{cellLength,1};
                SleepData.(modelName).NREM.data.hippocampus.alphaBandPower{size(SleepData.(modelName).NREM.data.hippocampus.alphaBandPower,1) + 1,1} = cellHip_AlphaPower{cellLength,1};
                SleepData.(modelName).NREM.data.hippocampus.betaBandPower{size(SleepData.(modelName).NREM.data.hippocampus.betaBandPower,1) + 1,1} = cellHip_BetaPower{cellLength,1};
                SleepData.(modelName).NREM.data.hippocampus.gammaBandPower{size(SleepData.(modelName).NREM.data.hippocampus.gammaBandPower,1) + 1,1} = cellHip_GammaPower{cellLength,1};
                SleepData.(modelName).NREM.data.hippocampus.muaPower{size(SleepData.(modelName).NREM.data.hippocampus.muaPower,1) + 1,1} = cellHip_MUAPower{cellLength,1};
                % CBV
                SleepData.(modelName).NREM.data.CBV.LH{size(SleepData.(modelName).NREM.data.CBV.LH,1) + 1,1} = cellLH_CBV{cellLength,1};
                SleepData.(modelName).NREM.data.CBV.RH{size(SleepData.(modelName).NREM.data.CBV.RH,1) + 1,1} = cellRH_CBV{cellLength,1};
                SleepData.(modelName).NREM.data.CBV.fLH{size(SleepData.(modelName).NREM.data.CBV.fLH,1) + 1,1} = cellfLH_CBV{cellLength,1};
                SleepData.(modelName).NREM.data.CBV.fRH{size(SleepData.(modelName).NREM.data.CBV.fRH,1) + 1,1} = cellfRH_CBV{cellLength,1};
                % HbT
                SleepData.(modelName).NREM.data.HbT.LH{size(SleepData.(modelName).NREM.data.HbT.LH,1) + 1,1} = cellLH_HbT{cellLength,1};
                SleepData.(modelName).NREM.data.HbT.RH{size(SleepData.(modelName).NREM.data.HbT.RH,1) + 1,1} = cellRH_HbT{cellLength,1};
                SleepData.(modelName).NREM.data.HbT.fLH{size(SleepData.(modelName).NREM.data.HbT.fLH,1) + 1,1} = cellfLH_HbT{cellLength,1};
                SleepData.(modelName).NREM.data.HbT.fRH{size(SleepData.(modelName).NREM.data.HbT.fRH,1) + 1,1} = cellfRH_HbT{cellLength,1};
                % GCaMP
                SleepData.(modelName).NREM.data.GCaMP.LH{size(SleepData.(modelName).NREM.data.GCaMP.LH,1) + 1,1} = cellLH_GCaMP{cellLength,1};
                SleepData.(modelName).NREM.data.GCaMP.RH{size(SleepData.(modelName).NREM.data.GCaMP.RH,1) + 1,1} = cellRH_GCaMP{cellLength,1};
                SleepData.(modelName).NREM.data.GCaMP.fLH{size(SleepData.(modelName).NREM.data.GCaMP.fLH,1) + 1,1} = cellfLH_GCaMP{cellLength,1};
                SleepData.(modelName).NREM.data.GCaMP.fRH{size(SleepData.(modelName).NREM.data.GCaMP.fRH,1) + 1,1} = cellfRH_GCaMP{cellLength,1};
                % HbR
                SleepData.(modelName).NREM.data.HbR.LH{size(SleepData.(modelName).NREM.data.HbR.LH,1) + 1,1} = cellLH_HbR{cellLength,1};
                SleepData.(modelName).NREM.data.HbR.RH{size(SleepData.(modelName).NREM.data.HbR.RH,1) + 1,1} = cellRH_HbR{cellLength,1};
                SleepData.(modelName).NREM.data.HbR.fLH{size(SleepData.(modelName).NREM.data.HbR.fLH,1) + 1,1} = cellfLH_HbR{cellLength,1};
                SleepData.(modelName).NREM.data.HbR.fRH{size(SleepData.(modelName).NREM.data.HbR.fRH,1) + 1,1} = cellfRH_HbR{cellLength,1};
                % HbO
                SleepData.(modelName).NREM.data.HbO.LH{size(SleepData.(modelName).NREM.data.HbO.LH,1) + 1,1} = cellLH_HbO{cellLength,1};
                SleepData.(modelName).NREM.data.HbO.RH{size(SleepData.(modelName).NREM.data.HbO.RH,1) + 1,1} = cellRH_HbO{cellLength,1};
                SleepData.(modelName).NREM.data.HbO.fLH{size(SleepData.(modelName).NREM.data.HbO.fLH,1) + 1,1} = cellfLH_HbO{cellLength,1};
                SleepData.(modelName).NREM.data.HbO.fRH{size(SleepData.(modelName).NREM.data.HbO.fRH,1) + 1,1} = cellfRH_HbO{cellLength,1};
                % whiskers
                SleepData.(modelName).NREM.data.whiskerAcceleration{size(SleepData.(modelName).NREM.data.whiskerAcceleration,1) + 1,1} = cellWhiskerAcceleration{cellLength,1};
                % file IDs & bin times
                SleepData.(modelName).NREM.FileIDs{size(SleepData.(modelName).NREM.FileIDs,1) + 1,1} = fileID;
                SleepData.(modelName).NREM.BinTimes{size(SleepData.(modelName).NREM.BinTimes,1) + 1,1} = cellBinTimes{cellLength,1};
            end
        end
    end
    disp(['Adding NREM sleeping epochs from ProcData file ' num2str(aa) ' of ' num2str(size(procDataFileIDs,1)) '...']); disp(' ')
end
% identify REM sleep epochs and place in SleepEventData.mat structure
sleepBins = REMsleepTime/5;
for aa = 1:size(procDataFileIDs,1) % loop through the list of ProcData files
    clearvars -except aa procDataFileIDs sleepBins NREMsleepTime REMsleepTime modelName SleepData startingDirectory
    procDataFileID = procDataFileIDs(aa,:); % pull character string associated with the current file
    load(procDataFileID); % load in procDataFile associated with character string
    [~,~,fileID] = GetFileInfo_IOS(procDataFileID); % gather file info
    remLogical = ProcData.sleep.logicals.(modelName).remLogical; % logical - ones denote potential sleep epoches (5 second bins)
    targetTime = ones(1,sleepBins); % target time
    sleepIndex = find(conv(remLogical,targetTime) >= sleepBins) - (sleepBins - 1); % find the periods of time where there are at least 11 more
    % 5 second epochs following. This is not the full list.
    if isempty(sleepIndex) % if sleepIndex is empty,skip this file
        % skip file
    else
        sleepCriteria = (0:(sleepBins - 1)); % this will be used to fix the issue in sleepIndex
        fixedSleepIndex = unique(sleepIndex + sleepCriteria); % sleep Index now has the proper time stamps from sleep logical
        for indexCount = 1:length(fixedSleepIndex) % loop through the length of sleep Index,and pull out associated data
            % cortex
            RH_deltaPower{indexCount,1} = ProcData.sleep.parameters.cortical_RH.deltaBandPower{fixedSleepIndex(indexCount),1};
            RH_thetaPower{indexCount,1} = ProcData.sleep.parameters.cortical_RH.thetaBandPower{fixedSleepIndex(indexCount),1};
            RH_alphaPower{indexCount,1} = ProcData.sleep.parameters.cortical_RH.alphaBandPower{fixedSleepIndex(indexCount),1};
            RH_betaPower{indexCount,1} = ProcData.sleep.parameters.cortical_RH.betaBandPower{fixedSleepIndex(indexCount),1};
            RH_gammaPower{indexCount,1} = ProcData.sleep.parameters.cortical_RH.gammaBandPower{fixedSleepIndex(indexCount),1};
            RH_muaPower{indexCount,1} = ProcData.sleep.parameters.cortical_RH.muaPower{fixedSleepIndex(indexCount),1};
            % hippocampus
            Hip_deltaPower{indexCount,1} = ProcData.sleep.parameters.hippocampus.deltaBandPower{fixedSleepIndex(indexCount),1};
            Hip_thetaPower{indexCount,1} = ProcData.sleep.parameters.hippocampus.thetaBandPower{fixedSleepIndex(indexCount),1};
            Hip_alphaPower{indexCount,1} = ProcData.sleep.parameters.hippocampus.alphaBandPower{fixedSleepIndex(indexCount),1};
            Hip_betaPower{indexCount,1} = ProcData.sleep.parameters.hippocampus.betaBandPower{fixedSleepIndex(indexCount),1};
            Hip_gammaPower{indexCount,1} = ProcData.sleep.parameters.hippocampus.gammaBandPower{fixedSleepIndex(indexCount),1};
            Hip_muaPower{indexCount,1} = ProcData.sleep.parameters.hippocampus.muaPower{fixedSleepIndex(indexCount),1};
            % CBV
            LH_CBV{indexCount,1} = ProcData.sleep.parameters.CBV.LH{fixedSleepIndex(indexCount),1};
            RH_CBV{indexCount,1} = ProcData.sleep.parameters.CBV.RH{fixedSleepIndex(indexCount),1};
            fLH_CBV{indexCount,1} = ProcData.sleep.parameters.CBV.fLH{fixedSleepIndex(indexCount),1};
            fRH_CBV{indexCount,1} = ProcData.sleep.parameters.CBV.fRH{fixedSleepIndex(indexCount),1};
            % HbT
            LH_HbT{indexCount,1} = ProcData.sleep.parameters.HbT.LH{fixedSleepIndex(indexCount),1};
            RH_HbT{indexCount,1} = ProcData.sleep.parameters.HbT.RH{fixedSleepIndex(indexCount),1};
            fLH_HbT{indexCount,1} = ProcData.sleep.parameters.HbT.fLH{fixedSleepIndex(indexCount),1};
            fRH_HbT{indexCount,1} = ProcData.sleep.parameters.HbT.fRH{fixedSleepIndex(indexCount),1};
            % GCaMP
            LH_GCaMP{indexCount,1} = ProcData.sleep.parameters.GCaMP.LH{fixedSleepIndex(indexCount),1};
            RH_GCaMP{indexCount,1} = ProcData.sleep.parameters.GCaMP.RH{fixedSleepIndex(indexCount),1};
            fLH_GCaMP{indexCount,1} = ProcData.sleep.parameters.GCaMP.fLH{fixedSleepIndex(indexCount),1};
            fRH_GCaMP{indexCount,1} = ProcData.sleep.parameters.GCaMP.fRH{fixedSleepIndex(indexCount),1};
            % HbR
            LH_HbR{indexCount,1} = ProcData.sleep.parameters.HbR.LH{fixedSleepIndex(indexCount),1};
            RH_HbR{indexCount,1} = ProcData.sleep.parameters.HbR.RH{fixedSleepIndex(indexCount),1};
            fLH_HbR{indexCount,1} = ProcData.sleep.parameters.HbR.fLH{fixedSleepIndex(indexCount),1};
            fRH_HbR{indexCount,1} = ProcData.sleep.parameters.HbR.fRH{fixedSleepIndex(indexCount),1};
            % HbO
            LH_HbO{indexCount,1} = ProcData.sleep.parameters.HbO.LH{fixedSleepIndex(indexCount),1};
            RH_HbO{indexCount,1} = ProcData.sleep.parameters.HbO.RH{fixedSleepIndex(indexCount),1};
            fLH_HbO{indexCount,1} = ProcData.sleep.parameters.HbO.fLH{fixedSleepIndex(indexCount),1};
            fRH_HbO{indexCount,1} = ProcData.sleep.parameters.HbO.fRH{fixedSleepIndex(indexCount),1};
            % whiskers
            whiskerAcceleration{indexCount,1} = ProcData.sleep.parameters.whiskerAcceleration{fixedSleepIndex(indexCount),1};
            binTimes{indexCount,1} = 5*fixedSleepIndex(indexCount);
        end
        % find if there are numerous sleep periods
        indexBreaks = find(fixedSleepIndex(2:end) - fixedSleepIndex(1:end - 1) > 1);
        % if there is only one period of sleep in this file and not multiple
        if isempty(indexBreaks)
            % RH delta
            matRH_DeltaPower = cell2mat(RH_deltaPower);
            arrayRH_DeltaPower = reshape(matRH_DeltaPower',[1,size(matRH_DeltaPower,2)*size(matRH_DeltaPower,1)]);
            cellRH_DeltaPower = {arrayRH_DeltaPower};
            % RH theta
            matRH_ThetaPower = cell2mat(RH_thetaPower);
            arrayRH_ThetaPower = reshape(matRH_ThetaPower',[1,size(matRH_ThetaPower,2)*size(matRH_ThetaPower,1)]);
            cellRH_ThetaPower = {arrayRH_ThetaPower};
            % RH alpha
            matRH_AlphaPower = cell2mat(RH_alphaPower);
            arrayRH_AlphaPower = reshape(matRH_AlphaPower',[1,size(matRH_AlphaPower,2)*size(matRH_AlphaPower,1)]);
            cellRH_AlphaPower = {arrayRH_AlphaPower};
            % RH beta
            matRH_BetaPower = cell2mat(RH_betaPower);
            arrayRH_BetaPower = reshape(matRH_BetaPower',[1,size(matRH_BetaPower,2)*size(matRH_BetaPower,1)]);
            cellRH_BetaPower = {arrayRH_BetaPower};
            % RH gamma
            matRH_GammaPower = cell2mat(RH_gammaPower);
            arrayRH_GammaPower = reshape(matRH_GammaPower',[1,size(matRH_GammaPower,2)*size(matRH_GammaPower,1)]);
            cellRH_GammaPower = {arrayRH_GammaPower};
            % RH MUA
            matRH_MUAPower = cell2mat(RH_muaPower);
            arrayRH_MUAPower = reshape(matRH_MUAPower',[1,size(matRH_MUAPower,2)*size(matRH_MUAPower,1)]);
            cellRH_MUAPower = {arrayRH_MUAPower};
            % hip delta
            matHip_DeltaPower = cell2mat(Hip_deltaPower);
            arrayHip_DeltaPower = reshape(matHip_DeltaPower',[1,size(matHip_DeltaPower,2)*size(matHip_DeltaPower,1)]);
            cellHip_DeltaPower = {arrayHip_DeltaPower};
            % hip theta
            matHip_ThetaPower = cell2mat(Hip_thetaPower);
            arrayHip_ThetaPower = reshape(matHip_ThetaPower',[1,size(matHip_ThetaPower,2)*size(matHip_ThetaPower,1)]);
            cellHip_ThetaPower = {arrayHip_ThetaPower};
            % hip alpha
            matHip_AlphaPower = cell2mat(Hip_alphaPower);
            arrayHip_AlphaPower = reshape(matHip_AlphaPower',[1,size(matHip_AlphaPower,2)*size(matHip_AlphaPower,1)]);
            cellHip_AlphaPower = {arrayHip_AlphaPower};
            % hip beta
            matHip_BetaPower = cell2mat(Hip_betaPower);
            arrayHip_BetaPower = reshape(matHip_BetaPower',[1,size(matHip_BetaPower,2)*size(matHip_BetaPower,1)]);
            cellHip_BetaPower = {arrayHip_BetaPower};
            % hip gamma
            matHip_GammaPower = cell2mat(Hip_gammaPower);
            arrayHip_GammaPower = reshape(matHip_GammaPower',[1,size(matHip_GammaPower,2)*size(matHip_GammaPower,1)]);
            cellHip_GammaPower = {arrayHip_GammaPower};
            % hip MUA
            matHip_MUAPower = cell2mat(Hip_muaPower);
            arrayHip_MUAPower = reshape(matHip_MUAPower',[1,size(matHip_MUAPower,2)*size(matHip_MUAPower,1)]);
            cellHip_MUAPower = {arrayHip_MUAPower};
            % whisker acceleration
            for x = 1:length(whiskerAcceleration)
                targetPoints = size(whiskerAcceleration{1,1},2);
                if size(whiskerAcceleration{x,1},2) ~= targetPoints
                    maxLength = size(whiskerAcceleration{x,1},2);
                    difference = targetPoints - size(whiskerAcceleration{x,1},2);
                    for y = 1:difference
                        whiskerAcceleration{x,1}(maxLength + y) = 0;
                    end
                end
            end
            matWhiskerAcceleration = cell2mat(whiskerAcceleration);
            arrayWhiskerAcceleration = reshape(matWhiskerAcceleration',[1,size(matWhiskerAcceleration,2)*size(matWhiskerAcceleration,1)]);
            cellWhiskerAcceleration = {arrayWhiskerAcceleration};
            % LH CBV
            matLH_CBV = cell2mat(LH_CBV);
            arrayLH_CBV = reshape(matLH_CBV',[1,size(matLH_CBV,2)*size(matLH_CBV,1)]);
            cellLH_CBV = {arrayLH_CBV};
            % RH CBV
            matRH_CBV = cell2mat(RH_CBV);
            arrayRH_CBV = reshape(matRH_CBV',[1,size(matRH_CBV,2)*size(matRH_CBV,1)]);
            cellRH_CBV = {arrayRH_CBV};
            % f LH CBV
            matfLH_CBV = cell2mat(fLH_CBV);
            arrayfLH_CBV = reshape(matfLH_CBV',[1,size(matfLH_CBV,2)*size(matfLH_CBV,1)]);
            cellfLH_CBV = {arrayfLH_CBV};
            % f RH CBV
            matfRH_CBV = cell2mat(fRH_CBV);
            arrayfRH_CBV = reshape(matfRH_CBV',[1,size(matfRH_CBV,2)*size(matfRH_CBV,1)]);
            cellfRH_CBV = {arrayfRH_CBV};
            % LH HbT
            matLH_HbT = cell2mat(LH_HbT);
            arrayLH_HbT = reshape(matLH_HbT',[1,size(matLH_HbT,2)*size(matLH_HbT,1)]);
            cellLH_HbT = {arrayLH_HbT};
            % RH HbT
            matRH_HbT = cell2mat(RH_HbT);
            arrayRH_HbT = reshape(matRH_HbT',[1,size(matRH_HbT,2)*size(matRH_HbT,1)]);
            cellRH_HbT = {arrayRH_HbT};
            % f LH HbT
            matfLH_HbT = cell2mat(fLH_HbT);
            arrayfLH_HbT = reshape(matfLH_HbT',[1,size(matfLH_HbT,2)*size(matfLH_HbT,1)]);
            cellfLH_HbT = {arrayfLH_HbT};
            % f RH HbT
            matfRH_HbT = cell2mat(fRH_HbT);
            arrayfRH_HbT = reshape(matfRH_HbT',[1,size(matfRH_HbT,2)*size(matfRH_HbT,1)]);
            cellfRH_HbT = {arrayfRH_HbT};
            % LH GCaMP
            matLH_GCaMP = cell2mat(LH_GCaMP);
            arrayLH_GCaMP = reshape(matLH_GCaMP',[1,size(matLH_GCaMP,2)*size(matLH_GCaMP,1)]);
            cellLH_GCaMP = {arrayLH_GCaMP};
            % RH GCaMP
            matRH_GCaMP = cell2mat(RH_GCaMP);
            arrayRH_GCaMP = reshape(matRH_GCaMP',[1,size(matRH_GCaMP,2)*size(matRH_GCaMP,1)]);
            cellRH_GCaMP = {arrayRH_GCaMP};
            % f LH GCaMP
            matfLH_GCaMP = cell2mat(fLH_GCaMP);
            arrayfLH_GCaMP = reshape(matfLH_GCaMP',[1,size(matfLH_GCaMP,2)*size(matfLH_GCaMP,1)]);
            cellfLH_GCaMP = {arrayfLH_GCaMP};
            % f RH GCaMP
            matfRH_GCaMP = cell2mat(fRH_GCaMP);
            arrayfRH_GCaMP = reshape(matfRH_GCaMP',[1,size(matfRH_GCaMP,2)*size(matfRH_GCaMP,1)]);
            cellfRH_GCaMP = {arrayfRH_GCaMP};
            % LH HbR
            matLH_HbR = cell2mat(LH_HbR);
            arrayLH_HbR = reshape(matLH_HbR',[1,size(matLH_HbR,2)*size(matLH_HbR,1)]);
            cellLH_HbR = {arrayLH_HbR};
            % RH HbR
            matRH_HbR = cell2mat(RH_HbR);
            arrayRH_HbR = reshape(matRH_HbR',[1,size(matRH_HbR,2)*size(matRH_HbR,1)]);
            cellRH_HbR = {arrayRH_HbR};
            % f LH HbR
            matfLH_HbR = cell2mat(fLH_HbR);
            arrayfLH_HbR = reshape(matfLH_HbR',[1,size(matfLH_HbR,2)*size(matfLH_HbR,1)]);
            cellfLH_HbR = {arrayfLH_HbR};
            % f RH HbR
            matfRH_HbR = cell2mat(fRH_HbR);
            arrayfRH_HbR = reshape(matfRH_HbR',[1,size(matfRH_HbR,2)*size(matfRH_HbR,1)]);
            cellfRH_HbR = {arrayfRH_HbR};
            % LH HbO
            matLH_HbO = cell2mat(LH_HbO);
            arrayLH_HbO = reshape(matLH_HbO',[1,size(matLH_HbO,2)*size(matLH_HbO,1)]);
            cellLH_HbO = {arrayLH_HbO};
            % RH HbO
            matRH_HbO = cell2mat(RH_HbO);
            arrayRH_HbO = reshape(matRH_HbO',[1,size(matRH_HbO,2)*size(matRH_HbO,1)]);
            cellRH_HbO = {arrayRH_HbO};
            % f LH HbO
            matfLH_HbO = cell2mat(fLH_HbO);
            arrayfLH_HbO = reshape(matfLH_HbO',[1,size(matfLH_HbO,2)*size(matfLH_HbO,1)]);
            cellfLH_HbO = {arrayfLH_HbO};
            % f RH HbO
            matfRH_HbO = cell2mat(fRH_HbO);
            arrayfRH_HbO = reshape(matfRH_HbO',[1,size(matfRH_HbO,2)*size(matfRH_HbO,1)]);
            cellfRH_HbO = {arrayfRH_HbO};
            % bin times
            matBinTimes = cell2mat(binTimes);
            arrayBinTimes = reshape(matBinTimes',[1,size(matBinTimes,2)*size(matBinTimes,1)]);
            cellBinTimes = {arrayBinTimes};
        else
            count = length(fixedSleepIndex);
            holdIndex = zeros(1,(length(indexBreaks) + 1));
            for indexCounter = 1:length(indexBreaks) + 1
                if indexCounter == 1
                    holdIndex(indexCounter) = indexBreaks(indexCounter);
                elseif indexCounter == length(indexBreaks) + 1
                    holdIndex(indexCounter) = count - indexBreaks(indexCounter - 1);
                else
                    holdIndex(indexCounter)= indexBreaks(indexCounter) - indexBreaks(indexCounter - 1);
                end
            end
            % go through each matrix counter
            splitCounter = 1:length(RH_deltaPower);
            convertedMat2Cell = mat2cell(splitCounter',holdIndex);
            for matCounter = 1:length(convertedMat2Cell)
                % cortex
                mat2CellRH_DeltaPower{matCounter,1} = RH_deltaPower(convertedMat2Cell{matCounter,1});
                mat2CellRH_ThetaPower{matCounter,1} = RH_thetaPower(convertedMat2Cell{matCounter,1});
                mat2CellRH_AlphaPower{matCounter,1} = RH_alphaPower(convertedMat2Cell{matCounter,1});
                mat2CellRH_BetaPower{matCounter,1} = RH_betaPower(convertedMat2Cell{matCounter,1});
                mat2CellRH_GammaPower{matCounter,1} = RH_gammaPower(convertedMat2Cell{matCounter,1});
                mat2CellRH_MUAPower{matCounter,1} = RH_muaPower(convertedMat2Cell{matCounter,1});
                % hippocampus
                mat2CellHip_DeltaPower{matCounter,1} = Hip_deltaPower(convertedMat2Cell{matCounter,1});
                mat2CellHip_ThetaPower{matCounter,1} = Hip_thetaPower(convertedMat2Cell{matCounter,1});
                mat2CellHip_AlphaPower{matCounter,1} = Hip_alphaPower(convertedMat2Cell{matCounter,1});
                mat2CellHip_BetaPower{matCounter,1} = Hip_betaPower(convertedMat2Cell{matCounter,1});
                mat2CellHip_GammaPower{matCounter,1} = Hip_gammaPower(convertedMat2Cell{matCounter,1});
                mat2CellHip_MUAPower{matCounter,1} = Hip_muaPower(convertedMat2Cell{matCounter,1});
                % CBV
                mat2CellLH_CBV{matCounter,1} = LH_CBV(convertedMat2Cell{matCounter,1});
                mat2CellRH_CBV{matCounter,1} = RH_CBV(convertedMat2Cell{matCounter,1});
                mat2CellfLH_CBV{matCounter,1} = fLH_CBV(convertedMat2Cell{matCounter,1});
                mat2CellfRH_CBV{matCounter,1} = fRH_CBV(convertedMat2Cell{matCounter,1});
                % HbT
                mat2CellLH_HbT{matCounter,1} = LH_HbT(convertedMat2Cell{matCounter,1});
                mat2CellRH_HbT{matCounter,1} = RH_HbT(convertedMat2Cell{matCounter,1});
                mat2CellfLH_HbT{matCounter,1} = fLH_HbT(convertedMat2Cell{matCounter,1});
                mat2CellfRH_HbT{matCounter,1} = fRH_HbT(convertedMat2Cell{matCounter,1});
                % GCaMP
                mat2CellLH_GCaMP{matCounter,1} = LH_GCaMP(convertedMat2Cell{matCounter,1});
                mat2CellRH_GCaMP{matCounter,1} = RH_GCaMP(convertedMat2Cell{matCounter,1});
                mat2CellfLH_GCaMP{matCounter,1} = fLH_GCaMP(convertedMat2Cell{matCounter,1});
                mat2CellfRH_GCaMP{matCounter,1} = fRH_GCaMP(convertedMat2Cell{matCounter,1});
                % HbR
                mat2CellLH_HbR{matCounter,1} = LH_HbR(convertedMat2Cell{matCounter,1});
                mat2CellRH_HbR{matCounter,1} = RH_HbR(convertedMat2Cell{matCounter,1});
                mat2CellfLH_HbR{matCounter,1} = fLH_HbR(convertedMat2Cell{matCounter,1});
                mat2CellfRH_HbR{matCounter,1} = fRH_HbR(convertedMat2Cell{matCounter,1});
                % HbO
                mat2CellLH_HbO{matCounter,1} = LH_HbO(convertedMat2Cell{matCounter,1});
                mat2CellRH_HbO{matCounter,1} = RH_HbO(convertedMat2Cell{matCounter,1});
                mat2CellfLH_HbO{matCounter,1} = fLH_HbO(convertedMat2Cell{matCounter,1});
                mat2CellfRH_HbO{matCounter,1} = fRH_HbO(convertedMat2Cell{matCounter,1});
                % whiskers
                mat2CellWhiskerAcceleration{matCounter,1} = whiskerAcceleration(convertedMat2Cell{matCounter,1});
                % bin counts
                mat2CellBinTimes{matCounter,1} = binTimes(convertedMat2Cell{matCounter,1});
            end
            % go through each cell counter
            for cellCounter = 1:length(mat2CellRH_DeltaPower)
                % RH delta
                matRH_DeltaPower = cell2mat(mat2CellRH_DeltaPower{cellCounter,1});
                arrayRH_DeltaPower = reshape(matRH_DeltaPower',[1,size(matRH_DeltaPower,2)*size(matRH_DeltaPower,1)]);
                cellRH_DeltaPower{cellCounter,1} = arrayRH_DeltaPower;
                % RH theta
                matRH_ThetaPower = cell2mat(mat2CellRH_ThetaPower{cellCounter,1});
                arrayRH_ThetaPower = reshape(matRH_ThetaPower',[1,size(matRH_ThetaPower,2)*size(matRH_ThetaPower,1)]);
                cellRH_ThetaPower{cellCounter,1} = arrayRH_ThetaPower;
                % RH alpha
                matRH_AlphaPower = cell2mat(mat2CellRH_AlphaPower{cellCounter,1});
                arrayRH_AlphaPower = reshape(matRH_AlphaPower',[1,size(matRH_AlphaPower,2)*size(matRH_AlphaPower,1)]);
                cellRH_AlphaPower{cellCounter,1} = arrayRH_AlphaPower;
                % RH beta
                matRH_BetaPower = cell2mat(mat2CellRH_BetaPower{cellCounter,1});
                arrayRH_BetaPower = reshape(matRH_BetaPower',[1,size(matRH_BetaPower,2)*size(matRH_BetaPower,1)]);
                cellRH_BetaPower{cellCounter,1} = arrayRH_BetaPower;
                % RH gamma
                matRH_GammaPower = cell2mat(mat2CellRH_GammaPower{cellCounter,1});
                arrayRH_GammaPower = reshape(matRH_GammaPower',[1,size(matRH_GammaPower,2)*size(matRH_GammaPower,1)]);
                cellRH_GammaPower{cellCounter,1} = arrayRH_GammaPower;
                % RH MUA
                matRH_MUAPower = cell2mat(mat2CellRH_MUAPower{cellCounter,1});
                arrayRH_MUAPower = reshape(matRH_MUAPower',[1,size(matRH_MUAPower,2)*size(matRH_MUAPower,1)]);
                cellRH_MUAPower{cellCounter,1} = arrayRH_MUAPower;
                % hip delta
                matHip_DeltaPower = cell2mat(mat2CellHip_DeltaPower{cellCounter,1});
                arrayHip_DeltaPower = reshape(matHip_DeltaPower',[1,size(matHip_DeltaPower,2)*size(matHip_DeltaPower,1)]);
                cellHip_DeltaPower{cellCounter,1} = arrayHip_DeltaPower;
                % hip theta
                matHip_ThetaPower = cell2mat(mat2CellHip_ThetaPower{cellCounter,1});
                arrayHip_ThetaPower = reshape(matHip_ThetaPower',[1,size(matHip_ThetaPower,2)*size(matHip_ThetaPower,1)]);
                cellHip_ThetaPower{cellCounter,1} = arrayHip_ThetaPower;
                % hip alpha
                matHip_AlphaPower = cell2mat(mat2CellHip_AlphaPower{cellCounter,1});
                arrayHip_AlphaPower = reshape(matHip_AlphaPower',[1,size(matHip_AlphaPower,2)*size(matHip_AlphaPower,1)]);
                cellHip_AlphaPower{cellCounter,1} = arrayHip_AlphaPower;
                % hip beta
                matHip_BetaPower = cell2mat(mat2CellHip_BetaPower{cellCounter,1});
                arrayHip_BetaPower = reshape(matHip_BetaPower',[1,size(matHip_BetaPower,2)*size(matHip_BetaPower,1)]);
                cellHip_BetaPower{cellCounter,1} = arrayHip_BetaPower;
                % hip gamma
                matHip_GammaPower = cell2mat(mat2CellHip_GammaPower{cellCounter,1});
                arrayHip_GammaPower = reshape(matHip_GammaPower',[1,size(matHip_GammaPower,2)*size(matHip_GammaPower,1)]);
                cellHip_GammaPower{cellCounter,1} = arrayHip_GammaPower;
                % hip MUA
                matHip_MUAPower = cell2mat(mat2CellHip_MUAPower{cellCounter,1});
                arrayHip_MUAPower = reshape(matHip_MUAPower',[1,size(matHip_MUAPower,2)*size(matHip_MUAPower,1)]);
                cellHip_MUAPower{cellCounter,1} = arrayHip_MUAPower;
                % LH CBV
                matLH_CBV = cell2mat(mat2CellLH_CBV{cellCounter,1});
                arrayLH_CBV = reshape(matLH_CBV',[1,size(matLH_CBV,2)*size(matLH_CBV,1)]);
                cellLH_CBV{cellCounter,1} = arrayLH_CBV;
                % RH CBV
                matRH_CBV = cell2mat(mat2CellRH_CBV{cellCounter,1});
                arrayRH_CBV = reshape(matRH_CBV',[1,size(matRH_CBV,2)*size(matRH_CBV,1)]);
                cellRH_CBV{cellCounter,1} = arrayRH_CBV;
                % f LH CBV
                matfLH_CBV = cell2mat(mat2CellfLH_CBV{cellCounter,1});
                arrayfLH_CBV = reshape(matfLH_CBV',[1,size(matfLH_CBV,2)*size(matfLH_CBV,1)]);
                cellfLH_CBV{cellCounter,1} = arrayfLH_CBV;
                % f RH CBV
                matfRH_CBV = cell2mat(mat2CellfRH_CBV{cellCounter,1});
                arrayfRH_CBV = reshape(matfRH_CBV',[1,size(matfRH_CBV,2)*size(matfRH_CBV,1)]);
                cellfRH_CBV{cellCounter,1} = arrayfRH_CBV;
                % LH HbT
                matLH_HbT = cell2mat(mat2CellLH_HbT{cellCounter,1});
                arrayLH_HbT = reshape(matLH_HbT',[1,size(matLH_HbT,2)*size(matLH_HbT,1)]);
                cellLH_HbT{cellCounter,1} = arrayLH_HbT;
                % RH HbT
                matRH_HbT = cell2mat(mat2CellRH_HbT{cellCounter,1});
                arrayRH_HbT = reshape(matRH_HbT',[1,size(matRH_HbT,2)*size(matRH_HbT,1)]);
                cellRH_HbT{cellCounter,1} = arrayRH_HbT;
                % f LH HbT
                matfLH_HbT = cell2mat(mat2CellfLH_HbT{cellCounter,1});
                arrayfLH_HbT = reshape(matfLH_HbT',[1,size(matfLH_HbT,2)*size(matfLH_HbT,1)]);
                cellfLH_HbT{cellCounter,1} = arrayfLH_HbT;
                % f RH HbT
                matfRH_HbT = cell2mat(mat2CellfRH_HbT{cellCounter,1});
                arrayfRH_HbT = reshape(matfRH_HbT',[1,size(matfRH_HbT,2)*size(matfRH_HbT,1)]);
                cellfRH_HbT{cellCounter,1} = arrayfRH_HbT;
                % LH GCaMP
                matLH_GCaMP = cell2mat(mat2CellLH_GCaMP{cellCounter,1});
                arrayLH_GCaMP = reshape(matLH_GCaMP',[1,size(matLH_GCaMP,2)*size(matLH_GCaMP,1)]);
                cellLH_GCaMP{cellCounter,1} = arrayLH_GCaMP;
                % RH GCaMP
                matRH_GCaMP = cell2mat(mat2CellRH_GCaMP{cellCounter,1});
                arrayRH_GCaMP = reshape(matRH_GCaMP',[1,size(matRH_GCaMP,2)*size(matRH_GCaMP,1)]);
                cellRH_GCaMP{cellCounter,1} = arrayRH_GCaMP;
                % f LH GCaMP
                matfLH_GCaMP = cell2mat(mat2CellfLH_GCaMP{cellCounter,1});
                arrayfLH_GCaMP = reshape(matfLH_GCaMP',[1,size(matfLH_GCaMP,2)*size(matfLH_GCaMP,1)]);
                cellfLH_GCaMP{cellCounter,1} = arrayfLH_GCaMP;
                % f RH GCaMP
                matfRH_GCaMP = cell2mat(mat2CellfRH_GCaMP{cellCounter,1});
                arrayfRH_GCaMP = reshape(matfRH_GCaMP',[1,size(matfRH_GCaMP,2)*size(matfRH_GCaMP,1)]);
                cellfRH_GCaMP{cellCounter,1} = arrayfRH_GCaMP;
                % LH HbR
                matLH_HbR = cell2mat(mat2CellLH_HbR{cellCounter,1});
                arrayLH_HbR = reshape(matLH_HbR',[1,size(matLH_HbR,2)*size(matLH_HbR,1)]);
                cellLH_HbR{cellCounter,1} = arrayLH_HbR;
                % RH HbR
                matRH_HbR = cell2mat(mat2CellRH_HbR{cellCounter,1});
                arrayRH_HbR = reshape(matRH_HbR',[1,size(matRH_HbR,2)*size(matRH_HbR,1)]);
                cellRH_HbR{cellCounter,1} = arrayRH_HbR;
                % f LH HbR
                matfLH_HbR = cell2mat(mat2CellfLH_HbR{cellCounter,1});
                arrayfLH_HbR = reshape(matfLH_HbR',[1,size(matfLH_HbR,2)*size(matfLH_HbR,1)]);
                cellfLH_HbR{cellCounter,1} = arrayfLH_HbR;
                % f RH HbR
                matfRH_HbR = cell2mat(mat2CellfRH_HbR{cellCounter,1});
                arrayfRH_HbR = reshape(matfRH_HbR',[1,size(matfRH_HbR,2)*size(matfRH_HbR,1)]);
                cellfRH_HbR{cellCounter,1} = arrayfRH_HbR;
                % LH HbO
                matLH_HbO = cell2mat(mat2CellLH_HbO{cellCounter,1});
                arrayLH_HbO = reshape(matLH_HbO',[1,size(matLH_HbO,2)*size(matLH_HbO,1)]);
                cellLH_HbO{cellCounter,1} = arrayLH_HbO;
                % RH HbO
                matRH_HbO = cell2mat(mat2CellRH_HbO{cellCounter,1});
                arrayRH_HbO = reshape(matRH_HbO',[1,size(matRH_HbO,2)*size(matRH_HbO,1)]);
                cellRH_HbO{cellCounter,1} = arrayRH_HbO;
                % f LH HbO
                matfLH_HbO = cell2mat(mat2CellfLH_HbO{cellCounter,1});
                arrayfLH_HbO = reshape(matfLH_HbO',[1,size(matfLH_HbO,2)*size(matfLH_HbO,1)]);
                cellfLH_HbO{cellCounter,1} = arrayfLH_HbO;
                % f RH HbO
                matfRH_HbO = cell2mat(mat2CellfRH_HbO{cellCounter,1});
                arrayfRH_HbO = reshape(matfRH_HbO',[1,size(matfRH_HbO,2)*size(matfRH_HbO,1)]);
                cellfRH_HbO{cellCounter,1} = arrayfRH_HbO;
                % whisker acceleration
                for x = 1:size(mat2CellWhiskerAcceleration{cellCounter,1},1)
                    targetPoints = size(mat2CellWhiskerAcceleration{cellCounter,1}{1,1},2);
                    if size(mat2CellWhiskerAcceleration{cellCounter,1}{x,1},2) ~= targetPoints
                        maxLength = size(mat2CellWhiskerAcceleration{cellCounter,1}{x,1},2);
                        difference = targetPoints - size(mat2CellWhiskerAcceleration{cellCounter,1}{x,1},2);
                        for y = 1:difference
                            mat2CellWhiskerAcceleration{cellCounter,1}{x,1}(maxLength + y) = 0;
                        end
                    end
                end
                matWhiskerAcceleration = cell2mat(mat2CellWhiskerAcceleration{cellCounter,1});
                arrayWhiskerAcceleration = reshape(matWhiskerAcceleration',[1,size(matWhiskerAcceleration,2)*size(matWhiskerAcceleration,1)]);
                cellWhiskerAcceleration{cellCounter,1} = arrayWhiskerAcceleration;
                % bin times
                matBinTimes = cell2mat(mat2CellBinTimes{cellCounter,1});
                arrayBinTimes = reshape(matBinTimes',[1,size(matBinTimes,2)*size(matBinTimes,1)]);
                cellBinTimes{cellCounter,1} = arrayBinTimes;
            end
        end
        % save the data in the SleepEventData struct
        if isfield(SleepData.(modelName),'REM') == false % if the structure is empty we need a special case to format the struct properly
            for cellLength = 1:size(cellRH_DeltaPower,2) % loop through however many sleep epochs this file has
                % cortex
                SleepData.(modelName).REM.data.cortical_RH.deltaBandPower{cellLength,1} = cellRH_DeltaPower{1,1};
                SleepData.(modelName).REM.data.cortical_RH.thetaBandPower{cellLength,1} = cellRH_ThetaPower{1,1};
                SleepData.(modelName).REM.data.cortical_RH.alphaBandPower{cellLength,1} = cellRH_AlphaPower{1,1};
                SleepData.(modelName).REM.data.cortical_RH.betaBandPower{cellLength,1} = cellRH_BetaPower{1,1};
                SleepData.(modelName).REM.data.cortical_RH.gammaBandPower{cellLength,1} = cellRH_GammaPower{1,1};
                SleepData.(modelName).REM.data.cortical_RH.muaPower{cellLength,1} = cellRH_MUAPower{1,1};
                % hippocampus
                SleepData.(modelName).REM.data.hippocampus.deltaBandPower{cellLength,1} = cellHip_DeltaPower{1,1};
                SleepData.(modelName).REM.data.hippocampus.thetaBandPower{cellLength,1} = cellHip_ThetaPower{1,1};
                SleepData.(modelName).REM.data.hippocampus.alphaBandPower{cellLength,1} = cellHip_AlphaPower{1,1};
                SleepData.(modelName).REM.data.hippocampus.betaBandPower{cellLength,1} = cellHip_BetaPower{1,1};
                SleepData.(modelName).REM.data.hippocampus.gammaBandPower{cellLength,1} = cellHip_GammaPower{1,1};
                SleepData.(modelName).REM.data.hippocampus.muaPower{cellLength,1} = cellHip_MUAPower{1,1};
                % CBV
                SleepData.(modelName).REM.data.CBV.LH{cellLength,1} = cellLH_CBV{1,1};
                SleepData.(modelName).REM.data.CBV.RH{cellLength,1} = cellRH_CBV{1,1};
                SleepData.(modelName).REM.data.CBV.fLH{cellLength,1} = cellfLH_CBV{1,1};
                SleepData.(modelName).REM.data.CBV.fRH{cellLength,1} = cellfRH_CBV{1,1};
                % HbT
                SleepData.(modelName).REM.data.HbT.LH{cellLength,1} = cellLH_HbT{1,1};
                SleepData.(modelName).REM.data.HbT.RH{cellLength,1} = cellRH_HbT{1,1};
                SleepData.(modelName).REM.data.HbT.fLH{cellLength,1} = cellfLH_HbT{1,1};
                SleepData.(modelName).REM.data.HbT.fRH{cellLength,1} = cellfRH_HbT{1,1};
                % GCaMP
                SleepData.(modelName).REM.data.GCaMP.LH{cellLength,1} = cellLH_GCaMP{1,1};
                SleepData.(modelName).REM.data.GCaMP.RH{cellLength,1} = cellRH_GCaMP{1,1};
                SleepData.(modelName).REM.data.GCaMP.fLH{cellLength,1} = cellfLH_GCaMP{1,1};
                SleepData.(modelName).REM.data.GCaMP.fRH{cellLength,1} = cellfRH_GCaMP{1,1};
                % HbR
                SleepData.(modelName).REM.data.HbR.LH{cellLength,1} = cellLH_HbR{1,1};
                SleepData.(modelName).REM.data.HbR.RH{cellLength,1} = cellRH_HbR{1,1};
                SleepData.(modelName).REM.data.HbR.fLH{cellLength,1} = cellfLH_HbR{1,1};
                SleepData.(modelName).REM.data.HbR.fRH{cellLength,1} = cellfRH_HbR{1,1};
                % HbO
                SleepData.(modelName).REM.data.HbO.LH{cellLength,1} = cellLH_HbO{1,1};
                SleepData.(modelName).REM.data.HbO.RH{cellLength,1} = cellRH_HbO{1,1};
                SleepData.(modelName).REM.data.HbO.fLH{cellLength,1} = cellfLH_HbO{1,1};
                SleepData.(modelName).REM.data.HbO.fRH{cellLength,1} = cellfRH_HbO{1,1};
                % whiskers
                SleepData.(modelName).REM.data.whiskerAcceleration{cellLength,1} = cellWhiskerAcceleration{1,1};
                % file IDs & bin times
                SleepData.(modelName).REM.FileIDs{cellLength,1} = fileID;
                SleepData.(modelName).REM.BinTimes{cellLength,1} = cellBinTimes{1,1};
            end
        else % if the struct is not empty,add each new iteration after previous data
            for cellLength = 1:size(cellRH_DeltaPower,1) % loop through however many sleep epochs this file has
                % cortex
                SleepData.(modelName).REM.data.cortical_RH.deltaBandPower{size(SleepData.(modelName).REM.data.cortical_RH.deltaBandPower,1) + 1,1} = cellRH_DeltaPower{cellLength,1};
                SleepData.(modelName).REM.data.cortical_RH.thetaBandPower{size(SleepData.(modelName).REM.data.cortical_RH.thetaBandPower,1) + 1,1} = cellRH_ThetaPower{cellLength,1};
                SleepData.(modelName).REM.data.cortical_RH.alphaBandPower{size(SleepData.(modelName).REM.data.cortical_RH.alphaBandPower,1) + 1,1} = cellRH_AlphaPower{cellLength,1};
                SleepData.(modelName).REM.data.cortical_RH.betaBandPower{size(SleepData.(modelName).REM.data.cortical_RH.betaBandPower,1) + 1,1} = cellRH_BetaPower{cellLength,1};
                SleepData.(modelName).REM.data.cortical_RH.gammaBandPower{size(SleepData.(modelName).REM.data.cortical_RH.gammaBandPower,1) + 1,1} = cellRH_GammaPower{cellLength,1};
                SleepData.(modelName).REM.data.cortical_RH.muaPower{size(SleepData.(modelName).REM.data.cortical_RH.muaPower,1) + 1,1} = cellRH_MUAPower{cellLength,1};
                % hippocampus
                SleepData.(modelName).REM.data.hippocampus.deltaBandPower{size(SleepData.(modelName).REM.data.hippocampus.deltaBandPower,1) + 1,1} = cellHip_DeltaPower{cellLength,1};
                SleepData.(modelName).REM.data.hippocampus.thetaBandPower{size(SleepData.(modelName).REM.data.hippocampus.thetaBandPower,1) + 1,1} = cellHip_ThetaPower{cellLength,1};
                SleepData.(modelName).REM.data.hippocampus.alphaBandPower{size(SleepData.(modelName).REM.data.hippocampus.alphaBandPower,1) + 1,1} = cellHip_AlphaPower{cellLength,1};
                SleepData.(modelName).REM.data.hippocampus.betaBandPower{size(SleepData.(modelName).REM.data.hippocampus.betaBandPower,1) + 1,1} = cellHip_BetaPower{cellLength,1};
                SleepData.(modelName).REM.data.hippocampus.gammaBandPower{size(SleepData.(modelName).REM.data.hippocampus.gammaBandPower,1) + 1,1} = cellHip_GammaPower{cellLength,1};
                SleepData.(modelName).REM.data.hippocampus.muaPower{size(SleepData.(modelName).REM.data.hippocampus.muaPower,1) + 1,1} = cellHip_MUAPower{cellLength,1};
                % CBV
                SleepData.(modelName).REM.data.CBV.LH{size(SleepData.(modelName).REM.data.CBV.LH,1) + 1,1} = cellLH_CBV{cellLength,1};
                SleepData.(modelName).REM.data.CBV.RH{size(SleepData.(modelName).REM.data.CBV.RH,1) + 1,1} = cellRH_CBV{cellLength,1};
                SleepData.(modelName).REM.data.CBV.fLH{size(SleepData.(modelName).REM.data.CBV.fLH,1) + 1,1} = cellfLH_CBV{cellLength,1};
                SleepData.(modelName).REM.data.CBV.fRH{size(SleepData.(modelName).REM.data.CBV.fRH,1) + 1,1} = cellfRH_CBV{cellLength,1};
                % HbT
                SleepData.(modelName).REM.data.HbT.LH{size(SleepData.(modelName).REM.data.HbT.LH,1) + 1,1} = cellLH_HbT{cellLength,1};
                SleepData.(modelName).REM.data.HbT.RH{size(SleepData.(modelName).REM.data.HbT.RH,1) + 1,1} = cellRH_HbT{cellLength,1};
                SleepData.(modelName).REM.data.HbT.fLH{size(SleepData.(modelName).REM.data.HbT.fLH,1) + 1,1} = cellfLH_HbT{cellLength,1};
                SleepData.(modelName).REM.data.HbT.fRH{size(SleepData.(modelName).REM.data.HbT.fRH,1) + 1,1} = cellfRH_HbT{cellLength,1};
                % GCaMP
                SleepData.(modelName).REM.data.GCaMP.LH{size(SleepData.(modelName).REM.data.GCaMP.LH,1) + 1,1} = cellLH_GCaMP{cellLength,1};
                SleepData.(modelName).REM.data.GCaMP.RH{size(SleepData.(modelName).REM.data.GCaMP.RH,1) + 1,1} = cellRH_GCaMP{cellLength,1};
                SleepData.(modelName).REM.data.GCaMP.fLH{size(SleepData.(modelName).REM.data.GCaMP.fLH,1) + 1,1} = cellfLH_GCaMP{cellLength,1};
                SleepData.(modelName).REM.data.GCaMP.fRH{size(SleepData.(modelName).REM.data.GCaMP.fRH,1) + 1,1} = cellfRH_GCaMP{cellLength,1};
                % HbR
                SleepData.(modelName).REM.data.HbR.LH{size(SleepData.(modelName).REM.data.HbR.LH,1) + 1,1} = cellLH_HbR{cellLength,1};
                SleepData.(modelName).REM.data.HbR.RH{size(SleepData.(modelName).REM.data.HbR.RH,1) + 1,1} = cellRH_HbR{cellLength,1};
                SleepData.(modelName).REM.data.HbR.fLH{size(SleepData.(modelName).REM.data.HbR.fLH,1) + 1,1} = cellfLH_HbR{cellLength,1};
                SleepData.(modelName).REM.data.HbR.fRH{size(SleepData.(modelName).REM.data.HbR.fRH,1) + 1,1} = cellfRH_HbR{cellLength,1};
                % HbO
                SleepData.(modelName).REM.data.HbO.LH{size(SleepData.(modelName).REM.data.HbO.LH,1) + 1,1} = cellLH_HbO{cellLength,1};
                SleepData.(modelName).REM.data.HbO.RH{size(SleepData.(modelName).REM.data.HbO.RH,1) + 1,1} = cellRH_HbO{cellLength,1};
                SleepData.(modelName).REM.data.HbO.fLH{size(SleepData.(modelName).REM.data.HbO.fLH,1) + 1,1} = cellfLH_HbO{cellLength,1};
                SleepData.(modelName).REM.data.HbO.fRH{size(SleepData.(modelName).REM.data.HbO.fRH,1) + 1,1} = cellfRH_HbO{cellLength,1};
                % whiskers
                SleepData.(modelName).REM.data.whiskerAcceleration{size(SleepData.(modelName).REM.data.whiskerAcceleration,1) + 1,1} = cellWhiskerAcceleration{cellLength,1};
                % file IDs & bin times
                SleepData.(modelName).REM.FileIDs{size(SleepData.(modelName).REM.FileIDs,1) + 1,1} = fileID;
                SleepData.(modelName).REM.BinTimes{size(SleepData.(modelName).REM.BinTimes,1) + 1,1} = cellBinTimes{cellLength,1};
            end
        end
    end
    disp(['Adding REM sleeping epochs from ProcData file ' num2str(aa) ' of ' num2str(size(procDataFileIDs,1)) '...']); disp(' ')
end
% save structure
disp([modelName ' model data added to SleepData structure.']); disp(' ')

end
