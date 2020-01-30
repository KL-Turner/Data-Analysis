function [SleepData] = CreateSleepData(procDataFiles, sleepTime)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: This function uses the sleep logicals in each ProcData file to find periods where there are 60 seconds of 
%            consecutive ones within the sleep logical (12 or more). If a ProcData file's sleep logical contains one or
%            more of these 60 second periods, each of those bins is gathered from the data and put into the SleepEventData.mat
%            struct along with the file's name. 
%________________________________________________________________________________________________________________________
%
%   Inputs: The function loops through each ProcData file within the current folder - no inputs to the function itself
%           This was done as it was easier to add to the SleepEventData struct instead of loading it and then adding to it
%           with each ProcData loop.
%
%   Outputs: SleepEventData.mat struct
%________________________________________________________________________________________________________________________

%% BLOCK PURPOSE: Create sleep scored data structure.
% Identify sleep epochs and place in SleepEventData.mat structure
SleepData = [];
sleepBins = sleepTime / 5;
for fileNumber = 1:size(procDataFiles, 1)           % Loop through the list of ProcData files
    procDataFile = procDataFiles(fileNumber, :);    % Pull character string associated with the current file
    load(procDataFile);                             % Load in procDataFile associated with character string
    [animal, ~, ~, fileID] = GetFileInfo(procDataFile);     % Gather file info
    
    clear LH_deltaPower RH_deltaPower LH_thetaPower RH_thetaPower LH_gammaPower RH_gammaPower WhiskerAcceleration HeartRate LH_CBV RH_CBV LH_ElectrodeCBV RH_ElectrodeCBV BinTimes
    clear cellLH_DeltaPower cellRH_DeltaPower cellLH_ThetaPower cellRH_ThetaPower cellLH_GammaPower cellRH_GammaPower cellWhiskerAcceleration cellHeartRate cellLH_CBV cellRH_CBV cellLH_ElectrodeCBV cellRH_ElectrodeCBV cellBinTimes
    clear mat2CellLH_DeltaPower mat2CellRH_DeltaPower mat2CellLH_ThetaPower mat2CellRH_TheatPower mat2CellLH_GammaPower mat2CellRH_GammaPower mat2CellWhiskerAcceleration mat2CellHeartRate mat2CellLH_CBV mat2CellRH_CBV mat2cellLH_ElectrodeCBV mat2cellRH_ElectrodeCBV mat2CellBinTimes
    clear matLH_DeltaPower matRH_DeltaPower matLH_ThetaPower matRH_ThetaPower matLH_GammaPower matRH_GammaPower matWhiskerAcceleration matHeartRate matLH_CBV matRH_CBV matLH_ElectrodeCBV matRH_ElectrodeCBV matBinTimes
    
    sleepLogical = ProcData.Sleep.Logicals.SleepLogical;    % Logical - ones denote potential sleep epoches (5 second bins)
    targetTime = ones(1, sleepBins);   % Target time 
    sleepIndex = find(conv(sleepLogical, targetTime) >= sleepBins) - (sleepBins - 1);   % Find the periods of time where there are at least 11 more
    % 5 second epochs following. This is not the full list.
    if isempty(sleepIndex)  % If sleepIndex is empty, skip this file
        % Skip file
    else
        sleepCriteria = (0:(sleepBins - 1));     % This will be used to fix the issue in sleepIndex
        fixedSleepIndex = unique(sleepIndex + sleepCriteria);   % Sleep Index now has the proper time stamps from sleep logical  
        for indexCount = 1:length(fixedSleepIndex)    % Loop through the length of Sleep Index, and pull out associated data
            LH_deltaPower{indexCount, 1} = ProcData.Sleep.Parameters.DeltaBand_Power.LH{fixedSleepIndex(indexCount), 1};
            RH_deltaPower{indexCount, 1} = ProcData.Sleep.Parameters.DeltaBand_Power.RH{fixedSleepIndex(indexCount), 1};
            LH_thetaPower{indexCount, 1} = ProcData.Sleep.Parameters.ThetaBand_Power.LH{fixedSleepIndex(indexCount), 1};
            RH_thetaPower{indexCount, 1} = ProcData.Sleep.Parameters.ThetaBand_Power.RH{fixedSleepIndex(indexCount), 1};
            LH_gammaPower{indexCount, 1} = ProcData.Sleep.Parameters.GammaBand_Power.LH{fixedSleepIndex(indexCount), 1};
            RH_gammaPower{indexCount, 1} = ProcData.Sleep.Parameters.GammaBand_Power.RH{fixedSleepIndex(indexCount), 1};
            WhiskerAcceleration{indexCount, 1} = ProcData.Sleep.Parameters.WhiskerAcceleration{fixedSleepIndex(indexCount), 1};
            HeartRate{indexCount, 1} = ProcData.Sleep.Parameters.HeartRate{fixedSleepIndex(indexCount), 1};
            LH_CBV{indexCount, 1} = ProcData.Sleep.Parameters.CBV.LH{fixedSleepIndex(indexCount), 1};
            RH_CBV{indexCount, 1} = ProcData.Sleep.Parameters.CBV.RH{fixedSleepIndex(indexCount), 1};
            LH_ElectrodeCBV{indexCount, 1} = ProcData.Sleep.Parameters.CBV.LH_Electrode{fixedSleepIndex(indexCount), 1};
            RH_ElectrodeCBV{indexCount, 1} = ProcData.Sleep.Parameters.CBV.RH_Electrode{fixedSleepIndex(indexCount), 1};
            BinTimes{indexCount, 1} = 5*fixedSleepIndex(indexCount);
        end
        
        indexBreaks = find(fixedSleepIndex(2:end) - fixedSleepIndex(1:end - 1) > 1);    % Find if there are numerous sleep periods
        
        if isempty(indexBreaks)   % If there is only one period of sleep in this file and not multiple
            matLH_DeltaPower = cell2mat(LH_deltaPower);
            arrayLH_DeltaPower = reshape(matLH_DeltaPower', [1, size(matLH_DeltaPower, 2)*size(matLH_DeltaPower, 1)]);
            cellLH_DeltaPower = {arrayLH_DeltaPower};
            
            matRH_DeltaPower = cell2mat(RH_deltaPower);
            arrayRH_DeltaPower = reshape(matRH_DeltaPower', [1, size(matRH_DeltaPower, 2)*size(matRH_DeltaPower, 1)]);
            cellRH_DeltaPower = {arrayRH_DeltaPower};
            
            matLH_ThetaPower = cell2mat(LH_thetaPower);
            arrayLH_ThetaPower = reshape(matLH_ThetaPower', [1, size(matLH_ThetaPower, 2)*size(matLH_ThetaPower, 1)]);
            cellLH_ThetaPower = {arrayLH_ThetaPower};
            
            matRH_ThetaPower = cell2mat(RH_thetaPower);
            arrayRH_ThetaPower = reshape(matRH_ThetaPower', [1, size(matRH_ThetaPower, 2)*size(matRH_ThetaPower, 1)]);
            cellRH_ThetaPower = {arrayRH_ThetaPower};
            
            matLH_GammaPower = cell2mat(LH_gammaPower);
            arrayLH_GammaPower = reshape(matLH_GammaPower', [1, size(matLH_GammaPower, 2)*size(matLH_GammaPower, 1)]);
            cellLH_GammaPower = {arrayLH_GammaPower};
            
            matRH_GammaPower = cell2mat(RH_gammaPower);
            arrayRH_GammaPower = reshape(matRH_GammaPower', [1, size(matRH_GammaPower, 2)*size(matRH_GammaPower, 1)]);
            cellRH_GammaPower = {arrayRH_GammaPower};
            
            for x = 1:length(WhiskerAcceleration)
                targetPoints = size(WhiskerAcceleration{1, 1}, 2);
                if size(WhiskerAcceleration{x, 1}, 2) ~= targetPoints
                    maxLength = size(WhiskerAcceleration{x, 1}, 2);
                    difference = targetPoints - size(WhiskerAcceleration{x, 1}, 2);
                    for y = 1:difference
                        WhiskerAcceleration{x, 1}(maxLength + y) = 0;
                    end
                end
            end
            
            matWhiskerAcceleration = cell2mat(WhiskerAcceleration);
            arrayWhiskerAcceleration = reshape(matWhiskerAcceleration', [1, size(matWhiskerAcceleration, 2)*size(matWhiskerAcceleration, 1)]);
            cellWhiskerAcceleration = {arrayWhiskerAcceleration};
            
            for x = 1:length(HeartRate)
                targetPoints = size(HeartRate{1, 1}, 2);
                if size(HeartRate{x, 1}, 2) ~= targetPoints
                    maxLength = size(HeartRate{x, 1}, 2);
                    difference = targetPoints - size(HeartRate{x, 1}, 2);
                    for y = 1:difference
                        HeartRate{x, 1}(maxLength + y) = mean(HeartRate{x, 1});
                    end
                end
            end
            
            matHeartRate = cell2mat(HeartRate);
            arrayHeartRate = reshape(matHeartRate', [1, size(matHeartRate, 2)*size(matHeartRate, 1)]);
            cellHeartRate = {arrayHeartRate};
            
            matLH_CBV = cell2mat(LH_CBV);
            arrayLH_CBV = reshape(matLH_CBV', [1, size(matLH_CBV, 2)*size(matLH_CBV, 1)]);
            cellLH_CBV = {arrayLH_CBV};
            
            matRH_CBV = cell2mat(RH_CBV);
            arrayRH_CBV = reshape(matRH_CBV', [1, size(matRH_CBV, 2)*size(matRH_CBV, 1)]);
            cellRH_CBV = {arrayRH_CBV};
            
            matLH_ElectrodeCBV = cell2mat(LH_ElectrodeCBV);
            arrayLH_ElectrodeCBV = reshape(matLH_ElectrodeCBV', [1, size(matLH_ElectrodeCBV, 2)*size(matLH_ElectrodeCBV, 1)]);
            cellLH_ElectrodeCBV = {arrayLH_ElectrodeCBV};
            
            matRH_ElectrodeCBV = cell2mat(RH_ElectrodeCBV);
            arrayRH_ElectrodeCBV = reshape(matRH_ElectrodeCBV', [1, size(matRH_ElectrodeCBV, 2)*size(matRH_ElectrodeCBV, 1)]);
            cellRH_ElectrodeCBV = {arrayRH_ElectrodeCBV};
            
            matBinTimes = cell2mat(BinTimes);
            arrayBinTimes = reshape(matBinTimes', [1, size(matBinTimes, 2)*size(matBinTimes, 1)]);
            cellBinTimes = {arrayBinTimes};
        else
            count = length(fixedSleepIndex);
            holdIndex = zeros(1, (length(indexBreaks) + 1));
            
            for indexCounter = 1:length(indexBreaks) + 1
                if indexCounter == 1
                    holdIndex(indexCounter) = indexBreaks(indexCounter);
                elseif indexCounter == length(indexBreaks) + 1
                    holdIndex(indexCounter) = count - indexBreaks(indexCounter - 1);
                else
                    holdIndex(indexCounter)= indexBreaks(indexCounter) - indexBreaks(indexCounter - 1);
                end
            end
            
            splitCounter = 1:length(LH_deltaPower);
            convertedMat2Cell = mat2cell(splitCounter', holdIndex);
            
            for matCounter = 1:length(convertedMat2Cell)
                mat2CellLH_DeltaPower{matCounter, 1} = LH_deltaPower(convertedMat2Cell{matCounter, 1});
                mat2CellRH_DeltaPower{matCounter, 1} = RH_deltaPower(convertedMat2Cell{matCounter, 1});
                mat2CellLH_ThetaPower{matCounter, 1} = LH_thetaPower(convertedMat2Cell{matCounter, 1});
                mat2CellRH_ThetaPower{matCounter, 1} = RH_thetaPower(convertedMat2Cell{matCounter, 1});
                mat2CellLH_GammaPower{matCounter, 1} = LH_gammaPower(convertedMat2Cell{matCounter, 1});
                mat2CellRH_GammaPower{matCounter, 1} = RH_gammaPower(convertedMat2Cell{matCounter, 1});
                mat2CellLH_CBV{matCounter, 1} = LH_CBV(convertedMat2Cell{matCounter, 1});
                mat2CellRH_CBV{matCounter, 1} = RH_CBV(convertedMat2Cell{matCounter, 1});
                mat2CellLH_ElectrodeCBV{matCounter, 1} = LH_ElectrodeCBV(convertedMat2Cell{matCounter, 1});
                mat2CellRH_ElectrodeCBV{matCounter, 1} = RH_ElectrodeCBV(convertedMat2Cell{matCounter, 1});
                mat2CellWhiskerAcceleration{matCounter, 1} = WhiskerAcceleration(convertedMat2Cell{matCounter, 1});
                mat2CellHeartRate{matCounter, 1} = HeartRate(convertedMat2Cell{matCounter, 1});
                mat2CellBinTimes{matCounter, 1} = BinTimes(convertedMat2Cell{matCounter, 1});
            end
            
            for cellCounter = 1:length(mat2CellLH_DeltaPower)
                matLH_DeltaPower = cell2mat(mat2CellLH_DeltaPower{cellCounter, 1});
                arrayLH_DeltaPower = reshape(matLH_DeltaPower', [1, size(matLH_DeltaPower, 2)*size(matLH_DeltaPower, 1)]);
                cellLH_DeltaPower{cellCounter, 1} = arrayLH_DeltaPower;
                
                matRH_DeltaPower = cell2mat(mat2CellRH_DeltaPower{cellCounter, 1});
                arrayRH_DeltaPower = reshape(matRH_DeltaPower', [1, size(matRH_DeltaPower, 2)*size(matRH_DeltaPower, 1)]);
                cellRH_DeltaPower{cellCounter, 1} = arrayRH_DeltaPower;
                
                matLH_ThetaPower = cell2mat(mat2CellLH_ThetaPower{cellCounter, 1});
                arrayLH_ThetaPower = reshape(matLH_ThetaPower', [1, size(matLH_ThetaPower, 2)*size(matLH_ThetaPower, 1)]);
                cellLH_ThetaPower{cellCounter, 1} = arrayLH_ThetaPower;
                
                matRH_ThetaPower = cell2mat(mat2CellRH_ThetaPower{cellCounter, 1});
                arrayRH_ThetaPower = reshape(matRH_ThetaPower', [1, size(matRH_ThetaPower, 2)*size(matRH_ThetaPower, 1)]);
                cellRH_ThetaPower{cellCounter, 1} = arrayRH_ThetaPower;
                
                matLH_GammaPower = cell2mat(mat2CellLH_GammaPower{cellCounter, 1});
                arrayLH_GammaPower = reshape(matLH_GammaPower', [1, size(matLH_GammaPower, 2)*size(matLH_GammaPower, 1)]);
                cellLH_GammaPower{cellCounter, 1} = arrayLH_GammaPower;
                
                matRH_GammaPower = cell2mat(mat2CellRH_GammaPower{cellCounter, 1});
                arrayRH_GammaPower = reshape(matRH_GammaPower', [1, size(matRH_GammaPower, 2)*size(matRH_GammaPower, 1)]);
                cellRH_GammaPower{cellCounter, 1} = arrayRH_GammaPower;
                
                matLH_CBV = cell2mat(mat2CellLH_CBV{cellCounter, 1});
                arrayLH_CBV = reshape(matLH_CBV', [1, size(matLH_CBV, 2)*size(matLH_CBV, 1)]);
                cellLH_CBV{cellCounter, 1} = arrayLH_CBV;
                
                matRH_CBV = cell2mat(mat2CellRH_CBV{cellCounter, 1});
                arrayRH_CBV = reshape(matRH_CBV', [1, size(matRH_CBV, 2)*size(matRH_CBV, 1)]);
                cellRH_CBV{cellCounter, 1} = arrayRH_CBV;
                
                matLH_ElectrodeCBV = cell2mat(mat2CellLH_ElectrodeCBV{cellCounter, 1});
                arrayLH_ElectrodeCBV = reshape(matLH_ElectrodeCBV', [1, size(matLH_ElectrodeCBV, 2)*size(matLH_ElectrodeCBV, 1)]);
                cellLH_ElectrodeCBV{cellCounter, 1} = arrayLH_ElectrodeCBV;
                
                matRH_ElectrodeCBV = cell2mat(mat2CellRH_ElectrodeCBV{cellCounter, 1});
                arrayRH_ElectrodeCBV = reshape(matRH_ElectrodeCBV', [1, size(matRH_ElectrodeCBV, 2)*size(matRH_ElectrodeCBV, 1)]);
                cellRH_ElectrodeCBV{cellCounter, 1} = arrayRH_ElectrodeCBV;
                
                for x = 1:size(mat2CellWhiskerAcceleration{cellCounter, 1}, 1)
                    targetPoints = size(mat2CellWhiskerAcceleration{cellCounter, 1}{1, 1}, 2);
                    if size(mat2CellWhiskerAcceleration{cellCounter, 1}{x, 1}, 2) ~= targetPoints
                        maxLength = size(mat2CellWhiskerAcceleration{cellCounter, 1}{x, 1}, 2);
                        difference = targetPoints - size(mat2CellWhiskerAcceleration{cellCounter, 1}{x, 1}, 2);
                        for y = 1:difference
                            mat2CellWhiskerAcceleration{cellCounter, 1}{x, 1}(maxLength + y) = 0;
                        end
                    end
                end
                
                matWhiskerAcceleration = cell2mat(mat2CellWhiskerAcceleration{cellCounter, 1});
                arrayWhiskerAcceleration = reshape(matWhiskerAcceleration', [1, size(matWhiskerAcceleration, 2)*size(matWhiskerAcceleration, 1)]);
                cellWhiskerAcceleration{cellCounter, 1} = arrayWhiskerAcceleration;
                
                for x = 1:size(mat2CellHeartRate{cellCounter, 1}, 1)
                    targetPoints = size(mat2CellHeartRate{cellCounter, 1}{1, 1}, 2);
                    if size(mat2CellHeartRate{cellCounter, 1}{x, 1}, 2) ~= targetPoints
                        maxLength = size(mat2CellHeartRate{cellCounter, 1}{x, 1}, 2);
                        difference = targetPoints - size(mat2CellHeartRate{cellCounter, 1}{x, 1}, 2);
                        for y = 1:difference
                            mat2CellHeartRate{cellCounter, 1}{x, 1}(maxLength + y) = mean(mat2CellHeartRate{cellCounter, 1}{x, 1});
                        end
                    end
                end
                
                matHeartRate = cell2mat(mat2CellHeartRate{cellCounter, 1});
                arrayHeartRate = reshape(matHeartRate', [1, size(matHeartRate, 2)*size(matHeartRate, 1)]);
                cellHeartRate{cellCounter, 1} = arrayHeartRate;
                
                matBinTimes = cell2mat(mat2CellBinTimes{cellCounter, 1});
                arrayBinTimes = reshape(matBinTimes', [1, size(matBinTimes, 2)*size(matBinTimes, 1)]);
                cellBinTimes{cellCounter, 1} = arrayBinTimes;
            end
        end
        
        %% BLOCK PURPOSE: Save the data in the SleepEventData struct
        if isempty(SleepData)  % If the structure is empty, we need a special case to format the struct properly
            for cellLength = 1:size(cellLH_DeltaPower, 2)   % Loop through however many sleep epochs this file has
                SleepData.NREM.Data.DeltaBand_Power.LH{cellLength, 1} = cellLH_DeltaPower{1, 1};
                SleepData.NREM.Data.DeltaBand_Power.RH{cellLength, 1} = cellRH_DeltaPower{1, 1};
                SleepData.NREM.Data.ThetaBand_Power.LH{cellLength, 1} = cellLH_ThetaPower{1, 1};
                SleepData.NREM.Data.ThetaBand_Power.RH{cellLength, 1} = cellRH_ThetaPower{1, 1};
                SleepData.NREM.Data.GammaBand_Power.LH{cellLength, 1} = cellLH_GammaPower{1, 1};
                SleepData.NREM.Data.GammaBand_Power.RH{cellLength, 1} = cellRH_GammaPower{1, 1};
                SleepData.NREM.Data.CBV.LH{cellLength, 1} = cellLH_CBV{1, 1};
                SleepData.NREM.Data.CBV.RH{cellLength, 1} = cellRH_CBV{1, 1};
                SleepData.NREM.Data.CBV.LH_Electrode{cellLength, 1} = cellLH_ElectrodeCBV{1, 1};
                SleepData.NREM.Data.CBV.RH_Electrode{cellLength, 1} = cellRH_ElectrodeCBV{1, 1};
                SleepData.NREM.Data.WhiskerAcceleration{cellLength, 1} = cellWhiskerAcceleration{1, 1};
                SleepData.NREM.Data.HeartRate{cellLength, 1} = cellHeartRate{1, 1};
                SleepData.NREM.FileIDs{cellLength, 1} = fileID;
                SleepData.NREM.BinTimes{cellLength, 1} = cellBinTimes{1, 1};
            end
        else    % If the struct is not empty, add each new iteration after previous data
            for cellLength = 1:size(cellLH_DeltaPower, 1)   % Loop through however many sleep epochs this file has
                SleepData.NREM.Data.DeltaBand_Power.LH{size(SleepData.NREM.Data.DeltaBand_Power.LH, 1) + 1, 1} = cellLH_DeltaPower{cellLength, 1};
                SleepData.NREM.Data.DeltaBand_Power.RH{size(SleepData.NREM.Data.DeltaBand_Power.RH, 1) + 1, 1} = cellRH_DeltaPower{cellLength, 1};
                SleepData.NREM.Data.ThetaBand_Power.LH{size(SleepData.NREM.Data.ThetaBand_Power.LH, 1) + 1, 1} = cellLH_ThetaPower{cellLength, 1};
                SleepData.NREM.Data.ThetaBand_Power.RH{size(SleepData.NREM.Data.ThetaBand_Power.RH, 1) + 1, 1} = cellRH_ThetaPower{cellLength, 1};
                SleepData.NREM.Data.GammaBand_Power.LH{size(SleepData.NREM.Data.GammaBand_Power.LH, 1) + 1, 1} = cellLH_GammaPower{cellLength, 1};
                SleepData.NREM.Data.GammaBand_Power.RH{size(SleepData.NREM.Data.GammaBand_Power.RH, 1) + 1, 1} = cellRH_GammaPower{cellLength, 1};
                SleepData.NREM.Data.CBV.LH{size(SleepData.NREM.Data.CBV.LH, 1) + 1, 1} = cellLH_CBV{cellLength, 1};
                SleepData.NREM.Data.CBV.RH{size(SleepData.NREM.Data.CBV.RH, 1) + 1, 1} = cellRH_CBV{cellLength, 1};
                SleepData.NREM.Data.CBV.LH_Electrode{size(SleepData.NREM.Data.CBV.LH_Electrode, 1) + 1, 1} = cellLH_ElectrodeCBV{cellLength, 1};
                SleepData.NREM.Data.CBV.RH_Electrode{size(SleepData.NREM.Data.CBV.RH_Electrode, 1) + 1, 1} = cellRH_ElectrodeCBV{cellLength, 1};
                SleepData.NREM.Data.WhiskerAcceleration{size(SleepData.NREM.Data.WhiskerAcceleration, 1) + 1, 1} = cellWhiskerAcceleration{cellLength, 1};
                SleepData.NREM.Data.HeartRate{size(SleepData.NREM.Data.HeartRate, 1) + 1, 1} = cellHeartRate{cellLength, 1};
                SleepData.NREM.FileIDs{size(SleepData.NREM.FileIDs, 1) + 1, 1} = fileID;
                SleepData.NREM.BinTimes{size(SleepData.NREM.BinTimes, 1) + 1, 1} = cellBinTimes{cellLength, 1};
            end
        end        
    end
        
    disp(['Adding NREM sleeping epochs from ProcData file ' num2str(fileNumber) ' of ' num2str(size(procDataFiles, 1)) '...']); disp(' ')
end

save([animal '_SleepData.mat'], 'SleepData'); % Save under animal name in current folder

disp('SleepData structure completed.'); disp(' ')

end
