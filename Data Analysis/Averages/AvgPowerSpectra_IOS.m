%_______________________________________________________________________________________________
% Written by Kevin L. Turner, Jr. 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%_______________________________________________________________________________________________
%
%   Purpose: 
%_______________________________________________________________________________________________
%
%   Inputs: 
%
%   Outputs: 
%________________________________________________________________________________________________

clear
clc

%% 
animalIDs = {'T101', 'T102', 'T103', 'T105', 'T108', 'T109'};
driveLetters = {'E', 'E', 'E', 'F', 'F', 'F'};
behavFields = {'Rest', 'NREM', 'REM'};
powerspec_dataTypes = {'CBV', 'deltaBandPower', 'thetaBandPower', 'alphaBandPower', 'betaBandPower', 'gammaBandPower'};
setA = {'LH', 'RH'};
setB = {'LH', 'RH', 'Hip'};
setC = {'LH', 'Hip'};
graphColors = {'sapphire', 'coral red', 'vegas gold', 'electric purple', 'jungle green', 'deep carrot orange', 'battleship grey'}; 

for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\' animalID '\Combined Imaging\'];
    cd(dataPath)
    load([animalID '_AnalysisResults.mat']);
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        for c = 1:length(powerspec_dataTypes)
            powerspec_dataType = powerspec_dataTypes{1,c};
            if strcmp(powerspec_dataType, 'CBV') == true
                subDataTypes = setA;
                for d = 1:length(subDataTypes)
                    subDataType = subDataTypes{1,d};
                    data.(behavField).(powerspec_dataType).(subDataType).S(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerspec_dataType).(subDataType).S;
                    data.(behavField).(powerspec_dataType).(subDataType).f(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerspec_dataType).(subDataType).f;
                end
            else
                subDataTypes = setB;
                for d = 1:length(subDataTypes)
                    subDataType = subDataTypes{1,d};
                    data.(behavField).(powerspec_dataType).(subDataType).S(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerspec_dataType).(subDataType).S;
                    data.(behavField).(powerspec_dataType).(subDataType).f(:,a) = AnalysisResults.PowerSpectra.(behavField).(powerspec_dataType).(subDataType).f;
                end
            end
        end
    end
end

for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    for f = 1:length(powerspec_dataTypes)
        powerspec_dataType = powerspec_dataTypes{1,f};
        LH_data_S = data.(behavField).(powerspec_dataType).LH.S;
        LH_data_f = data.(behavField).(powerspec_dataType).RH.f;
        RH_data_S = data.(behavField).(powerspec_dataType).RH.S;
        RH_data_f = data.(behavField).(powerspec_dataType).RH.f;
        data.(behavField).(powerspec_dataType).LH.cortCat_S = horzcat(LH_data_S, RH_data_S);
        data.(behavField).(powerspec_dataType).LH.cortCat_f = horzcat(LH_data_f, RH_data_f);
    end
end

for g = 1:length(behavFields)
    behavField = behavFields{1,g};
    for h = 1:length(powerspec_dataTypes)
        powerspec_dataType = powerspec_dataTypes{1,h};
        for j = 1:length(setC)
            subDataType = setC{1,j};
            if strcmp(subDataType, 'LH') == true
                data.(behavField).(powerspec_dataType).(subDataType).meanS = mean(data.(behavField).(powerspec_dataType).(subDataType).cortCat_S,2);
                data.(behavField).(powerspec_dataType).(subDataType).meanf = mean(data.(behavField).(powerspec_dataType).(subDataType).cortCat_f,2);
                data.(behavField).(powerspec_dataType).(subDataType).stdS = (std(data.(behavField).(powerspec_dataType).(subDataType).cortCat_S,0,2));
            elseif strcmp(subDataType, 'Hip') == true && strcmp(powerspec_dataType, 'CBV') == false
                data.(behavField).(powerspec_dataType).(subDataType).meanS = mean(data.(behavField).(powerspec_dataType).(subDataType).S,2);
                data.(behavField).(powerspec_dataType).(subDataType).meanf = mean(data.(behavField).(powerspec_dataType).(subDataType).f,2);
                data.(behavField).(powerspec_dataType).(subDataType).stdS = (std(data.(behavField).(powerspec_dataType).(subDataType).S,0,2));
            end
        end
    end
end

%% Figures
setD = {'deltaBandPower', 'thetaBandPower', 'alphaBandPower', 'betaBandPower', 'gammaBandPower', 'CBV'};
setE = {'deltaBandPower', 'thetaBandPower', 'alphaBandPower', 'betaBandPower', 'gammaBandPower'};

for k = 1:length(behavFields)
    behavField = behavFields{1,k};
    corticalFig = figure;
    for m = 1:length(setD)
        dataType = setD{1,m};
        graphColor = graphColors{1,m};
        xData = data.(behavField).(dataType).LH.meanf;
        yData = data.(behavField).(dataType).LH.meanS;
        stDev = data.(behavField).(dataType).LH.stdS;
        p(m) = loglog(xData, yData, 'color', colors_IOS(graphColor), 'LineWidth', 2);
        hold on
%         loglog(xData, yData + stDev, 'color', colors_IOS(graphColor))
%         loglog(xData, yData - stDev, 'color', colors_IOS(graphColor))
    end
    
    title(['Cortical ' behavField ' Power Spectra'])
    ylabel('Power')
    xlabel('Frequency (Hz)')
    legend(p(1:length(setD)), setD)
    xlim([0.1 1])
    axis square
   
    dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\Power Spectra\';
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end 
    savefig(corticalFig, [dirpath behavField '_Cortex_AveragePowerSpectra']);
end
    
for k = 1:length(behavFields)
    behavField = behavFields{1,k};
    hipFig = figure;
    for m = 1:length(setE)
        dataType = setE{1,m};
        graphColor = graphColors{1,m};
        xData = data.(behavField).(dataType).Hip.meanf;
        yData = data.(behavField).(dataType).Hip.meanS;
        stDev = data.(behavField).(dataType).Hip.stdS;
        p(m) = loglog(xData, yData, 'color', colors_IOS(graphColor), 'LineWidth', 2);
        hold on
%         loglog(xData, yData + stDev, 'color', colors_IOS(graphColor))
%         loglog(xData, yData - stDev, 'color', colors_IOS(graphColor))
    end
    
    title(['Hippocampal ' behavField ' Power Spectra'])
    ylabel('Power')
    xlabel('Frequency (Hz)')
    legend(p(1:length(setD)), setD)
    xlim([0.1 1])
    axis square
   
    dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\Power Spectra\';
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end 
    savefig(hipFig, [dirpath behavField '_Hippocampus_AveragePowerSpectra']);
end
