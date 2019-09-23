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

%% Rest and All data
animalIDs = {'T101', 'T102', 'T103', 'T105', 'T108', 'T109'};
driveLetters = {'E', 'E', 'E', 'F', 'F', 'F'};
behavFields = {'Rest', 'NREM', 'REM'};
coherr_dataTypes = {'CBV', 'deltaBandPower', 'thetaBandPower', 'alphaBandPower', 'betaBandPower', 'gammaBandPower'};
graphColors = {'sapphire', 'coral red', 'vegas gold', 'electric purple', 'jungle green', 'deep carrot orange', 'battleship grey'}; 

for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\' animalID '\Combined Imaging\'];
    cd(dataPath)
    load([animalID '_AnalysisResults.mat']);
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        for c = 1:length(coherr_dataTypes)
            coherr_dataType = coherr_dataTypes{1,c};
            data.(behavField).(coherr_dataType).C(:,a) = AnalysisResults.Coherence.(behavField).(coherr_dataType).C;
            data.(behavField).(coherr_dataType).f(:,a) = AnalysisResults.Coherence.(behavField).(coherr_dataType).f;
        end
    end
end

for d = 1:length(behavFields)
    behavField = behavFields{1,d};
    for e = 1:length(coherr_dataTypes)
        coherr_dataType = coherr_dataTypes{1,e};
        data.(behavField).(coherr_dataType).meanC = mean(data.(behavField).(coherr_dataType).C,2);
        data.(behavField).(coherr_dataType).meanf = mean(data.(behavField).(coherr_dataType).f,2);
        data.(behavField).(coherr_dataType).stdC = (std(data.(behavField).(coherr_dataType).C,0,2));
    end
end

%% Figures
for f = 1:length(behavFields)
    behavField = behavFields{1,f};
    behavFig = figure;
    for g = 1:length(coherr_dataTypes)
        coherr_dataType = coherr_dataTypes{1,g};
        graphColor = graphColors{1,g};
        xData = data.(behavField).(coherr_dataType).meanf;
        yData = data.(behavField).(coherr_dataType).meanC;
        stDev = data.(behavField).(coherr_dataType).stdC;
        p(g) = plot(xData, yData, 'color', colors_IOS(graphColor), 'LineWidth', 2);
        hold on
        plot(xData, yData + stDev, 'color', colors_IOS(graphColor))
        plot(xData, yData - stDev, 'color', colors_IOS(graphColor))
    end
    
    title(['L/R ' behavField ' Coherence'])
    ylabel('Coherence')
    xlabel('Frequency (Hz)')
    legend(p(1:length(coherr_dataTypes)), coherr_dataTypes)
    xlim([0 1])
    ylim([0 1])
    axis square
   
    dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\Coherence\';
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end 
    savefig(behavFig, [dirpath behavField '_AverageCoherence']);
    
end