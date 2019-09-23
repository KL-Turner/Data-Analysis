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

%% 
animalIDs = {'T101', 'T102', 'T103', 'T105', 'T108', 'T109'};
driveLetters = {'E', 'E', 'E', 'F', 'F', 'F'};
behavFields = {'Rest', 'NREM', 'REM'};

for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    driveLetter = driveLetters{1,a};
    dataPath = [driveLetter ':\' animalID '\Combined Imaging\'];
    cd(dataPath)
    load([animalID '_AnalysisResults.mat']);
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};    
        data.(behavField).LH_XC_Vals(:,:,a) = AnalysisResults.XCorr.(behavField).LH.XC_Vals;
        data.(behavField).LH_Lags(:,:,a) = AnalysisResults.XCorr.(behavField).LH.lags;
        data.(behavField).LH_F(:,:,a) = AnalysisResults.XCorr.(behavField).LH.F; 
        data.(behavField).RH_XC_Vals(:,:,a) = AnalysisResults.XCorr.(behavField).RH.XC_Vals;
        data.(behavField).RH_Lags(:,:,a) = AnalysisResults.XCorr.(behavField).RH.lags;
        data.(behavField).RH_F(:,:,a) = AnalysisResults.XCorr.(behavField).RH.F;
    end
end

for c = 1:length(behavFields)
    behavField = behavFields{1,c};
    data.(behavField).catXC_Vals = cat(3,data.(behavField).LH_XC_Vals, data.(behavField).RH_XC_Vals);
    data.(behavField).catLags = cat(3,data.(behavField).LH_Lags, data.(behavField).RH_Lags);
    data.(behavField).catF = cat(3,data.(behavField).LH_F, data.(behavField).RH_F);
end

for d = 1:length(behavFields)
    behavField = behavFields{1,d};
    data.(behavField).meanXC_Vals = mean(data.(behavField).catXC_Vals, 3);
    data.(behavField).meanLags = mean(data.(behavField).catLags, 3);
    data.(behavField).meanF = mean(data.(behavField).catF, 3);
end

for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    figName = [behavField 'fig'];
    figName = figure;
    maxLag = data.(behavField).meanLags(end);
    imagesc(data.(behavField).meanLags, data.(behavField).meanF, data.(behavField).meanXC_Vals);
    colormap parula
    colorbar
    title([behavField ' Cross Correlation'])
    ylabel('Frequency (Hz)')
    xlabel('Time (sec)')
    xticks([-maxLag -maxLag/2 0 maxLag/2 maxLag])
    xticklabels({'-5', '-2.5', '0', '2.5' '5'})
    axis square
    axis xy
    ylim([1 100])
    
    dirpath = 'C:\Users\klt8\Documents\Analysis Average Figures\Cross Correlation\';
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(figName, [dirpath behavField '_XCorr']);
    
end
