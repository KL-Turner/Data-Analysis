function [SpecData] = NormalizeSpectrograms_IOS(specDataFiles, neuralDataTypes, RestingBaselines)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Normalizes each spectrogram by the resting baseline for that day.
%________________________________________________________________________________________________________________________
%
%   Inputs: List of spectrogram files, and the RestingBaselines struct that contains the time indeces for each rest file.
%
%   Outputs: A normalized 'S' field for each Spectrogram.
%
%   Last Revised: March 22nd, 2019    
%________________________________________________________________________________________________________________________

for a = 1:size(specDataFiles,1)
    disp(['Normalizing spectrogram file ' num2str(a) ' of ' num2str(size(specDataFiles,1)) '...']); disp(' ')
    load(specDataFiles(a,:), '-mat');
    [animalID, fileDate, ~] = GetFileInfo_IOS(specDataFiles(a,:));
    strDay = ConvertDate_IOS(fileDate);
    for b = 1:length(neuralDataTypes)
        neuralDataType = neuralDataTypes{1,b};
        baseLine1 = RestingBaselines.Spectrograms.(neuralDataType).oneSec.(strDay);
        baseLine5 = RestingBaselines.Spectrograms.(neuralDataType).fiveSec.(strDay);
        
        S1 = SpecData.(neuralDataType).oneSec.S;
        S5 = SpecData.(neuralDataType).fiveSec.S;
        
        holdMatrix1 = baseLine1.*ones(size(S1));
        holdMatrix5 = baseLine5.*ones(size(S5));
        
        normS1 = (S1 - holdMatrix1)./holdMatrix1;
        normS5 = (S5 - holdMatrix5)./holdMatrix5;
        
        SpecData.(neuralDataType).oneSec.normS = normS1;
        SpecData.(neuralDataType).fiveSec.normS = normS5;
    end
    save(specDataFiles(a,:), 'SpecData')
end

AllSpecStructFileID = [animalID '_AllSpecStruct.mat'];
if ~exist(AllSpecStructFileID)
    % Character list of all SpecData files
    specDataFileStruct = dir('*_SpecData.mat');
    specDataFiles = {specDataFileStruct.name}';
    specDataFileIDs = char(specDataFiles);
    for c = 1: size(specDataFileIDs,1)
        specDataFileID = specDataFileIDs(c,:);
        disp(['Adding spectrogram data to full struct for file number ' num2str(c) ' of ' num2str(size(specDataFiles, 1)) '...']); disp(' ')
        load(specDataFileID)
        for d = 1:length(neuralDataTypes)
            neuralDataType = neuralDataTypes{1,d};
            AllSpecData.(neuralDataType).fileIDs{c,1} = specDataFileID;
            AllSpecData.(neuralDataType).fiveSec.S{c,1} =  SpecData.(neuralDataType).fiveSec.S;
            AllSpecData.(neuralDataType).fiveSec.normS{c,1} = SpecData.(neuralDataType).fiveSec.normS;
            AllSpecData.(neuralDataType).fiveSec.T{c,1} =  SpecData.(neuralDataType).fiveSec.T;
            AllSpecData.(neuralDataType).fiveSec.F{c,1} =  SpecData.(neuralDataType).fiveSec.F;
            AllSpecData.(neuralDataType).fiveSec.params =  SpecData.(neuralDataType).fiveSec.params;
            AllSpecData.(neuralDataType).fiveSec.movingwin =  SpecData.(neuralDataType).fiveSec.movingwin;
            
            AllSpecData.(neuralDataType).oneSec.S{c,1} = SpecData.(neuralDataType).oneSec.S;
            AllSpecData.(neuralDataType).oneSec.normS{c,1} = SpecData.(neuralDataType).oneSec.normS;
            AllSpecData.(neuralDataType).oneSec.T{c,1} = SpecData.(neuralDataType).oneSec.T;
            AllSpecData.(neuralDataType).oneSec.F{c,1} = SpecData.(neuralDataType).oneSec.F;
            AllSpecData.(neuralDataType).oneSec.params = SpecData.(neuralDataType).oneSec.params;
            AllSpecData.(neuralDataType).oneSec.movingwin = SpecData.(neuralDataType).oneSec.movingwin;
        end
    end
    disp('Saving structure...'); disp(' ')
    save(AllSpecStructFileID, 'AllSpecData', '-v7.3')
end

end
