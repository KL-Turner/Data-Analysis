function [SpectrogramData] = CreateTrialSpectrograms(animal, neuralDataType, rawDataFileIDs, SpectrogramData)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: //
%________________________________________________________________________________________________________________________
%
%   Inputs: //
%
%   Outputs: //
%________________________________________________________________________________________________________________________

for fileNumber = 1:length(rawDataFileIDs)
    rawDataFileID = rawDataFileIDs(fileNumber, :);
    load(rawDataFileID);
    [~, ~, ~, fileID] = GetFileInfo(rawDataFileID);
    RawNeuro = RawData.Data.(['Neural_' neuralDataType]);
    
    % Spectrogram parameters
    params.tapers = [5 9];
    params.Fs = RawData.Notes.analogSamplingRate;
    params.fpass = [1 100];
    movingwin1 = [1 1/5];
    movingwin5 = [5 1/5];

    disp(['Creating ' neuralDataType ' spectrogram for file number ' num2str(fileNumber) ' of ' num2str(size(rawDataFileIDs, 1)) '...']); disp(' ')
    
    [Neural_S1, Neural_T1, Neural_F1] = mtspecgramc(RawNeuro, movingwin1, params);
    [Neural_S5, Neural_T5, Neural_F5] = mtspecgramc(RawNeuro, movingwin5, params);
    
    SpectrogramData.(neuralDataType).FiveSec.S{fileNumber, 1} = Neural_S5';
    SpectrogramData.(neuralDataType).FiveSec.T{fileNumber, 1} = Neural_T5;
    SpectrogramData.(neuralDataType).FiveSec.F{fileNumber, 1} = Neural_F5;
    SpectrogramData.(neuralDataType).FileIDs{fileNumber, 1} = fileID;
    SpectrogramData.Notes.params = params;
    SpectrogramData.Notes.movingwin5 = movingwin5;
    
    SpectrogramData.(neuralDataType).OneSec.S{fileNumber, 1} = Neural_S1';
    SpectrogramData.(neuralDataType).OneSec.T{fileNumber, 1} = Neural_T1;
    SpectrogramData.(neuralDataType).OneSec.F{fileNumber, 1} = Neural_F1;
    SpectrogramData.Notes.movingwin1 = movingwin1; 
end

save([animal '_SpectrogramData.mat'], 'SpectrogramData', '-v7.3');

end
