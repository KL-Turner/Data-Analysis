function [SpectrogramData] = CreateTrialSpectrograms_2P(animalID, mergedDataFiles, SpectrogramData)
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

for fileNumber = 1:size(mergedDataFiles, 1)
    mergedDataFileID = mergedDataFiles(fileNumber, :);
    load(mergedDataFileID);
    duration = MergedData.Notes.trialDuration_Sec;
    analogFs = MergedData.Notes.MScan.MScan_analogSamplingRate;
    expectedLength = duration*analogFs;
    [~, ~, fileID, vesselID] = GetFileInfo_2P(mergedDataFileID);
    RawNeuro = detrend(MergedData.Data.Raw_NeuralData(1:expectedLength), 'constant');
    
    w0 = 60/(analogFs/2);  bw = w0/35;
    [num,den] = iirnotch(w0, bw);
    RawNeuro2 = filtfilt(num, den, RawNeuro);

    % Spectrogram parameters
    params.tapers = [1 1];
    params.Fs = MergedData.Notes.MScan.MScan_analogSamplingRate;
    params.fpass = [1 100];
    movingwin1 = [1 1/10];
    movingwin5 = [5 1/10];

    disp(['Creating spectrogram for file number ' num2str(fileNumber) ' of ' num2str(size(mergedDataFiles, 1)) '...']); disp(' ')
    
    [Neural_S1, Neural_T1, Neural_F1] = mtspecgramc(RawNeuro2, movingwin1, params);
    [Neural_S5, Neural_T5, Neural_F5] = mtspecgramc(RawNeuro2, movingwin5, params);
    
    SpectrogramData.FiveSec.S{fileNumber, 1} = Neural_S5';
    SpectrogramData.FiveSec.T{fileNumber, 1} = Neural_T5;
    SpectrogramData.FiveSec.F{fileNumber, 1} = Neural_F5;
    SpectrogramData.VesselIDs{fileNumber, 1} = vesselID;
    SpectrogramData.FileIDs{fileNumber, 1} = fileID;
    SpectrogramData.Notes.params = params;
    SpectrogramData.Notes.movingwin5 = movingwin5;
    
    SpectrogramData.OneSec.S{fileNumber, 1} = Neural_S1';
    SpectrogramData.OneSec.T{fileNumber, 1} = Neural_T1;
    SpectrogramData.OneSec.F{fileNumber, 1} = Neural_F1;
    SpectrogramData.Notes.movingwin1 = movingwin1; 
end

for a = 1:length(SpectrogramData.FiveSec.F)
    

save([animalID '_SpectrogramData.mat'], 'SpectrogramData', '-v7.3');

end
