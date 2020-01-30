function CreateTrialSpectrograms_IOS(rawDataFiles, neuralDataTypes)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyzes the raw neural data from each RawData.mat file and calculates two different spectrograms.
%________________________________________________________________________________________________________________________
%
%   Inputs: List of RawData.mat files.
%
%   Outputs: Saves a SpecData.mat file of the same name as RawData containing the analysis.
%
%   Last Revised: February 21st, 2019
%________________________________________________________________________________________________________________________

for a = 1:size(rawDataFiles, 1)
    rawDataFile = rawDataFiles(a, :);
    [animalID, ~, fileID] = GetFileInfo_IOS(rawDataFile);
    specDataFileID = [animalID '_' fileID '_SpecData.mat'];
    
    if ~exist(specDataFileID, 'file') == true
        load(rawDataFile);
        duration = RawData.notes.trialDuration_sec;
        analogFs = RawData.notes.analogSamplingRate;
        expectedLength = duration*analogFs;
        
        for b = 1:length(neuralDataTypes)
            neuralDataType = neuralDataTypes{1,b};
            try
                rawNeuro = detrend(RawData.data.(neuralDataType)(1:expectedLength), 'constant');
            catch
                sampleDiff = expectedLength - length(RawData.data.(neuralDataType));
                rawNeuro = detrend(horzcat(RawData.data.(neuralDataType), RawData.data.(neuralDataType)(end)*ones(1, sampleDiff)), 'constant');
            end
            
%             w0 = 60/(analogFs/2);
%             bw = w0/35;
%             [num,den] = iirnotch(w0, bw);
%             rawNeuro2 = filtfilt(num, den, rawNeuro);
            
            % Spectrogram parameters
            params1.tapers = [1 1];
            params1.Fs = analogFs;
            params1.fpass = [1 100];
            movingwin1 = [1 1/10];
            
            params5.tapers = [5 9];
            params5.Fs = analogFs;
            params5.fpass = [1 100];
            movingwin5 = [5 1/5];
            
            disp(['Creating ' neuralDataType ' spectrogram for file number ' num2str(a) ' of ' num2str(size(rawDataFiles, 1)) '...']); disp(' ')
            
            [S1, T1, F1] = mtspecgramc_IOS(rawNeuro, movingwin1, params1);
            [S5, T5, F5] = mtspecgramc_IOS(rawNeuro, movingwin5, params5);
            
            SpecData.(neuralDataType).fiveSec.S = S5';
            SpecData.(neuralDataType).fiveSec.T = T5;
            SpecData.(neuralDataType).fiveSec.F = F5;
            SpecData.(neuralDataType).fiveSec.params = params5;
            SpecData.(neuralDataType).fiveSec.movingwin = movingwin5;
            
            SpecData.(neuralDataType).oneSec.S = S1';
            SpecData.(neuralDataType).oneSec.T = T1;
            SpecData.(neuralDataType).oneSec.F = F1;
            SpecData.(neuralDataType).oneSec.params = params1;
            SpecData.(neuralDataType).oneSec.movingwin = movingwin1;
            
            save(specDataFileID, 'SpecData');
        end
    else
        disp([specDataFileID ' already exists. Continuing...']); disp(' ')
    end
end

end
