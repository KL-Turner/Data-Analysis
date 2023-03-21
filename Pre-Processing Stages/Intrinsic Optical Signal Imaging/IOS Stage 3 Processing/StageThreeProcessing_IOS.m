%________________________________________________________________________________________________________________________
% Written by Kevin L. Turnery
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: - Categorize behavioral (rest,whisk,stim) data using previously processed data structures, add 'flags'
%          - Create a temporary RestData structure that contains periods of rest - use this for initial figures
%          - Analyze neural data and create different spectrograms for each file's electrodes
%          - Uses periods when animal is not being stimulated or moving to establish an initial baseline
%          - Manually select awake files for a slightly different baseline not based on hard time vals
%          - Use the best baseline to convert reflectance changes to total hemoglobin
%          - Re-create the RestData structure now that we can deltaHbT
%          - Create an EventData structure looking at the different data types after whisking or stimulation
%          - Apply the resting baseline to each data type to create a percentage change
%          - Use the time indeces of the resting baseline file to apply a percentage change to the spectrograms
%          - Use the time indeces of the resting baseline file to create a reflectance pixel-based baseline
%          - Generate a summary figure for all of the analyzed and processed data
%________________________________________________________________________________________________________________________

zap;
% character list of all RawData files
rawDataFileStruct = dir('*_RawData.mat');
rawDataFiles = {rawDataFileStruct.name}';
rawDataFileIDs = char(rawDataFiles);
% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% categorize data
CategorizeData_IOS(procDataFileIDs)
% create RestData data structure
[RestData] = ExtractRestingData_IOS(procDataFileIDs,1);
% analyze the spectrogram for each session.
CreateTrialSpectrograms_IOS(rawDataFileIDs);
% create Baselines data structure
[RestingBaselines] = CalculateRestingBaselines_IOS(RestData,60);
% create Baselines for spectrogram data
[RestingBaselines] = CalculateSpectrogramBaselines_IOS(RestingBaselines,'setDuration');
% normalize spectrogram by baseline
NormalizeSpectrograms_IOS(RestingBaselines);
% manually select files for custom baseline calculation
[RestingBaselines] = CalculateManualRestingBaselinesTimeIndeces_IOS(procDataFileIDs,RestData,RestingBaselines,'reflectance');
% add delta HbT field to each processed data file
UpdateTotalHemoglobin_IOS(procDataFileIDs,RestingBaselines,'manualSelection')
% correct GCaMP attenuation
CorrectGCaMPattenuation_IOS(procDataFileIDs,RestingBaselines)
% re-create the RestData structure now that HbT (and/or corrected GCaMP) is available
[RestData] = ExtractRestingData_IOS(procDataFileIDs,2);
% create the EventData structure for CBV and neural data
[EventData] = ExtractEventTriggeredData_IOS(procDataFileIDs);
% normalize RestData structures by the resting baseline
[RestData] = NormRestDataStruct_IOS(RestData,RestingBaselines,'manualSelection');
% normalize EventData structures by the resting baseline
[EventData] = NormEventDataStruct_IOS(EventData,RestingBaselines,'manualSelection');
% find spectrogram baselines for each day
[RestingBaselines] = CalculateSpectrogramBaselines_IOS(RestingBaselines,'manualSelection');
% normalize spectrogram by baseline
NormalizeSpectrograms_IOS(RestingBaselines);
% create a structure with all spectrograms for convenient analysis further downstream
CreateAllSpecDataStruct_IOS()
%
% generate single trial figures
GenerateTrialFigures_IOS(procDataFileIDs,RestingBaselines);
