%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: 1) Categorize data using previously processed MergedData data structures, add 'flags'  
%            2) Create RestData structure that contains periods of rest.
%            3) Create EventData structure that contains periods after stimuli and whisks.
%            4) Uses periods when animal is not being stimulated or moving to establish a 
%               baseline for a given session of imaging.
%            5) Normalizes the different data structures.
%________________________________________________________________________________________________________________________
%
%   Inputs: 1) Select all _MergedData files from all days. Follow the command window prompts.
%           2) Select one single _RawData file for the animal information. The MergedData files are
%              already in the list and will be used to run the function.
%           3) No inputs. MergedData files already loaded.
%           4) No inputs. RestData.mat is already loaded.
%           5) No inputs. RestData.mat and EventData.mat are already loaded.
%
%   Outputs: 1) Additions to the MergedData structure including flags and scores.
%            2) A RestData.mat structure with periods of rest.
%            3) A EventData.mat structure with event-related information.
%            4) Baselines.mat containing the baselines for individual resting periods.
%            5) Creates NormData in the rest/event structures.
%
%   Last Revised: October 5th, 2018
%________________________________________________________________________________________________________________________

combDirectory = dir('*_MergedData.mat');
MergedDataFiles = {combDirectory.name}';
MergedDataFiles = char(MergedDataFiles);
animalID = 'DK12';
dataTypes = {'Vessel_Diameter', 'DeltaBand_Power', 'ThetaBand_Power', 'AlphaBand_Power', 'BetaBand_Power', 'GammaBand_Power', 'MUA_Power'};

%% BLOCK PURPOSE: [1] Categorize data 
disp('Analyzing Block [1] Categorizing data.'); disp(' ')
for fileNumber = 1:size(MergedDataFiles, 1)
    fileName = MergedDataFiles(fileNumber, :);
    disp(['Analyzing file ' num2str(fileNumber) ' of ' num2str(size(MergedDataFiles, 1)) '...']); disp(' ')
    CategorizeData2(fileName)
end

%% BLOCK PURPOSE: [2] Create RestData data structure
disp('Analyzing Block [2] Create RestData struct for CBV and neural data.'); disp(' ')
[RestData] = ExtractRestingData2(MergedDataFiles, dataTypes);
    
%% BLOCK PURPOSE: [3] Create EventData data structure
disp('Analyzing Block [3] Create EventData struct for CBV and neural data.'); disp(' ')
[EventData] = ExtractEventTriggeredData2(MergedDataFiles, dataTypes);

%% BLOCK PURPOSE: [4] Create Baselines data structure
disp('Analyzing Block [4] Create Baselines struct for CBV and neural data.'); disp(' ')
targetMinutes = 15;
disp(['Calculating the resting baselines for the first ' num2str(targetMinutes) ' minutes of each unique day...']);
[RestingBaselines] = CalculateRestingBaselines2(animalID, MergedDataFiles, targetMinutes, RestData);

% %% BLOCK PURPOSE: [5] Normalize RestData and behavioral data
% disp('Analyzing Block [5] Normalize EventData struct using Baselines for CBV and neural data.'); disp(' ')
% [RestData] = NormBehavioralDataStruct(RestData, RestingBaselines);
% [EventData] = NormBehavioralDataStruct(EventData, RestingBaselines);
% 
% save([animal '_RestData.mat'], 'RestData')
% save([animal '_EventData.mat'], 'EventData')
% 
% disp('Stage Three Processing - Complete.'); disp(' ')
