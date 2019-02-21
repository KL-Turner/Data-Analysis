function [RestingBaselines] = CalculateRestingBaselines_2P(animalID, mergedDataFiles, targetMinutes, RestData)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: This function finds the resting baseline for all fields of the RestData.mat structure, for each unique day
%________________________________________________________________________________________________________________________
%
%   Inputs: animal name (str) for saving purposes, targetMinutes (such as 30, 60, etc) which the code will interpret as 
%           the number of minutes/files at the beginning of each unique imaging day to use for baseline calculation. 
%           RestData.mat, which should have field names such as CBV, Delta Power, Gamma Power, etc. with ALL resting events
%
%   Outputs: SleepRestEventData.mat struct
%
%   Last Revision: October 4th, 2018
%________________________________________________________________________________________________________________________

% The RestData.mat struct has all resting events, regardless of duration. We want to set the threshold for rest as anything
% that is greater than 10 seconds.
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {10};

for f = 1:size(mergedDataFiles)
    mergedDataFile = mergedDataFiles(f, :);
    [animalID, fileDate, ~, vesselID] = GetFileInfo_2P(mergedDataFile);
    load(mergedDataFile)
    strDays{f, :} = ConvertDate(fileDate);
    vesselIDs{f,:} = vesselID;
end

uniqueVessels = unique(vesselIDs);   % Total number of unique days in this folder under a single animal
uniqueDays = unique(strDays);   % Total number of unique vessels in this folder under a single animal
dataTypes = fieldnames(RestData);   % List of fields in RestData that corresponds to the difference baselines we want to find

% Find the fieldnames of RestData and loop through each field. Each fieldname should be a different dataType of interest.
% These will typically be CBV, Delta, Theta, Gamma, and MUA
for u = 1:length(uniqueVessels)
    uV = uniqueVessels{u};
    
    for v = 1:length(uniqueDays)
        uD = uniqueDays{v};
        
        for dT = 1:length(dataTypes)
            dataType = char(dataTypes(dT));   % Load each loop iteration's fieldname as a character string
            [restLogical] = FilterEvents(RestData.(dataType), RestCriteria);   % RestData output is a logical
            allRestFiles = RestData.(dataType).fileIDs(restLogical, :);
            allRestDurations = RestData.(dataType).durations(restLogical, :);
            allRestEventTimes = RestData.(dataType).eventTimes(restLogical, :);
            allRestVesselIDs = RestData.(dataType).vesselIDs(restLogical, :);
            allRestData = RestData.(dataType).data(restLogical, :);

            vesselFilterLogical = zeros(length(allRestVesselIDs), 1);   % Find all the vessels that correspond to this loop's vessel (A1, A2, A3, etc)
            for vF = 1:length(allRestVesselIDs)
                vessel = char(allRestVesselIDs(vF, 1));
                if strcmp(uV, vessel)
                    vesselFilterLogical(vF, 1) = 1;
                end
            end

            vesselFilterLogical = logical(vesselFilterLogical);
            singleVesselFiles = allRestFiles(vesselFiltLogical);
            singleVesselDurations = allRestDurations(vesselFiltLogical);
            singleVesselEventTimes = allRestEventTimes(vesselFiltLogical);
            singleVesselIDs = allRestVesselIDs(vesselFiltLogical);
            singleVesselData = allRestData(vesselFiltLogical);
            
            allSingleVesselRestFiles = allRestFiles(vesselFilterLogical);
            restDurations = allRestDurations(vesselFilterLogical);
            restEventTimes = allRestEventTimes(vesselFilterLogical);
            restData = allRestData(vesselFilterLogical);

            % find the means of each unique day
            for x = 1:length(restData)
                tempData_means(1,x) = mean(restData{x});
            end

            RestingBaselines.(uV).(dataType).baseLine = mean(tempData_means);
            RestingBaselines.(uV).(dataType).restData = restData;
            RestingBaselines.(uV).(dataType).fileIDs = allSingleVesselRestFiles;
            RestingBaselines.(uV).(dataType).eventTimes = restEventTimes;
            RestingBaselines.(uV).(dataType).durations = restDurations;
        end
end
    save([animalID '_RestingBaselines.mat'], 'RestingBaselines');
    
end
