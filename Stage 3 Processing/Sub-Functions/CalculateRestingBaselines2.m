function [RestingBaselines] = CalculateRestingBaselines2(animalID, mergedDataFiles, targetMinutes, RestData)
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
RestCriteria.Value = {5};

for f = 1:size(mergedDataFiles)
    mergedDataFile = mergedDataFiles(f, :);
    [animalID, fileDate, fileID, imageID] = GetFileInfo2(mergedDataFile);
    load(mergedDataFile)
    strDay = ConvertDate(fileDate);
    vesselIDs{f,:} = MergedData.Notes.MScan.vesselID;
end

uVessels = unique(vesselIDs);

% Find the fieldnames of RestData and loop through each field. Each fieldname should be a different dataType of interest.
% These will typically be CBV, Delta, Theta, Gamma, and MUA
for u = 1:length(uVessels)
    uV = uVessels{u};
    for f = 1:length(vesselIDs)
        vesselID = vesselIDs{f};
        if strcmp(uV, vesselID)
            fileLogical(f,1) = 1;
        else
            fileLogical(f,1) = 0;
        end
    end
    
    fileLogical = logical(fileLogical);
    fileBreaks = strfind(mergedDataFiles(1, :), '_');
    uV_Files = mergedDataFiles(fileLogical, fileBreaks(1) + 1:fileBreaks(5) - 1);
    
    dataTypes = fieldnames(RestData);
    for dT = 1:length(dataTypes)
        dataType = char(dataTypes(dT));   % Load each loop iteration's fieldname as a character string
        
        % Use the RestCriteria we specified earlier to find all resting events that are greater than the criteria
        [restLogical] = FilterEvents(RestData.(dataType), RestCriteria);   % Output is a logical
        allRestDurations = RestData.(dataType).durations(restLogical, :);
        allRestEventTimes = RestData.(dataType).eventTimes(restLogical, :);
        allRestData = RestData.(dataType).data(restLogical, :);   % Pull out data from all those resting files that meet criteria
     
        allRestFiles = RestData.(dataType).fileIDs(restLogical, :);   % Overall logical for all resting file names that meet criteria
        vesselFiltLogical = zeros(length(allRestFiles), 1);
        for a = 1:length(allRestFiles)
            restFile = char(allRestFiles(a,1));
            for b = 1:size(uV_Files)
                uV_File = uV_Files(b, :);
                if strcmp(restFile, uV_File)
                    vesselFiltLogical(a, 1) = 1;
                end
            end
        end
        
        vesselFiltLogical = logical(vesselFiltLogical);
        restFiles = allRestFiles(vesselFiltLogical);
        restDurations = allRestDurations(vesselFiltLogical);
        restEventTimes = allRestEventTimes(vesselFiltLogical);
        restData = allRestData(vesselFiltLogical);
            
        % find the means of each unique day
        for x = 1:length(restData)
            tempData_means(1,x) = mean(restData{x});
        end
        
        RestingBaselines.(uV).(dataType).baseLine = mean(tempData_means);
        RestingBaselines.(uV).(dataType).restData = restData;
        RestingBaselines.(uV).(dataType).fileIDs = restFiles;
        RestingBaselines.(uV).(dataType).eventTimes = restEventTimes;
        RestingBaselines.(uV).(dataType).durations = restDurations;
    end
end
    save([animalID '_RestingBaselines.mat'], 'RestingBaselines');
    
end
