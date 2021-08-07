clear; clc;
currentFolder = pwd;
fileparts = strsplit(currentFolder,filesep);
rootFolder = fullfile(fileparts{1:end});
animalIDs = {'T115','T116','T117','T118','T125','T156'};
modelType = 'Manual';
baselineType = 'manualSelection';
leadTime = 10;
lagTime = 10;
vesselSamplingRate = 5;
emgSamplingRate = 30;
[z,p,k] = butter(4,1/(vesselSamplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    Data.(animalID) = [];
    dataLocation = [rootFolder '/' animalID '/2P Data/'];
    cd(dataLocation)
    % find and load SleepData.mat strut
    sleepDataFileStruct = dir('*_SleepData.mat');
    sleepDataFile = {sleepDataFileStruct.name}';
    sleepDataFileID = char(sleepDataFile);
    load(sleepDataFileID)
    % load the baseline structure
    baselinesFileStruct = dir('*_RestingBaselines.mat');
    baselinesFile = {baselinesFileStruct.name}';
    baselinesFileID = char(baselinesFile);
    load(baselinesFileID)
    % NREM file IDs
    sleepFileIDs = unique(SleepData.Manual.NREM.FileIDs);
    for bb = 1:length(sleepFileIDs)
        sleepFileID = sleepFileIDs(bb,1);
        sleepFileID = sleepFileID{1,1};
        directory = pwd;      % Full path of the directory to be searched in
        filesAndFolders = dir(directory);     % Returns all the files and folders in the directory
        filesInDir = filesAndFolders(~([filesAndFolders.isdir]));  % Returns only the files in the directory
        fileNames = {filesInDir.name};
        stringToBeFound = sleepFileID;
        numOfFiles = length(filesInDir);
        fileFound = false;
        while fileFound == false
            for cc = 1:numOfFiles
                fileName = fileNames{1,cc};
                found = strfind(fileName,stringToBeFound);
                if isempty(found) == false
                    fileFound = true;
                    mergedDataFileID = fileName;
                    [~,~,fileDate,~,~,vesselID] = GetFileInfo2_2P(fileName);
                    strDay = ConvertDate_2P_eLife2020(fileDate);
                    load(mergedDataFileID,'-mat')
                    % vessel diameter
                    vesselDiameter = MergedData.data.vesselDiameter.data;
                    normVesselDiameter = (vesselDiameter - RestingBaselines.(baselineType).vesselDiameter.data.(vesselID).(strDay))./(RestingBaselines.(baselineType).vesselDiameter.data.(vesselID).(strDay));
                    filtVesselDiameter = filtfilt(sos,g,normVesselDiameter)*100;
                    % EMG
                    EMG = MergedData.data.EMG.data;
                    if isfield(Data.(animalID),vesselID) == false
                        Data.(animalID).(vesselID).timePoint = [];
                        Data.(animalID).(vesselID).vesselData = [];
                        Data.(animalID).(vesselID).emgData = [];
                        Data.(animalID).(vesselID).fileID ={};
                    end
                    break
                end
            end
        end
        [figHandle] = GenerateSingleFigures_2P_eLife2020(mergedDataFileID,baselineType,'n',RestingBaselines);
        useFile = input('Use file for micro-arousal transitions? (y/n): ','s');
        if strcmp(useFile,'y') == true
            keepSearching = 'y';
            while strcmp(keepSearching,'y') == true
                timePoint = input('Input time of micro-arousal (s): ');
                vesselVector = filtVesselDiameter((timePoint - leadTime)*vesselSamplingRate:(timePoint + lagTime)*vesselSamplingRate);
                emgVector = EMG((timePoint - leadTime)*emgSamplingRate:(timePoint + lagTime)*emgSamplingRate);
                Data.(animalID).(vesselID).timePoint = cat(1,Data.(animalID).(vesselID).timePoint,timePoint);
                Data.(animalID).(vesselID).vesselData = cat(1,Data.(animalID).(vesselID).vesselData,vesselVector);
                Data.(animalID).(vesselID).emgData = cat(1,Data.(animalID).(vesselID).emgData,emgVector);
                Data.(animalID).(vesselID).fileID = cat(1,Data.(animalID).(vesselID).fileID,mergedDataFileID);
                keepSearching = input('Keep searching? (y/n)','s');
            end
            close(figHandle)
        else
            close(figHandle)
        end
    end
end