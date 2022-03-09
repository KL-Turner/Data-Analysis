Fig2_Supplement
% Power Spectra (HbT, Gamma, Pupil + Statistics)
% XCorr Time-to-Peak

%% variables for loops
resultsStruct = 'Results_PowerSpectrum';
load(resultsStruct);
animalIDs = fieldnames(Results_PowerSpectrum);
behavFields = {'Rest','NREM','REM','Awake','Asleep','All'};
dataTypes = {'mmArea','mmDiameter','zArea','zDiameter'};
% pre-allocate data structure
for aa = 1:length(behavFields)
    behavField = behavFields{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        data.(behavField).(dataType).S = [];
        data.(behavField).(dataType).f = [];
    end
end
% power spectra during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(Results_PowerSpectrum.(animalID).(behavField).(dataType).S) == false
                data.(behavField).(dataType).S = cat(2,data.(behavField).(dataType).S,Results_PowerSpectrum.(animalID).(behavField).(dataType).S);
                data.(behavField).(dataType).f = cat(1,data.(behavField).(dataType).f,Results_PowerSpectrum.(animalID).(behavField).(dataType).f);
            end
        end
    end
end
% take mean/StD of S/f
for aa = 1:length(behavFields)
    behavField = behavFields{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        data.(behavField).(dataType).meanS = mean(data.(behavField).(dataType).S,2);
        data.(behavField).(dataType).stdS = std(data.(behavField).(dataType).S,0,2);
        data.(behavField).(dataType).meanf = mean(data.(behavField).(dataType).f,1);
    end
end