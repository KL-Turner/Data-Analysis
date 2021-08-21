clear; clc;
microArousalFileID = 'MicroArousalData.mat';
load(microArousalFileID)
modelType = 'Manual';
baselineType = 'manualSelection';
leadTime = 20;
lagTime = 20;
vesselSamplingRate = 5;
dsFs = 30;
currentFolder = pwd;
fileparts = strsplit(currentFolder,filesep);
rootFolder = fullfile(fileparts{1:end});
animalIDs = {'T115','T116','T117','T118'};
[z,p,k] = butter(4,1/(vesselSamplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
[z1,p1,k1] = butter(4,0.1/(dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    dataLocation = [rootFolder '/' animalID '/2P Data/'];
    cd(dataLocation)
    baselinesFileStruct = dir('*_RestingBaselines.mat');
    baselinesFile = {baselinesFileStruct.name}';
    baselinesFileID = char(baselinesFile);
    load(baselinesFileID)
    vIDs = fieldnames(Data.(animalIDs{1,aa}));
    for bb = 1:length(vIDs)
        vID = vIDs{bb,1};
        for cc = 1:size(Data.(animalIDs{1,aa}).(vID).timePoint,1)
            timePoint = Data.(animalID).(vID).timePoint(cc,1);
            fileID = Data.(animalID).(vID).fileID(cc,1);
            fileID = fileID{1,1};
            [~,~,fileDate,~,~,~] = GetFileInfo2_2P(fileID);
            strDay = ConvertDate_2P_eLife2020(fileDate);
            load(fileID,'-mat')
            vesselDiameter = MergedData.data.vesselDiameter.data;
            normVesselDiameter = (vesselDiameter - RestingBaselines.(baselineType).vesselDiameter.data.(vID).(strDay))./(RestingBaselines.(baselineType).vesselDiameter.data.(vID).(strDay));
            filtVesselDiameter = filtfilt(sos,g,normVesselDiameter)*100;
            EMG = MergedData.data.EMG.data;
            deltaPower = MergedData.data.corticalNeural.deltaBandPower;
            normDeltaPower = (deltaPower - RestingBaselines.(baselineType).corticalNeural.deltaBandPower.(strDay))./(RestingBaselines.(baselineType).corticalNeural.deltaBandPower.(strDay));
            filtDeltaPower = filtfilt(sos1,g1,normDeltaPower)*100;
            vesselVec = filtVesselDiameter(floor((timePoint - leadTime)*vesselSamplingRate):floor((timePoint + lagTime)*vesselSamplingRate));
            emgVec = EMG(floor((timePoint - leadTime)*dsFs):floor((timePoint + lagTime)*dsFs));
            deltaVec = filtDeltaPower(floor((timePoint - leadTime)*dsFs):floor((timePoint + lagTime)*dsFs));
            Data.(animalID).(vID).vesselData_long(cc,:) = vesselVec(1:(leadTime + lagTime)*vesselSamplingRate);
            Data.(animalID).(vID).emgData_long(cc,:) = emgVec(1:(leadTime + lagTime)*dsFs);
            clear Data.(animalID).(vID).delta_long
            Data.(animalID).(vID).deltaData_long(cc,:) = deltaVec(1:(leadTime + lagTime)*dsFs);
        end
    end
end
cd(currentFolder)
save('MicroArousalData.mat','Data');
