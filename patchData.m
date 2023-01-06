% find and load RestData.mat struct
restDataFileStruct = dir('*_RestData.mat');
restDataFile = {restDataFileStruct.name}';
restDataFileID = char(restDataFile);
load(restDataFileID,'-mat')
RestData.CBV_HbT.LH = RestData.CBV_HbT.adjLH;
RestData.CBV_HbT.RH = RestData.CBV_HbT.adjRH;
RestData.CBV.LH = RestData.CBV.adjLH;
RestData.CBV.RH = RestData.CBV.adjRH;
RestData.CBV = rmfield(RestData.CBV,'adjLH');
RestData.CBV = rmfield(RestData.CBV,'adjRH');
RestData.CBV_HbT = rmfield(RestData.CBV_HbT,'adjLH');
RestData.CBV_HbT = rmfield(RestData.CBV_HbT,'adjRH');
save(restDataFileID,'RestData')

% find and load EventData.mat struct
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID,'-mat')
EventData.CBV_HbT.LH = EventData.CBV_HbT.adjLH;
EventData.CBV_HbT.RH = EventData.CBV_HbT.adjRH;
EventData.CBV.LH = EventData.CBV.adjLH;
EventData.CBV.RH = EventData.CBV.adjRH;
EventData.CBV = rmfield(EventData.CBV,'adjLH');
EventData.CBV = rmfield(EventData.CBV,'adjRH');
EventData.CBV_HbT = rmfield(EventData.CBV_HbT,'adjLH');
EventData.CBV_HbT = rmfield(EventData.CBV_HbT,'adjRH');
save(eventDataFileID,'EventData','-v7.3')

% find and load RestingBaselines.mat strut
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
RestingBaselines.manualSelection.CBV.LH = RestingBaselines.manualSelection.CBV.adjLH;
RestingBaselines.manualSelection.CBV.RH = RestingBaselines.manualSelection.CBV.adjRH;
RestingBaselines.manualSelection.CBV = rmfield(RestingBaselines.manualSelection.CBV,'adjLH');
RestingBaselines.manualSelection.CBV = rmfield(RestingBaselines.manualSelection.CBV,'adjRH');
save(baselineDataFileID,'RestingBaselines')

procDataFileIDs = ls('*ProcData.mat');
for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    load(procDataFileID);
    ProcData.data.CBV.LH = ProcData.data.CBV.adjLH;
    ProcData.data.CBV.RH = ProcData.data.CBV.adjRH;

    ProcData.data.CBV_HbT.LH = ProcData.data.CBV_HbT.adjLH;
    ProcData.data.CBV_HbT.RH = ProcData.data.CBV_HbT.adjRH;

    ProcData.data.CBV = rmfield(ProcData.data.CBV,'adjLH');
    ProcData.data.CBV = rmfield(ProcData.data.CBV,'adjRH');

    ProcData.data.CBV_HbT = rmfield(ProcData.data.CBV_HbT,'adjLH');
    ProcData.data.CBV_HbT = rmfield(ProcData.data.CBV_HbT,'adjRH');

    try
        ProcData.data.CBV = rmfield(ProcData.data.CBV,'LH_Cement');
        ProcData.data.CBV = rmfield(ProcData.data.CBV,'RH_Cement');
    catch
    end

    save(procDataFileID,'ProcData')
end
