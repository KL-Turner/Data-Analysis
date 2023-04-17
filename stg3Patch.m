zap
% character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% edit fieldname
for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    load(procDataFileID);
    disp(num2str(aa))
    ProcData.data.CBV.fLH = ProcData.data.CBV.frontalLH;
    ProcData.data.CBV.fRH = ProcData.data.CBV.frontalRH;
    ProcData.data.CBV = rmfield(ProcData.data.CBV,'Cement');
    ProcData.data.CBV = rmfield(ProcData.data.CBV,'frontalLH');
    ProcData.data.CBV = rmfield(ProcData.data.CBV,'frontalRH');
    ProcData.data = rmfield(ProcData.data,'CBV_HbT');
    ProcData.data = rmfield(ProcData.data,'GCaMP7s');
    ProcData.data = rmfield(ProcData.data,'Deoxy');
    ProcData.data = rmfield(ProcData.data,'stimulationsOriginal');
    save(procDataFileID,'ProcData')
end
