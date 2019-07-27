
% Character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);

AddSleepParameters_SVM(procDataFileIDs, RestingBaselines)
CreateTrainingDataSet_SVM(procDataFileIDs, RestingBaselines)