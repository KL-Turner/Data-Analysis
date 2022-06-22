function AnalizeFiberPhotometryData_FP
%READ ME
%This function will return a completed dataset from raw .csv files
%collected using Doric fiberphotometry data acquisition systems
%
%%Filenaming Format:
%%'AnimalName_Session1.csv'
%
%%Data Storage Format:
%%'I:\BilateralFiberPhotoMetry\SK_Bilaterl_FBR\GCaMP7s\T235\091721\Filename.csv'%

%Embeded Functions
%Calibrate_Correction
%RawExtraction

Home='I:\BilateralFiberPhotoMetry';
GFPDir='I:\BilateralFiberPhotoMetry\SK_Bilaterl_FBR\GCaMP7s';
OpticalChannelNames={'RHControl','RHGCaMP7s','RHBloodVolume','LHControl','LHGCaMP7s','LHBloodVolume','LH560Raw'};
popDir='I:\BilateralFiberPhotoMetry\PopulationData';
%% Get correction for Hemodynamic attenuation of GFP signal
cd(GFPDir);

if exist('Hemodynamic_Correction_Constant.mat','file')
    load('Hemodynamic_Correction_Constant.mat');
else
    AnimalFolders=dir;
    AnimalFolders(~[AnimalFolders.isdir])=[];
    tf=ismember({AnimalFolders.name},{'.','..'});
    AnimalFolders(tf)=[];
    for foldNum=1:size(AnimalFolders)
        cd([AnimalFolders(foldNum).folder '\' AnimalFolders(foldNum).name]);
        TheFiles= ls('*.csv');
        for fileNum=1:size(TheFiles,1)
            if ~strcmpi(TheFiles(fileNum).folder,'I:\BilateralFiberPhotoMetry\SK_Bilaterl_FBR\T235\091121')
                cd(TheFiles(fileNum).folder);
                filename=TheFiles(fileNum).name;
                [coeffVals_RH,coeffVals_LH]=Calibrate_Correction(filename);
                FitData(fileNum).coeffVals_RH=coeffVals_RH;
                FitData(fileNum).coeffVals_LH=coeffVals_LH;
                Slope_RH(fileNum)=FitData(fileNum).coeffVals_RH(1,1);
                Slope_LH(fileNum)=FitData(fileNum).coeffVals_LH(1,1);
            end
            if fileNum==size(TheFiles,1)
                Slope_RH(Slope_RH>=0)=NaN;
                CorrectionConst.RH=nanmean(Slope_RH);
                Slope_LH(Slope_LH>=0)=NaN;
                CorrectionConst.LH=nanmean(Slope_LH);
            end
        end
end
cd(GFPDir);
save('Hemodynamic_Correction_Constant','CorrectionConst');
end

cd(Home);
subfolders=dir;
subfolders(~[subfolders.isdir])=[];
tf=ismember({subfolders.name},{'.','..','Ignore','PopulationData'});
subfolders(tf)=[];
%% Process raw .csv files
for folderNum=1:size(subfolders,1)
    cd([subfolders(folderNum).folder '\' subfolders(folderNum).name]);
    FileList=dir(fullfile([subfolders(folderNum).folder '\' subfolders(folderNum).name],'**','*.csv'));
    for filNum=1:size(FileList,1)
        cd(FileList(filNum).folder);
        filename=FileList(filNum).name;
        correctFlag='y'; %If GCaMP signals need to be corrected for hemodynamic attenuation =='y' 
        RawExtraction(filename,OpticalChannelNames,correctFlag);
    end
end
close all;

