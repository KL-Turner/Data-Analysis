%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: 1) Track pupil diameter and detect periods of blinking
%            2) Patch any NaN values or droppped camera frames via interpolation
%            3) Manually check the first 5 and last 5 frames of each session to verify eye integrity/discharge
%            4) Manually check each blink for false positives
%            5) Extract resting pupil area and add it to RestData.mat structure
%            6) Determine baseline pupil area during rest and add it to RestingBaselines.mat
%            7) Extract whisking/stimulus triggered pupil area and add it to EventData.mat
%            8) Normalize RestData.mat and EventData.mat structures using resting baseline
%________________________________________________________________________________________________________________________

%% BLOCK PURPOSE: [0] Load the script's necessary variables and data structures.
% Clear the workspace variables and command windyow.
zap;
disp('Analyzing Block [0] Preparing the workspace and loading variables.'); disp(' ')
% Character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
[animalID,~,~] = GetFileInfo_IOS(procDataFileIDs(1,:));
%% BLOCK PURPOSE: [1] Track pupil diameter/area
disp('Analyzing Block [1] Tracking pupil diameter and blink detection.'); disp(' ')
TrackPupilDiameter_IOS(procDataFileIDs)
%% BLOCK PURPOSE: [2] Patch pupil area
disp('Analyzing Block [2]Patching NaN values and interpolating dropped frames.'); disp(' ')
for aa = 1:size(procDataFileIDs,1)
    disp(['Patching pupil area of file ' num2str(aa) '/' num2str(size(procDataFileIDs,1))]); disp(' ')
    PatchPupilArea_IOS(procDataFileIDs(aa,:))
end
%% BLOCK PURPOSE: [3] Check eye quality
disp('Analyzing Block [3] Manually check eye quality.'); disp(' ')
for bb = 1:size(procDataFileIDs,1)
    disp(['Manually checking eye of file ' num2str(bb) '/' num2str(size(procDataFileIDs,1))]); disp(' ')
    CheckPupilVideoFrames_IOS(procDataFileIDs(bb,:))
end
%% BLOCK PURPOSE: [4] Verify blinks
disp('Analyzing Block [4] Manually check blinks.'); disp(' ')
for cc = 1:size(procDataFileIDs)
    disp(['Manually checking blinks of file ' num2str(cc) '/' num2str(size(procDataFileIDs,1))]); disp(' ')
    CheckPupilBlinks_IOS(procDataFileIDs(cc,:))
end
%% BLOCK PURPOSE: [5] Add pupil area to RestData.mat
disp('Analyzing Block [5] Adding pupil area to RestData.mat'); disp(' ')
[RestData] = ExtractPupilRestingData_IOS(procDataFileIDs);
%% BLOCK PURPOSE: [6] Add pupil baseline to Restingbaselines.mat
disp('Analyzing Block [6] Adding pupil baseline to RestingBaselines.mat'); disp(' ')
[RestingBaselines] = AddPupilRestingBaseline_IOS();
%% BLOCK PURPOSE: [7] Add pupil area to EventData.mat
disp('Analyzing Block [7] Add pupil whisk/stim data to EventData.mat'); disp(' ')
[EventData] = ExtractPupilEventTriggeredData_IOS(procDataFileIDs);
%% BLOCK PURPOSE: [8] Normalize Rest/Event data structures
disp('Analyzing Block [8] Normalizing Rest/Event data structures.'); disp(' ')
[RestData] = NormBehavioralDataStruct_IOS(RestData,RestingBaselines,'manualSelection');
save([animalID '_RestData.mat'],'RestData','-v7.3')
[EventData] = NormBehavioralDataStruct_IOS(EventData,RestingBaselines,'manualSelection');
save([animalID '_EventData.mat'],'EventData','-v7.3')
%% fin
disp('Pupil Core Processing - Complete.'); disp(' ')
