function [RestingBaselines] = CalculatePixelWiselRestingBaselines_IOS(procDataFileIDs,RestingBaselines,runDecision)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Manually designate files with event times that correspond to appropriate rest
%________________________________________________________________________________________________________________________

if strcmp(runDecision,'y')
    [animalID,~,~] = GetFileInfo_IOS(procDataFileIDs(1,:));
    restFileList = unique(RestingBaselines.manualSelection.baselineFileInfo.fileIDs); % obtain the list of unique fileIDs
    redFrames = cell(size(restFileList,1),1);
    greenFrames = cell(size(restFileList,1),1);
    blueFrames = cell(size(restFileList,1),1);
    % obtain the information from all the resting files
    for bb = 1:length(restFileList)
        fileID = restFileList{bb,1}; % FileID of currently loaded file
        load([animalID '_' fileID '_ProcData.mat']);
        samplingRate = ProcData.notes.CBVCamSamplingRate;
        strDay = ConvertDate_IOS(fileID(1:6));
        imageWidth = ProcData.notes.CBVCamPixelWidth;
        imageHeight = ProcData.notes.CBVCamPixelHeight;
        % load in neural data from current file
        disp(['Reading frames from file ' num2str(bb)]); disp(' ');
        [frames] = ReadDalsaBinary_IOS(animalID,[fileID '_WindowCam.bin']);
        if ProcData.notes.greenFrames == 1
            greenFrames = frames(1:3:end - 1);
            blueFrames = frames(2:3:end);
            redFrames = frames(3:3:end);
        elseif ProcData.notes.greenFrames == 2
            greenFrames = frames(2:3:end);
            blueFrames = frames(3:3:end);
            redFrames = frames(1:3:end - 1);
        elseif ProcData.notes.greenFrames == 3
            greenFrames = frames(3:3:end);
            blueFrames = frames(1:3:end - 1);
            redFrames = frames(2:3:end);
        end
        redImageStack = reshape(cell2mat(redFrames),imageWidth,imageHeight,length(redFrames));
        greenImageStack = reshape(cell2mat(greenFrames),imageWidth,imageHeight,length(greenFrames));
        blueImageStack = reshape(cell2mat(blueFrames),imageWidth,imageHeight,length(blueFrames));
        restRedFileData = [];
        restGreenFileData = [];
        restBlueFileData = [];
        for cc = 1:length(RestingBaselines.manualSelection.baselineFileInfo.fileIDs)
            restFileID = RestingBaselines.manualSelection.baselineFileInfo.fileIDs{cc,1};
            if strcmp(fileID,restFileID)
                restDuration = round(RestingBaselines.manualSelection.baselineFileInfo.durations(cc,1),1);
                startTime = round(RestingBaselines.manualSelection.baselineFileInfo.eventTimes(cc,1),1);
                % conditions and indexing
                startTimeIndex = startTime*samplingRate;
                restDurationIndex = restDuration*samplingRate - 1;
                restRedEventData = redImageStack(:,:,(startTimeIndex:(startTimeIndex + restDurationIndex)));
                restGreenEventData = greenImageStack(:,:,(startTimeIndex:(startTimeIndex + restDurationIndex)));
                restBlueEventData = blueImageStack(:,:,(startTimeIndex:(startTimeIndex + restDurationIndex)));
                if sum(sum(isnan(restRedEventData))) == 0
                    restRedFileData = cat(3,restRedFileData,restRedEventData);
                end
                if sum(sum(isnan(restGreenEventData))) == 0
                    restGreenFileData = cat(3,restGreenFileData,restGreenEventData);
                end
                if sum(sum(isnan(restBlueEventData))) == 0
                    restBlueFileData = cat(3,restBlueFileData,restBlueEventData);
                end
            end
        end
        trialRestData.([strDay '_' fileID]).red = restRedFileData;
        trialRestData.([strDay '_' fileID]).green = restGreenFileData;
        trialRestData.([strDay '_' fileID]).blue = restBlueFileData;
    end
    fields = fieldnames(trialRestData);
    uniqueDays = GetUniqueDays_IOS(RestingBaselines.manualSelection.baselineFileInfo.fileIDs);
    for f = 1:length(uniqueDays)
        g = 1;
        stringDay = ConvertDate_IOS(uniqueDays{f});
        frameAvgs.red.(stringDay) = [];
        frameAvgs.green.(stringDay) = [];
        frameAvgs.blue.(stringDay) = [];
        for field = 1:length(fields)
            if strcmp(fields{field}(7:12),uniqueDays{f})
                frameAvgs.red.(stringDay) = cat(3,frameAvgs.red.(stringDay),trialRestData.(fields{field}).red);
                frameAvgs.green.(stringDay) = cat(3,frameAvgs.green.(stringDay),trialRestData.(fields{field}).green);
                frameAvgs.blue.(stringDay)= cat(3,frameAvgs.blue.(stringDay),trialRestData.(fields{field}).blue);
                g = g + 1;
            end
        end
    end
    dayFields = fieldnames(frameAvgs.red);
    for h = 1:length(dayFields)
        disp(['Adding pixel-wise baseline to baseline file for ' dayFields{h} '...']); disp(' ')
        RestingBaselines.Pixel.red.(dayFields{h}) = mean(frameAvgs.red.(dayFields{h}),3);
        RestingBaselines.Pixel.green.(dayFields{h}) = mean(frameAvgs.green.(dayFields{h}),3);
        RestingBaselines.Pixel.blue.(dayFields{h}) = mean(frameAvgs.blue.(dayFields{h}),3);
    end
    save([animalID '_RestingBaselines.mat'],'RestingBaselines');
end
