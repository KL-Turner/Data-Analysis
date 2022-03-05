clear; clc; close all;
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
qq = 1;
msFrameStack = [];
for aa = 1:size(procDataFileIDs)
    procDataFileID = procDataFileIDs(aa,:);
    load(procDataFileID)
    %% extract and condense stimulation times
    stimulationTimes = sort(cat(2,ProcData.data.stimulations.LPadSol,ProcData.data.stimulations.RPadSol),'ascend');
    condensedStimulationTimes = [];
    cc = 1;
    for bb = 1:length(stimulationTimes)
        if bb == 1
            condensedStimulationTimes(1,bb) = stimulationTimes(1,bb); %#ok<*SAGROW>
            cc = cc + 1;
        else
            timeDifference = stimulationTimes(1,bb) - stimulationTimes(1,bb - 1);
            if timeDifference > 1 % remove stimulations that are closer than 1 second to the previous
                condensedStimulationTimes(1,cc) = stimulationTimes(1,bb);
                cc = cc + 1;
            end
        end
    end
    %% extra movie frames
    [~,~,fileID] = GetFileInfo_IOS(procDataFileID);
    windowCamFileID = [fileID '_WindowCam.bin'];
    Fs = 30;
    imageHeight = 256;
    imageWidth = 256;
    trialDuration = 900;
    pixelsPerFrame = imageWidth*imageHeight;
    % open the file, get file size, back to the begining
    fid = fopen(windowCamFileID);
    fseek(fid,0,'eof');
    fileSize = ftell(fid);
    fseek(fid,0,'bof');
    % identify the number of frames to read
    nFramesToRead = floor(fileSize/(pixelsPerFrame*2));
    skippedPixels = pixelsPerFrame*2;
    % loop over each frame
    imageStack = NaN*ones(imageHeight,imageWidth,nFramesToRead);
    for n = 1:nFramesToRead
        fseek(fid,(n - 1)*skippedPixels,'bof');
        z = fread(fid,pixelsPerFrame,'*int16','b');
        img = reshape(z,imageHeight,imageWidth);
        imageStack(:,:,n) = rot90(img',2);
    end
    blueImageStack = imageStack(:,:,2:2:end);
    %% pull leading/lagging frames from each stimulation
    for mm = 1:length(condensedStimulationTimes)
        stimTime = condensedStimulationTimes(1,mm);
        startIndex = round((stimTime - 5)*15);
        endIndex = startIndex + 10*15;
        frames = blueImageStack(:,:,startIndex:endIndex);
        normFrame = mean(frames(:,:,1:70),3);
        for nn = 1:size(frames,3)
            frame = frames(:,:,nn);
            msFrames(:,:,nn) = ((frame - normFrame)./normFrame)*100;
        end
        msFrameStack(:,:,:,qq) = msFrames;
        qq = qq + 1;
    end
end
%% mean of stimuli
meanFrameStack = mean(msFrameStack,4);
figure
sliceViewer(meanFrameStack,'Colormap',jet)
