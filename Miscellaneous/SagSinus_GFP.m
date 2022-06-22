clear; clc; close all;
procDataFileIDs = ls('*_ProcData.mat');
rawDataFileIDs = ls('*_RawData.mat');
qq = 1;
puffFrameStack = [];
for aa = 1:size(procDataFileIDs) - 1
    puffDataFileID = procDataFileIDs(aa,:);
    load(puffDataFileID)
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
    %% movie frames
    [~,~,fileID] = GetFileInfo_IOS(puffDataFileID);
    windowCamFileID = [fileID '_WindowCam.bin'];
    Fs = ProcData.notes.dsFs;
    imageHeight = ProcData.notes.CBVCamPixelHeight;
    imageWidth = ProcData.notes.CBVCamPixelWidth;
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
    %% pull leading/lagging frames from each stimulation
    for mm = 1:length(condensedStimulationTimes)
        stimTime = condensedStimulationTimes(1,mm);
        startIndex = round((stimTime - 5)*Fs);
        endIndex = startIndex + 10*Fs;
        frames = imageStack(:,:,startIndex:endIndex);
        normFrame = mean(frames(:,:,1:140),3);
        for nn = 1:size(frames,3)
            frame = frames(:,:,nn);
            msFrames(:,:,nn) = ((frame - normFrame)./normFrame)*100;
        end
        puffFrameStack(:,:,:,qq) = msFrames;
        qq = qq + 1;
    end
end
% mean of stimuli
meanPuffFrameStack = mean(puffFrameStack,4);
figure
sliceViewer(meanPuffFrameStack,'Colormap',jet)
%% window image
windowImage = figure('Name','Window under 480 nm illumination');
imagesc(imageStack(:,:,1))
colormap gray
axis image
axis off
title('Window Image (CAG-EGFP, 480 nm illumination)')
savefig('windowImage.fig')
%% isoflurane image stack
[~,~,fileID] = GetFileInfo_IOS(procDataFileIDs(5,:));
load(procDataFileIDs(5,:))
load(rawDataFileIDs(5,:))
windowCamFileID = [fileID '_WindowCam.bin'];
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
% identify respiration
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1,k1] = butter(4,10/(ProcData.notes.analogSamplingRate/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
forceSensor = abs(filtfilt(sos1,g1,RawData.data.forceSensor.^2));
respiration = abs(filtfilt(sos1,g1,RawData.data.hippocampus.^2));
threshold = 0.015;
figure;
plot(forceSensor)
hold on;
yline(threshold)
% threshold for interpolation
respirationIndex = forceSensor > threshold;
[linkedDiffIndex] = LinkBinaryEvents_IOS(gt(respirationIndex,0),[ProcData.notes.analogSamplingRate,0]);
% identify edges
edgeFound = false;
xx = 1;
startEdge = [];
endEdge = [];
for aa = 1:length(linkedDiffIndex)
    if edgeFound == false
        if (linkedDiffIndex(1,aa) == 1) == true && (aa < length(linkedDiffIndex)) == true
            startEdge(xx,1) = aa;
            edgeFound = true;
        end
    elseif edgeFound == true
        if linkedDiffIndex(1,aa) == 0
            endEdge(xx,1) = aa;
            edgeFound = false;
            xx = xx + 1;
        elseif (length(linkedDiffIndex) == aa) == true && (linkedDiffIndex(1,aa) == 1) == true && edgeFound == true
            endEdge(xx,1) = aa;
        end
    end
end
for aa = 1:length(startEdge)
    breatheTimes(aa,1) = ((startEdge(aa,1) + endEdge(aa,1))/2)/ProcData.notes.analogSamplingRate;
end
isoFrameStack = [];
qq = 1;
%% pull leading/lagging frames from each stimulation
for mm = 1:length(breatheTimes)
    stimTime = breatheTimes(mm,1);
    if stimTime >= 6 && stimTime <= ProcData.notes.trialDuration_sec - 6
        startIndex = round((stimTime - 5)*Fs);
        endIndex = startIndex + 10*Fs;
        frames = imageStack(:,:,startIndex:endIndex);
        normFrame = mean(frames(:,:,1:140),3);
        for nn = 1:size(frames,3)
            frame = frames(:,:,nn);
            msFrames(:,:,nn) = ((frame - normFrame)./normFrame)*100;
        end
        isoFrameStack(:,:,:,qq) = msFrames;
        qq = qq + 1;
    end
end
% mean of stimuli
meanIsoFrameStack = mean(isoFrameStack,4);
figure
sliceViewer(meanIsoFrameStack,'Colormap',jet)
%% summary figure
summaryFig = figure('Name','0.5s Whisker Puff vs. Isoflurane');
%% puff stimulation
% t = 0
ax1 = subplot(2,4,1);
imagesc(meanPuffFrameStack(:,:,150))
colormap jet
caxis([-10,10])
axis image
axis off
title('puff @ t = 0')
% t = 1/3 sec
ax2 = subplot(2,4,2);
imagesc(meanPuffFrameStack(:,:,160))
colormap jet
caxis([-10,10])
axis image
axis off
title('t = 1/3 sec')
% t = 1 sec
ax3 = subplot(2,4,3);
imagesc(meanPuffFrameStack(:,:,180))
colormap jet
caxis([-10,10])
axis image
axis off
title('t = 1 sec')
% t = 2 sec
ax4 = subplot(2,4,4);
imagesc(meanPuffFrameStack(:,:,210))
colormap jet
caxis([-10,10])
c1 = colorbar;
ylabel(c1,'\DeltaF/F (%)','rotation',-90,'VerticalAlignment','bottom')
axis image
axis off
title('t = 2 sec')
% adjust and link axes
linkaxes([ax1,ax2,ax3,ax4],'xy')
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax4Pos(3:4) = ax1Pos(3:4);
set(ax4,'position',ax4Pos);
%% isoflurane
% t = 0
ax5 = subplot(2,4,5);
imagesc(isoFrameStack(:,:,75))
colormap jet
caxis([-10,10])
axis image
axis off
title('breathe @ t = 0')
% t = 1/3 sec
ax6 = subplot(2,4,6);
imagesc(isoFrameStack(:,:,80))
colormap jet
caxis([-10,10])
axis image
axis off
title('t = 1/3 sec')
% t = 1 sec
ax7 = subplot(2,4,7);
imagesc(isoFrameStack(:,:,90))
colormap jet
caxis([-10,10])
axis image
axis off
title('t = 1 sec')
% t = 2 sec
ax8 = subplot(2,4,8);
imagesc(isoFrameStack(:,:,105))
colormap jet
caxis([-10,10])
c2 = colorbar;
ylabel(c2,'\DeltaF/F (%)','rotation',-90,'VerticalAlignment','bottom')
axis image
axis off
title('t = 2 sec')
% adjust and link axes
linkaxes([ax5,ax6,ax7,ax8],'xy')
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
ax8Pos = get(ax8,'position');
ax8Pos(3:4) = ax5Pos(3:4);
set(ax8,'position',ax8Pos);
savefig('PuffvsIso.fig')
