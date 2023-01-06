clear; clc;
rawDataFileID = uigetfile('*RawData.mat','multiselect','off');
load(rawDataFileID)
[animalID,fileDate,fileID] = GetFileInfo_IOS(rawDataFileID);
[frames] = ReadDalsaBinary_IOS(animalID,[fileID '_WindowCam.bin']);
cbvFrames = frames(3:3:end);
gcampFrames = frames(1:3:end - 1);
deoxyFrames = frames(2:3:end);
figure;
subplot(1,3,1)
imagesc(cbvFrames{1,1})
axis image
colormap gray
title('green - CBV')
subplot(1,3,2)
imagesc(gcampFrames{1,1})
axis image
colormap gray
title('blue - GCaMP7s')
subplot(1,3,3)
imagesc(deoxyFrames{1,1})
axis image
colormap gray
title('red - deoxy')
drawnow
%
frameRate = RawData.notes.CBVCamSamplingRate/3;
imgWidth = RawData.notes.CBVCamPixelWidth;
imgHeight = RawData.notes.CBVCamPixelHeight;
%
% for aa = 1:length(cbvFrames)
%     LH_frames{1,aa} = flip(cbvFrames{1,aa}(:,1:imgWidth/2),2);
%     RH_frames{1,aa} = cbvFrames{1,aa}(:,imgWidth/2 + 1:end);
% end
for aa = 1:length(gcampFrames)
    LH_frames{1,aa} = flip(gcampFrames{1,aa}(:,1:imgWidth/2),2);
    RH_frames{1,aa} = gcampFrames{1,aa}(:,imgWidth/2 + 1:end);
end
%
slidingWindow = 10; % seconds
% correlationStack = nan(imgHeight,imgWidth/2,length(slidingWindow*frameRate + 1:length(LH_frames) - slidingWindow*frameRate));
correlationStack = [];
for aa = slidingWindow*frameRate + 1:length(LH_frames) - slidingWindow*frameRate
    LH_FrameSubset = LH_frames(aa - slidingWindow*frameRate:aa - 1);
    RH_FrameSubset = RH_frames(aa - slidingWindow*frameRate:aa - 1);
    for bb = 1:length(LH_FrameSubset)
        LH_Frame(bb,:) = reshape(LH_FrameSubset{1,bb},[1,imgWidth/2*imgHeight]);
        RH_frame(bb,:) = reshape(RH_FrameSubset{1,bb},[1,imgWidth/2*imgHeight]);
    end
    for cc = 1:size(LH_Frame,2)
        ccMatrix = corrcoef(LH_Frame(:,cc),RH_frame(:,cc));
        corR(1,cc) = ccMatrix(2,1);
    end
    corFrame = reshape(corR,[imgHeight,imgWidth/2]);
%     correlationStack(:,:,aa) = mat2gray(corFrame,[-1,1]);
    correlationStack = cat(3,correlationStack,mat2gray(corFrame,[-1,1]));
end
sliceViewer(correlationStack,'Colormap',jet)

% 
% 
% for aa = slidingWindow*frameRate + 1:length(gcampFrames) - slidingWindow*frameRate
%     gcampFrameSubset = gcampFrames(aa - slidingWindow*frameRate:aa - 1);
%     cbvFrameSubset = cbvFrames(aa - slidingWindow*frameRate:aa - 1);
%     for bb = 1:length(gcampFrameSubset)
%         gcampFrame(bb,:) = reshape(gcampFrameSubset{1,bb},[1,imgWidth*imgHeight]);
%         cbvFrame(bb,:) = reshape(cbvFrameSubset{1,bb},[1,imgWidth*imgHeight]);
%     end
%     for cc = 1:size(gcampFrame,2)
%         ccMatrix = corrcoef(gcampFrame(:,cc),cbvFrame(:,cc));
%         corR(1,cc) = ccMatrix(2,1);
%     end
%     corFrame = reshape(corR,[imgWidth,imgHeight]);
%     correlationStack = cat(3,correlationStack,mat2gray(corFrame,[0,1]));
% end
% sliceViewer(correlationStack,'Colormap',jet)