clear; clc;
% user inputs for file information
rawDataFileIDs = uigetfile('*_RawData.mat','MultiSelect','on');
for aa = 1:size(rawDataFileIDs,2)
    clearvars -except rawDataFileIDs aa ROIs Data
    load(rawDataFileIDs{1,aa})
    [~,~,fileID] = GetFileInfo_IOS(rawDataFileIDs{1,aa});
    windowCamFileID = [fileID '_WindowCam.bin'];
    Fs = RawData.notes.CBVCamSamplingRate;
    imageHeight = RawData.notes.CBVCamPixelHeight;
    imageWidth = RawData.notes.CBVCamPixelWidth;
    trialDuration = RawData.notes.trialDuration_sec;
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
    if aa == 1
        % create figure of the image frame
        roiFig = figure;
        imagesc(imageStack(:,:,3))
        colormap(gray)
        axis image
        xlabel('Caudal')
        caxis([0,2000])
        % draw ROI over the cement
        disp('Please select your region of interest for the left hemisphere'); disp(' ')
        [~,LH_xi,LH_yi] = roipoly;
        ROIs.LH_xi = LH_xi;
        ROIs.LH_yi = LH_yi;
        disp('Please select your region of interest for the right hemisphere'); disp(' ')
        [~,RH_xi,RH_yi] = roipoly;
        ROIs.RH_xi = RH_xi;
        ROIs.RH_yi = RH_yi;
        close(roiFig)
    end
    %% extract CBV/GCaMP frames
    cbvFrames = imageStack(:,:,3:2:end - 1);
    gcampFrames = imageStack(:,:,4:2:end);
    numFrames = size(cbvFrames,3);
    % LH - apply image mask to each frame of the stack
    LH_cbvRefl = zeros(1,numFrames);
    LH_gcampRefl = zeros(1,numFrames);
    for n = 1:numFrames
        LH_mask = roipoly(cbvFrames(:,:,1),ROIs.LH_xi,ROIs.LH_yi);
        LH_cbvMask = LH_mask.*double(cbvFrames(:,:,n));
        LH_gcampMask = LH_mask.*double(gcampFrames(:,:,n));
        LH_cbvRefl(n) = mean(nonzeros(LH_cbvMask));
        LH_gcampRefl(n) = mean(nonzeros(LH_gcampMask));
    end
    Data.LH_cbvRefl{aa,1} = LH_cbvRefl;
    Data.LH_gcampRefl{aa,1} = LH_gcampRefl;
    LH_cbvRefl_norm = (LH_cbvRefl - mean(LH_cbvRefl))./mean(LH_cbvRefl);
    LH_gcampRefl_norm = (LH_gcampRefl - mean(LH_gcampRefl))./mean(LH_gcampRefl);
    Data.LH_cbvRefl_norm{aa,1} = LH_cbvRefl_norm;
    Data.LH_gcampRefl_norm{aa,1} = LH_gcampRefl_norm;
    % RH - apply image mask to each frame of the stack
    RH_cbvRefl = zeros(1,numFrames);
    RH_gcampRefl = zeros(1,numFrames);
    for n = 1:numFrames
        RH_mask = roipoly(cbvFrames(:,:,1),ROIs.RH_xi,ROIs.RH_yi);
        RH_cbvMask = RH_mask.*double(cbvFrames(:,:,n));
        RH_gcampMask = RH_mask.*double(gcampFrames(:,:,n));
        RH_cbvRefl(n) = mean(nonzeros(RH_cbvMask));
        RH_gcampRefl(n) = mean(nonzeros(RH_gcampMask));
    end
    Data.RH_cbvRefl{aa,1} = RH_cbvRefl;
    Data.RH_gcampRefl{aa,1} = RH_gcampRefl;
    RH_cbvRefl_norm = (RH_cbvRefl - mean(RH_cbvRefl))./mean(RH_cbvRefl);
    RH_gcampRefl_norm = (RH_gcampRefl - mean(RH_gcampRefl))./mean(RH_gcampRefl);
    Data.RH_cbvRefl_norm{aa,1} = RH_cbvRefl_norm;
    Data.RH_gcampRefl_norm{aa,1} = RH_gcampRefl_norm;
end
catLH_cbvRefl = []; catRH_cbvRefl = [];
catLH_gcampRefl = []; catRH_gcampRefl = [];
for bb = 1:size(rawDataFileIDs,2)
    catLH_cbvRefl = cat(2,catLH_cbvRefl,Data.LH_cbvRefl{bb,1});
    catRH_cbvRefl = cat(2,catRH_cbvRefl,Data.RH_cbvRefl{bb,1});
    catLH_gcampRefl = cat(2,catLH_gcampRefl,Data.LH_gcampRefl{bb,1});
    catRH_gcampRefl = cat(2,catRH_gcampRefl,Data.RH_gcampRefl{bb,1});
end
%% generate figure
figure;
plot(detrend(catLH_cbvRefl,'constant'))
hold on
plot(detrend(catRH_cbvRefl,'constant'))
plot(detrend(catLH_gcampRefl,'constant'))
plot(detrend(catRH_gcampRefl,'constant'))
legend('LH CBV','RH CBV','LH GCaMP7s','RH GCaMP7s')
ylabel('Percent (%)')
%
figure;
plot(catLH_cbvRefl)
hold on
plot(catRH_cbvRefl)
plot(catLH_gcampRefl)
plot(catRH_gcampRefl)
legend('LH CBV','RH CBV','LH GCaMP7s','RH GCaMP7s')
ylabel('Percent (%)')
