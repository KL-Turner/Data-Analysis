%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate supplemental videos for Turner_Kederasetti_Gheres_Proctor_Costanzo_Drew_eLife2020
%________________________________________________________________________________________________________________________

clear; clc; close all;
%% information and data for first example
% dataLocation = [rootFolder '\Summary Figures and Structures\Supplemental Movies\'];
% cd(dataLocation)
fs1 = 30;   % cbv/pupil camera
fs2 = 150;  % whisker camera
colorRfcAwake = [(0/256),(64/256),(64/256)];
colorRfcNREM = [(0/256),(174/256),(239/256)];
colorRfcREM = [(190/256),(30/256),(45/256)];
%% Take frames from the CBV camera file for baseline/awake/rem examples
exampleRawDataFileID_A = 'T123_200301_14_48_14_RawData.mat';
load(exampleRawDataFileID_A,'-mat')
exampleCBVcamFileID_A = '200301_14_48_14_WindowCam.bin';
baseStartTime = 580;
baseEndTime = 600;
baseFrameIndex = baseStartTime*fs1:baseEndTime*fs1;
awakeStartTime = 400;
awakeEndTime = 600;
awakeFrameIndex = awakeStartTime*fs1:awakeEndTime*fs1;
awakeFrameWhiskIndex = awakeStartTime*fs2:awakeEndTime*fs2;
remStartTime = 100;
remEndTime = 300;
remFrameIndex = remStartTime*fs1:remEndTime*fs1;
remFrameWhiskIndex = remStartTime*fs2:remEndTime*fs2;
% Obtain subset of desired frames
cbvImageHeight = RawData.notes.CBVCamPixelHeight;
cbvImageWidth = RawData.notes.CBVCamPixelWidth;
baseCBVframes = GetCBVFrameSubset_IOS_eLife2020(exampleCBVcamFileID_A,cbvImageHeight,cbvImageWidth,baseFrameIndex);
awakeCBVframes = GetCBVFrameSubset_IOS_eLife2020(exampleCBVcamFileID_A,cbvImageHeight,cbvImageWidth,awakeFrameIndex);
remCBVframes = GetCBVFrameSubset_IOS_eLife2020(exampleCBVcamFileID_A,cbvImageHeight,cbvImageWidth,remFrameIndex);
%% Take frames from the CBV camera file for nrem example
exampleCBVcamFileID_B = '200301_15_03_33_WindowCam.bin';
nremStartTime = 300;
nremEndTime = 500;
nremFrameIndex = nremStartTime*fs1:nremEndTime*fs1;
nremFrameWhiskIndex = nremStartTime*fs2:nremEndTime*fs2;
% Obtain subset of desired frames
nremCBVframes = GetCBVFrameSubset_IOS_eLife2020(exampleCBVcamFileID_B,cbvImageHeight,cbvImageWidth,nremFrameIndex);
%% Use baseline frames to establish ROIs and resting baseline
% draw a rectangular ROI around the window to remove outside pixels
boxFig = figure;
imagesc(baseCBVframes(:,:,1))
colormap gray
caxis([0,2^12])
axis image
boxROI = drawrectangle;
boxPosition = round(boxROI.Vertices);
boxX = boxPosition(:,1);
boxY = boxPosition(:,2);
boxMask = poly2mask(boxX,boxY,cbvImageWidth,cbvImageHeight);
close(boxFig)
% take values from within a square ROI within the window
boxWidth = abs(boxPosition(1,1) - boxPosition(3,1));
boxHeight = abs(boxPosition(1,2) - boxPosition(2,2));
% baseline frames
for a = 1:size(baseCBVframes,3)
    baseCBVframe = baseCBVframes(:,:,a);
    baseBoxVals = baseCBVframe(boxMask);
    baseBoxFrames(:,:,a) = reshape(baseBoxVals,boxHeight,boxWidth); %#ok<*SAGROW>
end
% awake frames
for a = 1:size(awakeCBVframes,3)
    awakeCBVframe = awakeCBVframes(:,:,a);
    awakeBoxVals = awakeCBVframe(boxMask);
    awakeBoxFrames(:,:,a) = reshape(awakeBoxVals,boxHeight,boxWidth);
end
% nrem frames
for a = 1:size(nremCBVframes,3)
    nremCBVframe = nremCBVframes(:,:,a);
    nremBoxVals = nremCBVframe(boxMask);
    nremBoxFrames(:,:,a) = reshape(nremBoxVals,boxHeight,boxWidth);
end
% rem frames
for a = 1:size(remCBVframes,3)
    remCBVframe = remCBVframes(:,:,a);
    remBoxVals = remCBVframe(boxMask);
    remBoxFrames(:,:,a) = reshape(remBoxVals,boxHeight,boxWidth);
end
% set values from outside the window to NaN
windowFig = figure;
imagesc(baseBoxFrames(:,:,1))
colormap gray
caxis([0,2^12])
axis image
windowMask = roipoly;
close(windowFig)
% baseline frames
for w = 1:size(baseBoxFrames,3)
    baseWindowFrame = baseBoxFrames(:,:,w);
    baseWindowFrame(~windowMask) = NaN;
    baseWindowFrames(:,:,w) = baseWindowFrame;
end
% awake frames
for w = 1:size(awakeBoxFrames,3)
    awakeWindowFrame = awakeBoxFrames(:,:,w);
    awakeWindowFrame(~windowMask) = NaN;
    awakeWindowFrames(:,:,w) = awakeWindowFrame;
end
% nrem frames
for w = 1:size(nremBoxFrames,3)
    nremWindowFrame = nremBoxFrames(:,:,w);
    nremWindowFrame(~windowMask) = NaN;
    nremWindowFrames(:,:,w) = nremWindowFrame;
end
% rem frames
for w = 1:size(remBoxFrames,3)
    remWindowFrame = remBoxFrames(:,:,w);
    remWindowFrame(~windowMask) = NaN;
    remWindowFrames(:,:,w) = remWindowFrame;
end
%% Normalize each data set by the baseline frame
baselineFrame = mean(baseWindowFrames,3);
% awake data
for a = 1:size(awakeWindowFrames,3)
    awakeCBVImageStack(:,:,a) = ((awakeWindowFrames(:,:,a) - baselineFrame)./(baselineFrame)).*100;
end
% nrem data
for a = 1:size(nremWindowFrames,3)
    nremCBVImageStack(:,:,a) = ((nremWindowFrames(:,:,a) - baselineFrame)./(baselineFrame)).*100;
end
% rem data
for a = 1:size(remWindowFrames,3)
    remCBVImageStack(:,:,a) = ((remWindowFrames(:,:,a) - baselineFrame)./(baselineFrame)).*100;
end
%% Take frames from the Whisker camera file for baseline/awake/rem examples
exampleRawDataFileID_A = 'T123_200301_14_48_14_RawData.mat';
load(exampleRawDataFileID_A,'-mat')
exampleWhiskercamFileID_A = '200301_14_48_14_WhiskerCam.bin';
% Obtain subset of desired frames
whiskImageHeight = RawData.notes.whiskCamPixelHeight;
whiskImageWidth = RawData.notes.whiskCamPixelWidth;
whiskerPixelsPerFrame = whiskImageWidth*whiskImageHeight;
whiskerSkippedPixels = whiskerPixelsPerFrame; % (depends on cam setting) Multiply by two because there are 16 bits (2 bytes) per pixel
whiskFid = fopen(exampleWhiskercamFileID_A);
fseek(whiskFid,0,'eof');
% fileSize = ftell(whiskFid);
fseek(whiskFid,0,'bof');
% awake frames
awakeWhiskNFramesToRead = length(awakeFrameWhiskIndex);
awakeWhiskImageStack = zeros(whiskImageHeight,whiskImageWidth,awakeWhiskNFramesToRead);
for a = 1:awakeWhiskNFramesToRead
    fseek(whiskFid,awakeFrameWhiskIndex(a)*whiskerSkippedPixels,'bof');
    whiskZ = fread(whiskFid,whiskerPixelsPerFrame,'*uint8','b');
    whiskImg = reshape(whiskZ(1:whiskerPixelsPerFrame),whiskImageWidth,whiskImageHeight);
    awakeWhiskImageStack(:,:,a) = flip(imrotate(whiskImg,-90),2);
end
% rem frames
remWhiskNFramesToRead = length(remFrameWhiskIndex);
remWhiskImageStack = zeros(whiskImageHeight,whiskImageWidth,remWhiskNFramesToRead);
for a = 1:remWhiskNFramesToRead
    fseek(whiskFid,remFrameWhiskIndex(a)*whiskerSkippedPixels,'bof');
    whiskZ = fread(whiskFid,whiskerPixelsPerFrame,'*uint8','b');
    whiskImg = reshape(whiskZ(1:whiskerPixelsPerFrame),whiskImageWidth,whiskImageHeight);
    remWhiskImageStack(:,:,a) = flip(imrotate(whiskImg,-90),2);
end
fclose('all');
%% Take frames from the Whisker camera file for baseline/awake/rem examples
exampleWhiskercamFileID_B = '200301_15_03_33_WhiskerCam.bin';
% Obtain subset of desired frames
whiskFid = fopen(exampleWhiskercamFileID_B);
fseek(whiskFid,0,'eof');
fseek(whiskFid,0,'bof');
% nrem frames
nremWhiskNFramesToRead = length(nremFrameWhiskIndex);
nremWhiskImageStack = zeros(whiskImageHeight,whiskImageWidth,nremWhiskNFramesToRead);
for a = 1:nremWhiskNFramesToRead
    fseek(whiskFid,nremFrameWhiskIndex(a)*whiskerSkippedPixels,'bof');
    whiskZ = fread(whiskFid,whiskerPixelsPerFrame,'*uint8','b');
    whiskImg = reshape(whiskZ(1:whiskerPixelsPerFrame),whiskImageWidth,whiskImageHeight);
    nremWhiskImageStack(:,:,a) = flip(imrotate(whiskImg,-90),2);
end
fclose('all');
%% Downsample whisking image stack to match that of the cbv/whisker cameras
% awake frames
c = 1;
for b = 1:size(awakeWhiskImageStack,3)
    if rem(b,5) == 1
        dsAwakeWhiskImageStack(:,:,c) = awakeWhiskImageStack(:,:,b);
        c = c + 1;
    end
end
% nrem frames
c = 1;
for b = 1:size(nremWhiskImageStack,3)
    if rem(b,5) == 1
        dsNremWhiskImageStack(:,:,c) = nremWhiskImageStack(:,:,b);
        c = c + 1;
    end
end
% rem frames
c = 1;
for b = 1:size(remWhiskImageStack,3)
    if rem(b,5) == 1
        dsRemWhiskImageStack(:,:,c) = remWhiskImageStack(:,:,b);
        c = c + 1;
    end
end
%% Take frames from the Pupil camera file for baseline/awake/rem examples
exampleRawDataFileID_A = 'T123_200301_14_48_14_RawData.mat';
load(exampleRawDataFileID_A,'-mat')
examplePupilcamFileID_A = '200301_14_48_14_PupilCam.bin';
% Obtain subset of desired frames
pupilImageHeight = RawData.notes.pupilCamPixelHeight;
pupilImageWidth = RawData.notes.pupilCamPixelWidth;
pupilerPixelsPerFrame = pupilImageWidth*pupilImageHeight;
pupilerSkippedPixels = pupilerPixelsPerFrame; % (depends on cam setting) Multiply by two because there are 16 bits (2 bytes) per pixel
pupilFid = fopen(examplePupilcamFileID_A);
fseek(pupilFid,0,'eof');
% fileSize = ftell(pupilFid);
fseek(pupilFid,0,'bof');
% awake frames
awakePupilNFramesToRead = length(awakeFrameIndex);
awakePupilImageStack = zeros(pupilImageHeight,pupilImageWidth,awakePupilNFramesToRead);
for a = 1:awakePupilNFramesToRead
    fseek(pupilFid,awakeFrameIndex(a)*pupilerSkippedPixels,'bof');
    pupilZ = fread(pupilFid,pupilerPixelsPerFrame,'*uint8','b');
    pupilImg = reshape(pupilZ(1:pupilerPixelsPerFrame),pupilImageWidth,pupilImageHeight);
    awakePupilImageStack(:,:,a) = flip(imrotate(pupilImg,-90),2);
end
% rem frames
remPupilNFramesToRead = length(remFrameIndex);
remPupilImageStack = zeros(pupilImageHeight,pupilImageWidth,remPupilNFramesToRead);
for a = 1:remPupilNFramesToRead
    fseek(pupilFid,remFrameIndex(a)*pupilerSkippedPixels,'bof');
    pupilZ = fread(pupilFid,pupilerPixelsPerFrame,'*uint8','b');
    pupilImg = reshape(pupilZ(1:pupilerPixelsPerFrame),pupilImageWidth,pupilImageHeight);
    remPupilImageStack(:,:,a) = flip(imrotate(pupilImg,-90),2);
end
fclose('all');
%% Take frames from the Pupil camera file for baseline/awake/rem examples
examplePupilcamFileID_B = '200301_15_03_33_PupilCam.bin';
% Obtain subset of desired frames
pupilFid = fopen(examplePupilcamFileID_B);
fseek(pupilFid,0,'eof');
fseek(pupilFid,0,'bof');
% nrem frames
nremPupilNFramesToRead = length(nremFrameIndex);
nremPupilImageStack = zeros(pupilImageHeight,pupilImageWidth,nremPupilNFramesToRead);
for a = 1:nremPupilNFramesToRead
    fseek(pupilFid,nremFrameIndex(a)*pupilerSkippedPixels,'bof');
    pupilZ = fread(pupilFid,pupilerPixelsPerFrame,'*uint8','b');
    pupilImg = reshape(pupilZ(1:pupilerPixelsPerFrame),pupilImageWidth,pupilImageHeight);
    nremPupilImageStack(:,:,a) = flip(imrotate(pupilImg,-90),2);
end
fclose('all');
clearvars -except awakeCBVImageStack nremCBVImageStack remCBVImageStack...
    dsAwakeWhiskImageStack dsNremWhiskImageStack dsRemWhiskImageStack whiskImg...
    awakePupilImageStack nremPupilImageStack remPupilImageStack pupilImg...
    awakeStartTime awakeEndTime nremStartTime nremEndTime remStartTime remEndTime fs1
%% Generate supplemental video file for awake example
exampleProcDataFileID_A = 'T123_200301_14_48_14_ProcData.mat';
load(exampleProcDataFileID_A,'-mat')
exampleSpecDataFileID = 'T123_200301_14_48_14_SpecDataA.mat';
load(exampleSpecDataFileID,'-mat')
exampleBaselineFileID = 'T123_RestingBaselines.mat';
load(exampleBaselineFileID,'-mat')
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1,k1] = butter(4,10/(ProcData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(4,0.5/(ProcData.notes.dsFs/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
% whisker angle
filtWhiskerAngle = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
% heart rate
heartRate = ProcData.data.heartRate;
% CBV data
HbT = ProcData.data.CBV_HbT.adjBarrels;
filtHbT = filtfilt(sos2,g2,HbT);
% cortical and hippocampal spectrograms
cortical_LHnormS = SpecData.cortical_LH.normS.*100;
hippocampusNormS = SpecData.hippocampus.normS.*100;
T = SpecData.cortical_LH.T;
F = SpecData.cortical_LH.F;
try
    % movie file comparing processed with original data
    outputVideo = VideoWriter('SupplementalVideo1_Awake.mp4','MPEG-4');
    fps = 30;   % default fps from video acquisition
    speedUp = 2;   % speed up by factor of
    outputVideo.FrameRate = fps*speedUp;
    open(outputVideo);
    %%
    fig = figure('Position',get(0,'Screensize'));
    sgtitle({'Supplemental Video 1 (Awake)',' '})
    % Whisker angle and heart rate
    ax12 = subplot(6,4,[1,2,3,5,6,7]);
    p1 = plot((1:length(filtWhiskerAngle))/ProcData.notes.dsFs,-filtWhiskerAngle,'color','k','LineWidth',0.5);
    hold on;
    plot([400,600],[60,60],'color',colorRfcAwake,'LineWidth',5)
    ylabel({'Whisker','angle (deg)'})
    xlim([400,600])
    ylim([-20,61.5])
    yyaxis right
    p2 = plot((1:length(heartRate)),heartRate,'color',colors_eLife2020('deep carrot orange'),'LineWidth',0.5);
    ylabel('Heart rate (Hz)','rotation',-90,'VerticalAlignment','bottom')
    legend([p1,p2],'Whisker angle','Heart rate','Location','NorthWest','AutoUpdate','off')
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([400,450,500,550,600])
    xlim([400,600])
    ylim([5,10.25])
    ax12.TickLength = [0.01,0.01];
    ax12.YAxis(1).Color = 'k';
    ax12.YAxis(2).Color = colors_eLife2020('deep carrot orange');
    % CBV and behavioral indeces
    ax34 = subplot(6,4,[9,10,11,13,14,15]);
    p3 = plot((1:length(filtHbT))/ProcData.notes.CBVCamSamplingRate,filtHbT,'color','r','LineWidth',1);
    hold on
    p4 = plot([400,600],[100,100],'color',colorRfcAwake,'LineWidth',5);
    xlim([400,600])
    ylim([-45,100])
    ylabel('\Delta[HbT] (\muM)')
    legend([p3,p4],'\DeltaHbT','rfc-Awake','Location','NorthWest','AutoUpdate','off')
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([400,450,500,550,600])
    xlim([400,600])
    ax34.TickLength = [0.01,0.01];
    % Left cortical electrode spectrogram
    ax5 = subplot(6,4,[17,18,19]);
    semilog_imagesc_eLife2020(T,F,cortical_LHnormS,'y')
    axis xy
    hold on
    plot([400,600],[150,150],'color',colorRfcAwake,'LineWidth',5)
    c5 = colorbar;
    ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    yticks([10,100])
    yticklabels({'10','100'})
    caxis([-100,100])
    ylabel({'Cortical LFP','Freq (Hz)'})
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([400,450,500,550,600])
    xlim([400,600])
    ax5.TickLength = [0.01,0.01];
    % Hippocampal electrode spectrogram
    ax6 = subplot(6,4,[21,22,23]);
    semilog_imagesc_eLife2020(T,F,hippocampusNormS,'y')
    axis xy
    hold on
    plot([400,600],[150,150],'color',colorRfcAwake,'LineWidth',5)
    c6 = colorbar;
    ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    yticks([10,100])
    yticklabels({'10','100'})
    caxis([-100,100])
    xlabel('Time (s)')
    ylabel({'Hippocampal LFP','Freq (Hz)'})
    set(gca,'box','off')
    xticks([400,450,500,550,600])
    xticklabels({'0','50','100','150','200'})
    xlim([400,600])
    ax6.TickLength = [0.01,0.01];
    % Axes properties
    ax12Pos = get(ax12,'position');
    ax5Pos = get(ax5,'position');
    ax6Pos = get(ax6,'position');
    ax5Pos(3) = ax12Pos(3);
    ax6Pos(3) = ax12Pos(3);
    set(ax5,'position',ax5Pos);
    set(ax6,'position',ax6Pos);
    c6Pos = get(c6,'position');
    c5Pos = get(c5,'position');
    c5Pos(1) = c6Pos(1);
    set(c5,'position',c5Pos);
    %%
    for a = 1:size(awakePupilImageStack,3)
        % lines
        axes(ax12)
        x12 = xline(awakeStartTime + (a/fs1),'color',colorRfcAwake,'LineWidth',2);
        axes(ax34)
        x34 = xline(awakeStartTime + (a/fs1),'color',colorRfcAwake,'LineWidth',2);
        axes(ax5)
        x5 = xline(awakeStartTime + (a/fs1),'color',colorRfcAwake,'LineWidth',2);
        axes(ax6)
        x6 = xline(awakeStartTime + (a/fs1),'color',colorRfcAwake,'LineWidth',2);
        % whisker movie
        s1 = subplot(6,4,[4,8]);
        imagesc(dsAwakeWhiskImageStack(:,:,a))
        title('Whisker camera')
        colormap(gca,'gray')
        caxis([0,2^8])
        axis image
        axis off
        hold off
        % window movie
        s2 = subplot(6,4,[12,16]);
        imagesc(awakeCBVImageStack(:,:,a))
        title('IOS camera')
        colormap(gca,'gray')
        c0 = colorbar;
        ylabel(c0,'\DeltaR/R (%)','rotation',-90,'VerticalAlignment','bottom')
        caxis([-15,15])
        axis image
        axis off
        hold off
        % pupil movie
        s3 = subplot(6,4,[20,24]);
        imagesc(awakePupilImageStack(:,:,a))
        title('Eye camera')
        colormap(gca,'gray')
        caxis([0,2^8])
        axis image
        axis off
        hold off
        currentFrame = getframe(fig);
        writeVideo(outputVideo,currentFrame);
        delete(s1)
        delete(s2)
        delete(s3)
        delete(x12)
        delete(x34)
        delete(x5)
        delete(x6)
    end
    close(outputVideo)
    close(fig)
    sendmail('kevinlturnerjr@gmail.com','Video S1 Complete');
catch
    sendmail('kevinlturnerjr@gmail.com','Video S1 Error');
end
%% Generate supplemental video file for NREM example
exampleProcDataFileID_B = 'T123_200301_15_03_33_ProcData.mat';
load(exampleProcDataFileID_B,'-mat')
exampleSpecDataFileID = 'T123_200301_15_03_33_SpecDataA.mat';
load(exampleSpecDataFileID,'-mat')
exampleBaselineFileID = 'T123_RestingBaselines.mat';
load(exampleBaselineFileID,'-mat')
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1,k1] = butter(4,10/(ProcData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(4,0.5/(ProcData.notes.dsFs/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
% whisker angle
filtWhiskerAngle = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
% heart rate
heartRate = ProcData.data.heartRate;
% CBV data
HbT = ProcData.data.CBV_HbT.adjBarrels;
filtHbT = filtfilt(sos2,g2,HbT);
% cortical and hippocampal spectrograms
cortical_LHnormS = SpecData.cortical_LH.normS.*100;
hippocampusNormS = SpecData.hippocampus.normS.*100;
T = SpecData.cortical_LH.T;
F = SpecData.cortical_LH.F;
try
    % movie file comparing processed with original data
    outputVideo = VideoWriter('SupplementalVideo2_NREM.mp4','MPEG-4');
    fps = 30;   % default fps from video acquisition
    speedUp = 2;   % speed up by factor of
    outputVideo.FrameRate = fps*speedUp;
    open(outputVideo);
    %%
    fig = figure('Position',get(0,'Screensize'));
    sgtitle({'Supplemental Video 2 (NREM)',' '})
    % Whisker angle and heart rate
    ax12 = subplot(6,4,[1,2,3,5,6,7]);
    p1 = plot((1:length(filtWhiskerAngle))/ProcData.notes.dsFs,-filtWhiskerAngle,'color','k','LineWidth',0.5);
    hold on
    plot([300,336],[60,60],'color',colorRfcNREM,'LineWidth',5)
    plot([336,360],[60,60],'color',colorRfcAwake,'LineWidth',5)
    plot([360,500],[60,60],'color',colorRfcNREM,'LineWidth',5)
    ylabel({'Whisker','angle (deg)'})
    xlim([300,500])
    ylim([-20,61.5])
    yyaxis right
    p2 = plot((1:length(heartRate)),heartRate,'color',colors_eLife2020('deep carrot orange'),'LineWidth',0.5);
    ylabel('Heart rate (Hz)','rotation',-90,'VerticalAlignment','bottom')
    legend([p1,p2],'Whisker angle','Heart rate','Location','NorthWest','AutoUpdate','off')
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([300,350,400,450,500])
    xlim([300,500])
    ylim([5,10.25])
    ax12.TickLength = [0.01,0.01];
    ax12.YAxis(1).Color = 'k';
    ax12.YAxis(2).Color = colors_eLife2020('deep carrot orange');
    % CBV and behavioral indeces
    ax34 = subplot(6,4,[9,10,11,13,14,15]);
    p3 = plot((1:length(filtHbT))/ProcData.notes.CBVCamSamplingRate,filtHbT,'color','r','LineWidth',1);
    hold on
    p4 = plot([300,336],[100,100],'color',colorRfcNREM,'LineWidth',5);
    p5 = plot([336,360],[100,100],'color',colorRfcAwake,'LineWidth',5);
    plot([360,500],[100,100],'color',colorRfcNREM,'LineWidth',5)
    xlim([300,500])
    ylim([-45,100])
    ylabel('\Delta[HbT] (\muM)')
    legend([p3,p5,p4],'\DeltaHbT','rfc-Awake','rfc-NREM','Location','NorthWest','AutoUpdate','off')
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([300,350,400,450,500])
    ax34.TickLength = [0.01,0.01];
    % Left cortical electrode spectrogram
    ax5 = subplot(6,4,[17,18,19]);
    semilog_imagesc_eLife2020(T,F,cortical_LHnormS,'y')
    axis xy
    hold on
    plot([300,336],[150,150],'color',colorRfcNREM,'LineWidth',5)
    plot([336,360],[150,150],'color',colorRfcAwake,'LineWidth',5)
    plot([360,500],[150,150],'color',colorRfcNREM,'LineWidth',5)
    c5 = colorbar;
    ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    caxis([-100,100])
    ylabel({'Cortical LFP','Freq (Hz)'})
    yticks([10,100])
    yticklabels({'10','100'})
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([300,350,400,450,500])
    xlim([300,500])
    ax5.TickLength = [0.01,0.01];
    % Hippocampal electrode spectrogram
    ax6 = subplot(6,4,[21,22,23]);
    semilog_imagesc_eLife2020(T,F,hippocampusNormS,'y')
    hold on
    axis xy
    plot([300,336],[150,150],'color',colorRfcNREM,'LineWidth',5)
    plot([336,360],[150,150],'color',colorRfcAwake,'LineWidth',5)
    plot([360,500],[150,150],'color',colorRfcNREM,'LineWidth',5)
    c6 = colorbar;
    ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    caxis([-100,100])
    xlabel('Time (s)')
    ylabel({'Hippocampal LFP','Freq (Hz)'})
    yticks([10,100])
    yticklabels({'10','100'})
    set(gca,'box','off')
    xticks([300,350,400,450,500])
    xticklabels({'0','50','100','150','200'})
    xlim([300,500])
    ax6.TickLength = [0.01,0.01];
    % Axes properties
    ax12Pos = get(ax12,'position');
    ax5Pos = get(ax5,'position');
    ax6Pos = get(ax6,'position');
    ax5Pos(3) = ax12Pos(3);
    ax6Pos(3) = ax12Pos(3);
    set(ax5,'position',ax5Pos);
    set(ax6,'position',ax6Pos);
    c6Pos = get(c6,'position');
    c5Pos = get(c5,'position');
    c5Pos(1) = c6Pos(1);
    set(c5,'position',c5Pos);
    %%
    for a = 1:size(nremPupilImageStack,3)
        % lines
        axes(ax12)
        if nremStartTime + (a/fs1) <= 336
            x12 = xline(nremStartTime + (a/fs1),'color',colorRfcNREM,'LineWidth',2);
            axes(ax34)
            x34 = xline(nremStartTime + (a/fs1),'color',colorRfcNREM,'LineWidth',2);
            axes(ax5)
            x5 = xline(nremStartTime + (a/fs1),'color',colorRfcNREM,'LineWidth',2);
            axes(ax6)
            x6 = xline(nremStartTime + (a/fs1),'color',colorRfcNREM,'LineWidth',2);
        elseif nremStartTime + (a/fs1) > 336 && nremStartTime + (a/fs1) < 360
            x12 = xline(nremStartTime + (a/fs1),'color',colorRfcAwake,'LineWidth',2);
            axes(ax34)
            x34 = xline(nremStartTime + (a/fs1),'color',colorRfcAwake,'LineWidth',2);
            axes(ax5)
            x5 = xline(nremStartTime + (a/fs1),'color',colorRfcAwake,'LineWidth',2);
            axes(ax6)
            x6 = xline(nremStartTime + (a/fs1),'color',colorRfcAwake,'LineWidth',2);
        elseif nremStartTime + (a/fs1) >= 360
            x12 = xline(nremStartTime + (a/fs1),'color',colorRfcNREM,'LineWidth',2);
            axes(ax34)
            x34 = xline(nremStartTime + (a/fs1),'color',colorRfcNREM,'LineWidth',2);
            axes(ax5)
            x5 = xline(nremStartTime + (a/fs1),'color',colorRfcNREM,'LineWidth',2);
            axes(ax6)
            x6 = xline(nremStartTime + (a/fs1),'color',colorRfcNREM,'LineWidth',2);
        end
        % whisker movie
        s1 = subplot(6,4,[4,8]);
        imagesc(dsNremWhiskImageStack(:,:,a))
        title('Whisker camera')
        colormap(gca,'gray')
        caxis([0,2^8])
        axis image
        axis off
        hold off
        % window movie
        s2 = subplot(6,4,[12,16]);
        imagesc(nremCBVImageStack(:,:,a))
        title('IOS camera')
        colormap(gca,'gray')
        c0 = colorbar;
        ylabel(c0,'\DeltaR/R (%)','rotation',-90,'VerticalAlignment','bottom')
        caxis([-15,15])
        axis image
        axis off
        hold off
        % pupil movie
        s3 = subplot(6,4,[20,24]);
        imagesc(nremPupilImageStack(:,:,a))
        title('Eye camera')
        colormap(gca,'gray')
        caxis([0,2^8])
        axis image
        axis off
        hold off
        currentFrame = getframe(fig);
        writeVideo(outputVideo,currentFrame);
        delete(s1)
        delete(s2)
        delete(s3)
        delete(x12)
        delete(x34)
        delete(x5)
        delete(x6)
    end
    close(outputVideo)
    close(fig)
    sendmail('kevinlturnerjr@gmail.com','Video S2 Complete');
catch
    sendmail('kevinlturnerjr@gmail.com','Video S2 Error');
end
%% Generate supplemental video file for REM example
exampleProcDataFileID_A = 'T123_200301_14_48_14_ProcData.mat';
load(exampleProcDataFileID_A,'-mat')
exampleSpecDataFileID = 'T123_200301_14_48_14_SpecDataA.mat';
load(exampleSpecDataFileID,'-mat')
exampleBaselineFileID = 'T123_RestingBaselines.mat';
load(exampleBaselineFileID,'-mat')
[~,fileDate,~] = GetFileInfo_IOS_eLife2020(exampleProcDataFileID_A);
strDay = ConvertDate_IOS_eLife2020(fileDate);
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1,k1] = butter(4,10/(ProcData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(4,0.5/(ProcData.notes.dsFs/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
% whisker angle
filtWhiskerAngle = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
% heart rate
heartRate = ProcData.data.heartRate;
% CBV data
HbT = ProcData.data.CBV_HbT.adjBarrels;
filtHbT = filtfilt(sos2,g2,HbT);
% cortical and hippocampal spectrograms
cortical_LHnormS = SpecData.cortical_LH.normS.*100;
hippocampusNormS = SpecData.hippocampus.normS.*100;
T = SpecData.cortical_LH.T;
F = SpecData.cortical_LH.F;
try
    % movie file comparing processed with original data
    outputVideo = VideoWriter('SupplementalVideo3_REM.mp4','MPEG-4');
    fps = 30;   % default fps from video acquisition
    speedUp = 2;   % speed up by factor of
    outputVideo.FrameRate = fps*speedUp;
    open(outputVideo);
    %%
    fig = figure('Position',get(0,'Screensize'));
    sgtitle({'Supplemental Video 3 (REM)',' '})
    % Whisker angle and heart rate
    ax12 = subplot(6,4,[1,2,3,5,6,7]);
    p1 = plot((1:length(filtWhiskerAngle))/ProcData.notes.dsFs,-filtWhiskerAngle,'color','k','LineWidth',0.5);
    hold on
    plot([100,147],[60,60],'color',colorRfcNREM,'LineWidth',5)
    plot([147,300],[60,60],'color',colorRfcREM,'LineWidth',5)
    ylabel({'Whisker','angle (deg)'})
    xlim([100,300])
    ylim([-20,61.5])
    yyaxis right
    p2 = plot((1:length(heartRate)),heartRate,'color',colors_eLife2020('deep carrot orange'),'LineWidth',0.5);
    ylabel('Heart rate (Hz)','rotation',-90,'VerticalAlignment','bottom')
    legend([p1,p2],'Whisker angle','Heart rate','Location','NorthWest','AutoUpdate','off')
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([100,150,200,250,300])
    xlim([100,300])
    ylim([5,10.25])
    ax12.TickLength = [0.01,0.01];
    ax12.YAxis(1).Color = 'k';
    ax12.YAxis(2).Color = colors_eLife2020('deep carrot orange');
    % CBV and behavioral indeces
    ax34 = subplot(6,4,[9,10,11,13,14,15]);
    p3 = plot((1:length(filtHbT))/ProcData.notes.CBVCamSamplingRate,filtHbT,'color','r','LineWidth',1);
    hold on
    p4 = plot([100,147],[100,100],'color',colorRfcNREM,'LineWidth',5);
    p5 = plot([147,300],[100,100],'color',colorRfcREM,'LineWidth',5);
    xlim([100,300])
    ylim([-45,100])
    ylabel('\Delta[HbT] (\muM)')
    legend([p3,p4,p5],'\DeltaHbT','rfc-NREM','rfc-REM','Location','NorthWest','AutoUpdate','off')
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([100,150,200,250,300])
    ax34.TickLength = [0.01,0.01];
    % Left cortical electrode spectrogram
    ax5 = subplot(6,4,[17,18,19]);
    semilog_imagesc_eLife2020(T,F,cortical_LHnormS,'y')
    axis xy
    hold on
    plot([100,147],[150,150],'color',colorRfcNREM,'LineWidth',5)
    plot([147,300],[150,150],'color',colorRfcREM,'LineWidth',5)
    c5 = colorbar;
    ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    caxis([-100,100])
    ylabel({'Cortical LFP','Freq (Hz)'})
    yticks([10,100])
    yticklabels({'10','100'})
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([100,150,200,250,300])
    xlim([100,300])
    ax5.TickLength = [0.01,0.01];
    % Hippocampal electrode spectrogram
    ax6 = subplot(6,4,[21,22,23]);
    semilog_imagesc_eLife2020(T,F,hippocampusNormS,'y')
    axis xy
    hold on
    plot([100,147],[150,150],'color',colorRfcNREM,'LineWidth',5)
    plot([147,300],[150,150],'color',colorRfcREM,'LineWidth',5)
    c6 = colorbar;
    ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    caxis([-100,100])
    xlabel('Time (s)')
    ylabel({'Hippocampal LFP','Freq (Hz)'})
    yticks([10,100])
    yticklabels({'10','100'})
    set(gca,'box','off')
    xticks([100,150,200,250,300])
    xticklabels({'0','50','100','150','200'})
    xlim([100,300])
    ax6.TickLength = [0.01,0.01];
    % Axes properties
    ax12Pos = get(ax12,'position');
    ax5Pos = get(ax5,'position');
    ax6Pos = get(ax6,'position');
    ax5Pos(3) = ax12Pos(3);
    ax6Pos(3) = ax12Pos(3);
    set(ax5,'position',ax5Pos);
    set(ax6,'position',ax6Pos);
    c6Pos = get(c6,'position');
    c5Pos = get(c5,'position');
    c5Pos(1) = c6Pos(1);
    set(c5,'position',c5Pos);
    %%
    for a = 1:size(remPupilImageStack,3)
        % lines
        axes(ax12)
        if remStartTime + (a/fs1) <= 147
            x12 = xline(remStartTime + (a/fs1),'color',colorRfcNREM,'LineWidth',2);
            axes(ax34)
            x34 = xline(remStartTime + (a/fs1),'color',colorRfcNREM,'LineWidth',2);
            axes(ax5)
            x5 = xline(remStartTime + (a/fs1),'color',colorRfcNREM,'LineWidth',2);
            axes(ax6)
            x6 = xline(remStartTime + (a/fs1),'color',colorRfcNREM,'LineWidth',2);
        elseif remStartTime + (a/fs1) > 147
            x12 = xline(remStartTime + (a/fs1),'color',colorRfcREM,'LineWidth',2);
            axes(ax34)
            x34 = xline(remStartTime + (a/fs1),'color',colorRfcREM,'LineWidth',2);
            axes(ax5)
            x5 = xline(remStartTime + (a/fs1),'color',colorRfcREM,'LineWidth',2);
            axes(ax6)
            x6 = xline(remStartTime + (a/fs1),'color',colorRfcREM,'LineWidth',2);
        end
        % whisker movie
        s1 = subplot(6,4,[4,8]);
        imagesc(dsRemWhiskImageStack(:,:,a))
        title('Whisker camera')
        colormap(gca,'gray')
        caxis([0,2^8])
        axis image
        axis off
        hold off
        % window movie
        s2 = subplot(6,4,[12,16]);
        imagesc(remCBVImageStack(:,:,a))
        title('IOS camera')
        colormap(gca,'gray')
        c0 = colorbar;
        ylabel(c0,'\DeltaR/R (%)','rotation',-90,'VerticalAlignment','bottom')
        caxis([-15,15])
        axis image
        axis off
        hold off
        % pupil movie
        s3 = subplot(6,4,[20,24]);
        imagesc(remPupilImageStack(:,:,a))
        title('Eye camera')
        colormap(gca,'gray')
        caxis([0,2^8])
        axis image
        axis off
        hold off
        currentFrame = getframe(fig);
        writeVideo(outputVideo,currentFrame);
        delete(s1)
        delete(s2)
        delete(s3)
        delete(x12)
        delete(x34)
        delete(x5)
        delete(x6)
    end
    close(outputVideo)
    close(fig)
    sendmail('kevinlturnerjr@gmail.com','Video S3 Complete');
catch
    sendmail('kevinlturnerjr@gmail.com','Video S3 Error');
end