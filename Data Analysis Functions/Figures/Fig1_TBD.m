function [] = Fig1_TBD(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose:
%________________________________________________________________________________________________________________________

%% load data
exampleProcDataFileID = 'T141_201105_12_05_20_ProcData.mat';
load(exampleProcDataFileID,'-mat')
exampleRawDataFileID = 'T141_201105_12_05_20_RawData.mat';
load(exampleRawDataFileID,'-mat')
exampleSpecDataFileID = 'T141_201105_12_05_20_SpecDataA.mat';
load(exampleSpecDataFileID,'-mat')
exampleBaselineFileID = 'T141_RestingBaselines.mat';
load(exampleBaselineFileID,'-mat')
strDay = 'Nov05';
dsFs = ProcData.notes.dsFs;
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1,k1] = butter(4,10/(dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(4,0.5/(dsFs/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
% pupil area
filtPupilArea = medfilt1(ProcData.data.Pupil.pupilArea);
% whisker angle
filtWhiskerAngle = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
% EMG
normEMG = ProcData.data.EMG.emg - RestingBaselines.manualSelection.EMG.emg.(strDay).mean;
filtEMG = filtfilt(sos1,g1,normEMG);
% HbT
filtLH_HbT = filtfilt(sos2,g2,ProcData.data.CBV_HbT.adjLH);
filtRH_HbT = filtfilt(sos2,g2,ProcData.data.CBV_HbT.adjRH);
% hippocampal spectrogram
hippocampusNormS = SpecData.hippocampus.normS.*100;
T = SpecData.hippocampus.T;
F = SpecData.hippocampus.F;
% images
[~,~,fileID] = GetFileInfo_IOS(exampleProcDataFileID);
pupilCamFileID = [fileID '_PupilCam.bin'];
fid = fopen(pupilCamFileID); % reads the binary file in to the work space
fseek(fid,0,'eof'); % find the end of the video frame
fileSize = ftell(fid); % calculate file size
fseek(fid,0,'bof'); % find the begining of video frames
imageHeight = ProcData.notes.pupilCamPixelHeight; % how many pixels tall is the frame
imageWidth = ProcData.notes.pupilCamPixelWidth; % how many pixels wide is the frame
pixelsPerFrame = imageWidth*imageHeight;
skippedPixels = pixelsPerFrame;
nFramesToRead = floor(fileSize/(pixelsPerFrame));
imageStack = zeros(imageHeight,imageWidth,nFramesToRead);
for dd = 1:length(imageStack)
    fseek(fid,(dd - 1)*skippedPixels,'bof');
    z = fread(fid,pixelsPerFrame,'*uint8','b');
    img = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight);
    imageStack(:,:,dd) = flip(imrotate(img,-90),2);
end
fclose('all');
%% example pupil/eye images
Fig1A = figure;
subplot(1,7,1)
imagesc(imageStack(:,:,1200));
title(['t = ' num2str(round(1200/dsFs)) ' sec'])
colormap gray
axis image
axis off
subplot(1,7,2)
imagesc(imageStack(:,:,4200));
title(['t = ' num2str(round(4200/dsFs)) ' sec'])
colormap gray
axis image
axis off
subplot(1,7,3)
imagesc(imageStack(:,:,7866));
title(['t = ' num2str(round(7866/dsFs)) ' sec'])
colormap gray
axis image
axis off
subplot(1,7,4)
imagesc(imageStack(:,:,13200));
title(['t = ' num2str(round(13200/dsFs)) ' sec'])
colormap gray
axis image
axis off
subplot(1,7,5)
imagesc(imageStack(:,:,18510));
title(['t = ' num2str(round(18510/dsFs)) ' sec'])
colormap gray
axis image
axis off
subplot(1,7,6)
imagesc(imageStack(:,:,23458));
title(['t = ' num2str(round(23458/dsFs)) ' sec'])
colormap gray
axis image
axis off
subplot(1,7,7)
imagesc(imageStack(:,:,26332));
title(['t = ' num2str(round(26332/dsFs)) ' sec'])
colormap gray
axis image
axis off
% save figure
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(Fig1A,[dirpath 'Fig1A']);
end
%% example trial
Fig1B = figure;
sgtitle('Figure 1 - Turner et al. pupil tbd')
% pupil
ax12 = subplot(6,1,[1,2]);
p1 = plot((1:length(filtPupilArea))/dsFs,filtPupilArea,'color',colors('black'));
hold on;
s1 = scatter((1:length(ProcData.data.Pupil.blinkTimes))/dsFs,ProcData.data.Pupil.blinkTimes,'r');
x1 = xline(1200/dsFs,'g');
xline(4200/dsFs,'g')
xline(7866/dsFs,'g')
xline(13200/dsFs,'g')
xline(18510/dsFs,'g')
xline(23458/dsFs,'g')
xline(26332/dsFs,'g')
title('Pupil area')
xlabel('Time (sec)')
ylabel('\DeltaPupil area (pixels)')
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600,660,720,780,840,900])
legend([p1,s1,x1],'Pupil Area','Blinks','Rep Imgs','Location','Northwest')
% EMG and force sensor
ax3 = subplot(6,1,3);
p3 = plot((1:length(filtWhiskerAngle))/dsFs,-filtWhiskerAngle,'color',colors('black'),'LineWidth',0.5);
ylabel('Angle (deg)')
ylim([-10,50])
yyaxis right
p2 = plot((1:length(filtEMG))/dsFs,filtEMG,'color',colors('blue-violet'),'LineWidth',0.5);
ylabel('EMG pwr (a.u.)','rotation',-90,'VerticalAlignment','bottom')
ylim([-4,2.5])
legend([p2,p3],'Whisker Angle','EMG','Location','Northwest')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600,660,720,780,840,900])
ax3.TickLength = [0.01,0.01];
ax3.YAxis(1).Color = colors('black');
ax3.YAxis(2).Color = colors('blue-violet');
% HbT
ax45 =subplot(6,1,[4,5]);
p4 = plot((1:length(filtRH_HbT))/dsFs,filtRH_HbT,'color',colors('black'),'LineWidth',1);
hold on
p5 = plot((1:length(filtLH_HbT))/dsFs,filtLH_HbT,'color',colors('blue-green'),'LineWidth',1);
ylabel('\Delta[HbT] (\muM)')
legend([p4,p5,],'LH HbT','RH HbT','Location','Northwest')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600,660,720,780,840,900])
ylim([-35,135])
ax45.TickLength = [0.01,0.01];
% hippocampal electrode spectrogram
ax6 = subplot(6,1,6);
Semilog_ImageSC(T,F,hippocampusNormS,'y')
axis xy
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (min)')
ylabel({'Hipp LFP','Freq (Hz)'})
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600,660,720,780,840,900])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'})
ax6.TickLength = [0.01,0.01];
% axes properties
ax12Pos = get(ax12,'position');
ax6Pos = get(ax6,'position');
ax6Pos(3:4) = ax12Pos(3:4);
set(ax6,'position',ax6Pos);
% save figure
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(Fig1B,[dirpath 'Fig1B']);
    % set(Fig1B,'PaperPositionMode','auto');
    % print('-painters','-dpdf','-fillpage',[dirpath 'Fig1B'])
end

end
