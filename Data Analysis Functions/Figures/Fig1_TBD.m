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
dataStructure = 'Results_Fig1.mat';
if exist(dataStructure,'file') == 2
    load(dataStructure)
else
    filePath = [rootFolder delim 'Data' delim 'T141' delim 'Bilateral Imaging'];
    cd(filePath)
    exampleProcDataFileID = 'T141_201105_12_05_20_ProcData.mat';
    load(exampleProcDataFileID,'-mat')
    exampleSpecDataFileID = 'T141_201105_12_05_20_SpecDataA.mat';
    load(exampleSpecDataFileID,'-mat')
    exampleBaselineFileID = 'T141_RestingBaselines.mat';
    load(exampleBaselineFileID,'-mat')
    examplePupilData = 'T141_PupilData.mat';
    load(examplePupilData)
    [~,fileDate,fileID] = GetFileInfo_IOS(exampleProcDataFileID);
    pupilCamFileID = [fileID '_PupilCam.bin'];
    strDay = ConvertDate_IOS(fileDate);
    Results_Fig1.dsFs = ProcData.notes.dsFs;
    % setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
    [z1,p1,k1] = butter(4,10/(ProcData.notes.dsFs/2),'low');
    [sos1,g1] = zp2sos(z1,p1,k1);
    [z2,p2,k2] = butter(4,0.5/(ProcData.notes.dsFs/2),'low');
    [sos2,g2] = zp2sos(z2,p2,k2);
    % pupil area
    Results_Fig1.filtPupilArea = medfilt1(ProcData.data.Pupil.pupilArea);
    % blink times
    Results_Fig1.blinkTimes = ProcData.data.Pupil.blinkTimes;
    % whisker angle
    Results_Fig1.filtWhiskerAngle = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
    % EMG
    normEMG = ProcData.data.EMG.emg - RestingBaselines.manualSelection.EMG.emg.(strDay).mean;
    Results_Fig1.filtEMG = filtfilt(sos1,g1,normEMG);
    % HbT
    Results_Fig1.filtLH_HbT = filtfilt(sos2,g2,ProcData.data.CBV_HbT.adjLH);
    Results_Fig1.filtRH_HbT = filtfilt(sos2,g2,ProcData.data.CBV_HbT.adjRH);
    % hippocampal spectrogram
    Results_Fig1.hippocampusNormS = SpecData.hippocampus.normS.*100;
    Results_Fig1.T = SpecData.hippocampus.T;
    Results_Fig1.F = SpecData.hippocampus.F;
    % images
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
    % save images of interest
    Results_Fig1.images = cat(3,imageStack(:,:,1200),imageStack(:,:,4200),imageStack(:,:,7866),...
        imageStack(:,:,13200),imageStack(:,:,18510),imageStack(:,:,23458),imageStack(:,:,26332));
    % pupil tracking
    [data] = FuncRunPupilTracker(exampleProcDataFileID);
    Results_Fig1.workingImg = data.workingImg;
    Results_Fig1.x12 = data.x12;
    Results_Fig1.y12 = data.y12;
    Results_Fig1.threshImg = data.threshImg;
    Results_Fig1.pupilHistEdges = data.pupilHistEdges;
    Results_Fig1.normFit = data.normFit;
    Results_Fig1.intensityThresh = data.intensityThresh;
    Results_Fig1.saveRadonImg = data.saveRadonImg;
    Results_Fig1.overlay = data.overlay;
    cd([rootFolder delim])
    save('Results_Fig1.mat','Results_Fig1')
end
%% tracking algorithm images
% subplot for eye ROI
Fig1A = figure('Name','Figure Panel 1 - Turner et al. 2022');
subplot(2,2,1)
imagesc(Results_Fig1.workingImg)
hold on;
x1 = plot(Results_Fig1.x12,Results_Fig1.y12,'color','r','LineWidth',1');
title('ROI to measure changes in pupil area')
legend(x1,'eye ROI')
colormap gray
axis image
axis off
% subplot for ROI histrogram and threshold
subplot(2,2,2)
Results_Fig1.pupilHist = histogram(Results_Fig1.threshImg((Results_Fig1.threshImg ~= 0)),'BinEdges',Results_Fig1.pupilHistEdges,'Normalization','Probability');
hold on;
plot(Results_Fig1.pupilHist.BinEdges,Results_Fig1.normFit,'r','LineWidth',2);
xline(Results_Fig1.intensityThresh,'--m','LineWidth',1);
title('Histogram of image pixel intensities')
xlabel('Pixel intensities');
ylabel('Bin Counts');
legend({'Normalized Bin Counts','MLE fit of data','Pupil intensity threshold'},'Location','northwest');
xlim([0,256]);
axis square
% subplot for radon transform
subplot(2,2,3)
imagesc(Results_Fig1.saveRadonImg)
title('Radon transform back to image space')
colormap gray
caxis([-0.01,0.05])
axis image
axis off
% subplot for measured pupil area
subplot(2,2,4)
imagesc(Results_Fig1.overlay);
title('Calculated pupil area')
colormap gray
axis image
axis off
% save figure
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Figure Panels' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(Fig1A,[dirpath 'Fig1A']);
    set(Fig1A,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig1A'])
end
%% example pupil/eye images
Fig1B = figure('Name','Figure Panel 1 - Turner et al. 2022');
subplot(1,7,1)
imagesc(Results_Fig1.images(:,:,1));
title(['t = ' num2str(round(1200/Results_Fig1.dsFs)) ' sec'])
colormap gray
axis image
axis off
subplot(1,7,2)
imagesc(Results_Fig1.images(:,:,2));
title(['t = ' num2str(round(4200/Results_Fig1.dsFs)) ' sec'])
colormap gray
axis image
axis off
subplot(1,7,3)
imagesc(Results_Fig1.images(:,:,3));
title(['t = ' num2str(round(7866/Results_Fig1.dsFs)) ' sec'])
colormap gray
axis image
axis off
subplot(1,7,4)
imagesc(Results_Fig1.images(:,:,4));
title(['t = ' num2str(round(13200/Results_Fig1.dsFs)) ' sec'])
colormap gray
axis image
axis off
subplot(1,7,5)
imagesc(Results_Fig1.images(:,:,5));
title(['t = ' num2str(round(18510/Results_Fig1.dsFs)) ' sec'])
colormap gray
axis image
axis off
subplot(1,7,6)
imagesc(Results_Fig1.images(:,:,6));
title(['t = ' num2str(round(23458/Results_Fig1.dsFs)) ' sec'])
colormap gray
axis image
axis off
subplot(1,7,7)
imagesc(Results_Fig1.images(:,:,7));
title(['t = ' num2str(round(26332/Results_Fig1.dsFs)) ' sec'])
colormap gray
axis image
axis off
% save figure
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Figure Panels' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(Fig1B,[dirpath 'Fig1B']);
    set(Fig1B,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig1B'])
end
%% example trial
Fig1C =  figure('Name','Figure Panel 1 - Turner et al. 2022');
% pupil
ax12 = subplot(6,1,[1,2]);
p1 = plot((1:length(Results_Fig1.filtPupilArea))/Results_Fig1.dsFs,Results_Fig1.filtPupilArea,'color',colors('black'));
hold on;
s1 = scatter((1:length(Results_Fig1.blinkTimes))/Results_Fig1.dsFs,Results_Fig1.blinkTimes,'r');
x1 = xline(1200/Results_Fig1.dsFs,'g');
xline(4200/Results_Fig1.dsFs,'g')
xline(7866/Results_Fig1.dsFs,'g')
xline(13200/Results_Fig1.dsFs,'g')
xline(18510/Results_Fig1.dsFs,'g')
xline(23458/Results_Fig1.dsFs,'g')
xline(26332/Results_Fig1.dsFs,'g')
title('Pupil area')
xlabel('Time (sec)')
ylabel('\DeltaPupil area (pixels)')
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600,660,720,780,840,900])
legend([p1,s1,x1],'Pupil Area','Blinks','Rep Imgs','Location','Northwest')
% EMG and force sensor
ax3 = subplot(6,1,3);
p3 = plot((1:length(Results_Fig1.filtWhiskerAngle))/Results_Fig1.dsFs,-Results_Fig1.filtWhiskerAngle,'color',colors('black'),'LineWidth',0.5);
ylabel('Angle (deg)')
ylim([-10,50])
yyaxis right
p2 = plot((1:length(Results_Fig1.filtEMG))/Results_Fig1.dsFs,Results_Fig1.filtEMG,'color',colors('blue-violet'),'LineWidth',0.5);
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
p4 = plot((1:length(Results_Fig1.filtRH_HbT))/Results_Fig1.dsFs,Results_Fig1.filtRH_HbT,'color',colors('black'),'LineWidth',1);
hold on
p5 = plot((1:length(Results_Fig1.filtLH_HbT))/Results_Fig1.dsFs,Results_Fig1.filtLH_HbT,'color',colors('blue-green'),'LineWidth',1);
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
Semilog_ImageSC(Results_Fig1.T,Results_Fig1.F,Results_Fig1.hippocampusNormS,'y')
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
ax3Pos = get(ax3,'position');
ax6Pos = get(ax6,'position');
ax6Pos(3:4) = ax3Pos(3:4);
set(ax6,'position',ax6Pos);
% save figure
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Figure Panels' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(Fig1C,[dirpath 'Fig1C']);
    % set(Fig1C,'PaperPositionMode','auto');
    % print('-painters','-dpdf','-fillpage',[dirpath 'Fig1C'])
end

end
