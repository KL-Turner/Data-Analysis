function [AnalysisResults] = Fig1_S1_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate figure panel 1-S1 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% set-up and process data
if isfield(AnalysisResults,'CrossCorrROI') == false
    dataDir = [rootFolder '\Summary Figures and Structures\Cross Correlation ROI\'];
    cd(dataDir)
    % character list of RawData files
    rawDataFileStruct = dir('*_RawData.mat');
    rawDataFiles = {rawDataFileStruct.name}';
    rawDataFileIDs = char(rawDataFiles);
    % character list of ProcData files
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
    % character list of WindowCam files
    windowCamFileStruct = dir('*_WindowCam.bin');
    windowCamFiles = {windowCamFileStruct.name}';
    windowCamFileIDs = char(windowCamFiles);
    % animal ID
    [animalID,fileDate,~] = GetFileInfo_IOS_eLife2020(rawDataFileIDs(1,:));
    strDay = ConvertDate_IOS_eLife2020(fileDate);
    hem = 'Barrels';
    % take the cross correlation between the neural data and each pixel
    for aa = 1:size(rawDataFileIDs,1)
        rawDataFileID = rawDataFileIDs(aa,:);
        load(rawDataFileID,'-mat')
        procDataFileID = procDataFileIDs(aa,:);
        load(procDataFileID,'-mat')
        windowCamFileID = windowCamFileIDs(aa,:);
        % read camera frames from binary file
        [frames] = ReadDalsaBinary_IOS_eLife2020(animalID,windowCamFileID);
        if aa == 1
            % open figure for ROI drawing
            windowFig = figure;
            imagesc(frames{1})
            title([animalID ' CBV camera'])
            xlabel('Image size (pixels)')
            ylabel('Image size (pixels)')
            colormap gray
            colorbar
            axis image
            caxis([0,2^RawData.notes.CBVCamBitDepth])
            % draw ROI for the mask over the entire windows
            isok = false;
            while isok == false
                disp('Drag a rectangular ROI over the window '); disp(' ')
                [~,rect] = imcrop;
                hold on;
                ROIoutline = rectangle('Position',rect,'EdgeColor','r');
                checkMask = input('Is the ROI okay? (y/n): ','s'); disp(' ')
                if strcmp(checkMask,'y') == true
                    isok = true;
                    ROIs.([hem '_' strDay]).rect = rect;
                end
                delete(ROIoutline);
            end
            close(windowFig)
        end
        % extract the pixel values from the window ROIs
        imageMask = nan(size(frames{1}));
        rectMask = ROIs.([hem '_' strDay]).rect;
        rectMask = round(rectMask);
        imageMask(rectMask(2):(rectMask(2) + rectMask(4)),rectMask(1):(rectMask(1) + rectMask(3))) = 1;
        for c = 1:length(frames) - 1
            frame = frames{1,c};
            frameHold = double(frame).*imageMask;
            imageStack.(hem)(:,c) = frameHold(~isnan(frameHold));
        end
        % extract and process gamma band power
        [z,p,k] = butter(4,1/(ProcData.notes.dsFs/2),'low');
        [sos,g] = zp2sos(z,p,k);
        LH_gammaBandPower = detrend(filtfilt(sos,g,ProcData.data.cortical_LH.gammaBandPower - ProcData.data.cortical_LH.gammaBandPower(1)) + ProcData.data.cortical_LH.gammaBandPower(1),'constant');
        RH_gammaBandPower = detrend(filtfilt(sos,g,ProcData.data.cortical_RH.gammaBandPower - ProcData.data.cortical_RH.gammaBandPower(1)) + ProcData.data.cortical_RH.gammaBandPower(1),'constant');
        LH_muaPower = detrend(filtfilt(sos,g,ProcData.data.cortical_LH.muaPower - ProcData.data.cortical_LH.muaPower(1)) + ProcData.data.cortical_LH.muaPower(1),'constant');
        RH_muaPower = detrend(filtfilt(sos,g,ProcData.data.cortical_RH.muaPower - ProcData.data.cortical_RH.muaPower(1)) + ProcData.data.cortical_RH.muaPower(1),'constant');
        % cross correlation
        lagTime = 5;   % seconds
        maxLag = lagTime*RawData.notes.CBVCamSamplingRate;
        singleHem = ProcData.notes.hemisphere;
        if strcmp(singleHem,'LH') == true
            gammaBandArray = LH_gammaBandPower;
            muaArray = LH_muaPower;
        elseif strcmp(singleHem,'RH') == true
            gammaBandArray = RH_gammaBandPower;
            muaArray = RH_muaPower;
        end
        % extract pixel values from each numel index in matrix image
        for e = 1:size(imageStack.(hem),1)
            disp(['Analyzing matrix numel ' num2str(e) ' of ' num2str(size(imageStack.(hem),1))]); disp(' ')
            pixelArray = imageStack.(hem)(e,:);
            pixelArray = detrend(filtfilt(sos,g,pixelArray - pixelArray(1)) + pixelArray(1),'constant');
            % gamma band power cross-correlation
            [gammaXCorrVals,lags] = xcorr(pixelArray,gammaBandArray,maxLag,'coeff');
            zeroPoint = find(lags == 0);
            validGammaVals = gammaXCorrVals(zeroPoint:zeroPoint + 45);
            gammaMaxCorr = min(validGammaVals);
            if isnan(gammaMaxCorr) == true
                gammaCorrMatrix.(hem)(1,e) = 0;
            else
                gammaCorrMatrix.(hem)(1,e) = gammaMaxCorr;
            end
            % multi unit activity cross-correlation
            [muaXCorrVals,~] = xcorr(pixelArray,muaArray,maxLag,'coeff');
            muaValidVals = muaXCorrVals(zeroPoint:zeroPoint + 45);
            muaMaxCorr = min(muaValidVals);
            if isnan(muaMaxCorr) == true
                muaCorrMatrix.(hem)(1,e) = 0;
            else
                muaCorrMatrix.(hem)(1,e) = muaMaxCorr;
            end
        end
        % proper size of the ROI based on camera/lens magnification
        circRadius = 30;   % pixels ~= 1 mm diameter
        % place circle along the most correlation region of each hemisphere
        rectMask = ROIs.([hem '_' strDay]).rect;
        rectMask = round(rectMask);
        imgWidth = rectMask(3) + 1;
        imgHeight = rectMask(4) + 1;
        gammaCorrImg(:,:,aa) = reshape(gammaCorrMatrix.(hem),imgHeight,imgWidth); %#ok<*AGROW>
        muaCorrImg(:,:,aa) = reshape(muaCorrMatrix.(hem),imgHeight,imgWidth); %#ok<*AGROW>
    end
    % save results
    AnalysisResults.CrossCorrROI.gammaCorrImgAvg = mean(gammaCorrImg,3);
    AnalysisResults.CrossCorrROI.muaCorrImgAvg = mean(muaCorrImg,3);
    % generate image
    isok = false;
    while isok == false
        windowFig = figure;
        imagesc(AnalysisResults.CrossCorrROI.gammaCorrImgAvg)
        title('Peak pixel correlations')
        xlabel('Image size (pixels)')
        ylabel('Image size (pixels)')
        colormap parula
        colorbar
        axis image
        disp('Move the ROI over the most correlated region'); disp(' ')
        circ = drawcircle('Center',[0,0],'Radius',circRadius,'Color','r');
        checkCircle = input('Is the ROI okay? (y/n): ','s'); disp(' ')
        circPosition = round(circ.Center);
        if strcmp(checkCircle,'y') == true
            isok = true;
            rectBottomLeftCorner = [rectMask(1),rectMask(2) + rectMask(4)];
            rectTopLeftCorner = [rectMask(1),rectMask(2)];
            circPositionEdit = [rectBottomLeftCorner(1) + circPosition(1),rectTopLeftCorner(2) + circPosition(2)];
            ROIs.([hem '_' strDay]).circPosition = circPositionEdit;
            ROIs.([hem '_' strDay]).circRadius = circRadius;
        end
        delete(windowFig);
    end
    % save results
    AnalysisResults.CrossCorrROI.frame = frames{1};
    AnalysisResults.CrossCorrROI.ROIs = ROIs;
    % image mask to set background pixels to zero
    testFig = figure;
    imagesc(AnalysisResults.CrossCorrROI.gammaCorrImgAvg.*-1)
    title('Gamma-band power vs. pixel reflectance')
    xlabel('Image size (pixels)')
    ylabel('Image size (pixels)')
    caxis([.25,0.5])
    axis image
    set(gca,'Ticklength',[0,0])
    set(gca,'box','off')
    AnalysisResults.CrossCorrROI.imgMask = roipoly;
    close(testFig)
    save([rootFolder '\AnalysisResults.mat\'],'AnalysisResults')
end
strDay = 'Mar01';
%% Fig. 1-S1
summaryFigure = figure('Name','Fig1-S1 (a-c)');
sgtitle('Figure 1-S1 - Turner et al. 2020')
%% [1-S1a] original image with circular ROI
ax1 = subplot(1,3,1);
imagesc(AnalysisResults.CrossCorrROI.frame)
hold on;
drawcircle('Center',AnalysisResults.CrossCorrROI.ROIs.(['Barrels_' strDay]).circPosition,'Radius',AnalysisResults.CrossCorrROI.ROIs.(['Barrels_' strDay]).circRadius,'Color','r');
title({'[1-S1a] 1 mm OD ROI placement','on example window'})
xlabel('Image width (pixels)')
ylabel('Image height (pixels)')
set(gca,'YDir','reverse')
colormap(ax1,gray)
c1 = colorbar;
ylabel(c1,'Pixel intensity (12-bit)','rotation',-90,'VerticalAlignment','bottom')
caxis([0,2^12])
axis image
set(gca,'Ticklength',[0,0])
set(gca,'box','off')
%% [1-S1b] gamma cross correlation image
ax2 = subplot(1,3,2);
gammaImg = AnalysisResults.CrossCorrROI.gammaCorrImgAvg.*-1;
gammaImg(AnalysisResults.CrossCorrROI.imgMask == 0) = 0;
imagesc(gammaImg)
title({'[1-S1b] Gamma-band vs. \DeltaR','Cross correlation'})
xlabel('Image width (pixels)')
ylabel('Image height (pixels)')
colormap(ax2,parula)
c2 = colorbar;
ylabel(c2,{'Corr. coefficient';'Gamma-band vs. \DeltaR'},'rotation',-90,'VerticalAlignment','bottom')
caxis([.25,0.5])
axis image
set(gca,'Ticklength',[0,0])
set(gca,'box','off')
%% [1-S1c] MUA cross correlation image
ax3 = subplot(1,3,3);
muaImg = AnalysisResults.CrossCorrROI.muaCorrImgAvg*-1;
muaImg(AnalysisResults.CrossCorrROI.imgMask == 0) = 0;
imagesc(muaImg)
title({'[1-S1c] MUA vs. \DeltaR','Cross-correlation'})
xlabel('Image width (pixels)')
ylabel('Image height (pixels)')
colormap(ax3,parula)
c3 = colorbar;
ylabel(c3,{'Corr. coefficient';'MUA vs. \DeltaR'},'rotation',-90,'VerticalAlignment','bottom')
caxis([0.1,0.35])
axis image
set(gca,'Ticklength',[0,0])
set(gca,'box','off')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig1-S1']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Fig1-S1'])
end

end
