zap;
procDataFileID = uigetfile('*_ProcData.mat');
[animalID,fileDate,fileID] = GetFileInfo_IOS(procDataFileID);
load(procDataFileID)
cbvImageHeight = ProcData.notes.CBVCamPixelHeight;
cbvImageWidth = ProcData.notes.CBVCamPixelWidth;
[frames] = ReadDalsaBinary_IOS(animalID,[fileID '_WindowCam.bin']);
frames = frames(2:end);   % remove underexposed first frame
if strcmp(ProcData.notes.greenFrames,'odd') == true
    cbvFrames = frames(1:2:end);
    gcampFrames = frames(2:2:end);
elseif strcmp(ProcData.notes.greenFrames,'even') == true
    cbvFrames = frames(2:2:end);
    gcampFrames = frames(1:2:end);
end
%% LH mask
LH_boxFig = figure;
imagesc(frames{1})
colormap gray
caxis([0,1000])
axis image
[~,LH_xi,LH_yi] = roipoly;
ROIs.LH.xi = LH_xi;
ROIs.LH.yi = LH_yi;
LH_boxMask = roipoly(frames{1},ROIs.LH.xi,ROIs.LH.yi);
% LH_boxROI = drawrectangle;
% LH_boxPosition = round(LH_boxROI.Vertices);
% LH_boxX = LH_boxPosition(:,1);
% LH_boxY = LH_boxPosition(:,2);
% LH_boxMask = poly2mask(LH_boxX,LH_boxY,cbvImageWidth,cbvImageHeight);
close(LH_boxFig)
%% RH mask
RH_boxFig = figure;
imagesc(frames{1})
colormap gray
caxis([0,1000])
axis image
[~,RH_xi,RH_yi] = roipoly;
ROIs.RH.xi = RH_xi;
ROIs.RH.yi = RH_yi;
RH_boxMask = roipoly(frames{1},ROIs.RH.xi,ROIs.RH.yi);
% RH_boxROI = drawrectangle;
% RH_boxPosition = round(RH_boxROI.Vertices);
% RH_boxX = RH_boxPosition(:,1);
% RH_boxY = RH_boxPosition(:,2);
% RH_boxMask = poly2mask(RH_boxX,RH_boxY,cbvImageWidth,cbvImageHeight);
close(RH_boxFig)
%% apply mask to image stack
finalBinaryMask = LH_boxMask | RH_boxMask;
cbvWindowFrames = zeros(cbvImageWidth,cbvImageHeight,length(cbvFrames));
gcampWindowFrames = zeros(cbvImageWidth,cbvImageHeight,length(gcampFrames));
for aa = 1:size(cbvFrames,2)
    cbvFrame = cbvFrames{1,aa};
    gcampFrame = gcampFrames{1,aa};
    cbvWindowFrames(:,:,aa) = cbvFrame.*finalBinaryMask;
    gcampWindowFrames(:,:,aa) = gcampFrame.*finalBinaryMask;
end
%% normalize by baseline
baselineCBVframe = mean(cbvWindowFrames(:,:,50*30:55*30),3);
baselineGCaMPframe = mean(gcampWindowFrames(:,:,3000:3150),3);
normCBVframes = zeros(cbvImageWidth,cbvImageHeight,length(cbvFrames));
normGCaMPframes = zeros(cbvImageWidth,cbvImageHeight,length(gcampFrames));
for bb = 1:size(cbvWindowFrames,3)
    normCBVframes(:,:,bb) = ((cbvWindowFrames(:,:,bb) - baselineCBVframe)./baselineCBVframe).*100;
    normGCaMPframes(:,:,bb) = ((gcampWindowFrames(:,:,bb) - baselineGCaMPframe)./baselineGCaMPframe).*100;
end
%%
% handle = implay(normGCaMPframes,30);
% handle.Visual.ColorMap.Map = jet(256);
% handle.Visual.ColorMap.UserRange = 1; 
% handle.Visual.ColorMap.UserRangeMin = -20; 
% handle.Visual.ColorMap.UserRangeMax = 20;
%%
% handle = implay(normCBVframes,30);
% handle.Visual.ColorMap.Map = jet(256);
% handle.Visual.ColorMap.UserRange = 1; 
% handle.Visual.ColorMap.UserRangeMin = -10; 
% handle.Visual.ColorMap.UserRangeMax = 10;
%%
figure;
subplot(1,2,1)
sgtitle('GCaMP')
imagesc(normGCaMPframes(:,:,1500))
axis image
axis off
colormap summer
caxis([-20,20])
subplot(1,2,2)
imagesc(normGCaMPframes(:,:,2020))
axis image
axis off
colormap summer
caxis([-20,20])
figure;
sgtitle('Reflectance')
subplot(1,2,1)
imagesc(normCBVframes(:,:,1500))
axis image
axis off
colormap cool
caxis([-10,10])
colorbar
subplot(1,2,2)
imagesc(normCBVframes(:,:,2020))
axis image
axis off
colormap cool
colorbar
caxis([-10,10])