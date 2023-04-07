zap;
rawDataFileID = uigetfile('*RawData.mat');
load(rawDataFileID)
fileID = uigetfile('*.pcoraw');
stimulations = RawData.data.stimulations;
info = imfinfo(fileID);
numberOfPages = length(info);
for k = 1:numberOfPages
    disp(num2str(k))
    imageStack(:,:,k) = imread(fileID,k);
end
% place circle along the most relevant region of each hemisphere
% generate image
isok = false;
while isok == false
    windowFig = figure;
    imagesc(imageStack(:,:,1))
    xlabel('Image size (pixels)')
    ylabel('Image size (pixels)')
    colormap gray
    colorbar
    axis image
    disp('Move the ROI over the desired region'); disp(' ')
    drawnow
    circRadius = 60;
    circ = drawcircle('Center',[0,0],'Radius',circRadius,'Color','r');
    checkCircle = input('Is the ROI okay? (y/n): ','s'); disp(' ')
    circPosition = round(circ.Center);
    if strcmpi(checkCircle,'y') == true
        isok = true;
        ROIs.circPosition = circPosition;
        ROIs.circRadius = circRadius;
    end
    delete(windowFig);
end
% draw circular ROIs based on XCorr for LH/RH/Barrels, then free-hand for cement ROIs
maskFig = figure;
imagesc(imageStack(:,:,1));
axis image;
colormap gray
circROI = drawcircle('Center',ROIs.circPosition,'Radius',ROIs.circRadius);
ROImask = createMask(circROI,imageStack(:,:,1));
close(maskFig)
% apply image mask to each frame of the stack
for aa = 1:size(imageStack,3)
    mask = ROImask.*double(imageStack(:,:,aa));
    refl(aa,1) = mean(nonzeros(mask));
%     fluor(aa,1) = mean(nonzeros(mask));
end
normRefl = ((refl - mean(refl))./mean(refl))*100;
% normFluor = ((fluor - mean(fluor))./mean(fluor))*100;
figure;
x1 = plot((1:length(normRefl))/30,normRefl,'k');
ylabel('\DeltaR/R0 (%)')
yyaxis right
x2 = plot((1:length(stimulations))/20000,stimulations,'r');
ylabel('Volts')
xlabel('Time (s)')
legend([x1,x2],'LED','Trigger')