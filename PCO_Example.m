fileID = uigetfile('*.tif'); 
% fileID = uigetfile('*.pcoraw'); 
info = imfinfo(fileID);
numberOfPages = length(info);
for k = 1 : numberOfPages
    imageStack(:,:,k) = imread(fileID,k); %#ok<SAGROW>
end	
meanImg = mean(imageStack,3);
for bb = 1:length(imageStack)
    normImgStack(:,:,bb) = ((double(imageStack(:,:,bb)) - meanImg)./meanImg)*100; %#ok<SAGROW>
end
% first frame
figure;
imagesc(imageStack(:,:,1))
colormap gray
axis image
% movie
h = implay(normImgStack);
h.Visual.ColorMap.MapExpression = 'jet';

