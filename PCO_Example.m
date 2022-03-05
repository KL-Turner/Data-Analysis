fileID = uigetfile('*.tif');
info = imfinfo(fileID);
numberOfPages = length(info);
for k = 1 : numberOfPages
    imageStack(:,:,k) = imread(fileID,k); %#ok<SAGROW>
    % Now process thisPage somehow...
end	
meanImg = mean(imageStack,3);
for bb = 1:length(imageStack)
    normImgStack(:,:,bb) = ((double(imageStack(:,:,bb)) - meanImg)./meanImg)*100; %#ok<SAGROW>
end
h = implay(normImgStack);
h.Visual.ColorMap.MapExpression = 'jet';

