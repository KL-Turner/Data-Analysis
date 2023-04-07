clear; clc;
% load saved stacks in cur dir
blueImageStackID = 'blueImageStack.mat';
blueImageStack = load(blueImageStackID,'-mat');
shortStimStackID = 'shortStimulationAverage.mat'; 
shortStimStack = load(shortStimStackID,'-mat');
longStimStackID = 'longStimulationAverage.mat';
longStimStack = load(longStimStackID,'-mat');
%% short stim figure
summaryFig = figure;
% t = 0
subplot(2,4,1)
imagesc(shortStimStack.meanFrameStack(:,:,75))
colormap jet
caxis([-100,100])
axis image
axis off
title('t = 0')
% t = 1/3 sec
subplot(2,4,2)
imagesc(shortStimStack.meanFrameStack(:,:,80))
colormap jet
caxis([-100,100])
axis image
axis off
title('t = 1/3 sec')
% t = 1 sec
subplot(2,4,3)
imagesc(shortStimStack.meanFrameStack(:,:,90))
colormap jet
caxis([-100,100])
axis image
axis off
title('t = 1 sec')
% t = 2 sec
subplot(2,4,4)
imagesc(shortStimStack.meanFrameStack(:,:,105))
colormap jet
caxis([-100,100])
axis image
axis off
title('t = 2 sec')
%% long stim figure
% t = 0
subplot(2,4,5)
imagesc(longStimStack.meanFrameStack(:,:,75))
colormap jet
caxis([-100,100])
axis image
axis off
title('t = 0')
% t = 1/3 sec
subplot(2,4,6)
imagesc(longStimStack.meanFrameStack(:,:,80))
colormap jet
caxis([-100,100])
axis image
axis off
title('t = 1/3 sec')
% t = 1 sec
subplot(2,4,7)
imagesc(longStimStack.meanFrameStack(:,:,90))
colormap jet
caxis([-100,100])
axis image
axis off
title('t = 1 sec')
% t = 2 sec
subplot(2,4,8)
imagesc(longStimStack.meanFrameStack(:,:,105))
colormap jet
caxis([-100,100])
axis image
axis off
title('t = 2 sec')