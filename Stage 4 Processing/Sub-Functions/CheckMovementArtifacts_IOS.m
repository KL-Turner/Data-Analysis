function CheckMovementArtifacts(imageStack)
%_______________________________________________________________________________________________
% Written by Kevin L. Turner, Jr. 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%_______________________________________________________________________________________________
%% Create a figure to show the baseline frame in color and grey-scale
baselineFrame = mean(imageStack, 3);   % Average in the 3rd dimension (time) to create a baseline

% figure; 
% imagesc(baselineFrame);    % Show scaled baseline frame
% title('Color-scaled baseline frame (max 4096)');
% xlabel('width (pixels)');
% ylabel('height (pixels)');
% colorbar;
% axis square;
% 
% figure;
% imagesc(baselineFrame);
% colormap('gray');
% title('Grey-scale baseline frame');
% xlabel('width (pixels)');
% ylabel('height (pixels)');
% axis square;

%% Create implay movie for the desired timeframe
normalizedCBVFrames = zeros(size(imageStack));     % Pre-allocate space for normalized frames

for N = 1:size(imageStack, 3)
    normalizedCBVFrames(:,:, N) = imageStack(:,:, N)./ baselineFrame;     % Normalize by baseline
end

min = 0.95;
max = 1.05;

ImplayAutoColorMap(normalizedCBVFrames, min, max);
