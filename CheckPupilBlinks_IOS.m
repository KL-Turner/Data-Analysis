function [] = CheckPupilBlinks_IOS(procDataFileID)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Track pupil diameter and blinking
%________________________________________________________________________________________________________________________

load(procDataFileID)
if strcmp(ProcData.data.Pupil.frameCheck,'y') == true %#ok<NODEF>
    if isfield(ProcData.data.Pupil,'blinkCheckComplete') == false
        % load files and extract video information
        [animalID,~,fileID] = GetFileInfo_IOS(procDataFileID);
        pupilCamFileID = [fileID '_PupilCam.bin'];
        fid = fopen(pupilCamFileID); % reads the binary file in to the work space
        fseek(fid,0,'eof'); % find the end of the video frame
        fseek(fid,0,'bof'); % find the begining of video frames
        imageHeight = ProcData.notes.pupilCamPixelHeight;
        imageWidth = ProcData.notes.pupilCamPixelWidth;
        pixelsPerFrame = imageWidth*imageHeight;
        skippedPixels = pixelsPerFrame;
        nFramesToRead = ProcData.data.Pupil.blinkInds;
        imageStackA = zeros(imageHeight,imageWidth,length(nFramesToRead));
        imageStackB = zeros(imageHeight,imageWidth,length(nFramesToRead));
        imageStackC = zeros(imageHeight,imageWidth,length(nFramesToRead));
        imageStackD = zeros(imageHeight,imageWidth,length(nFramesToRead));
        imageStackE = zeros(imageHeight,imageWidth,length(nFramesToRead));
        % first 5 frames of video
        for dd = 1:length(nFramesToRead)
            % frame - 2
            fseek(fid,(nFramesToRead(dd) - 3)*skippedPixels,'bof');
            z = fread(fid,pixelsPerFrame,'*uint8','b');
            img = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight);
            imageStackA(:,:,dd) = flip(imrotate(img,-90),2);
            % frame - 1
            fseek(fid,(nFramesToRead(dd) - 2)*skippedPixels,'bof');
            z = fread(fid,pixelsPerFrame,'*uint8','b');
            img = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight);
            imageStackB(:,:,dd) = flip(imrotate(img,-90),2);
            % frame
            fseek(fid,(nFramesToRead(dd) - 1)*skippedPixels,'bof');
            z = fread(fid,pixelsPerFrame,'*uint8','b');
            img = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight);
            imageStackC(:,:,dd) = flip(imrotate(img,-90),2);
            % frame + 1
            fseek(fid,(nFramesToRead(dd))*skippedPixels,'bof');
            z = fread(fid,pixelsPerFrame,'*uint8','b');
            img = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight);
            imageStackD(:,:,dd) = flip(imrotate(img,-90),2);
            % frame + 2
            fseek(fid,(nFramesToRead(dd) + 1)*skippedPixels,'bof');
            z = fread(fid,pixelsPerFrame,'*uint8','b');
            img = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight);
            imageStackE(:,:,dd) = flip(imrotate(img,-90),2);
        end
        % request user input for this file
        for ee = 1:length(nFramesToRead)
            check = false;
            while check == false
                % figure showing 10 frames
                imageCheck = figure;
                subplot(1,5,1)
                imagesc(imageStackA(:,:,ee))
                colormap gray
                axis image
                axis off
                subplot(1,5,2)
                imagesc(imageStackB(:,:,ee))
                colormap gray
                axis image
                axis off
                subplot(1,5,3)
                imagesc(imageStackC(:,:,ee))
                colormap gray
                axis image
                axis off
                subplot(1,5,4)
                imagesc(imageStackD(:,:,ee))
                colormap gray
                axis image
                axis off
                subplot(1,5,5)
                imagesc(imageStackE(:,:,ee))
                colormap gray
                axis image
                axis off
                keepBlink = input(['(' num2str(ee) '/' num2str(size(imageStackA,3)) ') Is this a Blink [t = ' num2str(nFramesToRead(ee)) ']? (y/n): '],'s'); disp(' ')
                close(imageCheck)
                if strcmp(keepBlink,'y') == true || strcmp(keepBlink,'n') == true
                    ProcData.data.Pupil.blinkCheck{1,ee} = keepBlink;
                    check = true;
                end
            end
        end
        ProcData.data.Pupil.blinkCheckComplete = 'y';
        save(procDataFileID,'ProcData')
    else
        for bb = 1:length(ProcData.data.Pupil.blinkInds)
            ProcData.data.Pupil.blinkCheck{1,bb} = 'n';
        end
        ProcData.data.Pupil.blinkCheckComplete = 'y';
        save(procDataFileID,'ProcData')
    end
end

end
