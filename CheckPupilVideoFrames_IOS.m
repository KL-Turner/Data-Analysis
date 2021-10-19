%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Track pupil diameter and blinking
%________________________________________________________________________________________________________________________

zap;
% Character list of all ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% check each pupil file's first and last frames for eye clouding
for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    load(procDataFileID)
    if isfield(ProcData.data.Pupil,'frameCheck') == false
        % load files and extract video information
        [animalID,~,fileID] = GetFileInfo_IOS(procDataFileID);
        pupilCamFileID = [fileID '_PupilCam.bin'];
        fid = fopen(pupilCamFileID); % reads the binary file in to the work space
        fseek(fid,0,'eof'); % find the end of the video frame
        fileSize = ftell(fid); % calculate file size
        fseek(fid,0,'bof'); % find the begining of video frames
        imageHeight = ProcData.notes.pupilCamPixelHeight;
        imageWidth = ProcData.notes.pupilCamPixelWidth;
        pixelsPerFrame = imageWidth*imageHeight;
        skippedPixels = pixelsPerFrame;
        nFrames = floor(fileSize/(pixelsPerFrame));
        nFramesToRead_A = 1:5;
        nFramesToRead_B = (nFrames - 4):nFrames;
        imageStack = zeros(imageHeight,imageWidth,(length(nFramesToRead_A) + length(nFramesToRead_B)));
        % first 5 frames of video
        cc = 1;
        for dd = 1:length(nFramesToRead_A)
            fseek(fid,(nFramesToRead_A(dd) - 1)*skippedPixels,'bof');
            z = fread(fid,pixelsPerFrame,'*uint8','b');
            img = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight);
            imageStack(:,:,cc) = flip(imrotate(img,-90),2);
            cc = cc + 1;
        end
        % last 5 frames of video
        for ee = 1:length(nFramesToRead_B)
            fseek(fid,(nFramesToRead_B(ee) - 1)*skippedPixels,'bof');
            z = fread(fid,pixelsPerFrame,'*uint8','b');
            img = reshape(z(1:pixelsPerFrame),imageWidth,imageHeight);
            imageStack(:,:,cc) = flip(imrotate(img,-90),2);
            cc = cc + 1;
        end
        % figure showing 10 frames
        imageCheck = figure;
        sgtitle([animalID strrep(fileID,'_',' ')])
        for ff = 1:size(imageStack,3)
            subplot(2,5,ff)
            imagesc(imageStack(:,:,ff))
            colormap gray
            axis image
            axis off
        end
        % request user input for this file
        check = false;
        while check == false
            keepFigure = input(['(' num2str(aa) '/' num2str(size(procDataFileIDs,1)) ') Is eye data good for this session? (y/n): '],'s'); disp(' ')
            if strcmp(keepFigure,'y') == true || strcmp(keepFigure,'n') == true
                ProcData.data.Pupil.frameCheck = keepFigure;
                save(procDataFileID,'ProcData')
                % save the figure to directory.
                [pathstr,~,~] = fileparts(cd);
                dirpath = [pathstr '/Figures/Pupil Frame Check/'];
                if ~exist(dirpath,'dir')
                    mkdir(dirpath);
                end
                savefig(imageCheck,[dirpath animalID '_' fileID '_FrameCheck']);
                close(imageCheck)
                check = true;
            end
        end
    end
end