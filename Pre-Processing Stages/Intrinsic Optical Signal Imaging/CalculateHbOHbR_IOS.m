function [] = CalculateHbOHbR_IOS(procDataFileIDs,RestingBaselines)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner, adapted from code written by Qingguang Zhang
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
% pathlengths for each wavelength (cm)
X_Red = 0.3846483;    % 630 nm
X_Green = 0.0371713;  % 530 nm
X_Blue = 0.064898;    % 480 nm
% extinction coefficients (cm^-1*M^-1)
E_HbR_Red = 5148.8;   % 630 nm
E_HbO_Red = 610;
E_HbR_Green = 39500;  % 530 nm
E_HbO_Green = 39500;
E_HbR_Blue = 14550;   % 480 nm
E_HbO_Blue = 26629.2;
% load structure with ROI locations
ROIFileDir = dir('*_UpdatedROIs.mat');
ROIFileName = {ROIFileDir.name}';
ROIFileID = char(ROIFileName);
load(ROIFileID);
%% go through each file and do pixel-wise hbo-hbr analysis
for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    load(procDataFileID)
    imagingWavelengths = ProcData.notes.imagingWavelengths;
    if any(strcmp(imagingWavelengths,{'Red, Green, & Blue','Lime, Green, & Blue','Green & Blue','Lime & Blue'})) == true
        disp(['Adding HbO-HbR analysis to ProcData file (' num2str(aa) '/' num2str(size(procDataFileIDs,1)) ')...']); disp(' ')
        [animalID,fileDate,fileID] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(fileDate);
        windowCamFileID = [fileID '_WindowCam.bin'];
        % extract frames from window recording
        [frames] = ReadDalsaBinary_IOS(animalID,windowCamFileID);
        % split images R-G-B
        if ProcData.notes.greenFrames == 1
            greenFrames = frames(1:3:end - 1);
            blueFrames = frames(2:3:end);
            redFrames = frames(3:3:end);
        elseif ProcData.notes.greenFrames == 2
            greenFrames = frames(2:3:end);
            blueFrames = frames(3:3:end);
            redFrames = frames(1:3:end - 1);
        elseif ProcData.notes.greenFrames == 3
            greenFrames = frames(3:3:end);
            blueFrames = frames(1:3:end - 1);
            redFrames = frames(2:3:end);
        end
        % reshape cell aray to w x h x time image stack
        redImageStack = reshape(cell2mat(redFrames),ProcData.notes.CBVCamPixelWidth,ProcData.notes.CBVCamPixelHeight,length(redFrames));
        greenImageStack = reshape(cell2mat(greenFrames),ProcData.notes.CBVCamPixelWidth,ProcData.notes.CBVCamPixelHeight,length(greenFrames));
        blueImageStack = reshape(cell2mat(blueFrames),ProcData.notes.CBVCamPixelWidth,ProcData.notes.CBVCamPixelHeight,length(blueFrames));
        % normalize each stack by resting baseline frame
        normRedImageStack = (redImageStack)./RestingBaselines.Pixel.red.(strDay);
        normGreenImageStack = (greenImageStack)./RestingBaselines.Pixel.green.(strDay);
        normBlueImageStack = (blueImageStack)./RestingBaselines.Pixel.blue.(strDay);
        %% Calculate absorption coefficient for each pixel = -1/pathlength*ln(deltaR/R+1)
        Mu_Blue = -1/X_Blue.*log(normBlueImageStack); % cm^-1, natural logarithm
        Mu_Green = -1/X_Green.*log(normGreenImageStack); % cm^-1, natural logarithm
        Mu_Red = -1/X_Red.*log(normRedImageStack); % cm^-1, natural logarithm
        %% Calculate concentration of HbR & HbO for each pixel
        % Calculate concentrations (uM) using blue and green light
        C_HbR.BG = (E_HbO_Green.*Mu_Blue - E_HbO_Blue.*Mu_Green)/(E_HbO_Green.*E_HbR_Blue - E_HbO_Blue.*E_HbR_Green).*1e6; % in uM
        C_HbO.BG = (E_HbR_Green.*Mu_Blue - E_HbR_Blue.*Mu_Green)/(E_HbR_Green.*E_HbO_Blue - E_HbR_Blue.*E_HbO_Green).*1e6; % in uM
        D_HbO_HbR.BG = C_HbO.BG - C_HbR.BG; % difference between HbO and HbR
        % Calculate concentrations (uM) using green and red light
        C_HbR.RG = (E_HbO_Green.*Mu_Red - E_HbO_Red.*Mu_Green)/(E_HbO_Green.*E_HbR_Red - E_HbO_Red.*E_HbR_Green).*1e6; % in uM
        C_HbO.RG = (E_HbR_Green.*Mu_Red - E_HbR_Red.*Mu_Green)/(E_HbR_Green.*E_HbO_Red - E_HbR_Red.*E_HbO_Green).*1e6; % in uM
        D_HbO_HbR.RG = C_HbO.RG - C_HbR.RG; % difference between HbO and HbR
        % Calculate concentrations (uM) using blue and red light
        C_HbR.BR = (E_HbO_Blue.*Mu_Red - E_HbO_Red.*Mu_Blue)/(E_HbO_Blue.*E_HbR_Red - E_HbO_Red.*E_HbR_Blue).*1e6; % in uM
        C_HbO.BR = (E_HbR_Blue.*Mu_Red - E_HbR_Red.*Mu_Blue)/(E_HbR_Blue.*E_HbO_Red - E_HbR_Red.*E_HbO_Blue).*1e6; % in uM
        D_HbO_HbR.BR = C_HbO.BR - C_HbR.BR; % difference between HbO and HbR
        % HbO and HbR for each pixel
        out.C_HbR = C_HbR;
        out.C_HbO = C_HbO;
        out.D_HbO_HbR = D_HbO_HbR;
        %% extract analysis from each ROI
        hbNames = {'C_HbR','C_HbO','D_HbO_HbR'};
        wavelengthComps = {'BG','RG','BR'};
        roiNames = fieldnames(UpdatedROIs.(strDay));
        for bb = 1:length(roiNames)
            roiName = roiNames{bb,1};
            maskFig = figure;
            imagesc(frames{1});
            axis image;
            colormap gray
            circROI = drawcircle('Center',UpdatedROIs.(strDay).(roiName).circPosition,'Radius',UpdatedROIs.(strDay).(roiName).circRadius);
            mask = createMask(circROI,frames{1});
            close(maskFig)
            for cc = 1:length(hbNames)
                hbName = hbNames{1,cc};
                for dd = 1:length(wavelengthComps)
                    wavelengthComp = wavelengthComps{1,dd};
                    for ee = 1:size(out.(hbName).(wavelengthComp),3)
                        roiMask = mask.*double(out.(hbName).(wavelengthComp)(:,:,ee));
                        ProcData.data.(hbName).(wavelengthComp).(roiName)(1,ee) = mean(nonzeros(roiMask));
                    end
                end
            end
        end
       % save(procDataFileID,'ProcData')
    end
end