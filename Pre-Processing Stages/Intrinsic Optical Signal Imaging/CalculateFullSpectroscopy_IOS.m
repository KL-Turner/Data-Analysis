function [] = CalculateFullSpectroscopy_IOS(procDataFileIDs,RestingBaselines)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner, adapted from code written by Qingguang Zhang
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
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
        if isfield(ProcData.data,'HbO') == false
            [animalID,fileDate,fileID] = GetFileInfo_IOS(procDataFileID);
            strDay = ConvertDate_IOS(fileDate);
            windowCamFileID = [fileID '_WindowCam.bin'];
            conv2um = 1e6;
            %% extract frames from window recording
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
            normRedImageStack = ((redImageStack - RestingBaselines.Pixel.red.(strDay))./RestingBaselines.Pixel.red.(strDay)) + 1;
            normGreenImageStack = ((greenImageStack - RestingBaselines.Pixel.green.(strDay))./RestingBaselines.Pixel.green.(strDay)) + 1;
            normBlueImageStack = ((blueImageStack - RestingBaselines.Pixel.blue.(strDay))./RestingBaselines.Pixel.blue.(strDay)) + 1; %#ok<*NASGU>
            %% calculate absorption coefficient for each pixel = -1/pathlength*ln(deltaR/R+1)
            % pathlengths for each wavelength (cm)
            X_Red = 0.3846483;    % 630 nm
            X_Green = 0.0371713;  % 530 nm
            X_Blue = 0.064898;    % 480 nm
            % extinction coefficients (cm^-1*M^-1)
            E_HbO_Red = 610;      % 630 nm
            E_HbR_Red = 5148.8;
            E_HbO_Green = 39500;  % 530 nm
            E_HbR_Green = 39500;
            E_HbO_Blue = 14550;   % 480 nm
            E_HbR_Blue = 26629.2;
            Mu_Red = -1/X_Red*log(normRedImageStack); % cm^-1, natural logarithm
            Mu_Green = -1/X_Green*log(normGreenImageStack); % cm^-1, natural logarithm
            Mu_Blue = -1/X_Blue*log(normBlueImageStack); % cm^-1, natural logarithm
            %% calculate concentration of HbR & HbO for each pixel
            % Calculate concentrations (uM) using blue and green light
            HbR.BG = (E_HbO_Green*Mu_Blue - E_HbO_Blue*Mu_Green)./(E_HbO_Green*E_HbR_Blue - E_HbO_Blue*E_HbR_Green)*conv2um; % in uM
            HbO.BG = (E_HbR_Green*Mu_Blue - E_HbR_Blue*Mu_Green)./(E_HbR_Green*E_HbO_Blue - E_HbR_Blue*E_HbO_Green)*conv2um; % in uM
            % Calculate concentrations (uM) using green and red light
            HbR_ImageStack = (E_HbO_Green*Mu_Red - E_HbO_Red*Mu_Green)./(E_HbO_Green*E_HbR_Red - E_HbO_Red*E_HbR_Green)*conv2um; % in uM
            HbO_ImageStack = (E_HbR_Green*Mu_Red - E_HbR_Red*Mu_Green)./(E_HbR_Green*E_HbO_Red - E_HbR_Red*E_HbO_Green)*conv2um; % in uM
            % Calculate concentrations (uM) using blue and red light
            HbR.BR = (E_HbO_Blue*Mu_Red - E_HbO_Red*Mu_Blue)./(E_HbO_Blue*E_HbR_Red - E_HbO_Red*E_HbR_Blue)*conv2um; % in uM
            HbO.BR = (E_HbR_Blue*Mu_Red - E_HbR_Red*Mu_Blue)./(E_HbR_Blue*E_HbO_Red - E_HbR_Red*E_HbO_Blue)*conv2um; % in uM
            %% HbT analysis
            if any(strcmp(imagingWavelengths,{'Green','Green & Blue','Red, Green, & Blue'})) == true
                ledType = 'M530L3';
                bandfilterType = 'FB530-10';
                cutfilterType = 'MF525-39';
                [~,~,weightedcoeffHbT] = GetHbcoeffs_IOS(ledType,bandfilterType,cutfilterType);
            elseif any(strcmp(imagingWavelengths,{'Lime','Lime & Blue','Red, Lime, & Blue'})) == true
                ledType = 'M565L3';
                bandfilterType = 'FB570-10';
                cutfilterType = 'EO65160';
                [~,~,weightedcoeffHbT] = GetHbcoeffs_IOS(ledType,bandfilterType,cutfilterType);
            else
                weightedcoeffHbT = NaN;
            end
            HbT_ImageStack = (log(greenImageStack./RestingBaselines.Pixel.green.(strDay))).*weightedcoeffHbT*conv2um;
            %% correct hemodynamic attenuation of GCaMP fluorescence
            scale = RestingBaselines.Pixel.green.(strDay)./RestingBaselines.Pixel.blue.(strDay);
            correctedGCaMP_ImageStack = (blueImageStack./greenImageStack).*scale;
            %% extract pixel-wise ROIs
            ROInames = fieldnames(UpdatedROIs.(strDay));
            for bb = 1:length(ROInames)
                roiName = ROInames{bb,1};
                maskFig = figure;
                imagesc(greenFrames{1});
                axis image;
                colormap gray
                circROI = drawcircle('Center',UpdatedROIs.(strDay).(roiName).circPosition,'Radius',UpdatedROIs.(strDay).(roiName).circRadius);
                mask = createMask(circROI,greenFrames{1});
                close(maskFig)
                for cc = 1:size(HbT_ImageStack,3)
                    HbT_roiMask = mask.*double(HbT_ImageStack(:,:,cc));
                    ProcData.data.HbT.(roiName)(1,cc) = mean(nonzeros(HbT_roiMask));
                    HbO_roiMask = mask.*double(HbO_ImageStack(:,:,cc));
                    ProcData.data.HbO.(roiName)(1,cc) = mean(nonzeros(HbO_roiMask));
                    HbR_roiMask = mask.*double(HbR_ImageStack(:,:,cc));
                    ProcData.data.HbR.(roiName)(1,cc) = mean(nonzeros(HbR_roiMask));
                    GCaMP_roiMask = mask.*double(correctedGCaMP_ImageStack(:,:,cc));
                    ProcData.data.GCaMP.(roiName)(1,cc) = mean(nonzeros(GCaMP_roiMask));
                end
            end
            % save data
            save(procDataFileID,'ProcData')
        end
    end
end