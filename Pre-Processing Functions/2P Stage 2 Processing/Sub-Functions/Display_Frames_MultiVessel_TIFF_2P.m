function [imagefig,xstart,xstop,frames_hold,nframes,scan_type] = Display_Frames_MultiVessel_TIFF_2P(fname,mfname,tempData)
%this function displays the first frames of the TIFF linescan file and gets the
%user input to determine the x range to use for the radon transform
%Read header and take further action based on header information
% gives the user the opportunity to record vessel ID and depths
% removed depth and vessel ID input and moved it to body code DK 3/14/14
    nframes = tempData.notes.Frame_Count;
    tempData.xsize = str2double(tempData.notes.Frame_Width);
    tempData.ysize = str2double(tempData.notes.Frame_Height);
    frame_height = tempData.ysize;
 
    frames_hold = LoadTiffConcatenate_2P(fname,[1,5]);%
    imagefig = figure(2);
    colormap gray
    %plot the data
 
    hold off
    imagesc(double(frames_hold))
    hold on
    %axis image
    axis off
    title(fname)
    %plot the linescan trajectory
    if length(mfname)>0 %plot the
        try
            load(mfname)
            line_boundaries=find(abs(diff(scanData.pathObjSubNum)));
            scan_velocity=scanData.scanVelocity;
            for ss=1:length(line_boundaries)
                plot(line_boundaries(ss)*[1 1],[0 1000],'r')
            end
            subplot(212)
            hold off
            imagesc(scanData.im)
            axis image
            hold on
            pathImCoords(:,1) = scanData.path(:,1) * (size(scanData.im,2)-1)/(abs(diff(scanData.axisLimCol))) + 1 - (size(scanData.im,2)-1)/(abs(diff(scanData.axisLimCol)))*min(scanData.axisLimCol);
            pathImCoords(:,2) = scanData.path(:,2) * (size(scanData.im,1)-1)/(abs(diff(scanData.axisLimRow))) + 1 - (size(scanData.im,1)-1)/(abs(diff(scanData.axisLimRow)))*min(scanData.axisLimRow);
            plot(pathImCoords(:,1),pathImCoords(:,2),'r')

        end
    end
    nvessels = 1; %input('How many vessels?');
    scan_type = -ones(nvessels,1);
    vessel_ID = cell(nvessels,1);
    depths = -ones(nvessels,1);

    for vessel=1:nvessels
        disp('Select velocity calculation range.')
        [x,y]=ginput(2);
        xstart(vessel)=round(min(x));
        xstop(vessel)=round(max(x));
    %     depths(vessel)=input('depth?')
    %     vessel_ID{vessel}=input('Vessel ID: ','s')
        scan_type(vessel)=1; %input('Scan type? 0=diameter 1=velocity')
    end
end
