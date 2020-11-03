function mv_mpP = CapLinesScan_2P(mv_mpP,ImageID)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________

% get the image information
% draw the vessel boxes

%get file info
Info = imfinfo([ImageID '.TIF']);
%take file info and extrac magnification and frame rate
mv_mpP(1).Header.Filename=Info(1).Filename;
mv_mpP(1).Header.Frame_Width = num2str(Info(1).Width);
mv_mpP(1).Header.Frame_Height = num2str(Info(1).Height);
mv_mpP(1).Header.num_frames= length(Info);
mv_mpP(1).xsize = str2double(mv_mpP(1).Header.Frame_Width);
mv_mpP(1).ysize = str2double(mv_mpP(1).Header.Frame_Height);
mv_mpP(1).Header.Frame_Count=mv_mpP(1).Header.num_frames;
%Read header and take further action based on header information
text_hold = strread(Info(1).ImageDescription,'%s','delimiter','\n');
mag_start=strfind(text_hold{20},': ');
Magnification=text_hold{20}(mag_start+2:end-1);
mv_mpP(1).Header.Magnification = str2num(Magnification);
rotation_start=strfind(text_hold{19},': ');
mv_mpP(1).Header.Rotation=text_hold{19}(rotation_start+2:end);
framerate_start=strfind(text_hold{24},': ');
mv_mpP(1).Header.Frame_Rate=(text_hold{24}(framerate_start+2:end-3));
mv_mpP(1).frame_rate=1/str2num(mv_mpP(1).Header.Frame_Rate);
mv_mpP(1).startframe=1;
mv_mpP(1).endframe=mv_mpP(1).Header.num_frames;
mv_mpP(1).nframes_trial=NaN;

if (mv_mpP(1).Objective==1) %10X
    microns_per_pixel=1.2953;
end
if (mv_mpP(1).Objective==2) %small 20X
    microns_per_pixel=0.5595;
end
if (mv_mpP(1).Objective==3) %big 20X
    microns_per_pixel=0.64;
end
if (mv_mpP(1).Objective==4) %40X
    microns_per_pixel=0.3619;
end
% if (mv_mpP(1).Objective==5) %18X
%     microns_per_pixel=0.9174;
% end
if (mv_mpP(1).Objective==5) %16x
    microns_per_pixel=0.825;
end

mv_mpP(1).microns_per_pixel=microns_per_pixel;
Pixel_clock=1/((5/4)*mv_mpP(1).frame_rate*mv_mpP(1).xsize*mv_mpP(1).ysize);
time_per_line=str2num(mv_mpP(1).Header.Frame_Width)*(5/4)*Pixel_clock*(.05*1e-6);
mv_mpP(1).Header.time_per_line=1/(mv_mpP(1).frame_rate*str2num(mv_mpP(1).Header.Frame_Height));
Xfactor=microns_per_pixel/mv_mpP(1).Header.Magnification;
mv_mpP(1).Xfactor=Xfactor;
matfilename='';
[fig, start, stop,hold1,frames,scan_type] = Display_Frames_MultiVessel_TIFF(ImageID,matfilename,mv_mpP);

mv_mpP(1).imagefig = fig;
mv_mpP(1).xstart= start;
mv_mpP(1).xstop= stop;
mv_mpP(1).frames_hold= hold1;
mv_mpP(1).nframes= frames;
mv_mpP(1).startframe=1;%start at the beginning  %input('start frame: ');
mv_mpP(1).endframe=length(imfinfo([ImageID '.tif']));%start at the end input('end frame: ');
mv_mpP(1).the_decimate=1;%input('decimation factor: '); only for oversampled data
    

mv_mpP(1).VesselType=input('what is the type of this vessel?','s');
mv_mpP(1).Vessel.order=input('what is the order of this vessel?','s');

function [imagefig,xstart,xstop,frames_hold,nframes,scan_type] = Display_Frames_MultiVessel_TIFF(fname,mfname,mpP)
%this function displays the first frames of the TIFF linescan file and gets the
%user input to determine the x range to use for the radon transform
%Read header and take further action based on header information
% gives the user the opportunity to record vessel ID and depths
% removed depth and vessel ID input and moved it to body code DK 3/14/14
    nframes=mpP.Header.Frame_Count;
    mpP.xsize = str2double(mpP.Header.Frame_Width);
    mpP.ysize = str2double(mpP.Header.Frame_Height);
    frame_height=mpP.ysize;

    frames_hold=LoadTiffConcatenate(fname,[1 5]);%
    imagefig=figure(2);
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
    scan_type=-ones(nvessels,1);
    vessel_ID=cell(nvessels,1);
    depths=-ones(nvessels,1);

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

function [the_image] = LoadTiffConcatenate(the_tiff,the_frames)
% loads the min(frames):maxframes of the_tiff, can conatenates them into
% the_image, a double matrix.  If no frame range given, load all frames.

if isempty(the_frames)
    the_info=imfinfo(the_tiff);
    start_frame=1;%min(the_frames);
    end_frame=length(the_info);
    n_frames=end_frame-start_frame+1;
else
    start_frame=min(the_frames);
    end_frame=max(the_frames);
    n_frames=end_frame-start_frame+1;
end
tiff_file=Tiff([the_tiff '.tif']);
first_frame=tiff_file.read();
tiff_height=size(first_frame,1);
tiff_width=size(first_frame,2);
the_image=zeros(tiff_height*n_frames,tiff_width);
for n=1:n_frames
    tiff_file.setDirectory(n);
    the_image((1+(n-1)*tiff_height):(n*tiff_height),:)=double(tiff_file.read());
end
tiff_file.close();
end

end
