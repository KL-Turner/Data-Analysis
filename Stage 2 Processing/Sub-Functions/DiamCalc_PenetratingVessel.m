function mv_mpP=DiamCalc_PenetratingVessel(mv_mpP,ImageID)
% get the movie info
% draw the pa area

fname=[ImageID '.TIF'];
theimage=imread(fname,'TIFF','Index',1);
figure(2)
imagesc(double(theimage))
axis image
axis off
mv_mpP(1).nvessels=input('# of vessels on this slice:');

%get file info
Info = imfinfo(fname);
mv_mpP(1).Header.Filename=Info(1).Filename;
mv_mpP(1).Header.Frame_Width = num2str(Info(1).Width);
mv_mpP(1).Header.Frame_Height = num2str(Info(1).Height);
mv_mpP(1).Header.num_frames= length(Info);
mv_mpP(1).xsize = str2double(mv_mpP(1).Header.Frame_Width);
mv_mpP(1).ysize = str2double(mv_mpP(1).Header.Frame_Height);
mv_mpP(1).Header.Frame_Count=mv_mpP(1).Header.num_frames;

%take file info and extrac magnification and frame rate
text_hold=strread(Info(1).ImageDescription,'%s','delimiter','\n');
mag_start=strfind(text_hold{20},': ');
mv_mpP(1).Header.Magnification=text_hold{20}(mag_start+2:end);
rotation_start=strfind(text_hold{19},': ');
mv_mpP(1).Header.Rotation=text_hold{19}(rotation_start+2:end);
framerate_start=strfind(text_hold{24},': ');
mv_mpP(1).Header.Frame_Rate=(text_hold{24}(framerate_start+2:end-3));
mv_mpP(1).frame_rate=1/str2num(mv_mpP(1).Header.Frame_Rate);
%Read header and take further action based on header information
mv_mpP(1).startframe=1;
mv_mpP(1).endframe=mv_mpP(1).Header.num_frames;
mv_mpP(1).nframes_trial=NaN;%input('how many frames per trial?')

if (mv_mpP(1).Objective==1)
    microns_per_pixel=1.2953;
end
if (mv_mpP(1).Objective==2)
    microns_per_pixel=0.5595;
end
if (mv_mpP(1).Objective==3)
    microns_per_pixel=0.64;
end
if (mv_mpP(1).Objective==4)
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
Xfactor=microns_per_pixel/(str2num(mv_mpP(1).Header.Magnification(1:end-1)));

if (mv_mpP(1).nvessels>0)
    ystring=['y'];
    theinput=['n'];
    xsize=size(theimage,2);
    ysize=size(theimage,1);
    Y=repmat([1:ysize]',1,xsize);
    X=repmat([1:xsize],ysize,1);
    
    for vesselnumber=1:mv_mpP(1).nvessels
        
        mv_mpP(vesselnumber)=mv_mpP(1);
        area = imrect(gca);
        %area=impoly(gca,[1 1; 1 20;20 20;20 1]);
        while (strcmp(ystring,theinput)~=1)
            theinput=input('vessel box ok? y/n \n','s');
        end
        
        if    strcmp(ystring,theinput)
            api = iptgetapi(area);
            mv_mpP(vesselnumber).Vessel.box_position.xy=api.getPosition();
            mv_mpP(vesselnumber).Vessel.box_position.xy(3:4)=max(mv_mpP(vesselnumber).Vessel.box_position.xy(3:4));%constrain to be a square
            mv_mpP(vesselnumber).Vessel.xsize=xsize;
            mv_mpP(vesselnumber).Vessel.ysize=ysize;
            theinput=['n'];
        end
        
        mv_mpP(vesselnumber).Vessel.order=num2str(0);
        mv_mpP(vesselnumber).Xfactor=Xfactor;
        mv_mpP(vesselnumber).Vessel.order=num2str(0);
        mv_mpP(vesselnumber).VesselType=input('what is the type of this vessel?','s');
    end
end
close

end
