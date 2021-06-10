function [MScanData] = CapLinesScan_2P(tempData,ImageID)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________

%get file info
Info = imfinfo([ImageID '.TIF']);
%take file info and extrac magnification and frame rate
tempData.notes.Filename = Info(1).Filename;
tempData.notes.Frame_Width = num2str(Info(1).Width);
tempData.notes.Frame_Height = num2str(Info(1).Height);
tempData.notes.num_frames = length(Info(1));
tempData.notes.xsize = str2double(tempData.notes.Frame_Width);
tempData.notes.ysize = str2double(tempData.notes.Frame_Height);
tempData.notes.Frame_Count = tempData.notes.num_frames;
%Read header and take further action based on header information
text_hold = strread(Info(1).ImageDescription,'%s','delimiter','\n');
mag_start = strfind(text_hold{20},': ');
Magnification = text_hold{20}(mag_start + 2:end - 1);
tempData.notes.Magnification = str2num(Magnification);
rotation_start = strfind(text_hold{19},': ');
tempData.notes.Rotation = text_hold{19}(rotation_start+2:end);
framerate_start = strfind(text_hold{24},': ');
tempData.notes.Frame_Rate = (text_hold{24}(framerate_start + 2:end - 3));
tempData.notes.frame_rate = 1/str2num(tempData.notes.Frame_Rate);
tempData.notes.startframe = 1;
tempData.notes.endframe = tempData.notes.num_frames;
tempData.notes.nframes_trial = NaN;
% magnification uM per pixel
if (tempData.notes.objectiveID == 1) %10X
    microns_per_pixel=1.2953;
elseif (tempData.notes.objectiveID == 2) %small 20X
    microns_per_pixel=0.5595;
elseif (tempData.notes.objectiveID == 3) %big 20X
    microns_per_pixel=0.64;
elseif (tempData.notes.objectiveID == 4) %40X
    microns_per_pixel=0.3619;
elseif (tempData.notes.objectiveID == 5) %16x
    microns_per_pixel=0.825;
end
tempData.notes.microns_per_pixel = microns_per_pixel;
Pixel_clock = 1/((5/4)*tempData.notes.frame_rate*tempData.notes.xsize*tempData.notes.ysize);
time_per_line = str2num(tempData.notes.Frame_Width)*(5/4)*Pixel_clock*(.05*1e-6);
tempData.notes.time_per_line=1/(tempData.notes.frame_rate*str2num(tempData.notes.Frame_Height));
tempData.notes.LineRate = 1/time_per_line;
Xfactor = microns_per_pixel/tempData.notes.Magnification;
tempData.notes.Xfactor = Xfactor;
matfilename = '';
[~,start,stop,hold1,frames,~] = Display_Frames_MultiVessel_TIFF_2P(ImageID,matfilename,tempData);
tempData.notes.xstart = start;
tempData.notes.xstop = stop;
tempData.notes.frames_hold = hold1;
tempData.notes.nframes = frames;
tempData.notes.startframe = 1;%start at the beginning  %input('start frame: ');
tempData.notes.endframe = length(imfinfo([ImageID '.tif']));%start at the end input('end frame: ');
tempData.notes.the_decimate = 1;%input('decimation factor: '); only for oversampled data
MScanData = tempData;

end
