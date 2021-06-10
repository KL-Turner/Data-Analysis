function [the_image] = LoadTiffConcatenate_2P(the_tiff,the_frames)
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
tiff_file=Tiff(the_tiff);
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