function [MScanData] = GetVelocityLineScan_2P(MScanData,fileID)
%________________________________________________________________________________________________________________________
% Utilized in analysis by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Code unchanged with the exception of this block for record keeping. All rights belong to original author
%________________________________________________________________________________________________________________________

%% Takes a tiff file (movie) and an analog ascii file, and extracts the diameters
MScan_analogData = [fileID '.TXT'];
disp(['Loading MScan file: ' MScan_analogData '...']); disp(' ');
analogData = load(MScan_analogData,'-ascii');
MScanData.data.corticalNeural = analogData(:,2);
MScanData.data.EMG = analogData(:,3);
MScanData.data.hippocampalNeural = analogData(:,4);
MScanData.data.forceSensor = analogData(:,5);
MScanData.notes.analogSamplingRate = 20000;
%%
sampling_freq=MScanData.notes.frame_rate;
max_v=10000;%input('maximum velocity, in micrometers/sec: ')
windowsize=1/(sampling_freq);
angle_span=15;%number of degrees (+/-) around previous angle to look
angle_resolution=.1;%how accurate to determine the angle
the_scanmirrors=2;%we only use 6215 input('which mirrors? 1)6210(fast) 2)6215(slow)');
channel=1;
clear mpP mv_mpP thefigures save_filename analog_data
disp(['Processing File ' MScanData.notes.imageID])
thisfile=mfilename('fullpath');
MScanData.code=GetFunctionCode_2P(thisfile);
filename =[MScanData.notes.date '_' ([MScanData.notes.imageID])];
matfilename='';
        [MScanData,function_ID] = GetVelocityRadon_MPScan_L4_cc_TIFF_2P(filename,matfilename,windowsize,angle_span,angle_resolution,max_v,channel,the_scanmirrors,MScanData);
        try
            MScanData.radon_code = GetFunctionCode_2P(function_ID);        
        catch
            disp('code read fail!')
        end
    if the_scanmirrors==1
        MScanData.the_scan_mirrors=6210;
    elseif the_scanmirrors==2
        MScanData.the_scan_mirrors=6215;
    end
    MScanData = orderfields(MScanData);
    [MScanData,xcfig]=linescan_xcov_velocity_TIFF_2P(MScanData);%calculate the velocity using the cross correlation method of Kim et al, 2012
end
