function [mv_mpP] = GetVelocityLineScan_2P(MScanData,N_analog_channels)
%________________________________________________________________________________________________________________________
% Utilized in analysis by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Code unchanged with the exception of this block for record keeping. All rights belong to original author
%________________________________________________________________________________________________________________________

sampling_freq=MScanData.notes.frame_rate;
max_v=10000;%input('maximum velocity, in micrometers/sec: ')
windowsize=1/(sampling_freq);
angle_span=15;%number of degrees (+/-) around previous angle to look
angle_resolution=.1;%how accurate to determine the angle
the_scanmirrors=2;%we only use 6215 input('which mirrors? 1)6210(fast) 2)6215(slow)');
channel=1;
clear mpP mv_mpP thefigures save_filename analog_data
disp(['Processing File ' num2str(1)])
thisfile=mfilename('fullpath');
mpP.code=GetFunctionCode_2P(thisfile);
filename =[MScanData.notes.date '_' ([MScanData.notes.imageID])];
matfilename='';
        [mpP,function_ID]=GetVelocityRadon_MPScan_L4_cc_TIFF_2P(filename,matfilename,windowsize,angle_span, angle_resolution,max_v,channel,the_scanmirrors,MScanData);
        try
            mpP.radon_code=GetFunctionCode_2P(function_ID);
         
        catch
            disp('code read fail!')
        end
    if the_scanmirrors==1
        mpP.the_scan_mirrors=6210;
    elseif the_scanmirrors==2
        mpP.the_scan_mirrors=6215;
    end
    try
        mpP.ECoG.raw=mpP.Ch3;
    end
    %mpP.vessel=MScanData.notes.VesselNO;
    mpP.clock=clock;
    mpP.stimtime=NaN;
%     mpP.branchorder = MScanData.notes.Vessel.order;
    mpP.vesseltype = MScanData.notes.vesselType;
    mv_mpP=mpP;%put the substruct in the structure
    [mv_mpP,mpP] = MatchStructs_2P(mv_mpP,mpP);
    mv_mpP=orderfields(mv_mpP);
    mpP=orderfields(mpP);
    analog_data = load([MScanData.notes.date '_' MScanData.notes.imageID '.txt']);
    try
        mv_mpP.analog_data=analog_data;%put in the analog data
    catch
        mv_mpP.analog_data=[];
    end
    % Calculate locomotion velocity.
    maxvel=2*pi*0.06*10;
    maxVol=10;
    rate=maxvel/maxVol;
    mv_mpP.analog_data(:,3)=-mv_mpP.analog_data(:,2).*rate; 
    ntrials=1;%input
%     mv_mpP.shell=Header;
    [mv_mpP,xcfig]=linescan_xcov_velocity_TIFF_2P(mv_mpP);%calculate the velocity using the cross correlation method of Kim et al, 2012
end
