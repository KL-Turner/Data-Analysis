function [MScanData]=GetVelocityRadon_MPScan_L4_cc_TIFF_2P(fname,matfilename,window_time_size,angle_span, angle_resolution, max_v, channel,isinterleaved,MScanData);




theframes = [MScanData.notes.startframe,MScanData.notes.endframe];
the_x = [MScanData.notes.xstart,MScanData.notes.xstop];
the_decimate = MScanData.notes.the_decimate;
user_microns_per_pixel = MScanData.notes.microns_per_pixel;
time_per_line = MScanData.notes.time_per_line;
isinterleaved = 1;
Fs_Blood_Flow = MScanData.notes.LineRate;
windowsizea = round(Fs_Blood_Flow*window_time_size);%*.25 %for interlevaed lines
MScanData.notes.num_frames = theframes(end) - theframes(1);%str2double(mpContent.Header.notes.Frame_Count);
%adapt windowsize to interleaved/not interleaved
if isinterleaved==1
    windowsize=8*round(windowsizea/8);
else
    windowsize=4*round(windowsizea/4);
end
MScanData.Blood_flow.windowsize = windowsize;
% basic stuff
nframes = theframes(2) - theframes(1) + 1;
stepsize = .25*windowsize;
nlines = nframes*(MScanData.notes.ysize);
npoints = max(the_x) - min(the_x) + 1;
nsteps = floor(nlines/stepsize) - 3;
angles = (1:180);
angles_baseline = 1:15:180;
angles_adaptive = [-angle_span:1:angle_span]; %angle range to look over
angles_fine = -1.75:angle_resolution:1.75; % fine grained search
n_adaptive_angles=length(angles_adaptive);
theta_low = 5;%boundaries for for adaptivce theta
theta_high = n_adaptive_angles - 5;
spread_matrix = zeros(1,length(angles));
spread_matrix_adaptive = zeros(2,-nsteps,length(angles_adaptive));
baseline_var = zeros(2,nsteps,length(angles_baseline));
spread_matrix_fine = zeros(2,nsteps,length(angles_fine));
thetas = zeros(2,nsteps);%the angle of the lines
data_var = zeros(nsteps,1);
data_hold = zeros(windowsize,npoints);
%data caching
framerate = round(MScanData.notes.LineRate/MScanData.notes.ysize);
nframes_to_cache = round(MScanData.notes.LineRate/MScanData.notes.ysize);
data_temp = LoadTiffConcatenate_2P([fname '.tif'],theframes);%load froma tiff file
data_cache = data_temp(:,the_x(1):the_decimate:the_x(2) - 1);
MScanData.Blood_flow.mean_BG = (mean(data_cache));
cached_lines = 1:(nframes_to_cache*MScanData.notes.ysize);
use_lines = 1:windowsize;
if isinterleaved == 1
    data_hold=zeros(2,windowsize/2,length([the_x(1):the_decimate:the_x(2) - 1]));
    data_hold(1,:,:) = double(data_cache(1:2:windowsize,1:the_decimate:end));
    data_hold(2,:,:) = double(data_cache(2:2:windowsize,end:-the_decimate:1));
else
    data_hold = zeros(1,windowsize,length([the_x(1):the_decimate:the_x(2) - 1]));
    data_hold(1,:,:) = double(data_cache(1:windowsize,1:the_decimate:end));
end
line_counter = 4*stepsize;%this is the end of the window
data_cachemean_forward=mean(data_cache);
if isinterleaved == 1
    data_cachemean_backward = (fliplr(data_cachemean_forward));
end
if isinterleaved == 1
    for k = 1:nsteps
        frame_line_marker = line_counter;%find the new place in the frame
        data_hold(1,:,:) = double(data_cache(frame_line_marker-4*stepsize+1:2:frame_line_marker-1,1:the_decimate:end));
        data_hold(2,:,:) = double(data_cache(frame_line_marker+2-4*stepsize:2:frame_line_marker,end:-the_decimate:1));
        use_lines=use_lines+stepsize;
        %take out the local mean of the window
        the_t(k) = 1 + (k - 1)*stepsize+windowsize/2;
        npoints_decimated = size(data_hold,3);
        for n = 1:npoints_decimated
            data_hold_ms(1,:,n)=data_hold(1,:,n) - data_cachemean_forward(n);
            data_hold_ms(2,:,n)=data_hold(2,:,n) - data_cachemean_backward(n);
        end
        data_hold_ms = data_hold_ms - mean(data_hold_ms(:));
        if k == 1
            radon_hold(1,:,:) = (radon(squeeze(data_hold_ms(1,:,:)),angles));
            radon_hold(2,:,:) = radon(squeeze(data_hold_ms(2,:,:)),angles);
            spread_matrix(1,:) = var(squeeze(radon_hold(1,:,:)));
            spread_matrix(2,:) = var(squeeze(radon_hold(2,:,:)));
            [m(1,:),the_theta(1,:)] = max(spread_matrix(1,k));
            [m(2,:),the_theta(2,:)] = max(spread_matrix(2,k));
            thetas(1,k) = angles(the_theta(1,:));
            thetas(2,k) = angles(the_theta(2,:));
            radon_hold_adaptive(1,:,:) = radon(squeeze(data_hold_ms(1,:,:)),-thetas(1,k)+angles_adaptive);
            radon_hold_adaptive(2,:,:) = radon(squeeze(data_hold_ms(2,:,:)),-thetas(2,k)+angles_adaptive);
            spread_matrix_adaptive(1,k,1) = max(var(squeeze(radon_hold_adaptive(1,:,:))));
            spread_matrix_adaptive(2,k,1) = max(var(squeeze(radon_hold_adaptive(2,:,:))));
            radon_hold_fine(1,:,:) = radon(squeeze(data_hold_ms(1,:,:)-mean(data_hold_ms(:))),thetas(1,k)+angles_fine);
            radon_hold_fine(2,:,:) = radon(squeeze(data_hold_ms(2,:,:)-mean(data_hold_ms(:))),thetas(2,k)+angles_fine);
            spread_matrix_fine(1,k,:) = squeeze(var(squeeze(radon_hold_fine(1,:,:))));
            spread_matrix_fine(2,k,:) = squeeze(var(squeeze(radon_hold_fine(2,:,:))));
            [m the_theta(1,:)] = max(squeeze(spread_matrix_fine(1,k,:)));
            [m the_theta(2,:)] = max(squeeze(spread_matrix_fine(2,k,:)));
            thetas(1,k) = thetas(1,k)+angles_fine(the_theta(1,:));
            thetas(2,k) = thetas(2,k)+angles_fine(the_theta(2,:));
            baseline_var(1,k,:) = var(radon(squeeze(data_hold_ms(1,:,:)),angles_baseline));
            baseline_var(2,k,:)=var(radon(squeeze(data_hold_ms(2,:,:)),angles_baseline));
        else
            thetas(1,k)=thetas(1,k-1);
            thetas(2,k)=thetas(2,k-1);
            radon_hold_adaptive(1,:,:)=radon(squeeze(data_hold_ms(1,:,:)),thetas(1,k-1)+angles_adaptive);
            radon_hold_adaptive(2,:,:)=radon(squeeze(data_hold_ms(2,:,:)),thetas(2,k-1)+angles_adaptive);
            spread_matrix_adaptive(1,k,:)=var(squeeze(radon_hold_adaptive(1,:,:)));
            spread_matrix_adaptive(2,k,:)=var(squeeze(radon_hold_adaptive(2,:,:)));
            [m the_theta(1,:)]=max(spread_matrix_adaptive(1,k,:));
            [m the_theta(2,:)]=max(spread_matrix_adaptive(2,k,:));
            baseline_var(1,k,:)=var(radon(squeeze(data_hold_ms(1,:,:)),angles_baseline));%figure out the average radon of the section for calculating S/N
            baseline_var(2,k,:)=var(radon(squeeze(data_hold_ms(2,:,:)),angles_baseline));%figure out the average radon of the section for calculating S/N
            %if the peak of the variance is at an extreme, redo for angles centered
            %around the peak
            
            
            while ((the_theta(1,:)<=theta_low)||(the_theta(1,:)>=theta_high))
                thetas(1,k)=thetas(1,k)+angles_adaptive(the_theta(1,:));
                radon_hold_adaptive(1,:,:)=radon(squeeze(data_hold_ms(1,:,:)), thetas(1,k)+angles_adaptive);
                spread_matrix_adaptive(1,k,:)=var(squeeze(radon_hold_adaptive(1,:,:)));
                [m the_theta(1,:)]=max(squeeze(spread_matrix_adaptive(1,k,:)));
            end
            while ((the_theta(2,:)<=theta_low)||(the_theta(2,:)>=theta_high))
                thetas(2,k)=thetas(2,k)+angles_adaptive(the_theta(2,:));
                radon_hold_adaptive(2,:,:)=radon(squeeze(data_hold_ms(2,:,:)), thetas(2,k)+angles_adaptive);
                spread_matrix_adaptive(2,k,:)=var(squeeze(radon_hold_adaptive(2,:,:)));
                [m the_theta(2,:)]=max(squeeze(spread_matrix_adaptive(2,k,:)));
            end
            
            baseline_var(1,k,:)=var(radon(squeeze(data_hold_ms(1,:,:)),angles_baseline));%
            baseline_var(2,k,:)=var(radon(squeeze(data_hold_ms(2,:,:)),angles_baseline));%
            thetas(1,k)=thetas(1,k)+angles_adaptive(the_theta(1,:));
            radon_hold_fine(1,:,:)=radon(squeeze(data_hold_ms(1,:,:)),thetas(1,k)+angles_fine);
            spread_matrix_fine(1,k,:)=var(squeeze(radon_hold_fine(1,:,:)));
            [m the_theta(1,:)]=max(spread_matrix_fine(1,k,:));
            thetas(1,k)=thetas(1,k)+angles_fine(the_theta(1,:));
            thetas(2,k)=thetas(2,k)+angles_adaptive(the_theta(2,:));
            radon_hold_fine(2,:,:)=radon(squeeze(data_hold_ms(2,:,:)),thetas(2,k)+angles_fine);
            spread_matrix_fine(2,k,:)=var(squeeze(radon_hold_fine(2,:,:)));
            [m the_theta(2,:)]=max(spread_matrix_fine(2,k,:));
            thetas(2,k)=thetas(2,k)+angles_fine(the_theta(2,:));
        end
        line_counter=line_counter+stepsize;
    end
    
if((length(find(thetas<0))>=1)&&(length(find(thetas>0))>=1))
    disp('sign flip!') 
end
size(spread_matrix_adaptive);
MScanData.Blood_flow.thetas=thetas;
MScanData.Blood_flow.unscaled_velocity=(cotd(thetas-90));
MScanData.Blood_flow.spread_matrix_adaptive=spread_matrix_adaptive;
MScanData.Blood_flow.baseline_var=baseline_var;
Xfactor = MScanData.notes.Xfactor;
Tfactor = MScanData.notes.Tfactor;
MScanData.Blood_flow.the_decimate=the_decimate;
if isinterleaved ==1 %need to divide by 2 because of splitting due to inerleaving
    MScanData.Blood_flow.v(1,:)=the_decimate*Tfactor*Xfactor*(cotd(thetas(1,:)-90))/2;
    MScanData.Blood_flow.v(2,:)=the_decimate*Tfactor*Xfactor*(cotd(thetas(2,:)-90))/2;
else
    MScanData.Blood_flow.v(1,:)=the_decimate*Tfactor*Xfactor*(cotd(thetas(1,:)-90));
    MScanData.Blood_flow.v(2,:)=the_decimate*Tfactor*Xfactor*(cotd(thetas(1,:)-90));
end
MScanData.Blood_flow.velocity2(1,:)=VelocityCleanUp_2P(MScanData.Blood_flow.v(1,:),max_v);
if isinterleaved==1
    MScanData.Blood_flow.velocity2(2,:)=VelocityCleanUp_2P(MScanData.Blood_flow.v(2,:),max_v);
else
    MScanData.Blood_flow.velocity2(2,:)=-(MScanData.Blood_flow.velocity2(1,:));
end
v_hold=MScanData.Blood_flow.velocity2;
v_out=0.5*(v_hold(1,:)-v_hold(2,:));%take the average of the two directions
MScanData.Blood_flow.v_out=v_out;
MScanData.Blood_flow.the_t=the_t/(MScanData.notes.LineRate);
MScanData.Blood_flow.the_Fs=4/window_time_size;
disp('Number of frames: '); 
disp(MScanData.notes.num_frames);

the_date=date;
[fixed_v]=VelocityCleanUpSTDOutliers_2P(MScanData.Blood_flow.v_out, 3);% clear outliers > 3standard deviations
MScanData.Blood_flow.fixed_v=fixed_v;
params.Fs=1/(MScanData.Blood_flow.the_t(2)-MScanData.Blood_flow.the_t(1));
params.tapers=[3 5];%to deal with Chronux 2.0
[Shr,thr,fhr]=mtdspecgramc(diff(MScanData.Blood_flow.fixed_v),[2 .5],[pi/2],params);

MScanData.Blood_flow.velocity_specgram=Shr;
MScanData.Blood_flow.velocity_specgram_t=thr;
try %to load raw data
    MScanData.Blood_flow.Image = data_cache;
catch
    MScanData.Blood_flow.Image = zeros(MScanData.notes.num_frames*(MScanData.notes.Frame_Height),1);
end



separability=max(squeeze(mean(MScanData.Blood_flow.spread_matrix_adaptive,1))')./mean(squeeze(mean(MScanData.Blood_flow.baseline_var,1))');%
MScanData.Blood_flow.separability=separability; % separability should be >3
sep_threshold=3;
low_sep=find(separability<sep_threshold);
fileID = 'test';
params_ps.Fs=MScanData.Blood_flow.the_Fs;
params_ps.tapers=[3 5];
[MScanData.Blood_flow.powerSpec_s, MScanData.Blood_flow.powerSpec_f] = mtspectrumc(detrend(MScanData.Blood_flow.fixed_v),params_ps);
ps_plot = figure;
plot(MScanData.Blood_flow.powerSpec_f,MScanData.Blood_flow.powerSpec_s)
title(strcat(fileID,' Power Spectrum'));
xlabel('Frequency (Hz)');
ylabel('Power');
saveas(ps_plot,strcat(fileID,'_PowerSpectrum_',date));

psl_plot = figure;
semilogy(MScanData.Blood_flow.powerSpec_f,MScanData.Blood_flow.powerSpec_s);
title(strcat(fileID,' Power Spectrum (semilog)'));
xlabel('Frequency (Hz)');
ylabel('Power');
saveas(psl_plot,strcat(fileID,'_PowerSpectrumSemilog_',date));

v_plot = figure;
plot(MScanData.Blood_flow.the_t,MScanData.Blood_flow.fixed_v);
title(strcat(fileID,' Fixed Velocity'));
xlabel('Time (s)');
ylabel('Velocity (um/s)');
saveas(v_plot,strcat(fileID,'_FixedVelocity_',date));

theta_plot = figure; hold on
plot(MScanData.Blood_flow.the_t, MScanData.Blood_flow.thetas(1,:));
plot(MScanData.Blood_flow.the_t, MScanData.Blood_flow.thetas(2,:));
legend('Theta1','Theta2');
title(strcat(fileID,' Thetas'));
xlabel('Time (s)');
ylabel('Theta (degrees)');
saveas(theta_plot,strcat(fileID,'_Thetas_',date));

v2_plot = figure; hold on
plot(MScanData.Blood_flow.the_t, MScanData.Blood_flow.velocity2(1,:));
plot(MScanData.Blood_flow.the_t, MScanData.Blood_flow.velocity2(2,:));
legend('Forward','Backward');
title(strcat(fileID,' Velocities'));
xlabel('Time (s)');
ylabel('Velocity (um/s)');
saveas(v2_plot,strcat(fileID,'_Velocities_',date));

sep_plot =figure;
plot(MScanData.Blood_flow.the_t, MScanData.Blood_flow.separability);
title(strcat(fileID, ' Separability'));
xlabel('Time (s)');
ylabel('Separability');
saveas(sep_plot,strcat(fileID,'_Separability_',date));
    
end