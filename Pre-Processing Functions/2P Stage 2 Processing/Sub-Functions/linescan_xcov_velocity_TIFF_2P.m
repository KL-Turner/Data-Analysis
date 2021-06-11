function [MScanData,xcfig]=linescan_xcov_velocity_TIFF_2P(MScanData)
%this function uses the cross correlation method to measure velocity
%using the method of Kim et al, PLoS ONE 2012
%Patrick Drew 7/2012
%
[~,xc_image,xshift] = linescan_xcov_2P(MScanData);
[b,a]=butter(2,200/(2*MScanData.notes.Tfactor));%200Hz cutoff
xc_velocity=filtfilt(b,a,xshift*MScanData.notes.Xfactor*MScanData.notes.Tfactor);% calculate the velocity from the displacement and filter
MScanData.Blood_flow.xc_velocity=xc_velocity;%velocity obtained with the cross correlation method
params.Fs=MScanData.notes.Tfactor;
params.tapers=[20 39];
[S_xc,f_xc]=mtspectrumc(xshift-mean(xshift(:)),params);
MScanData.Blood_flow.S_xcor=S_xc;
MScanData.Blood_flow.f_xcor=f_xc;
MScanData.Blood_flow.params_xc=params;
MScanData.Blood_flow.xc_image=xc_image;

xcfig=figure(33);
subplot(1,4,1:3)
hold off
imagesc((1:length(MScanData.Blood_flow.xc_image))/MScanData.notes.Tfactor,(size(MScanData.Blood_flow.xc_image,1)/2:-1:-size(MScanData.Blood_flow.xc_image,1)/2)*MScanData.notes.Tfactor*MScanData.notes.Xfactor/1000,MScanData.Blood_flow.xc_image)
hold on
plot((1:length(MScanData.Blood_flow.xc_velocity))/MScanData.notes.Tfactor,MScanData.Blood_flow.xc_velocity/1000,'w')
xlabel('time, seconds')
ylabel('velocity, mm/sec')
title([MScanData.notes.imageID ' ' MScanData.notes.animalID ' ' MScanData.notes.vesselID]);
axis xy
subplot(1,4,4)
loglog(f_xc,S_xc)
saveas(gcf,fullfile(['Velocity_C_' MScanData.notes.animalID MScanData.notes.vesselID]),'fig');
figure
time = (1:length(MScanData.Blood_flow.velocity2(1,:)));
time = time*str2num(MScanData.notes.Frame_Rate);
plot(time,MScanData.Blood_flow.velocity2(1,:))
xlabel('time, seconds')
ylabel('velocity, mm/sec')
title([MScanData.notes.animalID ' ' MScanData.notes.vesselID ' XC-velocity']);
saveas(gcf,fullfile(['Velocity1_C_' MScanData.notes.imageID]),'fig');

end