function [mv_mpP,xcfig]=linescan_xcov_velocity_TIFF_2P(mv_mpP)
%this function uses the cross correlation method to measure velocity
%using the method of Kim et al, PLoS ONE 2012
%Patrick Drew 7/2012
%
[~,xc_image,xshift]=linescan_xcov_2P(mv_mpP);
[b,a]=butter(2,200/(2*mv_mpP.Tfactor));%200Hz cutoff
xc_velocity=filtfilt(b,a,xshift*mv_mpP.Xfactor*mv_mpP.Tfactor);% calculate the velocity from the displacement and filter
mv_mpP.Blood_flow.xc_velocity=xc_velocity;%velocity obtained with the cross correlation method
params.Fs=mv_mpP.Tfactor;
params.tapers=[20 39];
[S_xc,f_xc]=mtspectrumc(xshift-mean(xshift(:)),params);
mv_mpP.Blood_flow.S_xcor=S_xc;
mv_mpP.Blood_flow.f_xcor=f_xc;
mv_mpP.Blood_flow.params_xc=params;
mv_mpP.Blood_flow.xc_image=xc_image;

xcfig=figure(33);
subplot(1,4,1:3)
hold off
imagesc((1:length(mv_mpP.Blood_flow.xc_image))/mv_mpP.Tfactor,(size(mv_mpP.Blood_flow.xc_image,1)/2:-1:-size(mv_mpP.Blood_flow.xc_image,1)/2)*mv_mpP.Tfactor*mv_mpP.Xfactor/1000,mv_mpP.Blood_flow.xc_image)
hold on
plot((1:length(mv_mpP.Blood_flow.xc_velocity))/mv_mpP.Tfactor,mv_mpP.Blood_flow.xc_velocity/1000,'w')
xlabel('time, seconds')
ylabel('velocity, mm/sec')
title([mv_mpP.Header.ImageID ' ' mv_mpP.Header.animal ' ' mv_mpP.Header.VesselNO]);
axis xy
subplot(1,4,4)
loglog(f_xc,S_xc)
saveas(gcf,fullfile(['Velocity_C_' mv_mpP.Header.ImageID]),'fig');

figure
time=(1:length(mv_mpP.Blood_flow.velocity2(1,:)));
time=time*str2num(mv_mpP.Header.Header.Frame_Rate);
plot(time,mv_mpP.Blood_flow.velocity2(1,:))
xlabel('time, seconds')
ylabel('velocity, mm/sec')
title([mv_mpP.Header.ImageID ' ' mv_mpP.Header.animal ' ' mv_mpP.Header.VesselNO]);
saveas(gcf,fullfile(['Velocity1_C_' mv_mpP.Header.ImageID]),'fig');
end