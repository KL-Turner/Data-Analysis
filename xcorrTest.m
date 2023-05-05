clear; clc; close all
load('W37_RH_230302_11_14_14_003_A02_MergedDataBackup.mat')

samplingRate = 10;
[z1,p1,k1] = butter(4,1/(samplingRate/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);

[z2,p2,k2] = butter(4,1/(30/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);

diameter = sqrt(MergedData.data.pupil.data.patchedPupilAreaA_trim./pi)*2;
mmPerPixel = 0.018;
pupilDiameter = diameter*mmPerPixel;
pupilDiameter2 = filtfilt(sos2,g2,pupilDiameter);

dsForceSensor = resample(MergedData.data.forceSensorL,10,30);
binForceSensor = double([BinarizeForceSensor_IOS(dsForceSensor,0.2000),0]);

dsWhiskerAngle = resample(MergedData.data.whiskerAngle,10,30);
binWhiskerAngle = double([BinarizeWhiskers_IOS(dsWhiskerAngle,10,0,250),0,0]);

whiskerAngle = detrend(binWhiskerAngle,'constant');
forceSensor = detrend(binForceSensor,'constant');
grabNE = detrend(filtfilt(sos1,g1,MergedData.data.GRABNE),'constant');
grabNE2 = filtfilt(sos1,g1,MergedData.data.GRABNE);
pupilArea = detrend(filtfilt(sos1,g1,resample(pupilDiameter,10,30)),'constant');

figure
s1 = subplot(3,1,1);
plot((1:length(dsForceSensor))/samplingRate,dsForceSensor,'k')
hold on
xline(3396/samplingRate,'m')
xline(3845/samplingRate,'m')
xline(4884/samplingRate,'m')
xline(5239/samplingRate,'m')
xline(5961/samplingRate,'m')
ylabel('Force (V)')
xlabel('Time (sec)')
yyaxis right
plot((1:length(dsWhiskerAngle))/samplingRate,dsWhiskerAngle,'r')
ylabel('Whisker Angle (deg)')
s2 = subplot(3,1,2);
plot((1:length(grabNE))/samplingRate,((grabNE2 - mean(grabNE2))/mean(grabNE2))*100,'k')
hold on
xline(3396/samplingRate,'m')
xline(3845/samplingRate,'m')
xline(4884/samplingRate,'m')
xline(5239/samplingRate,'m')
xline(5961/samplingRate,'m')
ylabel('\DeltaF/F0 (%)')
xlabel('Time (sec)')
s3 = subplot(3,1,3);
plot((1:length(pupilDiameter2))/(samplingRate*3),((pupilDiameter2 - mean(pupilDiameter2))/mean(pupilDiameter2))*100,'k')
hold on
xline(3396/samplingRate,'m')
xline(3845/samplingRate,'m')
xline(4884/samplingRate,'m')
xline(5239/samplingRate,'m')
xline(5961/samplingRate,'m')
ylabel('\DeltaD/D (%)')
xlabel('Time (sec)')
linkaxes([s1,s2,s3],'x')

lagSec = 20;
maxLag = samplingRate*lagSec;
[r1,lags1] = xcorr(pupilArea,forceSensor,maxLag,'coeff');
[r2,lags2] = xcorr(pupilArea,grabNE,maxLag,'coeff');
[r3,lags3] = xcorr(grabNE,forceSensor,maxLag,'coeff');
[r4,lags4] = xcorr(whiskerAngle,forceSensor,maxLag,'coeff');
[r5,lags5] = xcorr(pupilArea,whiskerAngle,maxLag,'coeff');

figure;
plot(lags1/samplingRate,r1,'k')
hold on
[m1,idx] = max(r1);
xline(lags1(idx)/samplingRate,'r');
title('Pupil w/ lagged copies of binForce')
xlabel('time (sec)')
ylabel('corr coef')

figure;
plot(lags2/samplingRate,r2,'k')
hold on
[m2,idx] = max(r2);
xline(lags2(idx)/samplingRate,'r');
title('Pupil w/ lagged copies of NE')

figure;
plot(lags3/samplingRate,r3,'k')
hold on
[m3,idx] = max(r3);
xline(lags3(idx)/samplingRate,'r');
title('NE w/ lagged copies of binForce')
xlabel('time (sec)')
ylabel('corr coef')

figure;
plot(lags4/samplingRate,r4,'k')
hold on
[m4,idx] = max(r4);
xline(lags4(idx)/samplingRate,'r');
title('binWhiskers w/ lagged copies of binForce')
xlabel('time (sec)')
ylabel('corr coef')

figure;
plot(lags5/samplingRate,r5,'k')
hold on
[m5,idx] = max(r5);
xline(lags5(idx)/samplingRate,'r');
title('pupil area w/ lagged copies of binWhiskers')
xlabel('time (sec)')
ylabel('corr coef')
