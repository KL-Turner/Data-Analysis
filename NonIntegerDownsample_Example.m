lineScanFs = MScanData.data.bloodFlow.Fs; % 12.2070 Hz at 512x512 linerate of 1562.5
desiredFs = 10;
trialDuration = round(MScanData.notes.numberOfFrames/MScanData.notes.frameRate,1);

integerSineWave = dsp.SineWave(1,lineScanFs,1,'SampleRate',desiredFs,'SamplesPerFrame',length(MScanData.data.bloodFlow.fixedVelocity));
x = 0:1/desiredFs:trialDuration;
% check figure
figure
p1 = plot((1:length(MScanData.data.bloodFlow.fixedVelocity))/lineScanFs,MScanData.data.bloodFlow.fixedVelocity,'r');
hold on;
p2 = plot(x/y,'b');
xlabel('Time (s)')
ylabel('Velocity (mm/sec)')
legend([p1,p2],'Original signal','Downsampled')


y = integerSineWave();
plot(y)