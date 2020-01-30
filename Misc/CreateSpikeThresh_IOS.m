function [thresh] = CreateSpikeThresh(RawData,Thresholds,strday)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs: 
%
%   Last Revised: February 29th, 2019
%________________________________________________________________________________________________________________________

stdev = Thresholds.(['MUA_StDev_' strday]);
[b,a] = butter(2,300/(RawData.Fs.Analog/2),'high');
TestMU = filtfilt(b,a,RawData.Data.Neuro);
figure;
subplot(411); plot(RawData.Data.Sol); axis tight; 
ylabel('Solenoids');
subplot(412); plot(RawData.Data.WhiskerAngle); axis tight; 
ylabel('Whisker Angle');
subplot(413); plot((TestMU-mean(TestMU))/stdev); axis tight;
ylabel('Multi-Unit')
pass = 0;
while pass == 0;
    thresh = input(['No spike threshold found on ' strday ' please enter an acceptable threshold. ']);
    [TestSR,~] = MUspikeRATE_gauss(TestMU,thresh,RawData.Fs.Analog,stdev);
    subplot(414); plot(TestSR); axis tight;
    ok3 = input('Threshold ok? (y/n) ','s');
    if strcmp(ok3,'y')
        pass = 1;
    end
end
display('----');