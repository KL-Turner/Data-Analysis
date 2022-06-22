clear; clc; close all;
colorRfcAwake = [(0/256),(64/256),(64/256)];
colorRfcNREM = [(0/256),(174/256),(239/256)];
colorRfcREM = [(190/256),(30/256),(45/256)];
colorOrange = [(233/256),(105/256),(44/256)];
colorPink = [(256/256),(28/256),(207/256)];
colorSapphire = [(15/256),(82/256),(187/256)];
colorDarkRed = [(164/256),(0/256),(0/256)];
exampleProcDataFileID = uigetfile('*_ProcData.mat','MultiSelect','Off');
load(exampleProcDataFileID)
[animalID,fileDate,fileID] = GetFileInfo_Test(exampleProcDataFileID);
specDataFileID = [animalID '_' fileID '_SpecDataA.mat'];
load(specDataFileID);
baselineFileID = [animalID '_RestingBaselines.mat'];
load(baselineFileID)
strDay = ConvertDate_Test(fileDate);
dsFs = ProcData.notes.dsFs;
analogFs = ProcData.notes.analogSamplingRate;
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z,p,k] = butter(3,[300,3000]/(analogFs/2));
[sos,g] = zp2sos(z,p,k);
[z1,p1,k1] = butter(4,10/(dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(4,0.5/(dsFs/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
% whisker angle
filtWhiskerAngle = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
% force sensor
filtForceSensor = filtfilt(sos1,g1,abs(ProcData.data.forceSensor));
% EMG
normEMG = ProcData.data.EMG.emg - RestingBaselines.manualSelection.EMG.emg.(strDay);
filtEMG = filtfilt(sos1,g1,normEMG);
% heart rate
heartRate = ProcData.data.heartRate;
% HbT data
LH_HbT = ProcData.data.CBV_HbT.adjLH;
filtLH_HbT = filtfilt(sos2,g2,LH_HbT);
RH_HbT = ProcData.data.CBV_HbT.adjRH;
filtRH_HbT = filtfilt(sos2,g2,RH_HbT);
% cortical and hippocampal spectrograms
cortical_LHnormS = SpecData.cortical_LH.normS.*100;
cortical_RHnormS = SpecData.cortical_RH.normS.*100;
hippocampusNormS = SpecData.hippocampus.normS.*100;
T = SpecData.cortical_LH.T;
F = SpecData.cortical_LH.F;
%% IOS sleep exampl
summaryFigure = figure;
sgtitle('Figure 1 - Turner et al. 2020')
% EMG and force sensor
ax1 = subplot(7,1,1);
p1 = plot((1:length(filtEMG))/dsFs,filtEMG,'k','LineWidth',0.5);
ylabel({'EMG','power (a.u.)'})
ylim([-2,2.5])
yyaxis right
p2 = plot((1:length(filtForceSensor))/dsFs,filtForceSensor,'color',colorPink,'LineWidth',0.5);
ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p1,p2],'EMG','Pressure')
set(gca,'Xticklabel',[])
set(gca,'box','off')
ylim([-0.1,2.5])
ax1.TickLength = [0.01,0.01];
ax1.YAxis(1).Color = 'k';
ax1.YAxis(2).Color = colorPink;
% whisker angle and heart rate
ax2 = subplot(7,1,2);
p3 = plot((1:length(filtWhiskerAngle))/dsFs,-filtWhiskerAngle,'k','LineWidth',0.5);
ylabel({'Whisker','angle (deg)'})
xlim([0,600])
ylim([-10,50])
yyaxis right
p4 = plot((1:length(heartRate)),heartRate,'color',colorOrange,'LineWidth',0.5);
ylabel({'Heart rate','Freq (Hz)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p3,p4],'Whisker angle','Heart rate')
set(gca,'Xticklabel',[])
set(gca,'box','off')
ylim([5,10])
ax2.TickLength = [0.01,0.01];
ax2.YAxis(1).Color = 'k';
ax2.YAxis(2).Color = colorOrange;
% HbT and behavioral indeces
ax34 =subplot(7,1,[3,4]);
p6 = plot((1:length(filtRH_HbT))/dsFs,filtRH_HbT,'color',colorSapphire,'LineWidth',1);
hold on
p5 = plot((1:length(filtLH_HbT))/dsFs,filtLH_HbT,'color',colorDarkRed,'LineWidth',1);
x1 = xline(0,'color',colorRfcNREM,'LineWidth',2);
x2 = xline(105,'color',colorRfcREM,'LineWidth',2);
x3 = xline(285,'color',colorRfcAwake,'LineWidth',2);
ylabel('\Delta[HbT] (\muM)')
legend([p5,p6,x3,x1,x2],'Left hem','Right hem','Awake','NREM','REM')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
ylim([-35,135])
ax34.TickLength = [0.01,0.01];
% left cortical electrode spectrogram
ax5 = subplot(7,1,5);
semilog_imagesc(T,F,cortical_LHnormS,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'LH cort LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
ax5.TickLength = [0.01,0.01];
% right cortical electrode spectrogram
ax6 = subplot(7,1,6);
semilog_imagesc(T,F,cortical_RHnormS,'y')
axis xy
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'RH cort LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
ax6.TickLength = [0.01,0.01];
% hippocampal electrode spectrogram
ax7 = subplot(7,1,7);
semilog_imagesc(T,F,hippocampusNormS,'y')
axis xy
c7 = colorbar;
ylabel(c7,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (min)')
ylabel({'Hipp LFP','Freq (Hz)'})
set(gca,'box','off')
ax7.TickLength = [0.01,0.01];
% axes properties
ax1Pos = get(ax1,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
ax5Pos(3:4) = ax1Pos(3:4);
ax6Pos(3:4) = ax1Pos(3:4);
ax7Pos(3:4) = ax1Pos(3:4);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
set(ax7,'position',ax7Pos);
%% heart rate
[LH_Sr,tr,fr,LH_HR] = FindHeartRate_Test(ProcData.data.CBV.LH,ProcData.notes.CBVCamSamplingRate);
[RH_Sr,~,~,RH_HR] = FindHeartRate_Test(ProcData.data.CBV.RH,ProcData.notes.CBVCamSamplingRate);
% Average the two signals from the left and right windows
Sr = (LH_Sr + RH_Sr)/2;
HR = (LH_HR + RH_HR)/2;
% Smooth the signal with a 5 Hz low pass third-order butterworth filter
[B,A] = butter(3,5/(ProcData.notes.CBVCamSamplingRate/2),'low');
filtHeartRate = filtfilt(B,A,HR);   % Filtered heart rate signal
figure;
imagesc(tr,fr,Sr)
axis xy
caxis([0,5])
hold on;
plot(filtHeartRate,'color',colorOrange,'LineWidth',2)
title('Spectral estimation of heart rate')
xlabel('Time(sec')
ylabel('Freq (Hz)')
set(gca,'box','off')
%% HR calc
function [Sr,tr,fr,HR] = FindHeartRate_Test(r,Fr)
% mean subtract to remove slow drift
r = r - mean(r);
 % [time band width, number of tapers]
tapers_r = [2,3];
movingwin_r = [3.33,1];
% Frame rate
params_r.Fs = Fr;
params_r.fpass = [5,15];
params_r.tapers = tapers_r;
[Sr,tr,fr] = mtspecgramc(r,movingwin_r,params_r);
% Sr: spectrum; tr: time; fr: frequency
% largest elements along the frequency direction
[~,ridx] = max(Sr,[],2);
HR = fr(ridx);   % heart rate, in Hz

end
%% semilog axis
function [] = semilog_imagesc(x,y,C,logaxis)
surface(x,y,zeros(size(C)),(C),'LineStyle','none');%
q = gca;
q.Layer = 'top';
if strcmp(logaxis,'y') == 1
    set(gca,'YScale','log');
elseif strcmp(logaxis,'x') == 1
    set(gca,'XScale','log');
elseif strcmp(logaxis,'xy') == 1
    set(gca,'XScale','log');
    set(gca,'YScale','log');
end
axis xy
axis tight
end
%% animal info from file name
function [animalID,fileDate,fileID] = GetFileInfo_Test(fileName)
% Identify the extension
extInd = strfind(fileName(1,:),'.');
extension = fileName(1,extInd + 1:end);
% Identify the underscores
fileBreaks = strfind(fileName(1,:),'_');
switch extension
    case 'bin'
        animalID = [];
        fileDate = fileName(:,1:fileBreaks(1) - 1);
        fileID = fileName(:,1:fileBreaks(4) - 1);
    case 'mat'
        % Use the known format to parse
        animalID = fileName(:,1:fileBreaks(1) - 1);
        if numel(fileBreaks) > 3
            fileDate = fileName(:,fileBreaks(1) + 1:fileBreaks(2) - 1);
            fileID = fileName(:,fileBreaks(1) + 1:fileBreaks(5) - 1);
        else
            fileDate = [];
            fileID = [];
        end
end

end
%% numeric date to string
function [days] = ConvertDate_Test(dateTag)
days = cell(size(dateTag,1),1);
for f = 1:size(dateTag,1)
    days{f} = datestr([2000 + str2double(dateTag(f,1:2)),str2double(dateTag(f,3:4)),str2double(dateTag(f,5:6)),00,00,00],'mmmdd');
end
if length(days) == 1
    days = days{1};
end

end
