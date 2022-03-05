function AnalyzeFiberPhotometryData_FP_KG(csvFile)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Kyle W. Gheres and Md Shakhawat Hossain
%________________________________________________________________________________________________________________________
%% READ ME
%
% This script removes the first 0.25s of each trial to eliminate a chirp
% found at the start of each trial.
%
% This script does not remove behavior from data prior to correcting for
% slow metabolic signal decay
%
% This script lacks the correction of [HbT] signal attenuation which will
% need to be added later.
%
%This script does not contain any behavior metrics.
%
% This script loads in all .csv files from a single imaging day, removes
% any slow decrease in signal intensity from each trial, rescales all of
% the data on the same interval [0,1] and then z-score normalizes the data
%_________________________________________________________________________________________________________________
% csvFiles= ls('*.csv'); %this is temporary KWG

%% notes
correctionFlag = 'y'; % If GCaMP signal need to be corrected
opticalChannelNames = {'RH_control','RH_GCaMP7s','RH_bloodVolume','LH_control','LH_GCaMP7s','LH_bloodVolume','syncChannel'};
% identify the extension
extInd = strfind(csvFile(1,:),'.');
% identify the underscores
fileBreaks = strfind(csvFile(1,:),'_');
% file parameters
FiberData.notes.animalID = csvFile(1:fileBreaks(1) - 1);
FiberData.notes.date = csvFile(fileBreaks(1) + 1:fileBreaks(2) - 1);
FiberData.notes.sessionNumber = csvFile(fileBreaks(2) + 1:extInd - 1);
%% read.CSV
csvData = csvread(csvFile,2,0);
FiberData.notes.sessionDuration_sec = round(csvData(end,1));
% 2 = RH 405nm, 3 = RH 465nm, 5 = RH 560nm, 7 = LH 405nm, 8 = LH 465nm, 10 = RL 560nm
dataChannels = [2,3,5,7,8,10];
channelData = csvData(:,dataChannels);
FiberData.syncData = csvData(:,11);
%% parameters
FiberData.notes.channels = opticalChannelNames;
FiberData.notes.decimation = 10;
FiberData.notes.fitFreq = 0.01;
FiberData.notes.lowFreq = 1;
FiberData.notes.finalFreq = [0.01,1];
FiberData.notes.trialLength = 15.50; % minutes
FiberData.notes.idleTime = 4.50; % minutes
%% determine sampling rate
timeArray = csvData(:,1);
FiberData.notes.samplingRate = round(length(timeArray)/timeArray(end));
%% Correct [HbT] dependent GCaMP signal attenuation
if strcmpi(correctionFlag,'y')
    % EXAMPLE CODE ONLY REPLACE WITH TRUE CORRECTION
    % rescaledData(:,2)=rescaledData(:,2)-correctConst*rescaledData(:,3);
    % This removes the [HbT] dependent attenuation from the recorded GCaMP
    % rescaledData(:,5)=rescaledData(:,5)-correctConst*rescaledData(:,6);
    % This removes the [HbT] dependent attenuation from the recorded GCaMP
end
%% filter characteristics
[z1,p1,k1] = butter(3,FiberData.notes.fitFreq/(0.5*FiberData.notes.samplingRate),'low'); % design lowpass filter for hemodynamic correction
[sosFit,gFit] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(3,FiberData.notes.lowFreq/(0.5*FiberData.notes.samplingRate),'low'); % low pass for optical data to physiologically relevant range
[sosLow,gLow] = zp2sos(z2,p2,k2);
%% rescale data
expectedLength = FiberData.notes.samplingRate*60*60*5;
theSignals = channelData;
firstPoint = repmat(theSignals(1,:),length(theSignals),1);
theSignals = theSignals - firstPoint; % This sets the signal data to have a first value of 0 which eliminates any filter initiation transients KWG
theSignals = theSignals(1:expectedLength,:);
syncData = csvData(:,11);
FiberData.syncData = syncData(1:expectedLength);
fitSignal = filtfilt(sosFit,gFit,theSignals); % low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends KWG
[lowPassData,FiberData] = CorrectFiberDecay_FP(fitSignal,theSignals,sosLow,gLow,FiberData); % Fit and remove decay from imaging trial, then low pass filter data KWG
for q = 1:size(lowPassData,2)
    rescaledData(:,q) = rescale(lowPassData(:,q),0,1); % rescale all data between 0 to 1 ONLY USE THIS DATA IF YOUR ESTIMATION OF [HbT] ATTENUATION WAS PERFORMED ON RESCALED DATA KWG
end
%% Z-score data
avgData = nanmean(rescaledData,1);
stdData = std(rescaledData,0,1,'omitnan');
avgMatrix = repmat(avgData,length(rescaledData),1);
stdMatrix = repmat(stdData,length(rescaledData),1);
zScoredFiberData = (rescaledData - avgMatrix)./stdMatrix; %This is each channel z scored to its standard deviation. KWG
%% split the data into trials
fields = {'RH_405','RH_465','RH_560','LH_405','LH_465','LH_560'};
for bb = 1:length(fields)
    FiberData.OpticalData.(fields{1,bb}).raw = channelData(:,bb);
    FiberData.OpticalData.(fields{1,bb}).zScored = zScoredFiberData(:,bb);
    FiberData.OpticalData.(fields{1,bb}).rescaled = rescaledData(:,bb);
end
save([FiberData.notes.animalID '_' FiberData.notes.date '_FiberData.mat'],'FiberData','-v7.3');

end
