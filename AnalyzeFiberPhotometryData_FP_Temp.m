function AnalyzeFiberPhotometryData_FP_Temp(csvFiles)
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
clearvars -except csvFiles
for qq = 1:size(csvFiles,1)

    fileName = csvFiles(qq,:);
    correctionFlag = 'y'; % If GCaMP signal need to be corrected
    opticalChannelNames = {'RH_control','RH_GCaMP7s','RH_bloodVolume','LH_control','LH_GCaMP7s','LH_bloodVolume','syncChannel'};
    % identify the extension
    extInd = strfind(fileName(1,:),'.');
    % identify the underscores
    fileBreaks = strfind(fileName(1,:),'_');
    % file parameters
    FiberData(qq).notes.animalID = fileName(1:fileBreaks(1) - 1);
    FiberData(qq).notes.date = fileName(fileBreaks(1) + 1:fileBreaks(2) - 1);
    FiberData(qq).notes.sessionNumber = fileName(fileBreaks(2) + 1:extInd - 1);
    
    %% read.CSV
    csvData = csvread(fileName,2,0);
    FiberData(qq).notes.sessionDuration_sec = round(csvData(end,1));
    % 2 = RH 405nm, 3 = RH 465nm, 5 = RH 560nm, 7 = LH 405nm, 8 = LH 465nm, 10 = RL 560nm
    dataChannels = [2,3,5,7,8,10];
    channelData = csvData(:,dataChannels);
    syncData = csvData(:,11);
    
    %% parameters
    FiberData(qq).notes.channels = opticalChannelNames;
    FiberData(qq).notes.decimation = 10;
    FiberData(qq).notes.fitFreq = 0.01;
    FiberData(qq).notes.lowFreq = 1;
    FiberData(qq).notes.finalFreq = [0.01,1];
    FiberData(qq).notes.trialLength = 15.50; % minutes
    FiberData(qq).notes.idleTime = 4.50; % minutes
    
    %% determine sampling rate
    timeArray = csvData(:,1);
    dTimeArray = diff(timeArray);
    fileBreaks = dTimeArray > 1;
    fileBreaksIndex = find(fileBreaks == 1);
    FiberData(qq).notes.samplingRate = round(fileBreaksIndex(1)/timeArray(fileBreaksIndex(1))); %Redundant to line 36? KWG
    
    %% Create spacer for between trials
    time_between_trials=4.5; %time in minutes between trials KWG
    time_between_sessions=15; % time in minutes between sessions this should be replaced by actual recorded intertrial times at some point KWG
    num_samples=FiberData(qq).notes.samplingRate*60*time_between_trials;
    theSpacer((1:num_samples),1)=NaN; % this array of NaN will separate trials for estimation of fluorescence signal decay KWG
    spaceMat=repmat(theSpacer,1,size(channelData,2)); % Matrix that will be concatenated between trials KWG
    
    session_samples=FiberData(qq).notes.samplingRate*60*time_between_sessions;
    sessionSpacer((1:session_samples),1)=NaN;
    sessionMat=repmat(sessionSpacer,1,size(channelData,2));
    
    %% filter characteristics
    [z1,p1,k1] = butter(3,FiberData(qq).notes.fitFreq/(0.5*FiberData(qq).notes.samplingRate),'low'); % design lowpass filter for hemodynamic correction
    [sosFit,gFit] = zp2sos(z1,p1,k1);
    [z2,p2,k2] = butter(3,FiberData(qq).notes.lowFreq/(0.5*FiberData(qq).notes.samplingRate),'low'); % low pass for optical data to physiologically relevant range
    [sosLow,gLow] = zp2sos(z2,p2,k2);
    [z3,p3,k3] = butter(3,FiberData(qq).notes.finalFreq/(0.5*FiberData(qq).notes.samplingRate),'bandpass'); % low pass filter for locomotion data
    [sosFinal,gFinal] = zp2sos(z3,p3,k3);
    
    %% Find trial start/end indicies
    timeJump=find(diff(timeArray)>1); %Find indicies of imaging trial end
    offsetFrame=round(0.5*FiberData(qq).notes.samplingRate,0); %skip first 0.5 sec from each trial
%     trialDiff=diff(csvData(:,3));
%     stopInds=find(trialDiff<-0.1)-1;
%     startInds=find(trialDiff>0.1)+1;
    stopInds=[(timeJump-offsetFrame);(length(channelData)-offsetFrame)];
    startInds=[offsetFrame;(timeJump+offsetFrame)];
    
%     diffStop=find(diff(stopInds)<10); % Find any indicies that were triggered by signal bounce
%     diffStart=find(diff(startInds)<10);% Find any indicies that were triggered by signal bounce
%     
%     stopInds(diffStop)=[];% Remove any indicies that were triggered by signal bounce
%     startInds(diffStart)=[];% Remove any indicies that were triggered by signal bounce
    
    spacedData=[];
    spacedSync=[];
    for nn=1:length(stopInds)
        if nn==1
            syncSignals=syncData(offsetFrame:stopInds(nn));
            theSignals=channelData((offsetFrame:stopInds(nn)),:); % Keep same length as following trials
            firstPoint=repmat(theSignals(1,:),length(theSignals),1);
            theSignals=theSignals-firstPoint; %This sets the signal data to have a first value of 0 which eliminates any filter initiation transients KWG
            fitSignal=filtfilt(sosFit,gFit,theSignals);    % low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends KWG
            [lowPassData,FiberData]=correctDecay(fitSignal,theSignals,sosLow,gLow,FiberData, nn,qq); %Fit and remove decay from imaging trial, then low pass filter data KWG
            spacedData=[lowPassData;spaceMat]; %add intertrial spacing of NAN values for proper temporal alignment KWG
            spacedSync=[syncSignals;spaceMat(:,1)];
            trialwiseData{qq,nn}=lowPassData;
        elseif nn==length(stopInds)
            syncSignals=syncData(startInds(nn):stopInds(nn));
            theSignals=channelData((startInds(nn):stopInds(nn)),:);
            firstPoint=repmat(theSignals(1,:),length(theSignals),1);
            theSignals=theSignals-firstPoint; %This sets the signal data to have a first value of 0 which eliminates any filter initiation transients KWG
            fitSignal=filtfilt(sosFit,gFit,theSignals);    % low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends KWG
            [lowPassData,FiberData]=correctDecay(fitSignal,theSignals,sosLow,gLow,FiberData, nn,qq); %Fit and remove decay from imaging trial, then low pass filter data KWG
            spacedData=[spacedData;lowPassData]; %add intertrial spacing of NAN values for proper temporal alignment KWG
            spacedSync=[spacedSync;syncSignals];
            trialwiseData{qq,nn}=lowPassData;
        else
            syncSignals=syncData(startInds(nn):stopInds(nn));
            theSignals=channelData((startInds(nn):stopInds(nn)),:);
            firstPoint=repmat(theSignals(1,:),length(theSignals),1);
            theSignals=theSignals-firstPoint; %This sets the signal data to have a first value of 0 which eliminates any filter initiation transients KWG
            fitSignal=filtfilt(sosFit,gFit,theSignals);    % low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends KWG
            [lowPassData,FiberData]=correctDecay(fitSignal,theSignals,sosLow,gLow,FiberData, nn,qq); %Fit and remove decay from imaging trial, then low pass filter data KWG
            spacedData=[spacedData;lowPassData;spaceMat]; %This matrix contains all data with correct temporal separation of trials with between trial intervals filled with NaN KWG
            spacedSync=[spacedSync;syncSignals;spaceMat(:,1)];
            trialwiseData{qq,nn}=lowPassData;
        end
    end
    if qq~=size(csvFiles,1)
    spacedData=[spacedData((1:(length(spacedData)-1)),:);sessionMat]; %The last sample has a click get rid of it KWG
    spacedSync=[spacedSync((1:(length(spacedSync)-1)),:);sessionMat(:,1)];
    end
    sessionSignals{qq}=spacedData;
    sessionSync{qq}=spacedSync;
end
catData=[];
catSync=[];
for cc=1:size(sessionSignals,2)
    catData=[catData;sessionSignals{cc}]; %concatenate all of the imaging session data together KWG
    catSync=[catSync;sessionSync{cc}];
end
for q = 1:size(catData,2)
    rescaledData(:,q) = rescale(catData(:,q),0,1); % rescale all data between 0 to 1 ONLY USE THIS DATA IF YOUR ESTIMATION OF [HbT] ATTENUATION WAS PERFORMED ON RESCALED DATA KWG
end

%% Correct [HbT] dependent GCaMP signal attenuation
if strcmpi(correctionFlag,'y')
    % EXAMPLE CODE ONLY REPLACE WITH TRUE CORRECTION
    %rescaledData(:,2)=rescaledData(:,2)-correctConst*rescaledData(:,3);
    %This removes the [HbT] dependent attenuation from the recorded GCaMP
    %signals
    
    %rescaledData(:,5)=rescaledData(:,5)-correctConst*rescaledData(:,6);
    %This removes the [HbT] dependent attenuation from the recorded GCaMP
    %signals
end
%% Z-score data
avgData = nanmean(rescaledData,1);
stdData = std(rescaledData,0,1,'omitnan');
avgMatrix = repmat(avgData,length(rescaledData),1);
stdMatrix = repmat(stdData,length(rescaledData),1);
zScoredFiberData = (rescaledData - avgMatrix)./stdMatrix; %This is each channel z scored to its standard deviation. KWG

%% This might be redundant
nanLogical=double(isnan(zScoredFiberData(:,2)));
trialEnd=find(diff(nanLogical)==1); % Find where each imaging trial ends
trialEnd=[trialEnd;length(zScoredFiberData)];
trialStart=find(diff(nanLogical)==-1)+1; % Find where each imaginge trial following the first begins
trialStart=[1;trialStart];
%% split the data into trials
fields = {'RH_405','RH_465','RH_560','LH_405','LH_465','LH_560'};

trialCounter=1; % use this counter to track and save individual trial files.

for sessionNum=1:size(FiberData,2)
    for trialNum=1:size(FiberData(sessionNum).FitStruct,2)
        trialData.notes=FiberData(sessionNum).notes;
        trialData.FitStruct=FiberData(sessionNum).FitStruct(trialNum);
        for bb = 1:length(fields)
            trialData.Opticaldata.(fields{1,bb}).raw = catData((trialStart(trialCounter):trialEnd(trialCounter)),bb);
            trialData.Opticaldata.(fields{1,bb}).zScored = zScoredFiberData((trialStart(trialCounter):trialEnd(trialCounter)),bb);
            trialData.Opticaldata.(fields{1,bb}).rescaled=rescaledData((trialStart(trialCounter):trialEnd(trialCounter)),bb);
        end
        trialData.syncData=catSync(trialStart(trialCounter):trialEnd(trialCounter));
        if trialCounter < 10
            save([trialData.notes.animalID '_' trialData.notes.date '_Trial0'  num2str(trialCounter) '_FiberData.mat'],'trialData','-v7.3');
            trialCounter = trialCounter + 1;
        else
            save([trialData.notes.animalID '_' trialData.notes.date '_Trial'  num2str(trialCounter) '_FiberData.mat'],'trialData','-v7.3');
            trialCounter = trialCounter + 1;
        end
    end
end
end

function [lowPassData,FiberData]=correctDecay(filtData,channelData,sosLow,gLow,FiberData,trialNum,sessionNum)
    %% READ ME
    % This subfunction is used to remove slow trends from fiber photometry
    % Slow decay is fitted using a 2 exponent regression on 0.05Hz low pass
    % filtered data (filtData). This fit is then subtracted from the
    % raw data (channelData) prior to low-pass filtering @1Hz.
    %% INPUTS
    %filtData: optical signals low-pass filtered below 0.05Hz to regression
    %of metabolic decay/photobleaching
    %
    %channelData: raw optical channel data to have baseline drift
    %subtracted from
    %
    %sosLow/gLow: filter parameters for low-pass filtering detrended
    %signals
    %
    %FiberData: large data structure for loading fitting parameters and
    %results for storage and analysis.
    %
    %trialNum: which trial number within imaging session was fit on
    %function call
    %
    %sessionNum: which imaging session is being processed
    
    %% Remove photobleaching/metabolism baseline drift RH
    spacing = 1:1:length(filtData(:,3));
    % RH correct TRITC blood volume
    [fitVals] = fit(spacing',filtData(:,3),'exp2');
    coeffVals = coeffvalues(fitVals);
    RH_predictedCBV = (coeffVals(1)*exp((coeffVals(2).*spacing))) + (coeffVals(3)*exp((coeffVals(4).*spacing)));
    RH_correctedCBV = channelData(:,3) - RH_predictedCBV';
    FiberData(sessionNum).FitStruct(trialNum).RH_CBV = fitVals;
    % RH correct Ca2+ dependent GCaMP
    [fitVals] = fit(spacing',filtData(:,2),'exp2');
    coeffVals = coeffvalues(fitVals);
    RH_predictedGCaMP = (coeffVals(1)*exp((coeffVals(2).*spacing))) + (coeffVals(3)*exp((coeffVals(4).*spacing)));
    RH_correctedGCaMP = channelData(:,2) - RH_predictedGCaMP';
    FiberData(sessionNum).FitStruct(trialNum).RH_GCaMP = fitVals;
    % RH correct Ca2+ independent GCaMP
    [fitVals] = fit(spacing',filtData(:,1),'exp2');
    coeffVals = coeffvalues(fitVals);
    RH_predictedGFP = (coeffVals(1)*exp((coeffVals(2).*spacing))) + (coeffVals(3)*exp((coeffVals(4).*spacing)));
    RH_correctedGFP = channelData(:,1) - RH_predictedGFP';
    FiberData(sessionNum).FitStruct(trialNum).RH_GFP = fitVals;
    
    %% Remove photobleaching/metaloism baseline drift RH
    % LH correct TRITC blood volume
    [fitVals] = fit(spacing',filtData(:,6),'exp2');
    coeffVals = coeffvalues(fitVals);
    LH_predictedCBV = (coeffVals(1)*exp((coeffVals(2).*spacing))) + (coeffVals(3)*exp((coeffVals(4).*spacing)));
    LH_correctedCBV = channelData(:,6) - LH_predictedCBV';
    FiberData(sessionNum).FitStruct(trialNum).LH_CBV = fitVals;
    % LH correct Ca2+ dependent GCaMP
    [fitVals] = fit(spacing',filtData(:,5),'exp2');
    coeffVals = coeffvalues(fitVals);
    LH_predictedGCaMP = (coeffVals(1)*exp((coeffVals(2).*spacing))) + (coeffVals(3)*exp((coeffVals(4).*spacing)));
    LH_correctedGCaMP = channelData(:,5) - LH_predictedGCaMP';
    FiberData(sessionNum).FitStruct(trialNum).LH_GCaMP = fitVals;
    % LH correct Ca2+ independent GCaMP
    [fitVals] = fit(spacing',filtData(:,4),'exp2');
    coeffVals = coeffvalues(fitVals);
    LH_predictedGFP = (coeffVals(1)*exp((coeffVals(2).*spacing))) + (coeffVals(3)*exp((coeffVals(4).*spacing)));
    LH_correctedGFP = channelData(:,4) - LH_predictedGFP';
    FiberData(sessionNum).FitStruct(trialNum).LH_GFP = fitVals;
    
    %% detrend, low-pass filter, rescale data
    detrendData = [RH_correctedGFP,RH_correctedGCaMP,RH_correctedCBV,LH_correctedGFP,LH_correctedGCaMP,LH_correctedCBV]; %This is the optical data with photobleaching and metabolism dependent decreases removed- KWG
    lowPassData = filtfilt(sosLow,gLow,detrendData); %This is the data low-pass filtered below 1Hz for each of the sessions- KWG
end
