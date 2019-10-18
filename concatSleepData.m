params.targetMinutes = 30;   % minutes
params.minTime.Rest = 10;   % seconds

% list of unstim Procdata.mat files
procDataFileStruct = dir('*_Procdata.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);

% find and load RestData.mat struct
restDataFileStruct = dir('*_RestData.mat');
restDataFile = {restDataFileStruct.name}';
restDataFileID = char(restDataFile);
load(restDataFileID)

% find and load RestingBaselines.mat strut
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID)

% find and load SleepData.mat strut
sleepDataFileStruct = dir('*_SleepData.mat');
sleepDataFile = {sleepDataFileStruct.name}';
sleepDataFileID = char(sleepDataFile);
load(sleepDataFileID)

% identify animal's ID and pull important infortmat
fileBreaks = strfind(restDataFileID, '_');
animalID = restDataFileID(1:fileBreaks(1)-1);
trialDuration_min = RestData.CBV.LH.trialDuration_sec/60;   % min
manualFileIDs = unique(RestingBaselines.manualSelection.baselineFileInfo.fileIDs);
fileTarget = params.targetMinutes/trialDuration_min;
filterSets = {'manualSelection','setDuration','entireDuration'};
samplingRate = RestData.CBV.LH.CBVCamSamplingRate;

RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};

PuffCriteria.Fieldname = {'puffDistances'};
PuffCriteria.Comparison = {'gt'};
PuffCriteria.Value = {5};

%% Analyze coherence during periods of rest
% use the RestCriteria we specified earlier to find unstim resting events that are greater than the criteria
dataType = 'CBV';
[restLogical] = FilterEvents_IOS(RestData.(dataType).LH,RestCriteria);
[puffLogical] = FilterEvents_IOS(RestData.(dataType).LH,PuffCriteria);
combRestLogical = logical(restLogical.*puffLogical);
restFiles = RestData.(dataType).LH.fileIDs(combRestLogical,:);

% identify the unique days and the unique number of files from the list of unstim resting events
restUniqueDays = GetUniqueDays_IOS(restFiles);
restUniqueFiles = unique(restFiles);
restNumberOfFiles = length(unique(restFiles));

% decimate the file list to only include those files that occur within the desired number of target minutes
filterSet = 'setDuration';
clear restFiltLogical
for c = 1:length(restUniqueDays)
    restDay = restUniqueDays(c);
    d = 1;
    for e = 1:restNumberOfFiles
        restFile = restUniqueFiles(e);
        restFileID = restFile{1}(1:6);
        if strcmp(filterSet,'manualSelection') == true
            if strcmp(restDay,restFileID) && sum(strcmp(restFile,manualFileIDs)) == 1
                restFiltLogical{c,1}(e,1) = 1; %#ok<*AGROW>
                d = d + 1;
            else
                restFiltLogical{c,1}(e,1) = 0;
            end
        elseif strcmp(filterSet,'setDuration') == true
            if strcmp(restDay,restFileID) && d <= fileTarget
                restFiltLogical{c,1}(e,1) = 1;
                d = d + 1;
            else
                restFiltLogical{c,1}(e,1) = 0;
            end
        elseif strcmp(filterSet,'entireDuration') == true
            if strcmp(restDay,restFileID)
                restFiltLogical{c,1}(e,1) = 1;
                d = d + 1;
            else
                restFiltLogical{c,1}(e,1) = 0;
            end
        end
    end
end
restFinalLogical = any(sum(cell2mat(restFiltLogical'),2),2);

% extract unstim the resting events that correspond to the acceptable file list and the acceptable resting criteria
clear restFileFilter
filtRestFiles = restUniqueFiles(restFinalLogical,:);
for f = 1:length(restFiles)
    restLogic = strcmp(restFiles{f},filtRestFiles);
    restLogicSum = sum(restLogic);
    if restLogicSum == 1
        restFileFilter(f,1) = 1;
    else
        restFileFilter(f,1) = 0;
    end
end
restFinalFileFilter = logical(restFileFilter);

ExpectedTime=900;% trial time in seconds
for filenum=1:size(TheFiles,1)
    %% Load Data
    FinFo=whos('-file',TheFiles(filenum).name);
    fnOpto=zeros(1,size(FinFo,1));
    for subfield=1:size(FinFo,1)
        fnOpto(subfield)=strcmpi(FinFo(subfield).name,'Opto');
    end
    if ~strcmpi(FinFo(fnOpto==1).class,'struct')% Checks to see if correction for optogenetic stimulus artifacts occured, loads proper data fields
        load(TheFiles(filenum).name,'IOS','Spectrograms','Ephys','Behavior','AcquisitionParams');
        IOSData=IOS.forepaw.dHbT; %Because we are predicating our correlations on movement forepaw ROI will be most relevant
        clear IOS
    else
        load(TheFiles(filenum).name,'Opto','Spectrograms','Ephys','Behavior','AcquisitionParams');
        IOSData=Opto.IOS.forepaw.dHbT;%Because we are predicating our correlations on movement forepaw ROI will be most relevant
        clear Opto
    end
    if filenum==1
        HbTData=[];
        BallData=[];
        EMGData=[];
        AtoniaData=[];
        BinVelData=[];
        ThetaData=[];
        LowThetaData=[];
        HighThetaData=[];
        GammaData=[];
        DeltaData=[];
        expectedIOSlength=AcquisitionParams.StimParams.dal_fr*ExpectedTime;
        expectedANALOGlength=AcquisitionParams.StimParams.an_fs*ExpectedTime;
        
        % [z,p,k]=butter(3,[4 10]/(AcquisitionParams.StimParams.an_fs*0.5),'bandpass');
        % [sos_b,g_b]=zp2sos(z,p,k);
        %
        % [z,p,k]=butter(3,[4 6]/(AcquisitionParams.StimParams.an_fs*0.5),'bandpass');
        % [sos_l,g_l]=zp2sos(z,p,k);
        %
        % [z,p,k]=butter(3,[7 10]/(AcquisitionParams.StimParams.an_fs*0.5),'bandpass');
        % [sos_h,g_h]=zp2sos(z,p,k);
        %
        % [z,p,k]=butter(3,1/(AcquisitionParams.StimParams.an_fs*0.5),'low');
        % [sos_sm,g_sm]=zp2sos(z,p,k);
    end
    % AllTheta=resample(filtfilt(sos_sm,g_sm,filtfilt(sos_b,g_b,Ephys.RawNeuro).^2),AcquisitionParams.StimParams.dal_fr,AcquisitionParams.StimParams.an_fs);
    % LowTheta=resample(filtfilt(sos_sm,g_sm,filtfilt(sos_l,g_l,Ephys.RawNeuro).^2),AcquisitionParams.StimParams.dal_fr,AcquisitionParams.StimParams.an_fs);
    % HighTheta=resample(filtfilt(sos_sm,g_sm,filtfilt(sos_h,g_h,Ephys.RawNeuro).^2),AcquisitionParams.StimParams.dal_fr,AcquisitionParams.StimParams.an_fs);
    
    % if length(AllTheta)>=expectedIOSlength
    %     AllTheta=AllTheta(1:expectedIOSlength);
    %     LowTheta=LowTheta(1:expectedIOSlength);
    %     HighTheta=HighTheta(1:expectedIOSlength);
    % else
    %     AllTheta(length(AllTheta):expectedIOSlength)=0;
    %     LowTheta(length(LowTheta):expectedIOSlength)=0;
    %     HighTheta(length(HighTheta):expectedIOSlength)=0;
    % end
    
    %% Concatenate data from all trials
    TempVel=log10(abs(diff(Behavior.ballVelocity)));
    if length(TempVel)>=expectedIOSlength
        TempVel=TempVel(1:expectedIOSlength);
    else
        TempVel(length(TempVel):expectedIOSlength)=-5;
    end
    TempVel(TempVel<-4)=-5;
    
    if length(IOSData)>=expectedIOSlength
        IOSData=IOSData(1:expectedIOSlength);
    end
    
    if length(Behavior.LinkedBallVelocity)>=expectedIOSlength
        TempBinVel=Behavior.LinkedBallVelocity(1:expectedIOSlength);
    else
        TempBinVel=Behavior.LinkedBallVelocity;
        TempBinVel(length(TempBinVel):expectedIOSlength)=0;
    end
    
    if length(Ephys.downSampleEMG)>=expectedIOSlength
        TempEMG=Ephys.downSampleEMG(1:expectedIOSlength);
    else
        TempEMG=Ephys.downSampleEMG;
        TempEMG(length(TempEMG):expectedIOSlength)=0;
    end
    
    if length(Behavior.Flags.Atonia)>=expectedIOSlength
        TempAtonia=Behavior.Flags.Atonia(1:expectedIOSlength);
    else
        TempAtonia=Behavior.Flags.Atonia;
        TempAtonia(length(TempAtonia):expectedIOSlength)=0;
    end
    
    HbTData=[HbTData;IOSData];
    BallData=[BallData;TempVel];
    BinVelData=[BinVelData;TempBinVel];
    EMGData=[EMGData;log10(TempEMG)];
    AtoniaData=[AtoniaData;TempAtonia];
    % ThetaData=[ThetaData;AllTheta];
    % LowThetaData=[LowThetaData;LowTheta];
    % HighThetaData=[HighThetaData;HighTheta];
    GammaData=[GammaData;Ephys.normGammaBandPower];
    DeltaData=[DeltaData;Ephys.normDeltaBandPower];
end

%% Sliding window average of data
window_length=5; %time in seconds
winStep=0.2; %time in seconds
iosStep=winStep*AcquisitionParams.StimParams.dal_fr;
iosWin=window_length*AcquisitionParams.StimParams.dal_fr;
anStep=winStep*AcquisitionParams.StimParams.an_fs;
anWin=window_length*AcquisitionParams.StimParams.an_fs;
Win_Col(1:size(HbTData,1),:)=iosWin;
Step_Col(1:size(HbTData,1),:)=iosStep;
stepnum=(length(HbTData)-iosWin)/iosStep;

for thestep=1:stepnum
    winIOS(:,thestep)=mean(HbTData(:,(1+((thestep-1)*Step_Col):(((thestep-1)*Step_Col)+Win_Col))),2);
    winBall(:,thestep)=mean(BallData(:,(1+((thestep-1)*Step_Col):(((thestep-1)*Step_Col)+Win_Col))),2);
    winVel(:,thestep)=mean(BinVelData(:,(1+((thestep-1)*Step_Col):(((thestep-1)*Step_Col)+Win_Col))),2);
    winEMG(:,thestep)=mean(EMGData(:,(1+((thestep-1)*Step_Col):(((thestep-1)*Step_Col)+Win_Col))),2);
    winAtonia(:,thestep)=mean(AtoniaData(:,(1+((thestep-1)*Step_Col):(((thestep-1)*Step_Col)+Win_Col))),2);
%     winTheta(:,thestep)=mean(ThetaData(:,(1+((thestep-1)*Step_Col):(((thestep-1)*Step_Col)+Win_Col))),2);
%     winLoTheta(:,thestep)=mean(LowThetaData(:,(1+((thestep-1)*Step_Col):(((thestep-1)*Step_Col)+Win_Col))),2);
%     winHiTheta(:,thestep)=mean(HighThetaData(:,(1+((thestep-1)*Step_Col):(((thestep-1)*Step_Col)+Win_Col))),2);
    winDelta(:,thestep)=mean(DeltaData(:,(1+((thestep-1)*Step_Col):(((thestep-1)*Step_Col)+Win_Col))),2);
    winGamma(:,thestep)=mean(GammaData(:,(1+((thestep-1)*Step_Col):(((thestep-1)*Step_Col)+Win_Col))),2);
end

%% Linearize data
catIOS=reshape(winIOS',1,[]);
catBall=reshape(winBall',1,[]);
catVel=reshape(winVel',1,[]);
catEMG=reshape(winEMG',1,[]);
catAtonia=reshape(winAtonia',1,[]);
% catTheta=reshape(winTheta',1,[]);
catDelta=reshape(winDelta',1,[]);
catGamma=reshape(winGamma',1,[]);
% catLowTheta=reshape(winLoTheta',1,[]);
% catHiTheta=reshape(winHiTheta',1,[]);

%% Binarize data points based on behavior
ThreshAtonia=catAtonia;
ThreshAtonia(ThreshAtonia>=0.5)=1;
ThreshAtonia(ThreshAtonia<1)=0;

ThreshBinVel=catVel;
ThreshBinVel(ThreshBinVel>=0.5)=1;
ThreshBinVel(ThreshBinVel<1)=0;

AtoniaVals=intersect(find(ThreshAtonia==1),find(ThreshBinVel==0));
VelVals=intersect(find(ThreshAtonia==0),find(ThreshBinVel==1));

NonLoco=find(ThreshBinVel==0);
NonAtonia=find(ThreshAtonia==0);
StillVals=intersect(NonLoco,NonAtonia);

%% Create Colormap
cRange=round(catIOS,0);
Cnum=min(cRange):1:max(cRange);
[Cmap]=brewermap(numel(Cnum),'YlOrRd');

%% Save Data
IOSData=catIOS;
NormGamm=catGamma*100;
NormDelta=catDelta*100;
logEMG=catEMG;
logBall=catBall;
BinVals.Still=StillVals;
BinVals.Vel=VelVals;
BinVals.Atonia=AtoniaVals;

save('State_plot.mat','IOSData', 'NormGamm','NormDelta', 'logEMG','logBall','Cmap','BinVals');

%% Visualize data
%BeHaveScore=catBall.*real(catEMG); %larger values correspond to periods of low muscle tone and little to no movement, small values are high behavior periods
IOSData=IOSData+abs(min(catIOS));
pointSize=20;

%% Color reflects hemoglobin concentration/ shape behavior
figure(99);scatter3(NormDelta,real(logEMG),NormGamm,pointSize,IOSData,'filled');
colormap(jet);
% zlim([-100 100]);
% ylim([-4 1]);
% xlim([-100 100]);
% caxis([0 120]);
barObj=colorbar;
barObj.FontName='Arial';
barObj.FontSize=18;
barObj.Label.String='\Delta HbT \muM';
barObj.Label.FontName='Arial';
barObj.Label.FontSize=18;
zlabel('Gamma(40-100Hz) Power % change from baseline','FontName','Arial','FontSize',18);
ylabel('EMG Power','FontName','Arial','FontSize',18);
xlabel('Delta(1-4Hz) Power % change from baseline','FontName','Arial','FontSize',18);
title('P15: Blood volume increases can  be seen in two distinct behavioral states','FontName','Arial','FontSize',24);

figure(101);scatter(real(logEMG(BinVals.Vel)),NormGamm(BinVals.Vel),pointSize,IOSData(BinVals.Vel),'Marker','*'); hold on;
scatter(real(logEMG(BinVals.Still)),NormGamm(BinVals.Still),pointSize,IOSData(BinVals.Still),'Marker','+');
scatter(real(logEMG(BinVals.Atonia)),NormGamm(BinVals.Atonia),pointSize,IOSData(BinVals.Atonia),'Marker','o');
colormap(jet);
% ylim([-100 100]);
% xlim([-4 1]);   
% caxis([0 120]);
barObj=colorbar;
barObj.FontName='Arial';
barObj.FontSize=18;
barObj.Label.String='\Delta HbT \muM';
barObj.Label.FontName='Arial';
barObj.Label.FontSize=18;
ylabel('Gamma(40-100Hz) Power % change from baseline','FontName','Arial','FontSize',18);
xlabel('EMG Power','FontName','Arial','FontSize',18);
title('P15: Blood volume increases can  be seen in two distinct behavioral states','FontName','Arial','FontSize',24);


% %% Color represents HbT at time point
% figure(104);scatter3(BeHaveScore,catDelta*100,catGamma*100,15,VizIOS,'filled');
% colormap(Cmap);
% xlim([-2 16]);
% ylim([-100 100]);
% zlim([-50 150]);
% caxis([0 120]);
% barObj=colorbar;
% barObj.FontName='Arial';
% barObj.FontSize=18;
% barObj.Label.String='\Delta HbT \muM';
% barObj.Label.FontName='Arial';
% barObj.Label.FontSize=18;
% zlabel('Gamma(40-100Hz) Power % change from baseline','FontName','Arial','FontSize',18);
% xlabel('Behavior coefficient','FontName','Arial','FontSize',18);
% ylabel('Delta(1-4Hz) Power % change from baseline','FontName','Arial','FontSize',18);
% title('P15: Blood volume increases can  be seen in two distinct behavioral states','FontName','Arial','FontSize',24);
% 
% figure(106);scatter(BeHaveScore,catGamma*100,15,VizIOS,'filled');
% colormap(Cmap);
% xlim([-2 16]);
% ylim([-50 150]);
% caxis([0 120]);
% barObj=colorbar;
% barObj.FontName='Arial';
% barObj.FontSize=18;
% barObj.Label.String='\Delta HbT \muM';
% barObj.Label.FontName='Arial';
% barObj.Label.FontSize=18;
% ylabel('Gamma(40-100Hz) Power % change from baseline','FontName','Arial','FontSize',18);
% xlabel('Behavior coefficient','FontName','Arial','FontSize',18);
% title('P15: Blood volume increases can  be seen in two distinct behavioral states','FontName','Arial','FontSize',24);
% 
% %% Color is binarized behavior time
% figure(105);scatter3(BeHaveScore(AtoniaVals),catDelta(AtoniaVals)*100,catGamma(AtoniaVals)*100,15,'g','filled');
% hold on; scatter3(BeHaveScore(VelVals),catDelta(VelVals)*100,catGamma(VelVals)*100,15,'r','filled');
% scatter3(BeHaveScore(StillVals),catDelta(StillVals)*100,catGamma(StillVals)*100,15,'b','filled');
% xlim([-2 16]);
% ylim([-100 100]);
% zlim([-50 150]);
% zlabel('Gamma(40-100Hz) Power % change from baseline','FontName','Arial','FontSize',18);
% xlabel('Behavior coefficient','FontName','Arial','FontSize',18);
% ylabel('Delta(1-4Hz) Power % change from baseline','FontName','Arial','FontSize',18);
% title('P15: Binarized behavior can be clustered based on Behavior coefficient','FontName','Arial','FontSize',24);
% legend('Muscle Atonia','Binarized Locomotion','Quiescent','FontName','Arial','FontSize',18);
% 
% figure(107);scatter(BeHaveScore(AtoniaVals),catGamma(AtoniaVals)*100,15,'g','filled');
% hold on; scatter(BeHaveScore(VelVals),catGamma(VelVals)*100,15,'r','filled');
% scatter(BeHaveScore(StillVals),catGamma(StillVals)*100,15,'b','filled');
% xlim([-2 16]);
% ylim([-50 150]);
% ylabel('Gamma(40-100Hz) Power % change from baseline','FontName','Arial','FontSize',18);
% xlabel('Behavior coefficient','FontName','Arial','FontSize',18);
% title('P15: Binarized behavior can be clustered based on Behavior coefficient','FontName','Arial','FontSize',24);
% legend('Muscle Atonia','Binarized Locomotion','Quiescent','FontName','Arial','FontSize',18);
