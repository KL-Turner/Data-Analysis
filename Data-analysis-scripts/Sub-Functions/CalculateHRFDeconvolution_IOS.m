function [AnalysisResults] = CalculateHRFDeconvolution_IOS(neuralBand,hemisphere,behavior,AnalysisResults)

%% Load and Setup
HRFLims = [0 5];
HRFParams.offset = 3;
HRFParams.dur = 10;
Event_Inds.CalcStart = 1;
Event_Inds.TestStart = 2;
Event_Inds.Increment = 2;

disp(['CalculateHRFDeconvolution: ' hemisphere ' ' neuralBand ' ' behavior]); disp(' ')
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID)
fileBreaks = strfind(baselineDataFileID, '_');
animalID = baselineDataFileID(1:fileBreaks(1)-1);
manualFileIDs = unique(RestingBaselines.manualSelection.baselineFileInfo.fileIDs);

if strcmp(behavior,'Rest')
    restDataFileStruct = dir('*_RestData.mat');
    restDataFile = {restDataFileStruct.name}';
    restDataFileID = char(restDataFile);
    load(restDataFileID)
    BehData = RestData;
    clear RestData;
elseif strcmp(behavior,'Whisk') || strcmp(behavior,'Contra')
    eventDataFileStruct = dir('*_EventData.mat');
    eventDataFile = {eventDataFileStruct.name}';
    eventDataFileID = char(eventDataFile);
    load(eventDataFileID)
    BehData = EventData;
    clear EventData;
elseif strcmp(behavior,'NREM') == true || strcmp(behavior,'REM')
    sleepDataFileStruct = dir('*_SleepData.mat');
    sleepDataFile = {sleepDataFileStruct.name}';
    sleepDataFileID = char(sleepDataFile);
    load(sleepDataFileID)
end

%% Get the arrays for the calculation
if strcmp(behavior,'Contra') == true || strcmp(behavior,'Whisk') == true || strcmp(behavior,'Rest') == true
    [NeuralDataStruct,NeuralFiltArray] = SelectConvolutionBehavioralEvents_IOS(BehData.(['cortical_' hemisphere(4:end)]).(neuralBand),behavior,hemisphere);
    [HemoDataStruct,HemoFiltArray] = SelectConvolutionBehavioralEvents_IOS(BehData.CBV.(hemisphere),behavior,hemisphere);
    fileIDs = NeuralDataStruct.fileIDs;
    restUniqueDays = GetUniqueDays_IOS(fileIDs);
    restUniqueFiles = unique(fileIDs);
    restNumberOfFiles = length(unique(fileIDs));
    clear restFiltLogical
    for c = 1:length(restUniqueDays)
        restDay = restUniqueDays(c);
        d = 1;
        for e = 1:restNumberOfFiles
            restFile = restUniqueFiles(e);
            restFileID = restFile{1}(1:6);
            if strcmp(restDay,restFileID) && sum(strcmp(restFile,manualFileIDs)) == 1
                restFiltLogical{c,1}(e,1) = 1; %#ok<*AGROW>
                d = d + 1;
            else
                restFiltLogical{c,1}(e,1) = 0;
            end
        end
    end
    restFinalLogical = any(sum(cell2mat(restFiltLogical'),2),2);
    
    clear restFileFilter
    filtRestFiles = restUniqueFiles(restFinalLogical,:);
    for f = 1:length(fileIDs)
        restLogic = strcmp(fileIDs{f},filtRestFiles);
        restLogicSum = sum(restLogic);
        if restLogicSum == 1
            restFileFilter(f,1) = 1;
        else
            restFileFilter(f,1) = 0;
        end
    end
    restFinalFileFilter = logical(restFileFilter);
    filtArrayEdit1 = logical(NeuralFiltArray.*restFinalFileFilter);
    NormData1 = NeuralDataStruct.NormData(filtArrayEdit1,:);
    filtArrayEdit2 = logical(HemoFiltArray.*restFinalFileFilter);
    NormData2 = HemoDataStruct.NormData(filtArrayEdit2,:);
    % [B, A] = butter(3,1/(30/2),'low');
    % if strcmp(behavior,'Contra') == true || strcmp(behavior,'Whisk') == true
    %     for a = 1:size(NormData1,1)
    %         NormData1(a,:) = filtfilt(B,A,NormData1(a,:));
    %         NormData2(a,:) = filtfilt(B,A,NormData2(a,:));
    %     end
    % elseif strcmp(behavior,'Rest') == true
    %     for a = 1:size(NormData1,1)
    %         NormData1{a,:} = filtfilt(B,A,NormData1{a,:});
    %         NormData2{a,:} = filtfilt(B,A,NormData2{a,:});
    %     end
    % end
elseif strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true
    NormData1 = SleepData.(behavior).data.(['cortical_' hemisphere(4:end)]).(neuralBand);
    NormData2 = SleepData.(behavior).data.CBV.(hemisphere(4:end));
end

%% Separate events for HRF calculation from events used for later testing.
% Insert padding of zeros with size equal to the HRF between individual
% events.
zpad_inds = Event_Inds.TestStart:Event_Inds.Increment:size(NormData1,1);
calc_inds = Event_Inds.CalcStart:Event_Inds.Increment:size(NormData1,1);
NeuralDataStruct.samplingRate = 30;
HemoDataStruct.samplingRate = 30;
zpad1 = zeros(1,HRFParams.dur*NeuralDataStruct.samplingRate);
zpad2 = zeros(1,HRFParams.dur*HemoDataStruct.samplingRate);

if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true 
    NormData1(zpad_inds) = {zpad1};
    % Mean subtract the data
    Processed1 = cell(size(NormData1));
    for c = 1:length(NormData1)
        template = zeros(size(NormData1{c}));
        strt = 2*NeuralDataStruct.samplingRate;
        stp = size(template,2);
        template(:,strt:stp) = detrend(NormData1{c}(:,strt:stp)-mean(NormData1{c}(:,strt:stp)));
        Processed1{c} = template;
    end
    Data1 = [Processed1{:}];
    clear Processed1;
elseif strcmp(behavior,'Whisk')
    Data1_end = 4;
    strt = round((NeuralDataStruct.epoch.offset-1)*NeuralDataStruct.samplingRate);
    stp = strt + round(Data1_end*NeuralDataStruct.samplingRate);
    template = zeros(size(NormData1(calc_inds,:)));
    offset1 = mean(NormData1(calc_inds,1:strt),2)*ones(1,stp-strt+1);
    template(:,strt:stp) = NormData1(calc_inds,strt:stp)-offset1;
    Data1Pad = [template ones(length(calc_inds),1)*zpad1];
    Data1 = reshape(Data1Pad',1,numel(Data1Pad));
elseif strcmp(behavior,'Contra')
    Data1_end = 1.5;
    strt = round((NeuralDataStruct.epoch.offset)*NeuralDataStruct.samplingRate);
    stp = strt + (Data1_end*NeuralDataStruct.samplingRate);
    template = zeros(size(NormData1(calc_inds,:)));
    offset1 = mean(NormData1(calc_inds,1:strt),2)*ones(1,stp-strt+1);
    template(:,strt:stp) = NormData1(calc_inds,strt:stp) - offset1;
    Data1Pad = [template ones(length(calc_inds),1)*zpad1];
    Data1 = reshape(Data1Pad',1,numel(Data1Pad));
end

if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true 
    NormData2(zpad_inds) = {zpad2};
    Processed2 = cell(size(NormData2));
    for c = 1:length(NormData2)
        template = zeros(size(NormData2{c}));
        strt = 2*NeuralDataStruct.samplingRate;
        stp = size(template,2);
        offset = mean(NormData2{c})*ones(1,stp-strt+1);
        template(:,strt:stp) = detrend(NormData2{c}(:,strt:stp)-offset);
        Processed2{c} = template;
    end
    Data2 = [Processed2{:}];
    clear Processed2
elseif strcmp(behavior,'Whisk')
    Data2_end = 6;
    strt = round((HemoDataStruct.epoch.offset-1)*HemoDataStruct.samplingRate);
    stp = strt + round(Data2_end*HemoDataStruct.samplingRate);
    template = zeros(size(NormData2(calc_inds,:)));
    offset2 = mean(NormData2(calc_inds,1:strt),2)*ones(1,stp-strt+1);
    template(:,strt:stp) = NormData2(calc_inds,strt:stp)-offset2;
    Data2Pad = [template ones(length(calc_inds),1)*zpad2];
    Data2 = reshape(Data2Pad',1,numel(Data2Pad));
elseif strcmp(behavior,'Contra')
    Data2_end = 3;
    strt = round(HemoDataStruct.epoch.offset*HemoDataStruct.samplingRate);
    stp = strt + (Data2_end*HemoDataStruct.samplingRate);
    template = zeros(size(NormData2(calc_inds,:)));
    offset2 = mean(NormData2(calc_inds,1:strt),2)*ones(1,stp-strt+1);
    template(:,strt:stp) = NormData2(calc_inds,strt:stp)-offset2;
    Data2Pad = [template ones(length(calc_inds),1)*zpad2];
    Data2 = reshape(Data2Pad',1,numel(Data2Pad));
end

%% Calculate HRF based on deconvolution
samplingRate = HemoDataStruct.samplingRate;
IR_est=IR_analytic_IOS(Data1',Data2',HRFParams.offset*samplingRate,HRFParams.dur*samplingRate);

HRF = sgolayfilt(IR_est.IR',3,samplingRate+1);
timevec = (1:length(HRF))/samplingRate-HRFParams.offset;
num_events = length(calc_inds);
TimeLims = timevec>=HRFLims(1) & timevec<=HRFLims(2);
timeLimHRF = HRF(TimeLims);
timeLimVec = timevec(TimeLims);

AnalysisResults.HRFs.(neuralBand).(hemisphere).(behavior).IR = timeLimHRF;
AnalysisResults.HRFs.(neuralBand).(hemisphere).(behavior).IRtimeVec = timeLimVec;
AnalysisResults.HRFs.(neuralBand).(hemisphere).(behavior).HRFParams = HRFParams;
AnalysisResults.HRFs.(neuralBand).(hemisphere).(behavior).num_calc_events = num_events;
AnalysisResults.HRFs.(neuralBand).(hemisphere).(behavior).Event_Inds = Event_Inds;

kernelFig = figure;
sgtitle([animalID ' ' hemisphere ' ' neuralBand ' during ' behavior])
subplot(1,2,1)
plot(AnalysisResults.HRFs.(neuralBand).(hemisphere).(behavior).IRtimeVec,AnalysisResults.HRFs.(neuralBand).(hemisphere).(behavior).IR,'k')
title('Impulse response function')
ylabel('A.U')
xlabel('Time (s)')
axis square
axis tight

%% Calculate the gamma HRF
options = optimset('MaxFunEvals',2e3,'MaxIter',2e3,'TolFun',1e-7,'TolX',1e-7);
initvals = [1e-1,1,1];
HRFDur = 5; % seconds
[gam_params,~,~] = fminsearch(@(x)gammaconvolve_IOS(x,Data1,Data2,HemoDataStruct.samplingRate,HRFDur),initvals,options);
t = 0:1/HemoDataStruct.samplingRate:HRFDur;
a = ((gam_params(2)/gam_params(3))^2*8*log10(2));
beta = ((gam_params(3)^2)/gam_params(2)/8/log10(2));
gamma = gam_params(1)*(t/gam_params(2)).^a.*exp((t-gam_params(2))/(-1*beta));
AnalysisResults.HRFs.(neuralBand).(hemisphere).(behavior).gammaFunc = gamma;
AnalysisResults.HRFs.(neuralBand).(hemisphere).(behavior).gammaTimeVec = t;
subplot(1,2,2)
plot(AnalysisResults.HRFs.(neuralBand).(hemisphere).(behavior).gammaTimeVec,AnalysisResults.HRFs.(neuralBand).(hemisphere).(behavior).gammaFunc,'k')
title('Gamma function')
ylabel('A.U')
xlabel('Time (s)')
axis square
axis tight

% save figures
[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Combined Imaging/Figures/HRF Kernels/'];
if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end
savefig(kernelFig,[dirpath animalID '_' hemisphere '_' neuralBand '_' behavior '_HRFs']);
close(kernelFig)
%% save results struct
save([animalID '_AnalysisResults.mat'],'AnalysisResults');

end
