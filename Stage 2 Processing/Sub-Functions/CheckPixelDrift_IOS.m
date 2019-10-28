function [] = CheckPixelDrift_IOS(procDataFileIDs,imagingType)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Determine the slow exponential drift for each day of imaging and correct the drift if desired. The slow
%            drift is caused by the CCD sensor's sensitivity changing as the camera heats up over multiple hours.
%________________________________________________________________________________________________________________________

% establish the number of unique days based on file IDs
[~,fileDates,~] = GetFileInfo_IOS(procDataFileIDs);
[uniqueDays, ~, DayID] = GetUniqueDays_IOS(fileDates);
firstsFileOfDay = cell(1,length(uniqueDays));
for a = 1:length(uniqueDays)
    FileInd = DayID == a;
    dayFilenames = procDataFileIDs(FileInd,:);
    firstsFileOfDay(a) = {dayFilenames(1,:)};
end

% go through each day and concate the data to observe slow drift
for b = 1:length(firstsFileOfDay)
    indDayProcDataFileList = {};
    catLH_CBVdata = [];
    catRH_CBVdata = [];
    catBarrels_CBVdata = [];
    catCementData = [];
    fileName = firstsFileOfDay{1,b};
    [~, fileDate, ~] = GetFileInfo_IOS(fileName);
    p = 1;
    for c = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(c,:);
        if strfind(procDataFileID, fileDate) >= 1 
            indDayProcDataFileList{p,1} = procDataFileID;
            p = p + 1;
        end
    end
    % load the processed CBV/cement data from each file and concat it into one array
    for d = 1:length(indDayProcDataFileList)
        indDayProcDataFile = indDayProcDataFileList{d,1};
        load(indDayProcDataFile)
        samplingRate = ProcData.notes.CBVCamSamplingRate;
        trialDuration = ProcData.notes.trialDuration_sec;
        if strcmp(imagingType,'bilateral') == true
            LH_CBVdata = ProcData.data.CBV.LH;
            RH_CBVdata = ProcData.data.CBV.RH;
            cementData = ProcData.data.CBV.Cement;
            catLH_CBVdata = horzcat(catLH_CBVdata,LH_CBVdata);
            catRH_CBVdata = horzcat(catRH_CBVdata,RH_CBVdata);
            catCementData = horzcat(catCementData,cementData);
        elseif strcmp(imagingType,'single') == true
            barrels_CBVdata = ProcData.data.CBV.Barrels;
            cementData = ProcData.data.CBV.Cement;
            catBarrels_CBVdata = horzcat(catBarrels_CBVdata,barrels_CBVdata);
            catCementData = horzcat(catCementData,cementData);
        end
    end
    % establish whether a slow exponential trend exists for the data
    [B,A] = butter(3,0.01/(samplingRate/2),'low');
    if strcmp(imagingType,'bilateral') == true
        filtCatCementData = filtfilt(B,A,catCementData);
        x = ((1:length(filtCatCementData))/samplingRate)';
        % create a weight vector for the trend
        weightVec = ones(1,length(x));
        secondHalfMean = mean(filtCatCementData(floor(length(filtCatCementData/2)):end));
        for t = 1:length(weightVec)
            if filtCatCementData(t) > secondHalfMean
                weightVec(t) = 10;
            end
        end
        % compare weighted vs. unweighted fits
        modelFit1 = fit(x,filtCatCementData','exp2');
        modelFit1_Y = modelFit1(x);
        modelFit2 = fit(x,filtCatCementData','exp2','Weight',weightVec);
        modelFit2_Y = modelFit2(x);
        modelFit2_norm = (modelFit2_Y - min(modelFit2_Y))./min(modelFit2_Y);
        modelFit2_flip = 1 - modelFit2_norm;
        adjCatLH_CBVdata = catLH_CBVdata.*modelFit2_flip';
        rsAdjCatLH_CBVdata = reshape(adjCatLH_CBVdata,[length(indDayProcDataFileList) samplingRate*trialDuration]);
        adjCatRH_CBVdata = catRH_CBVdata.*modelFit2_flip';
        rsAdjCatRH_CBVdata = reshape(adjCatRH_CBVdata,[length(indDayProcDataFileList) samplingRate*trialDuration]);
        % comparison showing original LH data and the corrected data
        figure;
        subplot(2,2,1)
        plot(x,catLH_CBVdata,'k')
        title('Original LH data')
        xlabel('Time (sec)')
        ylabel('12-bit Pixel val')
        axis tight
        subplot(2,2,2)
        plot(x,filtCatCementData,'k')
        hold on
        plot(x,modelFit1_Y,'b')
        plot(x,modelFit2_Y,'r')
        title('Cement drift')
        xlabel('Time (sec)')
        ylabel('12-bit Pixel val')
        legend('Cement data','exp fit','weighted fit')
        axis tight
        subplot(2,2,3)
        plot(x,modelFit2_flip,'r')
        title('Correction profile')
        xlabel('Time (sec)')
        ylabel('Normalized val')
        axis tight
        subplot(2,2,4)
        plot(x,catLH_CBVdata,'k')
        hold on
        plot(x,adjCatLH_CBVdata,'r')
        title('Corrected LH data')
        xlabel('Time (sec)')
        ylabel('12-bit Pixel val')
        legend('Original data','Corrected data')
        axis tight
        % comparison showing original RH data and the corrected data
        figure;
        subplot(2,2,1)
        plot(x,catRH_CBVdata,'k')
        title('Original RH data')
        xlabel('Time (sec)')
        ylabel('12-bit Pixel val')
        axis tight
        subplot(2,2,2)
        plot(x,filtCatCementData,'k')
        hold on
        plot(x,modelFit1_Y,'b')
        plot(x,modelFit2_Y,'r')
        title('Cement drift')
        xlabel('Time (sec)')
        ylabel('12-bit Pixel val')
        legend('Cement data','exp fit','weighted fit')
        axis tight
        subplot(2,2,3)
        plot(x,modelFit2_flip,'r')
        title('Correction profile')
        xlabel('Time (sec)')
        ylabel('Normalized val')
        axis tight
        subplot(2,2,4)
        plot(x,catRH_CBVdata,'k')
        hold on
        plot(x,adjCatRH_CBVdata,'r')
        title('Corrected RH data')
        xlabel('Time (sec)')
        ylabel('12-bit Pixel val')
        legend('Original data','Corrected data')
        axis tight
    elseif strcmp(imagingType,'single') == true
        filtCatCementData = filtfilt(B,A,catCementData);
        x = ((1:length(filtCatCementData))/samplingRate)';
        % create a weight vector for the trend
        weightVec = ones(1,length(x));
        secondHalfMean = mean(filtCatCementData(floor(length(filtCatCementData/2)):end));
        for t = 1:length(weightVec)
            if filtCatCementData(t) > secondHalfMean
                weightVec(t) = 10;
            end
        end
        % compare weighted vs. unweighted fits
        modelFit1 = fit(x,filtCatCementData','exp2');
        modelFit1_Y = modelFit1(x);
        modelFit2 = fit(x,filtCatCementData','exp2','Weight',weightVec);
        modelFit2_Y = modelFit2(x);
        modelFit2_norm = (modelFit2_Y - min(modelFit2_Y))./min(modelFit2_Y);
        modelFit2_flip = 1 - modelFit2_norm;
        adjCatBarrels_CBVdata = catBarrels_CBVdata.*modelFit2_flip';
        rsAdjCatBarrels_CBVdata = reshape(adjCatBarrels_CBVdata,[length(indDayProcDataFileList) samplingRate*trialDuration_sec]);
        % comparison showing original LH data and the corrected data
        figure;
        subplot(2,2,1)
        plot(x,catBarrels_CBVdata,'k')
        title('Original LH data')
        xlabel('Time (sec)')
        ylabel('12-bit Pixel val')
        axis tight
        subplot(2,2,2)
        plot(x,filtCatCementData,'k')
        hold on
        plot(x,modelFit1_Y,'b')
        plot(x,modelFit2_Y,'r')
        title('Cement drift')
        xlabel('Time (sec)')
        ylabel('12-bit Pixel val')
        legend('Cement data','exp fit','weighted fit')
        axis tight
        subplot(2,2,3)
        plot(x,modelFit2_flip,'r')
        title('Correction profile')
        xlabel('Time (sec)')
        ylabel('Normalized val')
        axis tight
        subplot(2,2,4)
        plot(x,catBarrels_CBVdata,'k')
        hold on
        plot(x,adjCatBarrels_CBVdata,'r')
        title('Corrected LH data')
        xlabel('Time (sec)')
        ylabel('12-bit Pixel val')
        legend('Original data','Corrected data')
        axis tight
    end
    correctionDecision = 'n';
    while strcmp(correctionDecision,'n') == true
        applyCorrection = input('Apply correction profile to pixel values? (y/n): ','s'); disp(' ')
        if strcmp(applyCorrection,'y') == true || strcmp(applyCorrection,'n') == true
            correctionDecision = 'y';
        else
            disp('Invalid input. Must be ''y'' or ''n'''); disp(' ')
        end
    end
    % apply correction to each file
    if strcmp(correctionDecision,'y') == true
        for d = 1:length(indDayProcDataFileList)
            indDayProcDataFile = indDayProcDataFileList{d,1};
            load(indDayProcDataFile)
            if strcmp(imagingType,'bilateral') == true
                ProcData.data.CBV.correctedLH = rsAdjCatLH_CBVdata(:,d);
                ProcData.data.CBV.correctedRH = rsAdjCatRH_CBVdata(:,d);
            elseif strcmp(imagingType,'single') == true
                ProcData.data.CBV.correctedBarrels = rsAdjCatBarrels_CBVdata(:,d);
            end
            save(indDayProcDataFile,'ProcData')
        end
    elseif strcmp(correctionDecision,'n') == true
        for d = 1:length(indDayProcDataFileList)
            indDayProcDataFile = indDayProcDataFileList{d,1};
            load(indDayProcDataFile)
            if strcmp(imagingType,'bilateral') == true
                ProcData.data.CBV.correctedLH = ProcData.data.CBV.LH;
                ProcData.data.CBV.correctedRH = ProcData.data.CBV.RH;
            elseif strcmp(imagingType,'single') == true
                ProcData.data.CBV.correctedBarrels = ProcData.data.CBV.Barrels;
            end
            save(indDayProcDataFile,'ProcData')
        end
    end
end

end
