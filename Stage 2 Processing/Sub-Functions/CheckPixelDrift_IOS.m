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
[animalIDs,fileDates,~] = GetFileInfo_IOS(procDataFileIDs);
animalID = animalIDs(1,:);
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
    catLH_cementData = [];
    catRH_cementData = [];
    catBarrels_CBVdata = [];
    catCementData = [];
    fileName = firstsFileOfDay{1,b};
    [~, fileDate, ~] = GetFileInfo_IOS(fileName);
    strDay = ConvertDate_IOS(fileDate);
    p = 1;
    for c = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(c,:);
        if strfind(procDataFileID, fileDate) >= 1 
            indDayProcDataFileList{p,1} = procDataFileID;
            p = p + 1;
        end
    end
    %% load the processed CBV/cement data from each file and concat it into one array
    for d = 1:length(indDayProcDataFileList)
        indDayProcDataFile = indDayProcDataFileList{d,1};
        load(indDayProcDataFile)
        samplingRate = ProcData.notes.CBVCamSamplingRate;
        trialDuration = ProcData.notes.trialDuration_sec;
        if strcmp(imagingType,'bilateral') == true
            LH_CBVdata = ProcData.data.CBV.LH;
            RH_CBVdata = ProcData.data.CBV.RH;
            LH_cementData = ProcData.data.CBV.LH_Cement;
            RH_cementData = ProcData.data.CBV.RH_Cement;
            catLH_CBVdata = horzcat(catLH_CBVdata,LH_CBVdata);
            catRH_CBVdata = horzcat(catRH_CBVdata,RH_CBVdata);
            catLH_cementData = horzcat(catLH_cementData,LH_cementData);
            catRH_cementData = horzcat(catRH_cementData,RH_cementData);
        elseif strcmp(imagingType,'single') == true
            barrels_CBVdata = ProcData.data.CBV.Barrels;
            LH_cementData = ProcData.data.CBV.Cement;
            catBarrels_CBVdata = horzcat(catBarrels_CBVdata,barrels_CBVdata);
            catCementData = horzcat(catCementData,LH_cementData);
        end
    end
    %% establish whether a slow exponential trend exists for the data
    [B,A] = butter(3,0.01/(samplingRate/2),'low');
    if strcmp(imagingType,'bilateral') == true
        filtCatLH_cementData = filtfilt(B,A,catLH_cementData);
        filtCatRH_cementData = filtfilt(B,A,catRH_cementData);
        x = ((1:length(filtCatLH_cementData))/samplingRate)';
        % create a weight vector for the trend
        LH_weightVec = ones(1,length(x));
        RH_weightVec = ones(1,length(x));
        LH_secondHalfMean = mean(filtCatLH_cementData(floor(length(filtCatLH_cementData/2)):end));
        RH_secondHalfMean = mean(filtCatRH_cementData(floor(length(filtCatRH_cementData/2)):end));
        for t = 1:length(LH_weightVec)
            if filtCatLH_cementData(t) > LH_secondHalfMean
                LH_weightVec(t) = 10;
            end
            if filtCatRH_cementData(t) > RH_secondHalfMean
                RH_weightVec(t) = 10;
            end
        end
        % compare weighted vs. unweighted fits
        LH_modelFit1 = fit(x,filtCatLH_cementData','exp2');
        RH_modelFit1 = fit(x,filtCatRH_cementData','exp2');
        LH_modelFit1_Y = LH_modelFit1(x);
        RH_modelFit1_Y = RH_modelFit1(x);
        LH_modelFit2 = fit(x,filtCatLH_cementData','exp2','Weight',LH_weightVec);
        RH_modelFit2 = fit(x,filtCatRH_cementData','exp2','Weight',RH_weightVec);
        LH_modelFit2_Y = LH_modelFit2(x);
        RH_modelFit2_Y = RH_modelFit2(x);
        LH_modelFit2_norm = (LH_modelFit2_Y - min(LH_modelFit2_Y))./min(LH_modelFit2_Y);
        RH_modelFit2_norm = (RH_modelFit2_Y - min(RH_modelFit2_Y))./min(RH_modelFit2_Y);
        LH_modelFit2_flip = 1 - LH_modelFit2_norm;
        RH_modelFit2_flip = 1 - RH_modelFit2_norm;
        % apply exponential correction to original data
        adjCatLH_CBVdata = catLH_CBVdata.*LH_modelFit2_flip';
        rsAdjCatLH_CBVdata = reshape(adjCatLH_CBVdata,[length(indDayProcDataFileList) samplingRate*trialDuration]);
        adjCatRH_CBVdata = catRH_CBVdata.*RH_modelFit2_flip';
        rsAdjCatRH_CBVdata = reshape(adjCatRH_CBVdata,[length(indDayProcDataFileList) samplingRate*trialDuration]);
        % comparison showing original LH data and the corrected data
        fixPixels = figure;
        subplot(2,4,1)
        plot(x,catLH_CBVdata,'k')
        title('LH original data')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        axis tight
        subplot(2,4,2)
        plot(x,filtCatLH_cementData,'k')
        hold on
        plot(x,LH_modelFit1_Y,'b')
        plot(x,LH_modelFit2_Y,'r')
        title('LH cement drift')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        legend('cement data','exp fit','weighted fit')
        axis tight
        subplot(2,4,5)
        plot(x,LH_modelFit2_flip,'r')
        title('LH correction profile')
        xlabel('Time (sec)')
        ylabel('Normalized val')
        axis tight
        subplot(2,4,6)
        plot(x,catLH_CBVdata,'k')
        hold on
        plot(x,adjCatLH_CBVdata,'r')
        title('LH corrected data')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        legend('original data','Corrected data')
        axis tight
        subplot(2,4,3)
        plot(x,catRH_CBVdata,'k')
        title('RH original data')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        axis tight
        subplot(2,4,4)
        plot(x,filtCatRH_cementData,'k')
        hold on
        plot(x,RH_modelFit1_Y,'b')
        plot(x,RH_modelFit2_Y,'r')
        title('RH cement drift')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        legend('cement data','exp fit','weighted fit')
        axis tight
        subplot(2,4,7)
        plot(x,RH_modelFit2_flip,'r')
        title('RH correction profile')
        xlabel('Time (sec)')
        ylabel('Normalized val')
        axis tight
        subplot(2,4,8)
        plot(x,catRH_CBVdata,'k')
        hold on
        plot(x,adjCatRH_CBVdata,'r')
        title('RH corrected data')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        legend('original data','corrected data')
        axis tight
        correctionDecision = 'n';
        while strcmp(correctionDecision,'n') == true
            applyCorrection = input('Apply correction profile to pixel values? (y/n): ','s'); disp(' ')
            if strcmp(applyCorrection,'y') == true || strcmp(applyCorrection,'n') == true
                correctionDecision = 'y';
            else
                disp('Invalid input. Must be ''y'' or ''n'''); disp(' ')
            end
        end
        sgtitle([animalID ' ' strDay ' pixel correction applied: ' correctionDecision]) 
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
        rsAdjCatBarrels_CBVdata = reshape(adjCatBarrels_CBVdata,[length(indDayProcDataFileList) samplingRate*trialDuration]);
        % comparison showing original LH data and the corrected data
        fixPixels = figure;
        subplot(2,2,1)
        plot(x,catBarrels_CBVdata,'k')
        title('Barrels original data')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        axis tight
        subplot(2,2,2)
        plot(x,filtCatCementData,'k')
        hold on
        plot(x,modelFit1_Y,'b')
        plot(x,modelFit2_Y,'r')
        title('cement drift')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        legend('cement data','exp fit','weighted fit')
        axis tight
        subplot(2,2,3)
        plot(x,modelFit2_flip,'r')
        title('correction profile')
        xlabel('Time (sec)')
        ylabel('Normalized val')
        axis tight
        subplot(2,2,4)
        plot(x,catBarrels_CBVdata,'k')
        hold on
        plot(x,adjCatBarrels_CBVdata,'r')
        title('Corrected data')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        legend('original data','corrected data')
        axis tight
        correctionDecision = 'n';
        while strcmp(correctionDecision,'n') == true
            applyCorrection = input('Apply correction profile to pixel values? (y/n): ','s'); disp(' ')
            if strcmp(applyCorrection,'y') == true || strcmp(applyCorrection,'n') == true
                correctionDecision = 'y';
            else
                disp('Invalid input. Must be ''y'' or ''n'''); disp(' ')
            end
        end
        sgtitle([animalID ' ' strDay ' pixel correction applied: ' correctionDecision])
    end
    
    %% Save the file to directory.
    [pathstr, ~, ~] = fileparts(cd);
    if strcmp(imagingType,'bilateral') == true
        dirpath = [pathstr '/Combined Imaging/Figures/Pixel Drift Correction/'];
    elseif strcmp(imagingType,'single') == true
        dirpath = [pathstr '/Single Hemisphere/Figures/Pixel Drift Correction/'];
    end
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(fixPixels, [dirpath animalID '_' strDay '_PixelDriftCorrection']);
    close(fixPixels)
    
    %% apply corrected data to each file from reshaped matrix
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
