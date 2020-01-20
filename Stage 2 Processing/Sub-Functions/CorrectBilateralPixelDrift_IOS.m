function [] = CorrectBilateralPixelDrift_IOS(procDataFileIDs)
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
    catCement_cementData = [];
    fileName = firstsFileOfDay{1,b};
    [~,fileDate,~] = GetFileInfo_IOS(fileName);
    strDay = ConvertDate_IOS(fileDate);
    p = 1;
    for c = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(c,:);
        if strfind(procDataFileID, fileDate) >= 1
            indDayProcDataFileList{p,1} = procDataFileID; %#ok<AGROW>
            p = p + 1;
        end
    end
    % load the processed CBV/cement data from each file and concat it into one array
    for d = 1:length(indDayProcDataFileList)
        indDayProcDataFile = indDayProcDataFileList{d,1};
        load(indDayProcDataFile)
        samplingRate = ProcData.notes.CBVCamSamplingRate;
        trialDuration = ProcData.notes.trialDuration_sec;
        LH_CBVdata = ProcData.data.CBV.LH;
        RH_CBVdata = ProcData.data.CBV.RH;
        LH_cementData = ProcData.data.CBV.LH_Cement;
        RH_cementData = ProcData.data.CBV.RH_Cement;
        Cement_cementData = ProcData.data.CBV.Cement;
        catLH_CBVdata = horzcat(catLH_CBVdata,LH_CBVdata); %#ok<AGROW>
        catRH_CBVdata = horzcat(catRH_CBVdata,RH_CBVdata); %#ok<AGROW>
        catLH_cementData = horzcat(catLH_cementData,LH_cementData); %#ok<AGROW>
        catRH_cementData = horzcat(catRH_cementData,RH_cementData); %#ok<AGROW>
        catCement_cementData = horzcat(catCement_cementData,Cement_cementData); %#ok<AGROW>
    end
    % establish whether a slow exponential trend exists for the data
    [B,A] = butter(3,0.01/(samplingRate/2),'low');
    filtCatLH_cementData = filtfilt(B,A,catLH_cementData);
    filtCatRH_cementData = filtfilt(B,A,catRH_cementData);
    filtCatCement_cementData = filtfilt(B,A,catCement_cementData);
    x = ((1:length(filtCatLH_cementData))/samplingRate)';
    % create a weight vector for the trend
    LH_weightVec = ones(1,length(x));
    RH_weightVec = ones(1,length(x));
    Cement_weightVec = ones(1,length(x));
    LH_secondHalfMean = mean(filtCatLH_cementData(floor(length(filtCatLH_cementData/2)):end));
    RH_secondHalfMean = mean(filtCatRH_cementData(floor(length(filtCatRH_cementData/2)):end));
    Cement_secondHalfMean = mean(filtCatCement_cementData(floor(length(filtCatCement_cementData/2)):end));
    for t = 1:length(LH_weightVec)
        if filtCatLH_cementData(t) > LH_secondHalfMean
            LH_weightVec(t) = 10;
        end
        if filtCatRH_cementData(t) > RH_secondHalfMean
            RH_weightVec(t) = 10;
        end
        if filtCatCement_cementData(t) > Cement_secondHalfMean
            Cement_weightVec(t) = 10;
        end
    end
    % compare weighted models
    LH_modelFit = fit(x,filtCatLH_cementData','exp2','Weight',LH_weightVec);
    RH_modelFit = fit(x,filtCatRH_cementData','exp2','Weight',RH_weightVec);
    Cement_modelFit = fit(x,filtCatCement_cementData','exp2','Weight',Cement_weightVec);
    LH_modelFit_Y = LH_modelFit(x);
    RH_modelFit_Y = RH_modelFit(x);
    Cement_modelFit_Y = Cement_modelFit(x);
    LH_modelFit_norm = (LH_modelFit_Y - min(LH_modelFit_Y))./min(LH_modelFit_Y);
    RH_modelFit_norm = (RH_modelFit_Y - min(RH_modelFit_Y))./min(RH_modelFit_Y);
    Cement_modelFit_norm = (Cement_modelFit_Y - min(Cement_modelFit_Y))./min(Cement_modelFit_Y);
    LH_modelFit_flip = 1 - LH_modelFit_norm;
    RH_modelFit_flip = 1 - RH_modelFit_norm;
    Cement_modelFit_flip = 1 - Cement_modelFit_norm;
    % apply exponential correction to original data
    LH_adjCatA_CBVdata = catLH_CBVdata.*LH_modelFit_flip';
    LH_rsAdjCatA_CBVdata = reshape(LH_adjCatA_CBVdata',[samplingRate*trialDuration,length(indDayProcDataFileList)]);
    LH_adjCatB_CBVdata = catLH_CBVdata.*RH_modelFit_flip';
    LH_rsAdjCatB_CBVdata = reshape(LH_adjCatB_CBVdata,[samplingRate*trialDuration,length(indDayProcDataFileList)]);
    LH_adjCatC_CBVdata = catLH_CBVdata.*Cement_modelFit_flip';
    LH_rsAdjCatC_CBVdata = reshape(LH_adjCatC_CBVdata,[samplingRate*trialDuration,length(indDayProcDataFileList)]);
    RH_adjCatA_CBVdata = catRH_CBVdata.*LH_modelFit_flip';
    RH_rsAdjCatA_CBVdata = reshape(RH_adjCatA_CBVdata,[samplingRate*trialDuration,length(indDayProcDataFileList)]);
    RH_adjCatB_CBVdata = catRH_CBVdata.*RH_modelFit_flip';
    RH_rsAdjCatB_CBVdata = reshape(RH_adjCatB_CBVdata,[samplingRate*trialDuration,length(indDayProcDataFileList)]);
    RH_adjCatC_CBVdata = catRH_CBVdata.*Cement_modelFit_flip';
    RH_rsAdjCatC_CBVdata = reshape(RH_adjCatC_CBVdata,[samplingRate*trialDuration,length(indDayProcDataFileList)]);
    % evaluate fits
    LH_Atest = corrcoef(catLH_CBVdata,LH_modelFit_Y);
    LH_Atest_R = LH_Atest(1,2);
    LH_Btest = corrcoef(catLH_CBVdata,RH_modelFit_Y);
    LH_Btest_R = LH_Btest(1,2);
    LH_Ctest = corrcoef(catLH_CBVdata,Cement_modelFit_Y);
    LH_Ctest_R = LH_Ctest(1,2);
    RH_Atest = corrcoef(catRH_CBVdata,LH_modelFit_Y);
    RH_Atest_R = RH_Atest(1,2);
    RH_Btest = corrcoef(catRH_CBVdata,RH_modelFit_Y);
    RH_Btest_R = RH_Btest(1,2);
    RH_Ctest = corrcoef(catRH_CBVdata,Cement_modelFit_Y);
    RH_Ctest_R = RH_Ctest(1,2);
    % comparison showing original LH data and the corrected data
    fixPixels = figure;
    subplot(4,3,1)
    plot(x,filtCatLH_cementData,'k')
    hold on
    plot(x,LH_modelFit_Y,'r')
    title('Cement drift fit A')
    xlabel('Time (sec)')
    ylabel('12-bit pixel val')
    legend('cement data','exp fit')
    axis tight
    
    subplot(4,3,2)
    plot(x,filtCatRH_cementData,'k')
    hold on
    plot(x,RH_modelFit_Y,'b')
    title('Cement drift fit B')
    xlabel('Time (sec)')
    ylabel('12-bit pixel val')
    legend('cement data','exp fit')
    axis tight
    
    subplot(4,3,3)
    plot(x,filtCatCement_cementData,'k')
    hold on
    plot(x,Cement_modelFit_Y,'g')
    title('Cement drift fit C')
    xlabel('Time (sec)')
    ylabel('12-bit pixel val')
    legend('cement data','exp fit')
    axis tight
    
    subplot(4,3,4)
    plot(x,catLH_CBVdata,'k')
    title('LH original data')
    xlabel('Time (sec)')
    ylabel('12-bit pixel val')
    axis tight
    
    subplot(4,3,5)
    plot(x,catRH_CBVdata,'k')
    title('RH original data')
    xlabel('Time (sec)')
    ylabel('12-bit pixel val')
    axis tight
    
    subplot(4,3,6)
    plot(x,LH_modelFit_flip,'r')
    hold on
    plot(x,RH_modelFit_flip,'b')
    plot(x,Cement_modelFit_flip,'g')
    title('Correction profile')
    xlabel('Time (sec)')
    ylabel('Normalized val')
    legend('ROI A','ROI B','ROI C')
    axis tight
    
    subplot(4,3,7)
    plot(x,catLH_CBVdata,'k')
    hold on
    plot(x,LH_adjCatA_CBVdata,'r')
    title({'A corrected data';['CorrCoef: ' num2str(LH_Atest_R)]})
    xlabel('Time (sec)')
    ylabel('12-bit pixel val')
    axis tight
    
    subplot(4,3,8)
    plot(x,catLH_CBVdata,'k')
    hold on
    plot(x,LH_adjCatB_CBVdata,'b')
    title({'B corrected data';['CorrCoef: ' num2str(LH_Btest_R)]})
    xlabel('Time (sec)')
    ylabel('12-bit pixel val')
    axis tight
    
    subplot(4,3,9)
    plot(x,catLH_CBVdata,'k')
    hold on
    plot(x,LH_adjCatC_CBVdata,'g')
    title({'C corrected data';['CorrCoef: ' num2str(LH_Ctest_R)]})
    xlabel('Time (sec)')
    ylabel('12-bit pixel val')
    axis tight
    
    subplot(4,3,10)
    plot(x,catRH_CBVdata,'k')
    hold on
    plot(x,RH_adjCatA_CBVdata,'r')
    title({'A corrected data';['CorrCoef: ' num2str(RH_Atest_R)]})
    xlabel('Time (sec)')
    ylabel('12-bit pixel val')
    axis tight
    
    subplot(4,3,11)
    plot(x,catRH_CBVdata,'k')
    hold on
    plot(x,RH_adjCatB_CBVdata,'b')
    title({'B corrected data';['CorrCoef: ' num2str(RH_Btest_R)]})
    xlabel('Time (sec)')
    ylabel('12-bit pixel val')
    axis tight
    
    subplot(4,3,12)
    plot(x,catRH_CBVdata,'k')
    hold on
    plot(x,RH_adjCatC_CBVdata,'g')
    title({'C corrected data';['CorrCoef: ' num2str(RH_Ctest_R)]})
    xlabel('Time (sec)')
    ylabel('12-bit pixel val')
    axis tight
    
    % determine which correction profile to use for RH data
    correctionDecision = 'n';
    while strcmp(correctionDecision,'n') == true
        applyCorrection = input(['Apply correction profile to ' strDay ' pixel values? (O/A/B/C): '],'s'); disp(' ')
        if strcmp(applyCorrection,'O') == true || strcmp(applyCorrection,'A') == true || strcmp(applyCorrection,'B') == true || strcmp(applyCorrection,'C') == true
            correctionDecision = 'y';
        else
            disp('Invalid input. Must be ''O'', ''A'', ''B'', ''C'''); disp(' ')
        end
    end
    sgtitle([animalID ' ' strDay ' pixel correction applied: ' applyCorrection])
    % Save the file to directory.
    [pathstr,~,~] = fileparts(cd);
    dirpath = [pathstr '/Bilateral Imaging/Pixel Drift Correction/'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(fixPixels, [dirpath animalID '_' strDay '_PixelDriftCorrection']);
    close(fixPixels)
    % apply corrected data to each file from reshaped matrix
    for d = 1:length(indDayProcDataFileList)
        indDayProcDataFile = indDayProcDataFileList{d,1};
        load(indDayProcDataFile)
        % LH pixel correction
        if strcmp(applyCorrection,'O') == true
            ProcData.data.CBV.adjLH = ProcData.data.CBV.LH;
            ProcData.data.CBV.adjRH = ProcData.data.CBV.RH;
        elseif strcmp(applyCorrection,'A') == true
            ProcData.data.CBV.adjLH = LH_rsAdjCatA_CBVdata(:,d)';
            ProcData.data.CBV.adjRH = RH_rsAdjCatA_CBVdata(:,d)';
        elseif strcmp(applyCorrection,'B') == true
            ProcData.data.CBV.adjLH = LH_rsAdjCatB_CBVdata(:,d)';
            ProcData.data.CBV.adjRH = RH_rsAdjCatB_CBVdata(:,d)';
        elseif strcmp(applyCorrection,'C') == true
            ProcData.data.CBV.adjLH = LH_rsAdjCatC_CBVdata(:,d)';
            ProcData.data.CBV.adjRH = RH_rsAdjCatC_CBVdata(:,d)';
        end
        disp(['Saving pixel corrections to ' strDay ' ProcData file ' num2str(d) ' of ' num2str(length(indDayProcDataFileList))]); disp(' ')
        save(indDayProcDataFile,'ProcData')
    end
end

end

