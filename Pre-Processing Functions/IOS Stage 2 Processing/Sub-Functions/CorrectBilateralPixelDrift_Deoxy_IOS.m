function [] = CorrectBilateralPixelDrift_Deoxy_IOS(procDataFileIDs)
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
    catLH_Deoxydata = [];
    catRH_Deoxydata = [];
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
    % load the processed Deoxy/cement data from each file and concat it into one array
    for d = 1:length(indDayProcDataFileList)
        indDayProcDataFile = indDayProcDataFileList{d,1};
        load(indDayProcDataFile)
        samplingRate = ProcData.notes.CBVCamSamplingRate;
        trialDuration = ProcData.notes.trialDuration_sec;
        LH_Deoxydata = ProcData.data.Deoxy.LH;
        RH_Deoxydata = ProcData.data.Deoxy.RH;
        Cement_cementData = ProcData.data.Deoxy.Cement;
        catLH_Deoxydata = horzcat(catLH_Deoxydata,LH_Deoxydata); %#ok<AGROW>
        catRH_Deoxydata = horzcat(catRH_Deoxydata,RH_Deoxydata); %#ok<AGROW>
        catCement_cementData = horzcat(catCement_cementData,Cement_cementData); %#ok<AGROW>
    end
    % establish whether a slow exponential trend exists for the data
    [B,A] = butter(3,0.01/(samplingRate/2),'low');
    filtCatCement_cementData = filtfilt(B,A,catCement_cementData);
    x = ((1:length(filtCatCement_cementData))/samplingRate)';
    % create a weight vector for the trend
    Cement_weightVec = ones(1,length(x));
    Cement_secondHalfMean = mean(filtCatCement_cementData(floor(length(filtCatCement_cementData/2)):end));
    for t = 1:length(Cement_weightVec)
        if filtCatCement_cementData(t) > Cement_secondHalfMean
            Cement_weightVec(t) = 10;
        end
    end
    % compare weighted models
    Cement_modelFit = fit(x,filtCatCement_cementData','exp2','Weight',Cement_weightVec);
    Cement_modelFit_Y = Cement_modelFit(x);
    Cement_modelFit_norm = (Cement_modelFit_Y - min(Cement_modelFit_Y))./min(Cement_modelFit_Y);
    Cement_modelFit_flip = 1 - Cement_modelFit_norm;
    % apply exponential correction to original data
    LH_adjCat_Deoxydata = catLH_Deoxydata.*Cement_modelFit_flip';
    LH_rsAdjCat_Deoxydata = reshape(LH_adjCat_Deoxydata,[samplingRate*trialDuration,length(indDayProcDataFileList)]);
    RH_adjCat_Deoxydata = catRH_Deoxydata.*Cement_modelFit_flip';
    RH_rsAdjCat_Deoxydata = reshape(RH_adjCat_Deoxydata,[samplingRate*trialDuration,length(indDayProcDataFileList)]);
    % comparison showing original LH data and the corrected data
    fixPixels = figure;
    subplot(1,3,1)
    plot(x,filtCatCement_cementData,'k')
    hold on
    plot(x,Cement_modelFit_Y,'g')
    title('Cement ROI drift')
    xlabel('Time (sec)')
    ylabel('12-bit pixel val')
    legend('cement reflectance','exp fit')
    axis tight
    axis square
    subplot(1,3,2)
    p1 = plot(x,catLH_Deoxydata,'k');
    hold on
    p2 = plot(x,LH_adjCat_Deoxydata,'g');
    title('LH ROI Deoxy Reflectance')
    xlabel('Time (sec)')
    ylabel('12-bit pixel val')
    legend([p1,p2],'original data','proposed corrected')
    axis tight
    axis square
    subplot(1,3,3)
    p1 = plot(x,catRH_Deoxydata,'k');
    hold on
    p2 = plot(x,RH_adjCat_Deoxydata,'g');
    title('RH ROI Deoxy Reflectance')
    xlabel('Time (sec)')
    ylabel('12-bit pixel val')
    legend([p1,p2],'original data','proposed corrected')
    axis tight
    axis square
    % determine which correction profile to use for RH data
    correctionDecision = 'n';
    while strcmp(correctionDecision,'n') == true
        applyCorrection = input(['Apply correction profile to ' strDay ' pixel values? (y/n): '],'s'); disp(' ')
        if strcmp(applyCorrection,'y') == true || strcmp(applyCorrection,'n') == true 
            correctionDecision = 'y';
        else
            disp('Invalid input'); disp(' ')
        end
    end
    sgtitle({[animalID ' ' strDay ' green frames'];['pixel correction applied: ' applyCorrection]})
    savefig(fixPixels,[animalID '_' strDay '_Deoxy_PixelDriftCorrection']);
    close(fixPixels)
    % apply corrected data to each file from reshaped matrix
    for d = 1:length(indDayProcDataFileList)
        indDayProcDataFile = indDayProcDataFileList{d,1};
        load(indDayProcDataFile)
        % pixel correction
        if strcmp(applyCorrection,'n') == true
            ProcData.data.Deoxy.adjLH = ProcData.data.Deoxy.LH;
            ProcData.data.Deoxy.adjRH = ProcData.data.Deoxy.RH;
        elseif strcmp(applyCorrection,'y') == true
            ProcData.data.Deoxy.adjLH = LH_rsAdjCat_Deoxydata(:,d)';
            ProcData.data.Deoxy.adjRH = RH_rsAdjCat_Deoxydata(:,d)';
        end
        disp(['Saving pixel corrections to ' strDay ' ProcData file ' num2str(d) ' of ' num2str(length(indDayProcDataFileList))]); disp(' ')
        save(indDayProcDataFile,'ProcData')
    end
end

end
