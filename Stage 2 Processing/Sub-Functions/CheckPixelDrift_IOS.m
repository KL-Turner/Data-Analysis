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
    catCement_cementData = [];
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
            Cement_cementData = ProcData.data.CBV.Cement;
            catLH_CBVdata = horzcat(catLH_CBVdata,LH_CBVdata);
            catRH_CBVdata = horzcat(catRH_CBVdata,RH_CBVdata);
            catLH_cementData = horzcat(catLH_cementData,LH_cementData);
            catRH_cementData = horzcat(catRH_cementData,RH_cementData);
            catCement_cementData = horzcat(catCement_cementData,Cement_cementData);
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
        LH_rsAdjCatA_CBVdata = reshape(LH_adjCatA_CBVdata,[length(indDayProcDataFileList) samplingRate*trialDuration]);
        LH_adjCatB_CBVdata = catLH_CBVdata.*RH_modelFit_flip';
        LH_rsAdjCatB_CBVdata = reshape(LH_adjCatB_CBVdata,[length(indDayProcDataFileList) samplingRate*trialDuration]);
        LH_adjCatC_CBVdata = catLH_CBVdata.*Cement_modelFit_flip';
        LH_rsAdjCatC_CBVdata = reshape(LH_adjCatC_CBVdata,[length(indDayProcDataFileList) samplingRate*trialDuration]);
        
        RH_adjCatA_CBVdata = catRH_CBVdata.*LH_modelFit_flip';
        RH_rsAdjCatA_CBVdata = reshape(RH_adjCatA_CBVdata,[length(indDayProcDataFileList) samplingRate*trialDuration]);
        RH_adjCatB_CBVdata = catRH_CBVdata.*RH_modelFit_flip';
        RH_rsAdjCatB_CBVdata = reshape(RH_adjCatB_CBVdata,[length(indDayProcDataFileList) samplingRate*trialDuration]);
        RH_adjCatC_CBVdata = catRH_CBVdata.*Cement_modelFit_flip';
        RH_rsAdjCatC_CBVdata = reshape(RH_adjCatC_CBVdata,[length(indDayProcDataFileList) samplingRate*trialDuration]);
        
        %% comparison showing original LH data and the corrected data
        LH_fixPixels = figure;
        subplot(2,4,1)
        plot(x,catLH_CBVdata,'k')
        title('LH original data')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        axis tight
        
        subplot(2,4,2)
        plot(x,filtCatLH_cementData,'k')
        hold on
        plot(x,LH_modelFit_Y,'r')
        title('Cement drift fit A')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        legend('cement data','exp fit')
        axis tight
        
        subplot(2,4,3)
        plot(x,filtCatRH_cementData,'k')
        hold on
        plot(x,RH_modelFit_Y,'b')
        title('Cement drift fit B')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        legend('cement data','exp fit')
        axis tight
        
        subplot(2,4,4)
        plot(x,filtCatCement_cementData,'k')
        hold on
        plot(x,Cement_modelFit_Y,'g')
        title('Cement drift fit C')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        legend('cement data','exp fit')
        axis tight
        
        subplot(2,4,5)
        plot(x,LH_modelFit_flip,'r')
        hold on
        plot(x,RH_modelFit_flip,'b')
        plot(x,Cement_modelFit_flip,'g')
        title('Correction profile')
        xlabel('Time (sec)')
        ylabel('Normalized val')
        legend('ROI A','ROI B','ROI C')
        axis tight
        
        subplot(2,4,6)
        plot(x,catLH_CBVdata,'k')
        hold on
        plot(x,LH_adjCatA_CBVdata,'r')
        title('A corrected data')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        legend('original data','corrected data')
        axis tight
        
        subplot(2,4,7)
        plot(x,catLH_CBVdata,'k')
        hold on
        plot(x,LH_adjCatB_CBVdata,'b')
        title('B corrected data')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        legend('original data','corrected data')
        axis tight
        
        subplot(2,4,8)
        plot(x,catLH_CBVdata,'k')
        hold on
        plot(x,LH_adjCatC_CBVdata,'g')
        title('C corrected data')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        legend('original data','corrected data')
        axis tight
        
        %%
        correctionDecision = 'n';
        while strcmp(correctionDecision,'n') == true
            LH_applyCorrection = input('Apply correction profile to pixel values? (O/A/B/C): ','s'); disp(' ')
            if strcmp(LH_applyCorrection,'O') == true || strcmp(LH_applyCorrection,'A') == true || strcmp(LH_applyCorrection,'B') == true || strcmp(LH_applyCorrection,'C') == true
                correctionDecision = 'y';
            else
                disp('Invalid input. Must be ''O'', ''A'', ''B'', ''C'''); disp(' ')
            end
        end
        sgtitle([animalID ' ' strDay ' pixel correction applied: ' LH_applyCorrection])
        
        %% comparison showing original LH data and the corrected data
        RH_fixPixels = figure;
        subplot(2,4,1)
        plot(x,catRH_CBVdata,'k')
        title('RH original data')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        axis tight
        
        subplot(2,4,2)
        plot(x,filtCatLH_cementData,'k')
        hold on
        plot(x,LH_modelFit_Y,'r')
        title('Cement drift fit A')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        legend('cement data','exp fit')
        axis tight
        
        subplot(2,4,3)
        plot(x,filtCatRH_cementData,'k')
        hold on
        plot(x,RH_modelFit_Y,'b')
        title('Cement drift fit B')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        legend('cement data','exp fit')
        axis tight
        
        subplot(2,4,4)
        plot(x,filtCatCement_cementData,'k')
        hold on
        plot(x,Cement_modelFit_Y,'g')
        title('Cement drift fit C')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        legend('cement data','exp fit')
        axis tight
        
        subplot(2,4,5)
        plot(x,LH_modelFit_flip,'r')
        hold on
        plot(x,RH_modelFit_flip,'b')
        plot(x,Cement_modelFit_flip,'g')
        title('Correction profile')
        xlabel('Time (sec)')
        ylabel('Normalized val')
        legend('ROI A','ROI B','ROI C')
        axis tight
        
        subplot(2,4,6)
        plot(x,catRH_CBVdata,'k')
        hold on
        plot(x,RH_adjCatA_CBVdata,'r')
        title('A corrected data')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        legend('original data','corrected data')
        axis tight
        
        subplot(2,4,7)
        plot(x,catRH_CBVdata,'k')
        hold on
        plot(x,RH_adjCatB_CBVdata,'b')
        title('B corrected data')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        legend('original data','corrected data')
        axis tight
        
        subplot(2,4,8)
        plot(x,catRH_CBVdata,'k')
        hold on
        plot(x,RH_adjCatC_CBVdata,'g')
        title('C corrected data')
        xlabel('Time (sec)')
        ylabel('12-bit pixel val')
        legend('original data','corrected data')
        axis tight
        
        %%
        correctionDecision = 'n';
        while strcmp(correctionDecision,'n') == true
            RH_applyCorrection = input('Apply correction profile to pixel values? (O/A/B/C): ','s'); disp(' ')
            if strcmp(RH_applyCorrection,'O') == true || strcmp(RH_applyCorrection,'A') == true || strcmp(RH_applyCorrection,'B') == true || strcmp(RH_applyCorrection,'C') == true
                correctionDecision = 'y';
            else
                disp('Invalid input. Must be ''O'', ''A'', ''B'', ''C'''); disp(' ')
            end
        end
        sgtitle([animalID ' ' strDay ' pixel correction applied: ' RH_applyCorrection])
        
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
        savefig(LH_fixPixels, [dirpath animalID '_' strDay '_LH_PixelDriftCorrection']);
        close(LH_fixPixels)
        savefig(RH_fixPixels, [dirpath animalID '_' strDay '_RH_PixelDriftCorrection']);
        close(RH_fixPixels)
        
        %% apply corrected data to each file from reshaped matrix
        for d = 1:length(indDayProcDataFileList)
            indDayProcDataFile = indDayProcDataFileList{d,1};
            load(indDayProcDataFile)
            if strcmp(LH_applyCorrection,'O') == true
                ProcData.data.CBV.correctedLH = ProcData.data.CBV.LH;
            elseif strcmp(LH_applyCorrection,'A') == true
                ProcData.data.CBV.correctedLH = LH_rsAdjCatA_CBVdata(:,d);
            elseif strcmp(LH_applyCorrection,'B') == true
                ProcData.data.CBV.correctedLH = LH_rsAdjCatB_CBVdata(:,d);
            elseif strcmp(LH_applyCorrection,'C') == true
                ProcData.data.CBV.correctedLH = LH_rsAdjCatC_CBVdata(:,d);
            end
            if strcmp(RH_applyCorrection,'O') == true
                ProcData.data.CBV.correctedRH = ProcData.data.CBV.RH;
            elseif strcmp(RH_applyCorrection,'A') == true
                ProcData.data.CBV.correctedRH = RH_rsAdjCatA_CBVdata(:,d);
            elseif strcmp(RH_applyCorrection,'B') == true
                ProcData.data.CBV.correctedRH = RH_rsAdjCatB_CBVdata(:,d);
            elseif strcmp(RH_applyCorrection,'C') == true
                ProcData.data.CBV.correctedRH = RH_rsAdjCatC_CBVdata(:,d);
            end
            disp(['Saving pixel corrections to ProcData file' num2str(d) ' of ' num2str(length(indDayProcDataFileList))]); disp(' ')
            save(indDayProcDataFile,'ProcData')
        end
    end
end

%     elseif strcmp(imagingType,'single') == true
        %         filtCatCementData = filtfilt(B,A,catCementData);
        %         x = ((1:length(filtCatCementData))/samplingRate)';
        %         % create a weight vector for the trend
        %         weightVec = ones(1,length(x));
        %         secondHalfMean = mean(filtCatCementData(floor(length(filtCatCementData/2)):end));
        %         for t = 1:length(weightVec)
        %             if filtCatCementData(t) > secondHalfMean
        %                 weightVec(t) = 10;
        %             end
        %         end
        %         % compare weighted vs. unweighted fits
        %         modelFit1 = fit(x,filtCatCementData','exp2');
        %         modelFit1_Y = modelFit1(x);
        %         modelFit2 = fit(x,filtCatCementData','exp2','Weight',weightVec);
        %         modelFit2_Y = modelFit2(x);
        %         modelFit2_norm = (modelFit2_Y - min(modelFit2_Y))./min(modelFit2_Y);
        %         modelFit2_flip = 1 - modelFit2_norm;
        %         adjCatBarrels_CBVdata = catBarrels_CBVdata.*modelFit2_flip';
        %         rsAdjCatBarrels_CBVdata = reshape(adjCatBarrels_CBVdata,[length(indDayProcDataFileList) samplingRate*trialDuration]);
        %         % comparison showing original LH data and the corrected data
        %         LH_fixPixels = figure;
        %         subplot(2,2,1)
        %         plot(x,catBarrels_CBVdata,'k')
        %         title('Barrels original data')
        %         xlabel('Time (sec)')
        %         ylabel('12-bit pixel val')
        %         axis tight
        %         subplot(2,2,2)
        %         plot(x,filtCatCementData,'k')
        %         hold on
        %         plot(x,modelFit1_Y,'b')
        %         plot(x,modelFit2_Y,'r')
        %         title('cement drift')
        %         xlabel('Time (sec)')
        %         ylabel('12-bit pixel val')
        %         legend('cement data','exp fit','weighted fit')
        %         axis tight
        %         subplot(2,2,3)
        %         plot(x,modelFit2_flip,'r')
        %         title('correction profile')
        %         xlabel('Time (sec)')
        %         ylabel('Normalized val')
        %         axis tight
        %         subplot(2,2,4)
        %         plot(x,catBarrels_CBVdata,'k')
        %         hold on
        %         plot(x,adjCatBarrels_CBVdata,'r')
        %         title('Corrected data')
        %         xlabel('Time (sec)')
        %         ylabel('12-bit pixel val')
        %         legend('original data','corrected data')
        %         axis tight
        %         correctionDecision = 'n';
        %         while strcmp(correctionDecision,'n') == true
        %             RH_applyCorrection = input('Apply correction profile to pixel values? (y/n): ','s'); disp(' ')
        %             if strcmp(RH_applyCorrection,'y') == true || strcmp(RH_applyCorrection,'n') == true
        %                 correctionDecision = 'y';
        %             else
        %                 disp('Invalid input. Must be ''y'' or ''n'''); disp(' ')
        %             end
        %         end
        %         sgtitle([animalID ' ' strDay ' pixel correction applied: ' correctionDecision])
        %     end
        %
        %% Save the file to directory.
%         [pathstr, ~, ~] = fileparts(cd);
%         if strcmp(imagingType,'bilateral') == true
%             dirpath = [pathstr '/Combined Imaging/Figures/Pixel Drift Correction/'];
%         elseif strcmp(imagingType,'single') == true
%             dirpath = [pathstr '/Single Hemisphere/Figures/Pixel Drift Correction/'];
%         end
%         if ~exist(dirpath,'dir')
%             mkdir(dirpath);
%         end
%         savefig(LH_fixPixels, [dirpath animalID '_' strDay '_LH_PixelDriftCorrection']);
%         close(LH_fixPixels)
%         savefig(RH_fixPixels, [dirpath animalID '_' strDay '_RH_PixelDriftCorrection']);
%         close(RH_fixPixels)
%         
%         %% apply corrected data to each file from reshaped matrix
%             for d = 1:length(indDayProcDataFileList)
%                 indDayProcDataFile = indDayProcDataFileList{d,1};
%                 load(indDayProcDataFile)
%                 if strcmp(correctionDecision,'O') == true
%                     
%                 elseif strcmp(correctionDecision,'A') == true
%                     
%                 elseif strcmp(correctionDecision,'B') == true
%                     
%                 elseif strcmp(correctionDecision,'C') == true
%                     
%                 end
%             end
%     end
% end

end
