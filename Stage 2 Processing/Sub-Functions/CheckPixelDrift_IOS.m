function [] = CheckPixelDrift_IOS(procDataFileIDs,imagingType)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: 
%________________________________________________________________________________________________________________________
%
%   Inputs: 
%
%   Outputs: 
%
%   Last Revised: September 11th, 2019
%________________________________________________________________________________________________________________________

% Check how many different days are in the procDataFileID list
[~,fileDates,~] = GetFileInfo_IOS(procDataFileIDs);
[uniqueDays, ~, DayID] = GetUniqueDays_IOS(fileDates);
firstsFileOfDay = cell(1,length(uniqueDays));
for a = 1:length(uniqueDays)
    FileInd = DayID == a;
    dayFilenames = procDataFileIDs(FileInd,:);
    firstsFileOfDay(a) = {dayFilenames(1,:)};
end

for b = 1:length(firstsFileOfDay)
    indDayProcDataFileList = {};
    catLH_Data = [];
    catRH_Data = [];
    catCBV_Data = [];
    catLH_CementData = [];
    catRH_CementData = [];
    catCementData = [];
    fs = 30;
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
    
    for d = 1:length(indDayProcDataFileList)
        indDayProcDataFile = indDayProcDataFileList{d,1};
        load(indDayProcDataFile)
        if strcmp(imagingType,'bilateral') == true
            LH_Data = ProcData.data.CBV.LH;
            RH_Data = ProcData.data.CBV.RH;
            LH_cementData = ProcData.data.CBV.LH_Cement;
            RH_cementData = ProcData.data.CBV.RH_Cement;
            catLH_Data = horzcat(catLH_Data,LH_Data);
            catRH_Data = horzcat(catRH_Data,RH_Data);
            catLH_CementData = horzcat(catLH_CementData,LH_cementData);
            catRH_CementData = horzcat(catRH_CementData,RH_cementData);
        elseif strcmp(imagingType,'single') == true
            CBV_Data = ProcData.data.CBV.Barrels;
            cementData = ProcData.data.CBV.Cement;
            catCBV_Data = horzcat(catCBV_Data,CBV_Data);
            catCementData = horzcat(catCementData,cementData);
        end
    end
    
    
    [B, A] = butter(4,0.1/(30/2),'low');
    filtCatLH_CementData = filtfilt(B,A,catLH_CementData);
    filtCatRH_CementData = filtfilt(B,A,catRH_CementData);
    if strcmp(imagingType,'bilateral') == true
        correctedCatLH_Data = (catLH_Data - filtCatLH_CementData);
        correctedCatRH_Data = (catRH_Data - filtCatRH_CementData);
        LH_DC_reset = mean(correctedCatLH_Data);
        RH_DC_reset = mean(correctedCatRH_Data);
        LH_1kDiff = 1000 - LH_DC_reset;
        RH_1kDiff = 1000 - RH_DC_reset;
        LH_DC_shift = ones(1,length(correctedCatLH_Data))*LH_1kDiff;
        RH_DC_shift = ones(1,length(correctedCatRH_Data))*RH_1kDiff;
        shiftedCorrectedCatLH_Data = correctedCatLH_Data + LH_DC_shift;
        shiftedCorrectedCatRH_Data = correctedCatRH_Data + RH_DC_shift;
        
        figure;
        plot((1:length(catLH_Data))/fs,catLH_Data,'r')
        hold on
        plot((1:length(catRH_Data))/fs,catRH_Data,'b')
        plot((1:length(filtCatLH_CementData))/fs,filtCatLH_CementData,'k')
        plot((1:length(filtCatRH_CementData))/fs,filtCatRH_CementData,'g')
        plot((1:length(shiftedCorrectedCatLH_Data))/fs,shiftedCorrectedCatLH_Data,'m')
        plot((1:length(shiftedCorrectedCatRH_Data))/fs,shiftedCorrectedCatRH_Data,'c')
        title('Pixel Drift')
        xlabel('Time (sec)')
        ylabel('Mean pixel val')
        legend('Original LH','Original RH','LH Cement','RH Cement','Corrected LH','Corrected RH')
    elseif strcmp(imagingType,'single') == true
        correctedCatCBV_Data = (catCBV_Data - filtCatLH_CementData);
        CBV_DC_reset = mean(correctedCatCBV_Data);
        CBV_1kDiff = 1000 - CBV_DC_reset;
        CBV_DC_shift = ones(1,length(correctedCatCBV_Data))*CBV_1kDiff;
        shiftedCorrectedCatCBV_Data = correctedCatCBV_Data + CBV_DC_shift;
        
        figure;
        plot((1:length(catCBV_Data))/fs, catCBV_Data, 'r')
        hold on
        plot((1:length(filtCatLH_CementData))/fs,filtCatLH_CementData, 'g')
        plot((1:length(filtCatLH_CementData))/fs,shiftedCorrectedCatCBV_Data, 'm')
        title('Pixel Drift')
        xlabel('Time (sec)')
        ylabel('Mean pixel val')
        legend('Barrels','Cement','Barrels corrected')
        
    end
end

end
