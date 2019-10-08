function [] = CheckPixelDrift_IOS(procDataFileIDs)
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
[~, fileDates, ~] = GetFileInfo_IOS(procDataFileIDs);
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
    catCement_Data = [];
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
        LH_Data = ProcData.data.CBV.LH;
        RH_Data = ProcData.data.CBV.RH;
        cement_Data = ProcData.data.CBV.Cement; 
        catLH_Data = horzcat(catLH_Data, LH_Data);
        catRH_Data = horzcat(catRH_Data, RH_Data);
        catCement_Data = horzcat(catCement_Data, cement_Data);
    end
    
      
    [B, A] = butter(4, 0.1/(30/2), 'low');
    filtCatCement_Data = filtfilt(B, A, catCement_Data);    % Filtered heart rate signal
    correctedCatLH_Data = (catLH_Data - filtCatCement_Data);
    correctedCatRH_Data = (catRH_Data - filtCatCement_Data);
    LH_DC_reset = mean(correctedCatLH_Data);
    RH_DC_reset = mean(correctedCatRH_Data);
    LH_1kDiff = 1000 - LH_DC_reset;
    RH_1kDiff = 1000 - RH_DC_reset;
    LH_DC_shift = ones(1,length(correctedCatLH_Data))*LH_1kDiff;
    RH_DC_shift = ones(1,length(correctedCatRH_Data))*RH_1kDiff;
    shiftedCorrectedCatLH_Data = correctedCatLH_Data + LH_DC_shift;
    shiftedCorrectedCatRH_Data = correctedCatRH_Data + RH_DC_shift;

    figure;
    plot((1:length(catLH_Data))/fs, catLH_Data, 'r')
    hold on
    plot((1:length(catRH_Data))/fs, catRH_Data, 'b')
    plot((1:length(catCement_Data))/fs, catCement_Data, 'k')
    plot((1:length(catCement_Data))/fs, filtCatCement_Data, 'g')
    plot((1:length(catCement_Data))/fs, shiftedCorrectedCatLH_Data, 'm')
    plot((1:length(catCement_Data))/fs, shiftedCorrectedCatRH_Data, 'c')
    title('Pixel Drift')
    xlabel('Time (sec)')
    ylabel('Mean pixel val')
    legend('LH', 'RH', 'Cement', 'LH corrected', 'RH corrected')
    
end
