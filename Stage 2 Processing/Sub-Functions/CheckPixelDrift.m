function [] = CheckPixelDrift(procDataFileIDs)
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
        if strfind(procDataFileID, fileDate) == 6 
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
    
    LH_p = polyfit((1:length(catLH_Data))/fs,catLH_Data,3);
    RH_p = polyfit((1:length(catLH_Data))/fs,catRH_Data,3);
    
    LH_y = polyval(LH_p,(1:length(catLH_Data))/fs,LH_p);
    RH_y = polyval(RH_p,(1:length(catRH_Data))/fs,RH_p);
    
    figure;
    plot((1:length(catLH_Data))/fs, catLH_Data, 'r')
    hold on
    plot((1:length(catRH_Data))/fs, catRH_Data, 'b')
    plot((1:length(catCement_Data))/fs, catCement_Data, 'k')
    plot((1:length(catLH_Data))/fs, LH_y, 'k');
    plot((1:length(catRH_Data))/fs, RH_y, 'k');

    title('Pixel Drift')
    xlabel('Time (sec)')
    ylabel('Mean pixel val')
    legend('LH', 'RH', 'Cement', 'LH polyfit', 'RH polyfit')
    
    catLH_Data = catLH_Data - LH_y;
    catRH_Data = catRH_Data - RH_y;
    
end
