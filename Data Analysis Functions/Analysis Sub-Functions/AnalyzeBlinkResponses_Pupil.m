function [Results_BlinkResponses] = AnalyzeBlinkResponses_Pupil(animalID,rootFolder,delim,Results_BlinkResponses)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________

%% only run analysis for valid animal IDs
dataLocation = [rootFolder delim 'Data' delim animalID delim 'Bilateral Imaging'];
cd(dataLocation)
% find and load RestingBaselines.mat struct
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
% procdata file IDs
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
samplingRate = 30;    % lowpass filter
[z,p,k] = butter(4,10/(samplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
[z2,p2,k2] = butter(4,1/(samplingRate/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
specSamplingRate = 10;
data.pupilArea =[];
data.LH_HbT = [];
data.RH_HbT = [];
data.LH_cort = [];
data.RH_cort = [];
data.hip = [];
data.whisk = [];
data.EMG = [];
for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    [animalID,fileDate,fileID] = GetFileInfo_IOS(procDataFileID);
    strDay = ConvertDate_IOS(fileDate);
    load(procDataFileID)
    if strcmp(ProcData.data.Pupil.frameCheck,'y') == true
        specDataFileID = [animalID '_' fileID '_SpecDataB.mat'];
        load(specDataFileID)
        if isfield(ProcData.data.Pupil,'shiftedBlinks') == true
            blinks = ProcData.data.Pupil.shiftedBlinks;
        elseif isempty(ProcData.data.Pupil.blinkInds) == false
            blinks = ProcData.data.Pupil.blinkInds;
        else
            blinks = [];
        end
        % condense blinks
        if isempty(blinks) == false
            cc = 1;
            for bb = 1:length(blinks)
                if bb == 1
                    condensedBlinkTimes(1,bb) = blinks(1,bb);
                    cc = cc + 1;
                else
                    timeDifference = blinks(1,bb) - blinks(1,bb - 1);
                    if timeDifference > 6
                        condensedBlinkTimes(1,cc) = blinks(1,bb);
                        cc = cc + 1;
                    end
                end
            end
            % extract blink triggered data
            for dd = 1:length(condensedBlinkTimes)
                blink = condensedBlinkTimes(1,dd);
                if (blink/samplingRate) >= 3 && (blink/samplingRate) <= 889
                    LH_HbT = filtfilt(sos2,g2,ProcData.data.CBV_HbT.adjLH);
                    LH_hbtArray = LH_HbT((blink - (2*samplingRate)):(blink + (10*samplingRate)));
                    LH_hbtArray = LH_hbtArray - mean(LH_hbtArray(1:2*samplingRate));
                    data.LH_HbT = cat(1,data.LH_HbT,LH_hbtArray);
                    
                    RH_HbT = filtfilt(sos2,g2,ProcData.data.CBV_HbT.adjRH);
                    RH_hbtArray = RH_HbT((blink - (2*samplingRate)):(blink + (10*samplingRate)));
                    RH_hbtArray = RH_hbtArray - mean(RH_hbtArray(1:2*samplingRate));
                    data.RH_HbT = cat(1,data.RH_HbT,RH_hbtArray);
                    
                    LH_corticalS_Data = SpecData.cortical_LH.normS;
                    data.F = SpecData.cortical_LH.F;
                    T = round(SpecData.cortical_LH.T,1);
                    data.T = -2:(1/specSamplingRate):10;
                    startTimeIndex = find(T == round((blink/samplingRate) - 2),1);
                    durationIndex = startTimeIndex + 12*specSamplingRate;
                    LH_corticalS_Vals = LH_corticalS_Data(:,startTimeIndex:durationIndex);
                    data.LH_cort = cat(3,data.LH_cort,LH_corticalS_Vals);
                    
                    RH_corticalS_Data = SpecData.cortical_RH.normS;
                    RH_corticalS_Vals = RH_corticalS_Data(:,startTimeIndex:durationIndex);
                    data.RH_cort = cat(3,data.RH_cort,RH_corticalS_Vals);
                    
                    hipS_Data = SpecData.hippocampus.normS;
                    hipS_Vals = hipS_Data(:,startTimeIndex:durationIndex);
                    data.hip = cat(3,data.hip,hipS_Vals);
                    
                    EMG = filtfilt(sos,g,ProcData.data.EMG.emg - RestingBaselines.manualSelection.EMG.emg.(strDay));
                    emgArray = EMG((blink - (2*samplingRate)):(blink + (10*samplingRate)));
                    data.EMG = cat(1,data.EMG,emgArray);
                    
                    if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
                        try
                            pupilArea = filtfilt(sos,g,ProcData.data.Pupil.pupilArea);
                        catch
                            pupilArea = ProcData.data.Pupil.pupilArea;
                        end
                        areaArray = pupilArea((blink - (2*samplingRate)):(blink + (10*samplingRate)));
                        data.EMG = cat(1,data.pupilArea,areaArray);
                    end
                end
            end
        end
    end
end
Results_BlinkResponses.(animalID).pupilArea = mean(data.pupilArea,1);
Results_BlinkResponses.(animalID).LH_HbT = mean(data.LH_HbT,1);
Results_BlinkResponses.(animalID).RH_HbT = mean(data.RH_HbT,1);
Results_BlinkResponses.(animalID).LH_cort = mean(data.LH_cort,3);
Results_BlinkResponses.(animalID).RH_cort = mean(data.RH_cort,3);
Results_BlinkResponses.(animalID).hip = mean(data.hip,3);
Results_BlinkResponses.(animalID).EMG = mean(data.EMG,1);
Results_BlinkResponses.(animalID).T = data.T;
Results_BlinkResponses.(animalID).F = data.F;
% save data
cd([rootFolder delim])
save('Results_BlinkResponses.mat','Results_BlinkResponses')

end
