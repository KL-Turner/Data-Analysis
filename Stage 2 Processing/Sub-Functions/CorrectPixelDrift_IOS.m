function [] = CorrectPixelDrift_IOS(procDataFileIDs,imagingType)
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

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    load(procDataFileID)
    disp(['Correcting pixel drift in ProcData file ' num2str(a) ' of ' num2str(size(procDataFileIDs,1)) '...']); disp(' ')
    if strcmp(imagingType,'bilateral') == true
        LH_CBV = ProcData.data.CBV.LH;
        RH_CBV = ProcData.data.CBV.RH;
        CementDrift = ProcData.data.CBV.Cement;
        [B, A] = butter(3, 0.1/(30/2), 'low');
        filtCement_Data = filtfilt(B, A, CementDrift);
        correctedLH_CBV = (LH_CBV - filtCement_Data);
        correctedRH_CBV = (RH_CBV - filtCement_Data);
        LH_DC_reset = mean(correctedLH_CBV);
        RH_DC_reset = mean(correctedRH_CBV);
        LH_1kDiff = 1000 - LH_DC_reset;
        RH_1kDiff = 1000 - RH_DC_reset;
        LH_DC_shift = ones(1,length(correctedLH_CBV))*LH_1kDiff;
        RH_DC_shift = ones(1,length(correctedRH_CBV))*RH_1kDiff;
        shiftedCorrectedLH_CBV = correctedLH_CBV + LH_DC_shift;
        shiftedCorrectedRH_CBV = correctedRH_CBV + RH_DC_shift;
        ProcData.data.CBV.LH = shiftedCorrectedLH_CBV;
        ProcData.data.CBV.RH = shiftedCorrectedRH_CBV;
    elseif strcmp(imagingType,'single') == true
        CBV = ProcData.data.CBV.Barrels;
        CementDrift = ProcData.data.CBV.Cement;
        [B, A] = butter(3, 0.1/(30/2), 'low');
        filtCement_Data = filtfilt(B, A, CementDrift);
        correctedCBV = (CBV - filtCement_Data);
        CBV_DC_reset = mean(correctedCBV);
        CBV_1kDiff = 1000 - CBV_DC_reset;
        CBV_DC_shift = ones(1,length(correctedCBV))*CBV_1kDiff;
        shiftedCorrected_CBV = correctedCBV + CBV_DC_shift;
        ProcData.data.CBV.Barrels = shiftedCorrected_CBV;
    end
    save(procDataFileID, 'ProcData');
end
