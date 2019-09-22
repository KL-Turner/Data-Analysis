function [] = ApplySleepLogical_SVM(procDataFileIDs, SVMResults)
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
%   Last Revised: July 27th, 2019
%________________________________________________________________________________________________________________________

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    load(procDataFileID)
    fileID = procDataFileID(1:end-13);
    c = 1;
    for b = 1:length(SVMResults.fileIDs)
        svmFileCheck = SVMResults.fileIDs{b,1};
        svmFileID = svmFileCheck(1:end-14);
        if strcmp(fileID, svmFileID) == true
            behavState{c,1} = SVMResults.labels{b,1};
            c = c+1;
        end
    end
    
    for d = 1:length(behavState)
        if strcmp(behavState{d,1}, 'Not Sleep') == true
            awakeLogical(d,1) = 1;
            nremLogical(d,1) = 0;
            remLogical(d,1) = 0;
        elseif strcmp(behavState{d,1}, 'NREM Sleep') == true
            awakeLogical(d,1) = 0;
            nremLogical(d,1) = 1;
            remLogical(d,1) = 0;
        elseif strcmp(behavState{d,1}, 'REM Sleep') == true
            awakeLogical(d,1) = 0;
            nremLogical(d,1) = 0;
            remLogical(d,1) = 1;
        end
    end
    
    ProcData.sleep.logicals.awakeLogical = logical(awakeLogical);
    ProcData.sleep.logicals.nremLogical = logical(nremLogical);
    ProcData.sleep.logicals.remLogical = logical(remLogical);
    save(procDataFileID, 'ProcData')
    
end

end