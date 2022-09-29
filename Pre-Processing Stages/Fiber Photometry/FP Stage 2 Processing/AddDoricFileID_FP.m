function [] = AddDoricFileID_FP(doricDataFiles,labviewDataFiles)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Add the corresponding LabVIEW filename to each Doric file set
%________________________________________________________________________________________________________________________

for aa = 1:size(labviewDataFileIDs)
    labviewDataFile = labviewDataFiles(aa,:);
    doricDataFile = doricDataFiles(aa,:);
    
    for d = 1:size(doricDataFiles,1)
        mscanDataFile = doricDataFiles(d,:);
        fileBreaks = strfind(mscanDataFile,'_');
        if strcmp(mscanDataFile(fileBreaks(2) + 1:fileBreaks(3) - 1),imageID)
            load(mscanDataFile)
            if ~isfield(MScanData.notes,'labviewFileID')
                disp(['Adding LabVIEW file ID for ' mscanDataFile '...']); disp(' ')
                MScanData.notes.labviewFileID = labviewDataFile;
                save(mscanDataFile,'MScanData')
            else
                disp(['LabVIEW file already added for ' mscanDataFile '. Continuing...']); disp(' ')
            end
        end
    end
end

end
