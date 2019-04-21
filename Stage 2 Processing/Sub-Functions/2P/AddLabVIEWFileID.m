function [] = AddLabVIEWFileID(mscanDataFiles, labviewDataFiles, msExcelFile)

[~, ~, alldata] = xlsread(msExcelFile);
b = 1;
for a = 2:size(alldata,1)
    imageIDs{b,1} = alldata{a,3};
    b = b + 1;
end

uniqueImages = unique(imageIDs);

for c = 1:size(uniqueImages)
    labviewDataFile = labviewDataFiles(c,:);
    imageID = string(uniqueImages(c,1));
    
    for d = 1:size(mscanDataFiles,1)
        mscanDataFile = mscanDataFiles(d,:);
        fileBreaks = strfind(mscanDataFile, '_');
        if strcmp(mscanDataFile(fileBreaks(3) + 1:fileBreaks(4) - 1), imageID)
            load(mscanDataFile)
            
            if ~isfield(MScanData.notes, 'labviewFileID')
                disp(['Adding LabVIEW file ID for ' mscanDataFile '...']); disp(' ')
                MScanData.notes.labviewFileID = labviewDataFile;
                save(mscanDataFile, 'MScanData')
            else
                disp(['LabVIEW file already added for ' mscanDataFile '. Continuing...']); disp(' ')
            end
        end
    end
end


end