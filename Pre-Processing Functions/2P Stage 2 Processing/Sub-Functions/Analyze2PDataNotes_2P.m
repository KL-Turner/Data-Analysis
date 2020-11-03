function [] = Analyze2PDataNotes_2P(msExcelFile)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab
%________________________________________________________________________________________________________________________
%
%   Purpose: Extract the vessel notes and information from a MS Excel file, and create a struct titled '_MScanData.mat'.
%________________________________________________________________________________________________________________________

% Read the image info from the formated xls file
[~,~,alldata] = xlsread(msExcelFile);
for a = 2:size(alldata,1)   % Loop through all rows of the excel sheet except the first row
    clear MScanData
    %% notes
    tempData.notes.date = num2str(alldata{a,1});
    tempData.notes.animalID = alldata{a,2};
    tempData.notes.imageID = alldata{a,3};
    tempData.notes.movieType = alldata{a,4};
    tempData.notes.laserPower = alldata{a,5};
    tempData.notes.objectiveID = alldata{a,6};
    tempData.notes.frameRate = alldata{a,7};
    tempData.notes.numberOfFrames = alldata{a,8};
    tempData.notes.vesselType = alldata{a,9};
    tempData.notes.vesselDepth = alldata{a,10};
    tempData.notes.comments = alldata{a,11};
    tempData.notes.vesselID = alldata{a,12};
    tempData.notes.drug = alldata{a,13};
    currentFileID = ([tempData.notes.animalID '_' tempData.notes.date '_' tempData.notes.imageID '_' tempData.notes.vesselID '_MScanData.mat']);
    if ~exist(currentFileID,'file')   % Only run analysis if the current file doesn't exist yet
        % Vessel diameter calculation for surface vessels
        if strcmp(tempData.notes.movieType,'MS') == true
            [MScanData] = DiamCalcSurfaceVessel_2P(tempData,[tempData.notes.date '_' tempData.notes.imageID]);     
        % Vessel diameter (area) calculation for penetrating vessels
        elseif strcmp(tempData.notes.movieType,'MP') == true
            [MScanData] = AreaCalcPenVessel_2P(tempData,[tempData.notes.date '_' tempData.notes.imageID]);
        % Line scan calculation for capillaries
        elseif strcmp(tempData.notes.movieType,'C') == true
            [MScanData] = CapLinesScan_2P(tempData,[tempData.notes.date '_' tempData.notes.imageID]);
        end     
        % Checklist for analysis steps - debugging purposes
        MScanData.notes.checklist.analyzeDiam = false;
        MScanData.notes.checklist.processData = false;
        MScanData.notes.checklist.offsetCorrect = false;   
        % Save the RawData file for the current movie type
        disp(['File Created. Saving MScanData File ' num2str(a - 1) '...']); disp(' ')
        save([MScanData.notes.animalID '_' MScanData.notes.date '_' MScanData.notes.imageID '_' MScanData.notes.vesselID '_MScanData'],'MScanData')
        close all
    else
        disp([currentFileID ' already exists in the current directory. Continuing...']); disp(' ')
    end
end

end
