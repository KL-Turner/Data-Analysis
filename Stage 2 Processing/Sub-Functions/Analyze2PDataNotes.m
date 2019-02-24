function Analyze2PDataNotes(msExcel_File)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Patrick J. Drew: https://github.com/DrewLab
%________________________________________________________________________________________________________________________
%
%   Purpose: 
%________________________________________________________________________________________________________________________
%
%   Inputs: 
%
%   Outputs:
%
%   Last Revised: January 29th, 2019
%________________________________________________________________________________________________________________________

% Read the image info from the formated xls file, save as a RawData file with surface, penetrating arterioles and capillaries separate.
[~, ~, alldata] = xlsread(msExcel_File);
for row = 2:size(alldata, 1)   % Loop through all rows of the excel sheet
    clear MScanData
    %% Notes
    tempData.Notes.date = num2str(alldata{row, 1});
    tempData.Notes.animalID = alldata{row, 2};
    tempData.Notes.imageID = alldata{row, 3};
    tempData.Notes.movieType = alldata{row, 4};
    tempData.Notes.laserPower = alldata{row, 5};
    tempData.Notes.objectiveID = alldata{row, 6};
    tempData.Notes.frameRate = alldata{row, 7};
    tempData.Notes.numberOfFrames = alldata{row, 8};
    tempData.Notes.vesselType = alldata{row, 9};
    tempData.Notes.vesselDepth = alldata{row, 10};
    tempData.Notes.comments = alldata{row, 11};
    tempData.Notes.vesselID = alldata{row, 12};
    tempData.Notes.drug = alldata{row, 13};
    
    currentFileID = ([tempData.Notes.animalID '_' tempData.Notes.imageID '_' tempData.Notes.date '_MScanData.mat']);
    if ~exist(currentFileID) %#ok<EXIST>
        %% Vessel diameter calculation for movie surface vessels
        if strcmp(tempData.Notes.movieType, 'MS')
            MScanData = DiamCalcSurfaceVessel(tempData, [tempData.Notes.date '_' tempData.Notes.imageID]);
            
            % Vessel diameter calculation for movie penetrating
        elseif strcmp(tempData.Notes.movieType, 'MP')
            MScanData = DiamCalcPenetratingVessel(tempData, [tempData.Notes.date '_' tempData.Notes.imageID]);
            
            % Vessel diameter calculation for capillaries
        elseif strcmp(tempData.Notes.movieType, 'C')
            MScanData = CapillaryLineScan(tempData, [tempData.Notes.date '_' tempData.Notes.imageID]);
        end
        
        % Save the RawData file for the current movie type
        disp(['File Created. Saving MScanData File ' num2str(row - 1) '...']); disp(' ')
        save([MScanData.Notes.animalID '_' MScanData.Notes.date '_' MScanData.Notes.imageID '_MScanData'], 'MScanData')
        close all
    end
end

end
