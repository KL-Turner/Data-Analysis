function CombineLabVIEWMScanFiles_2P(~, mscanDataFiles)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Combine the MScan and LabVIEW data structures into one.
%________________________________________________________________________________________________________________________
%
%   Inputs: List of MScan and LabVIEW data files.
%
%   Outputs: Single MergedData structure with the important information from both.
%
%   Last Revised: February 29th, 2019
%________________________________________________________________________________________________________________________

for a = 1:size(mscanDataFiles,1)
    disp(['Combining the data from LabVIEW and MScan file(s) number ' num2str(a) ' of ' num2str(size(mscanDataFiles, 1)) '...']); disp(' ');
    mscanDataFile = mscanDataFiles(a,:);
    load(mscanDataFile);
    labviewDataFile = MScanData.notes.labviewFileID;
    load(labviewDataFile);

    [animalID, ~, ~, fileID] = GetFileInfo_2P(labviewDataFile);
    vesselID = MScanData.notes.vesselID;
    imageID = MScanData.notes.imageID;
    
    % Pull the notes and data from LabVIEW
    MergedData.notes.LabVIEW = LabVIEWData.notes;
    MergedData.data.whiskerAngle = LabVIEWData.data.dsWhiskerAngle_trim;
    MergedData.data.binWhiskerAngle = LabVIEWData.data.binWhiskerAngle_trim;
    MergedData.data.forceSensorL = LabVIEWData.data.dsForceSensorL_trim;
    MergedData.data.binForceSensorL = LabVIEWData.data.binForceSensorL_trim;
    
    % Pull the notes and data from MScan
    MergedData.notes.MScan = MScanData.notes;
    MergedData.data.rawNeuralData = MScanData.data.rawNeuralData_trim;
    MergedData.data.muaPower = MScanData.data.muaPower_trim';
    MergedData.data.gammaPower = MScanData.data.gammaPower_trim';
    MergedData.data.betaPower = MScanData.data.betaPower_trim';
    MergedData.data.alphaPower = MScanData.data.alphaPower_trim';
    MergedData.data.thetaPower = MScanData.data.thetaPower_trim';
    MergedData.data.deltaPower = MScanData.data.deltaPower_trim';
    MergedData.data.forceSensorM = MScanData.data.dsForceSensorM_trim;
    MergedData.data.EMG = MScanData.data.dsEMG_trim;
    MergedData.data.binForceSensorM = MScanData.data.binForceSensorM_trim;
    MergedData.data.vesselDiameter = MScanData.data.vesselDiameter_trim;
    
    % Most useful notes to be referenced in future analysis
    MergedData.notes.trialDuration_Sec = LabVIEWData.notes.trialDuration_Seconds_trim;
    MergedData.notes.p2Fs = MScanData.notes.frameRate;
    MergedData.notes.dsFs = MScanData.notes.downSampledFs;
    MergedData.notes.anFs = LabVIEWData.notes.analogSamplingRate_Hz;
    
    save([animalID '_' vesselID '_' fileID '_' imageID '_MergedData'], 'MergedData')
end

end
