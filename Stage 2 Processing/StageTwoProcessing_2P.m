%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner 
% Ph.D. Candidate, Department of Bioengineering 
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: 1) Generate bilateral ROIs for CBV analysis.
%            2) Create ProcData structure using threshholds for the observed data.
%________________________________________________________________________________________________________________________
%
%   Inputs: 1) Selects all _RawData files from all days. Draw all left ROIs first, then the right.
%           2) Selects all _RawData files from all days. Follow the command window prompts for each threshold value.
%
%   Outputs: 1) An animal_hemisphere_ROIs.mat file with the xi, yi coordinates for each day.
%            2) A ProcData.mat structure for each inputed rawdata.mat file, as well as a Thresholds.mat file that serves
%               as a record for each variable/day's set threshold.        
%
%   Last Revised: October 3rd, 2018    
%________________________________________________________________________________________________________________________

msExcel_File = uigetfile('*.xlsx');
Analyze_2P_DataNotes(msExcel_File);

mscanDirectory = dir('*_MScanData.mat');
mscanDataFiles = {mscanDirectory.name}';
mscanDataFiles = char(mscanDataFiles);

Analyze_2P_Diameter(mscanDirectory);

labviewDirectory = dir('*_LabVIEWData.mat');
labviewDataFiles = {labviewDirectory.name}';
labviewDataFiles = char(labviewDataFiles);

Combine_LabVIEW_MScan_Files(labviewDataFiles, mscanDataFiles)

combDirectory = dir('*_CombData.mat');
combDataFiles = {combDirectory.name}';
combDataFiles = char(combDataFiles);

Correct_LabVIEW_Offset(combDataFiles)

ProcessCombDataFile(combDataFiles)



