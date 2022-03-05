function [] = MainScript_Pupil()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
% Purpose:
%
% Scripts used to pre-process the original data are located in the folder "Pre-Processing Scripts".
% Functions that are used in both the analysis and pre-processing are located in the analysis folder.
%________________________________________________________________________________________________________________________

zap;
multiWaitbar('CloseAll');
%% verify code repository and data are in the current directory/added path
currentFolder = pwd;
addpath(genpath(currentFolder));
fileparts = strsplit(currentFolder,filesep);
if ismac
    rootFolder = fullfile(filesep,fileparts{1:end});
    delim = '/';
else
    rootFolder = fullfile(fileparts{1:end});
    delim = '\';
end
% add root folder to Matlab's working directory
addpath(genpath(rootFolder))
%% analysis subfunctions
runAnalysis = true;
if runAnalysis == true
    %     AnalyzeBehavioralArea_Pupil_Handler(rootFolder,delim,false)
    %     AnalyzeEvokedResponses_Pupil_Handler(rootFolder,delim,false)
    %     AnalyzeSleepModelAccuracy_Pupil_Handler(rootFolder,delim,false)
    %     AnalyzePupilSleepModelAccuracy_Pupil_Handler(rootFolder,delim,false)
    %     AnalyzePupilAreaSleepProbability_Pupil_Handler(rootFolder,delim,false)
    %     AnalyzeBlinkResponses_Pupil_Handler(rootFolder,delim,false)
    %     AnalyzePowerSpectrum_Pupil_Handler(rootFolder,delim,false)
    %     AnalyzeCoherence_Pupil_Handler(rootFolder,delim,true)
    %     AnalyzeCrossCorrelation_Pupil_Handler(rootFolder,delim,false)
    %     AnalyzeStimulusBlinks_Pupil_Handler(rootFolder,delim,false)
    %     AnalyzeBlinkPeriodogram_Pupil_Handler(rootFolder,delim,false)
    %     AnalyzePupilHbTRelationship_Pupil_Handler(rootFolder,delim,false)
    %     AnalyzePupilGammaRelationship_Pupil_Handler(rootFolder,delim,false)
    %     AnalyzeBlinkCoherogram_Pupil_Handler(rootFolder,delim,false)
    AnalyzeBlinkTransition_Pupil_Handler(rootFolder,delim,false)
    multiWaitbar('CloseAll');
end
%% generate figures
disp('Loading analysis results and generating figures...'); disp(' ')
saveFigs = true;
% ArousalStateDiameter_Pupil(rootFolder,saveFigs,delim);
% WhiskStimEvoked_Pupil(rootFolder,saveFigs,delim);
% PupilSleepModelAccuracy2_Pupil(rootFolder,saveFigs,delim);
% PupilSleepProbability_Pupil(rootFolder,saveFigs,delim)
% BlinkResponses_Pupil(rootFolder,saveFigs,delim)
% PowerSpectrum_Pupil(rootFolder,saveFigs,delim)
% Coherence_Pupil(rootFolder,saveFigs,delim)
% CrossCorrelation_Pupil(rootFolder,saveFigs,delim);
% BlinkPeriodogram_Pupil(rootFolder,saveFigs,delim);
% StimulusBlinks_Pupil(rootFolder,saveFigs,delim)
% PupilHbTRelationship(rootFolder,saveFigs,delim)
% PupilGammaRelationship(rootFolder,saveFigs,delim)
% BlinkCoherogram(rootFolder,saveFigs,delim)
BlinkTransition(rootFolder,saveFigs,delim)
%% main figures
% Fig1_TBD(rootFolder,saveFigs,delim)
% Fig2_TBD(rootFolder,saveFigs,delim)
% Fig3_TBD(rootFolder,saveFigs,delim)
% Fig4_TBD(rootFolder,saveFigs,delim)
% Fig5_TBD(rootFolder,saveFigs,delim)
%% supplement
% FigS1_TBD(rootFolder,saveFigs,delim)
% FigS2_TBD(rootFolder,saveFigs,delim)
% FigS3_TBD(rootFolder,saveFigs,delim)
% FigS4_TBD(rootFolder,saveFigs,delim)

end
