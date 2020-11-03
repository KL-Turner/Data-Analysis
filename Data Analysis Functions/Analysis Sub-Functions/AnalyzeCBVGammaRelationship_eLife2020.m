function [AnalysisResults] = AnalyzeCBVGammaRelationship_eLife2020(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the relationship between gamma-band power and hemodynamics [HbT] (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
IOS_animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
%% only run analysis for valid animal IDs
if any(strcmp(IOS_animalIDs,animalID))
    dataLocation = [rootFolder '\' animalID '\Bilateral Imaging\'];
    cd(dataLocation)
    % find and load manual baseline event information
    load('Forest_ScoringResults.mat')
    % find and load EventData.mat struct
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFile = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFile);
    % extract/concatenate data from each file
    catLH_Gamma = [];
    catRH_Gamma = [];
    catLH_HbT = [];
    catRH_HbT = [];
    for aa = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(aa,:);
        load(procDataFileID)
        for bb = 1:length(ProcData.sleep.parameters.cortical_LH.gammaBandPower)
            LH_Gamma = mean(ProcData.sleep.parameters.cortical_LH.specGammaBandPower{bb,1}{1,1});
            RH_Gamma = mean(ProcData.sleep.parameters.cortical_RH.specGammaBandPower{bb,1}{1,1});
            LH_HbT = mean(ProcData.sleep.parameters.CBV.hbtLH{bb,1});
            RH_HbT = mean(ProcData.sleep.parameters.CBV.hbtRH{bb,1});
            % group date based on arousal-state classification
            state = ScoringResults.labels{aa,1}{bb,1};
            if strcmp(state,'Not Sleep') == true
                if isfield(catLH_Gamma,'Awake') == false
                    catLH_Gamma.Awake = [];
                    catRH_Gamma.Awake = [];
                    catLH_HbT.Awake = [];
                    catRH_HbT.Awake = [];
                end
                catLH_Gamma.Awake = cat(1,catLH_Gamma.Awake,LH_Gamma);
                catRH_Gamma.Awake = cat(1,catRH_Gamma.Awake,RH_Gamma);
                catLH_HbT.Awake = cat(1,catLH_HbT.Awake,LH_HbT);
                catRH_HbT.Awake = cat(1,catRH_HbT.Awake,RH_HbT);
            elseif strcmp(state,'NREM Sleep') == true
                if isfield(catLH_Gamma,'NREM') == false
                    catLH_Gamma.NREM = [];
                    catRH_Gamma.NREM = [];
                    catLH_HbT.NREM = [];
                    catRH_HbT.NREM = [];
                end
                catLH_Gamma.NREM = cat(1,catLH_Gamma.NREM,LH_Gamma);
                catRH_Gamma.NREM = cat(1,catRH_Gamma.NREM,RH_Gamma);
                catLH_HbT.NREM = cat(1,catLH_HbT.NREM,LH_HbT);
                catRH_HbT.NREM = cat(1,catRH_HbT.NREM,RH_HbT);
            elseif strcmp(state,'REM Sleep') == true
                if isfield(catLH_Gamma,'REM') == false
                    catLH_Gamma.REM = [];
                    catRH_Gamma.REM = [];
                    catLH_HbT.REM = [];
                    catRH_HbT.REM = [];
                end
                catLH_Gamma.REM = cat(1,catLH_Gamma.REM,LH_Gamma);
                catRH_Gamma.REM = cat(1,catRH_Gamma.REM,RH_Gamma);
                catLH_HbT.REM = cat(1,catLH_HbT.REM,LH_HbT);
                catRH_HbT.REM = cat(1,catRH_HbT.REM,RH_HbT);
            end
        end
    end
    % save results
    AnalysisResults.(animalID).HbTvsGamma = [];
    AnalysisResults.(animalID).HbTvsGamma.Awake.Gamma = cat(1,catLH_Gamma.Awake,catRH_Gamma.Awake);
    AnalysisResults.(animalID).HbTvsGamma.NREM.Gamma = cat(1,catLH_Gamma.NREM,catRH_Gamma.NREM);
    AnalysisResults.(animalID).HbTvsGamma.REM.Gamma = cat(1,catLH_Gamma.REM,catRH_Gamma.REM);
    AnalysisResults.(animalID).HbTvsGamma.Awake.HbT = cat(1,catLH_HbT.Awake,catRH_HbT.Awake);
    AnalysisResults.(animalID).HbTvsGamma.NREM.HbT = cat(1,catLH_HbT.NREM,catRH_HbT.NREM);
    AnalysisResults.(animalID).HbTvsGamma.REM.HbT = cat(1,catLH_HbT.REM,catRH_HbT.REM);
    % save data
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end

end
