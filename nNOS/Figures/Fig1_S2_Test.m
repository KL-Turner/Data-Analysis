function [AnalysisResults] = Fig1_S2_Test(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S2 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
animalIDs = {'T135','T141','T142','T144','T151','T155','T156','T157','T159'};
solenoidNames = {'LPadSol','RPadSol','AudSol'};
% compDataTypes = {'Ipsi','Contra','Auditory'};
compDataTypes = {'Contra'};
dataTypes = {'adjLH','adjRH'};
treatments = {'C57BL6J','SSP_SAP'};
%% cd through each animal's directory and extract the appropriate analysis results
xx = 1;
zz = 1;
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    if strcmp(animalID,'T141') == true || strcmp(animalID,'T155') == true || strcmp(animalID,'T156') == true || strcmp(animalID,'T157') == true
        for bb = 1:length(dataTypes)
            dataType = dataTypes{1,bb};
            for dd = 1:length(solenoidNames)
                solenoidName = solenoidNames{1,dd};
                data.C57BL6J.(dataType).(solenoidName).count(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).count;
                data.C57BL6J.(dataType).(solenoidName).HbT(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).CBV_HbT.HbT;
                data.C57BL6J.(dataType).(solenoidName).CBV(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).CBV.CBV;
                data.C57BL6J.(dataType).(solenoidName).cortMUA(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).MUA.corticalData;
                data.C57BL6J.(dataType).(solenoidName).hipMUA(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).MUA.hippocampalData;
                data.C57BL6J.(dataType).(solenoidName).cortGam(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).Gam.corticalData;
                data.C57BL6J.(dataType).(solenoidName).hipGam(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).Gam.hippocampalData;
                data.C57BL6J.(dataType).(solenoidName).timeVector(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).timeVector;
                data.C57BL6J.(dataType).(solenoidName).cortS(:,:,xx) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).LFP.corticalS;
                data.C57BL6J.(dataType).(solenoidName).cortS_Gam(:,:,xx) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).LFP.corticalS(49:end,20:23);
                data.C57BL6J.(dataType).(solenoidName).hipS(:,:,xx) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).LFP.hippocampalS;
                data.C57BL6J.(dataType).(solenoidName).hipS_Gam(:,:,xx) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).LFP.hippocampalS(49:end,20:23);
                data.C57BL6J.(dataType).(solenoidName).T(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).LFP.T;
                data.C57BL6J.(dataType).(solenoidName).F(:,xx) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).LFP.F;
            end
        end
        xx = xx + 1;
    elseif strcmp(animalID,'T135') == true || strcmp(animalID,'T142') == true || strcmp(animalID,'T144') == true || strcmp(animalID,'T151') == true || strcmp(animalID,'T159') == true
        for bb = 1:length(dataTypes)
            dataType = dataTypes{1,bb};
            for dd = 1:length(solenoidNames)
                solenoidName = solenoidNames{1,dd};
                data.SSP_SAP.(dataType).(solenoidName).count(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).count;
                data.SSP_SAP.(dataType).(solenoidName).HbT(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).CBV_HbT.HbT;
                data.SSP_SAP.(dataType).(solenoidName).CBV(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).CBV.CBV;
                data.SSP_SAP.(dataType).(solenoidName).cortMUA(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).MUA.corticalData;
                data.SSP_SAP.(dataType).(solenoidName).hipMUA(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).MUA.hippocampalData;
                data.SSP_SAP.(dataType).(solenoidName).cortGam(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).Gam.corticalData;
                data.SSP_SAP.(dataType).(solenoidName).hipGam(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).Gam.hippocampalData;
                data.SSP_SAP.(dataType).(solenoidName).timeVector(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).timeVector;
                data.SSP_SAP.(dataType).(solenoidName).cortS(:,:,zz) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).LFP.corticalS;
                data.SSP_SAP.(dataType).(solenoidName).cortS_Gam(:,:,zz) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).LFP.corticalS(49:end,20:23);
                data.SSP_SAP.(dataType).(solenoidName).hipS(:,:,zz) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).LFP.hippocampalS;
                data.SSP_SAP.(dataType).(solenoidName).hipS_Gam(:,:,zz) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).LFP.hippocampalS(49:end,20:23);
                data.SSP_SAP.(dataType).(solenoidName).T(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).LFP.T;
                data.SSP_SAP.(dataType).(solenoidName).F(:,zz) = AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoidName).LFP.F;
            end
        end
        zz = zz + 1;
    end
end
%% concatenate the data from the contra and ipsi data
data.C57BL6J.adjLH.Contra.count = data.C57BL6J.adjLH.RPadSol.count;
data.C57BL6J.adjRH.Contra.count = data.C57BL6J.adjRH.LPadSol.count;

data.C57BL6J.adjLH.Contra.HbT = data.C57BL6J.adjLH.RPadSol.HbT;
data.C57BL6J.adjRH.Contra.HbT = data.C57BL6J.adjRH.LPadSol.HbT;

data.C57BL6J.adjLH.Contra.CBV = data.C57BL6J.adjLH.RPadSol.CBV;
data.C57BL6J.adjRH.Contra.CBV = data.C57BL6J.adjRH.LPadSol.CBV;

data.C57BL6J.adjLH.Contra.cortMUA = data.C57BL6J.adjLH.RPadSol.cortMUA;
data.C57BL6J.adjRH.Contra.cortMUA = data.C57BL6J.adjRH.LPadSol.cortMUA;

data.C57BL6J.adjLH.Contra.hipMUA = data.C57BL6J.adjLH.RPadSol.hipMUA;
data.C57BL6J.adjRH.Contra.hipMUA = data.C57BL6J.adjRH.LPadSol.hipMUA;

data.C57BL6J.adjLH.Contra.cortGam = data.C57BL6J.adjLH.RPadSol.cortGam;
data.C57BL6J.adjRH.Contra.cortGam = data.C57BL6J.adjRH.LPadSol.cortGam;

data.C57BL6J.adjLH.Contra.hipGam = data.C57BL6J.adjLH.RPadSol.hipGam;
data.C57BL6J.adjRH.Contra.hipGam = data.C57BL6J.adjRH.LPadSol.hipGam;

data.C57BL6J.adjLH.Contra.timeVector = data.C57BL6J.adjLH.RPadSol.timeVector;
data.C57BL6J.adjRH.Contra.timeVector = data.C57BL6J.adjRH.LPadSol.timeVector;

data.C57BL6J.adjLH.Contra.cortS = data.C57BL6J.adjLH.RPadSol.cortS;
data.C57BL6J.adjRH.Contra.cortS = data.C57BL6J.adjRH.LPadSol.cortS;

data.C57BL6J.adjLH.Contra.cortS_Gam = data.C57BL6J.adjLH.RPadSol.cortS_Gam;
data.C57BL6J.adjRH.Contra.cortS_Gam = data.C57BL6J.adjRH.LPadSol.cortS_Gam;

data.C57BL6J.adjLH.Contra.hipS = data.C57BL6J.adjLH.RPadSol.hipS;
data.C57BL6J.adjRH.Contra.hipS = data.C57BL6J.adjRH.LPadSol.hipS;

data.C57BL6J.adjLH.Contra.hipS_Gam = data.C57BL6J.adjLH.RPadSol.hipS_Gam;
data.C57BL6J.adjRH.Contra.hipS_Gam = data.C57BL6J.adjRH.LPadSol.hipS_Gam;

data.C57BL6J.adjLH.Contra.T = data.C57BL6J.adjLH.RPadSol.T;
data.C57BL6J.adjRH.Contra.T = data.C57BL6J.adjRH.LPadSol.T;

data.C57BL6J.adjLH.Contra.F = data.C57BL6J.adjLH.RPadSol.F;
data.C57BL6J.adjRH.Contra.F = data.C57BL6J.adjRH.LPadSol.F;

data.SSP_SAP.adjLH.Contra.count = data.SSP_SAP.adjLH.RPadSol.count;
data.SSP_SAP.adjRH.Contra.count = data.SSP_SAP.adjRH.LPadSol.count;

data.SSP_SAP.adjLH.Contra.HbT = data.SSP_SAP.adjLH.RPadSol.HbT;
data.SSP_SAP.adjRH.Contra.HbT = data.SSP_SAP.adjRH.LPadSol.HbT;

data.SSP_SAP.adjLH.Contra.CBV = data.SSP_SAP.adjLH.RPadSol.CBV;
data.SSP_SAP.adjRH.Contra.CBV = data.SSP_SAP.adjRH.LPadSol.CBV;

data.SSP_SAP.adjLH.Contra.cortMUA = data.SSP_SAP.adjLH.RPadSol.cortMUA;
data.SSP_SAP.adjRH.Contra.cortMUA = data.SSP_SAP.adjRH.LPadSol.cortMUA;

data.SSP_SAP.adjLH.Contra.hipMUA = data.SSP_SAP.adjLH.RPadSol.hipMUA;
data.SSP_SAP.adjRH.Contra.hipMUA = data.SSP_SAP.adjRH.LPadSol.hipMUA;

data.SSP_SAP.adjLH.Contra.cortGam = data.SSP_SAP.adjLH.RPadSol.cortGam;
data.SSP_SAP.adjRH.Contra.cortGam = data.SSP_SAP.adjRH.LPadSol.cortGam;

data.SSP_SAP.adjLH.Contra.hipGam = data.SSP_SAP.adjLH.RPadSol.hipGam;
data.SSP_SAP.adjRH.Contra.hipGam = data.SSP_SAP.adjRH.LPadSol.hipGam;

data.SSP_SAP.adjLH.Contra.timeVector = data.SSP_SAP.adjLH.RPadSol.timeVector;
data.SSP_SAP.adjRH.Contra.timeVector = data.SSP_SAP.adjRH.LPadSol.timeVector;

data.SSP_SAP.adjLH.Contra.cortS = data.SSP_SAP.adjLH.RPadSol.cortS;
data.SSP_SAP.adjRH.Contra.cortS = data.SSP_SAP.adjRH.LPadSol.cortS;

data.SSP_SAP.adjLH.Contra.cortS_Gam = data.SSP_SAP.adjLH.RPadSol.cortS_Gam;
data.SSP_SAP.adjRH.Contra.cortS_Gam = data.SSP_SAP.adjRH.LPadSol.cortS_Gam;

data.SSP_SAP.adjLH.Contra.hipS = data.SSP_SAP.adjLH.RPadSol.hipS;
data.SSP_SAP.adjRH.Contra.hipS = data.SSP_SAP.adjRH.LPadSol.hipS;

data.SSP_SAP.adjLH.Contra.hipS_Gam = data.SSP_SAP.adjLH.RPadSol.hipS_Gam;
data.SSP_SAP.adjRH.Contra.hipS_Gam = data.SSP_SAP.adjRH.LPadSol.hipS_Gam;

data.SSP_SAP.adjLH.Contra.T = data.SSP_SAP.adjLH.RPadSol.T;
data.SSP_SAP.adjRH.Contra.T = data.SSP_SAP.adjRH.LPadSol.T;

data.SSP_SAP.adjLH.Contra.F = data.SSP_SAP.adjLH.RPadSol.F;
data.SSP_SAP.adjRH.Contra.F = data.SSP_SAP.adjRH.LPadSol.F;

% data.Contra.HbT = cat(2,data.adjLH.RPadSol.HbT,data.adjRH.LPadSol.HbT);
% data.Contra.CBV = cat(2,data.adjLH.RPadSol.CBV,data.adjRH.LPadSol.CBV);
% data.Contra.cortMUA = cat(2,data.adjLH.RPadSol.cortMUA,data.adjRH.LPadSol.cortMUA);
% data.Contra.hipMUA = data.adjRH.RPadSol.hipMUA;
% data.Contra.cortGam = cat(2,data.adjLH.RPadSol.cortGam,data.adjRH.LPadSol.cortGam);
% data.Contra.hipGam = data.adjRH.RPadSol.hipGam;
% data.Contra.timeVector = cat(2,data.adjLH.RPadSol.timeVector,data.adjRH.LPadSol.timeVector);
% data.Contra.cortS = cat(3,data.adjLH.RPadSol.cortS,data.adjRH.LPadSol.cortS);
% data.Contra.cortS_Gam = cat(3,data.adjLH.RPadSol.cortS_Gam,data.adjRH.LPadSol.cortS_Gam);
% data.Contra.hipS = data.adjRH.RPadSol.hipS;
% data.Contra.hipS_Gam = data.adjRH.RPadSol.hipS_Gam;
% data.Contra.T = cat(2,data.adjLH.RPadSol.T,data.adjRH.LPadSol.T);
% data.Contra.F = cat(2,data.adjLH.RPadSol.F,data.adjRH.LPadSol.F);
% data.Ipsi.count = cat(2,data.adjLH.LPadSol.count,data.adjRH.RPadSol.count);
% data.Ipsi.HbT = cat(2,data.adjLH.LPadSol.HbT,data.adjRH.RPadSol.HbT);
% data.Ipsi.CBV = cat(2,data.adjLH.LPadSol.CBV,data.adjRH.RPadSol.CBV);
% data.Ipsi.cortMUA = cat(2,data.adjLH.LPadSol.cortMUA,data.adjRH.RPadSol.cortMUA);
% data.Ipsi.hipMUA = data.adjRH.LPadSol.hipMUA;
% data.Ipsi.cortGam = cat(2,data.adjLH.LPadSol.cortGam,data.adjRH.RPadSol.cortGam);
% data.Ipsi.hipGam = data.adjRH.LPadSol.hipGam;
% data.Ipsi.timeVector = cat(2,data.adjLH.LPadSol.timeVector,data.adjRH.RPadSol.timeVector);
% data.Ipsi.cortS = cat(3,data.adjLH.LPadSol.cortS,data.adjRH.RPadSol.cortS);
% data.Ipsi.cortS_Gam = cat(3,data.adjLH.LPadSol.cortS_Gam,data.adjRH.RPadSol.cortS_Gam);
% data.Ipsi.hipS = data.adjRH.LPadSol.hipS;
% data.Ipsi.hipS_Gam = data.adjRH.LPadSol.hipS_Gam;
% data.Ipsi.T = cat(2,data.adjLH.LPadSol.T,data.adjRH.RPadSol.T);
% data.Ipsi.F = cat(2,data.adjLH.LPadSol.F,data.adjRH.RPadSol.F);
% data.Auditory.count = cat(2,data.adjLH.AudSol.count,data.adjRH.AudSol.count);
% data.Auditory.HbT = cat(2,data.adjLH.AudSol.HbT,data.adjRH.AudSol.HbT);
% data.Auditory.CBV = cat(2,data.adjLH.AudSol.CBV,data.adjRH.AudSol.CBV);
% data.Auditory.cortMUA = cat(2,data.adjLH.AudSol.cortMUA,data.adjRH.AudSol.cortMUA);
% data.Auditory.hipMUA = data.adjRH.AudSol.hipMUA;
% data.Auditory.cortGam = cat(2,data.adjLH.AudSol.cortGam,data.adjRH.AudSol.cortGam);
% data.Auditory.hipGam = data.adjRH.AudSol.hipGam;
% data.Auditory.timeVector = cat(2,data.adjLH.AudSol.timeVector,data.adjRH.AudSol.timeVector);
% data.Auditory.cortS = cat(3,data.adjLH.AudSol.cortS,data.adjRH.AudSol.cortS);
% data.Auditory.cortS_Gam = cat(3,data.adjLH.AudSol.cortS_Gam,data.adjRH.AudSol.cortS_Gam);
% data.Auditory.hipS = data.adjRH.AudSol.hipS;
% data.Auditory.hipS_Gam = data.adjRH.AudSol.hipS_Gam;
% data.Auditory.T = cat(2,data.adjLH.AudSol.T,data.adjRH.AudSol.T);
% data.Auditory.F = cat(2,data.adjLH.AudSol.F,data.adjRH.AudSol.F);
%% take the averages of each field through the proper dimension
for ee = 1:length(treatments)
    treatment = treatments{1,ee};
    for gg = 1:length(dataTypes)
        dataType = dataTypes{1,gg};
        for ff = 1:length(compDataTypes)
            compDataType = compDataTypes{1,ff};
            data.(treatment).(dataType).(compDataType).mean_Count = mean(data.(treatment).(dataType).(compDataType).count,2);
            data.(treatment).(dataType).(compDataType).std_Count = std(data.(treatment).(dataType).(compDataType).count,0,2);
            data.(treatment).(dataType).(compDataType).mean_HbT = mean(data.(treatment).(dataType).(compDataType).HbT,2);
            data.(treatment).(dataType).(compDataType).std_HbT = std(data.(treatment).(dataType).(compDataType).HbT,0,2);
            data.(treatment).(dataType).(compDataType).mean_CBV = mean(data.(treatment).(dataType).(compDataType).CBV,2);
            data.(treatment).(dataType).(compDataType).std_CBV = std(data.(treatment).(dataType).(compDataType).CBV,0,2);
            data.(treatment).(dataType).(compDataType).mean_CortMUA = mean(data.(treatment).(dataType).(compDataType).cortMUA,2);
            data.(treatment).(dataType).(compDataType).std_CortMUA = std(data.(treatment).(dataType).(compDataType).cortMUA,0,2);
            data.(treatment).(dataType).(compDataType).mean_HipMUA = mean(data.(treatment).(dataType).(compDataType).hipMUA,2);
            data.(treatment).(dataType).(compDataType).std_HipMUA = std(data.(treatment).(dataType).(compDataType).hipMUA,0,2);
            data.(treatment).(dataType).(compDataType).mean_CortGam = mean(data.(treatment).(dataType).(compDataType).cortGam,2);
            data.(treatment).(dataType).(compDataType).std_CortGam = std(data.(treatment).(dataType).(compDataType).cortGam,0,2);
            data.(treatment).(dataType).(compDataType).mean_HipGam = mean(data.(treatment).(dataType).(compDataType).hipGam,2);
            data.(treatment).(dataType).(compDataType).std_HipGam = std(data.(treatment).(dataType).(compDataType).hipGam,0,2);
            data.(treatment).(dataType).(compDataType).mean_timeVector = mean(data.(treatment).(dataType).(compDataType).timeVector,2);
            data.(treatment).(dataType).(compDataType).mean_CortS = mean(data.(treatment).(dataType).(compDataType).cortS,3).*100;
            data.(treatment).(dataType).(compDataType).mean_CortS_Gam = mean(mean(mean(data.(treatment).(dataType).(compDataType).cortS_Gam.*100,2),1),3);
            data.(treatment).(dataType).(compDataType).std_CortS_Gam = std(mean(mean(data.(treatment).(dataType).(compDataType).cortS_Gam.*100,2),1),0,3);
            data.(treatment).(dataType).(compDataType).mean_HipS = mean(data.(treatment).(dataType).(compDataType).hipS,3).*100;
            data.(treatment).(dataType).(compDataType).mean_HipS_Gam = mean(mean(mean(data.(treatment).(dataType).(compDataType).hipS_Gam.*100,2),1),3);
            data.(treatment).(dataType).(compDataType).std_HipS_Gam = std(mean(mean(data.(treatment).(dataType).(compDataType).hipS_Gam.*100,2),1),0,3);
            data.(treatment).(dataType).(compDataType).mean_T = mean(data.(treatment).(dataType).(compDataType).T,2);
            data.(treatment).(dataType).(compDataType).mean_F = mean(data.(treatment).(dataType).(compDataType).F,2);
        end
    end
end
%% Fig. 1-S2
summaryFigure = figure('Name','Fig1-S2 (a-r)');
sgtitle('Figure 1-S2 (& 1d) - Turner et al. 2020')


ax1 = subplot(2,4,1);
imagesc(data.C57BL6J.adjLH.Contra.mean_T,data.C57BL6J.adjLH.Contra.mean_F,data.C57BL6J.adjLH.Contra.mean_CortS)
title('Control LH Contra stim cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-50,75])
axis square
axis xy
set(gca,'box','off')
%% [1-S2m] HbT contra stim
ax5 = subplot(2,4,5);
plot(data.C57BL6J.adjLH.Contra.mean_timeVector,data.C57BL6J.adjLH.Contra.mean_HbT,'color',colors('rich black'),'LineWidth',1)
hold on
plot(data.C57BL6J.adjLH.Contra.mean_timeVector,data.C57BL6J.adjLH.Contra.mean_HbT + data.C57BL6J.adjLH.Contra.std_HbT,'color',colors('battleship grey'),'LineWidth',0.5)
plot(data.C57BL6J.adjLH.Contra.mean_timeVector,data.C57BL6J.adjLH.Contra.mean_HbT - data.C57BL6J.adjLH.Contra.std_HbT,'color',colors('battleship grey'),'LineWidth',0.5)
title('Control LH Contra stim \Delta[HbT] (\muM)')
ylabel('\Delta[HbT] (\muM)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')

ax2 = subplot(2,4,2);
imagesc(data.C57BL6J.adjRH.Contra.mean_T,data.C57BL6J.adjRH.Contra.mean_F,data.C57BL6J.adjRH.Contra.mean_CortS)
title('Control RH Contra stim cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-50,75])
axis square
axis xy
set(gca,'box','off')
%% [1-S2m] HbT contra stim
ax6 = subplot(2,4,6);
plot(data.C57BL6J.adjRH.Contra.mean_timeVector,data.C57BL6J.adjRH.Contra.mean_HbT,'color',colors('rich black'),'LineWidth',1)
hold on
plot(data.C57BL6J.adjRH.Contra.mean_timeVector,data.C57BL6J.adjRH.Contra.mean_HbT + data.C57BL6J.adjRH.Contra.std_HbT,'color',colors('battleship grey'),'LineWidth',0.5)
plot(data.C57BL6J.adjRH.Contra.mean_timeVector,data.C57BL6J.adjRH.Contra.mean_HbT - data.C57BL6J.adjRH.Contra.std_HbT,'color',colors('battleship grey'),'LineWidth',0.5)
title('Control RH Contra stim \Delta[HbT] (\muM)')
ylabel('\Delta[HbT] (\muM)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')

ax3 = subplot(2,4,3);
imagesc(data.SSP_SAP.adjLH.Contra.mean_T,data.SSP_SAP.adjLH.Contra.mean_F,data.SSP_SAP.adjLH.Contra.mean_CortS)
title('SSP-SAP UnRx LH Contra stim cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-50,75])
axis square
axis xy
set(gca,'box','off')
%% [1-S2m] HbT contra stim
ax7 = subplot(2,4,7);
plot(data.SSP_SAP.adjLH.Contra.mean_timeVector,data.SSP_SAP.adjLH.Contra.mean_HbT,'color',colors('rich black'),'LineWidth',1)
hold on
plot(data.SSP_SAP.adjLH.Contra.mean_timeVector,data.SSP_SAP.adjLH.Contra.mean_HbT + data.SSP_SAP.adjLH.Contra.std_HbT,'color',colors('battleship grey'),'LineWidth',0.5)
plot(data.SSP_SAP.adjLH.Contra.mean_timeVector,data.SSP_SAP.adjLH.Contra.mean_HbT - data.SSP_SAP.adjLH.Contra.std_HbT,'color',colors('battleship grey'),'LineWidth',0.5)
title('SSP-SAP UnRx LH Contra stim \Delta[HbT] (\muM)')
ylabel('\Delta[HbT] (\muM)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')

ax4 = subplot(2,4,4);
imagesc(data.SSP_SAP.adjRH.Contra.mean_T,data.SSP_SAP.adjRH.Contra.mean_F,data.SSP_SAP.adjRH.Contra.mean_CortS)
title('SSP-SAP treated RH Contra stim cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-50,75])
axis square
axis xy
set(gca,'box','off')
%% [1-S2m] HbT contra stim
ax8 = subplot(2,4,8);
plot(data.SSP_SAP.adjRH.Contra.mean_timeVector,data.SSP_SAP.adjRH.Contra.mean_HbT,'color',colors('rich black'),'LineWidth',1)
hold on
plot(data.SSP_SAP.adjRH.Contra.mean_timeVector,data.SSP_SAP.adjRH.Contra.mean_HbT + data.SSP_SAP.adjRH.Contra.std_HbT,'color',colors('battleship grey'),'LineWidth',0.5)
plot(data.SSP_SAP.adjRH.Contra.mean_timeVector,data.SSP_SAP.adjRH.Contra.mean_HbT - data.SSP_SAP.adjRH.Contra.std_HbT,'color',colors('battleship grey'),'LineWidth',0.5)
title('SSP-SAP treated RH Contra stim \Delta[HbT] (\muM)')
ylabel('\Delta[HbT] (\muM)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')

% %% [1-S2a] cortical MUA contra stim
% ax1 = subplot(6,3,1);
% plot(data.Contra.mean_timeVector,data.Contra.mean_CortMUA,'color',colors('rich black'),'LineWidth',1);
% hold on
% plot(data.Contra.mean_timeVector,data.Contra.mean_CortMUA + data.Contra.std_CortMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.Contra.mean_timeVector,data.Contra.mean_CortMUA - data.Contra.std_CortMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1d,1-S2a] Contra stim cortical MUA')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% axis square
% set(gca,'box','off')
% ax1.TickLength = [0.03,0.03];
% %% [1-S2b] cortical MUA ispi stim
% ax2 = subplot(6,3,2);
% plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CortMUA,'color',colors('rich black'),'LineWidth',1)
% hold on
% plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CortMUA + data.Ipsi.std_CortMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CortMUA - data.Ipsi.std_CortMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S2b] Ipsi stim cortical MUA')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% axis square
% set(gca,'box','off')
% ax2.TickLength = [0.03,0.03];
% %% [1-S2c] cortical MUA auditory stim
% ax3 = subplot(6,3,3);
% plot(data.Auditory.mean_timeVector,data.Auditory.mean_CortMUA,'color',colors('rich black'),'LineWidth',1)
% hold on
% plot(data.Auditory.mean_timeVector,data.Auditory.mean_CortMUA + data.Auditory.std_CortMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.Auditory.mean_timeVector,data.Auditory.mean_CortMUA - data.Auditory.std_CortMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S2c] Aud stim cortical MUA')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% axis square
% set(gca,'box','off')
% ax3.TickLength = [0.03,0.03];
% %% [1-S2d] cortical LFP contra stim
% ax4 = subplot(6,3,4);
% imagesc(data.Contra.mean_T,data.Contra.mean_F,data.Contra.mean_CortS)
% title('[1d,1-S2d] Contra stim cortical LFP')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% c4 = colorbar;
% ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-50,75])
% axis square
% axis xy
% set(gca,'box','off')
% ax4.TickLength = [0.03,0.03];
% %% [1-S2e] cortical LFP ispi stim
% ax5 = subplot(6,3,5);
% imagesc(data.Ipsi.mean_T,data.Ipsi.mean_F,data.Ipsi.mean_CortS)
% title('[1-S2e] Ipsi stim cortical LFP')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% c5 = colorbar;
% ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-50,75])
% axis square
% axis xy
% set(gca,'box','off')
% ax5.TickLength = [0.03,0.03];
% %% [1-S2f] cortical LFP auditory stim
% ax6 = subplot(6,3,6);
% imagesc(data.Auditory.mean_T,data.Auditory.mean_F,data.Auditory.mean_CortS)
% title('[1-S2f] Aud stim cortical LFP')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% c6 = colorbar;
% ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-50,75])
% axis square
% axis xy
% set(gca,'box','off')
% ax6.TickLength = [0.03,0.03];
% %% [1-S2g] hippocampal MUA contra stim
% ax7 = subplot(6,3,7);
% plot(data.Contra.mean_timeVector,data.Contra.mean_HipMUA,'color',colors('rich black'),'LineWidth',1)
% hold on
% plot(data.Contra.mean_timeVector,data.Contra.mean_HipMUA + data.Contra.std_HipMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.Contra.mean_timeVector,data.Contra.mean_HipMUA - data.Contra.std_HipMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S2g] Contra stim hippocampal MUA')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% axis square
% set(gca,'box','off')
% ax7.TickLength = [0.03,0.03];
% %% [1-S2h] hippocampal MUA ispi stim
% ax8 = subplot(6,3,8);
% plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HipMUA,'color',colors('rich black'),'LineWidth',1)
% hold on
% plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HipMUA + data.Ipsi.std_HipMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HipMUA - data.Ipsi.std_HipMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S2h] Ipsi stim hippocampal MUA')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% axis square
% set(gca,'box','off')
% ax8.TickLength = [0.03,0.03];
% %% [1-S2i] hippocampal MUA auditory stim
% ax9 = subplot(6,3,9);
% plot(data.Auditory.mean_timeVector,data.Auditory.mean_HipMUA,'color',colors('rich black'),'LineWidth',1)
% hold on
% plot(data.Auditory.mean_timeVector,data.Auditory.mean_HipMUA + data.Auditory.std_HipMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.Auditory.mean_timeVector,data.Auditory.mean_HipMUA - data.Auditory.std_HipMUA,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S2i] Aud stim hippocampal MUA')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% axis square
% set(gca,'box','off')
% ax9.TickLength = [0.03,0.03];
% %% [1-S2j] hippocampal LFP contra stim
% ax10 = subplot(6,3,10);
% imagesc(data.Contra.mean_T,data.Contra.mean_F,data.Contra.mean_HipS)
% title('[1-S2j] Contra stim hippocampal LFP')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% c10 = colorbar;
% ylabel(c10,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-50,75])
% axis square
% axis xy
% set(gca,'box','off')
% ax10.TickLength = [0.03,0.03];
% %% [1-S2j] hippocampal LFP ispi stim
% ax11 = subplot(6,3,11);
% imagesc(data.Ipsi.mean_T,data.Ipsi.mean_F,data.Ipsi.mean_HipS)
% title('[1-S2j] Ipsi stim hippocampal LFP')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% c11 = colorbar;
% ylabel(c11,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-50,75])
% axis square
% axis xy
% set(gca,'box','off')
% ax11.TickLength = [0.03,0.03];
% %% [1-S2l] hippocampal LFP auditory stim
% ax12 = subplot(6,3,12);
% imagesc(data.Auditory.mean_T,data.Auditory.mean_F,data.Auditory.mean_HipS)
% title('[1-S2l] Aud stim hippocampal LFP')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% c12 = colorbar;
% ylabel(c12,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-50,75])
% axis square
% axis xy
% set(gca,'box','off')
% ax12.TickLength = [0.03,0.03];
% %% [1-S2m] HbT contra stim
% ax13 = subplot(6,3,13);
% plot(data.Contra.mean_timeVector,data.Contra.mean_HbT,'color',colors('rich black'),'LineWidth',1)
% hold on
% plot(data.Contra.mean_timeVector,data.Contra.mean_HbT + data.Contra.std_HbT,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.Contra.mean_timeVector,data.Contra.mean_HbT - data.Contra.std_HbT,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1d,1-S2m] Contra stim \Delta[HbT] (\muM)')
% ylabel('\Delta[HbT] (\muM)')
% xlabel('Peri-stimulus time (s)')
% axis square
% set(gca,'box','off')
% ax13.TickLength = [0.03,0.03];
% %% [1-S2n] HbT ispi stim
% ax14 = subplot(6,3,14);
% plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HbT,'color',colors('rich black'),'LineWidth',1)
% hold on
% plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HbT + data.Ipsi.std_HbT,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HbT - data.Ipsi.std_HbT,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S2n] Ipsi stim \Delta[HbT] (\muM)')
% ylabel('\Delta[HbT] (\muM)')
% xlabel('Peri-stimulus time (s)')
% axis square
% set(gca,'box','off')
% ax14.TickLength = [0.03,0.03];
% %% [1-S2o] HbT auditory stim
% ax15 = subplot(6,3,15);
% plot(data.Auditory.mean_timeVector,data.Auditory.mean_HbT,'color',colors('rich black'),'LineWidth',1)
% hold on
% plot(data.Auditory.mean_timeVector,data.Auditory.mean_HbT + data.Auditory.std_HbT,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.Auditory.mean_timeVector,data.Auditory.mean_HbT - data.Auditory.std_HbT,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S2o] Aud stim \Delta[HbT] (\muM)')
% ylabel('\Delta[HbT] (\muM)')
% xlabel('Peri-stimulus time (s)')
% axis square
% set(gca,'box','off')
% ax15.TickLength = [0.03,0.03];
% %% [1-S2p] refl contra stim
% ax16 = subplot(6,3,16);
% plot(data.Contra.mean_timeVector,data.Contra.mean_CBV,'color',colors('rich black'),'LineWidth',1)
% hold on
% plot(data.Contra.mean_timeVector,data.Contra.mean_CBV + data.Contra.std_CBV,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.Contra.mean_timeVector,data.Contra.mean_CBV - data.Contra.std_CBV,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S2p] Contra stim reflectance')
% ylabel('\DeltaR/R (%)')
% xlabel('Peri-stimulus time (s)')
% axis square
% set(gca,'box','off')
% ax16.TickLength = [0.03,0.03];
% %% [1-S2q] refl ispi stim
% ax17 = subplot(6,3,17);
% plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CBV,'color',colors('rich black'),'LineWidth',1)
% hold on
% plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CBV + data.Ipsi.std_CBV,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CBV - data.Ipsi.std_CBV,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S2q] Ipsi stim reflectance')
% ylabel('\DeltaR/R (%)')
% xlabel('Peri-stimulus time (s)')
% axis square
% set(gca,'box','off')
% ax17.TickLength = [0.03,0.03];
% %% [1-S2r] refl auditory stim
% ax18 = subplot(6,3,18);
% plot(data.Auditory.mean_timeVector,data.Auditory.mean_CBV,'color',colors('rich black'),'LineWidth',1)
% hold on
% plot(data.Auditory.mean_timeVector,data.Auditory.mean_CBV + data.Auditory.std_CBV,'color',colors('battleship grey'),'LineWidth',0.5)
% plot(data.Auditory.mean_timeVector,data.Auditory.mean_CBV - data.Auditory.std_CBV,'color',colors('battleship grey'),'LineWidth',0.5)
% title('[1-S2r] Aud stim reflectance')
% ylabel('\DeltaR/R (%)')
% xlabel('Peri-stimulus time (s)')
% axis square
% set(gca,'box','off')
% ax18.TickLength = [0.03,0.03];
% %% adjust and link axes
linkaxes([ax5,ax6,ax7,ax8],'xy')
% linkaxes([ax4,ax5,ax6,ax10,ax11,ax12],'xy')
% linkaxes([ax13,ax14,ax15],'xy')
% linkaxes([ax16,ax17,ax18],'xy')
% ax1Pos = get(ax1,'position');
% ax2Pos = get(ax2,'position');
% ax3Pos = get(ax3,'position');
% ax4Pos = get(ax4,'position');
% ax5Pos = get(ax5,'position');
% ax6Pos = get(ax6,'position');
% ax10Pos = get(ax10,'position');
% ax11Pos = get(ax11,'position');
% ax12Pos = get(ax12,'position');
% ax4Pos(3:4) = ax1Pos(3:4);
% ax5Pos(3:4) = ax2Pos(3:4);
% ax6Pos(3:4) = ax3Pos(3:4);
% ax10Pos(3:4) = ax1Pos(3:4);
% ax11Pos(3:4) = ax2Pos(3:4);
% ax12Pos(3:4) = ax3Pos(3:4);
% set(ax4,'position',ax4Pos);
% set(ax5,'position',ax5Pos);
% set(ax6,'position',ax6Pos);
% set(ax10,'position',ax10Pos);
% set(ax11,'position',ax11Pos);
% set(ax12,'position',ax12Pos);
% %% save figure(s)
% if strcmp(saveFigs,'y') == true
%     dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
%     if ~exist(dirpath,'dir')
%         mkdir(dirpath);
%     end
%     savefig(summaryFigure,[dirpath 'Fig1-S2']);
%     set(summaryFigure,'PaperPositionMode','auto');
%     print('-painters','-dpdf','-fillpage',[dirpath 'Fig1-S2'])
%     %% text diary
%     diaryFile = [dirpath 'Fig1-S2_Statistics.txt'];
%     if exist(diaryFile,'file') == 2
%         delete(diaryFile)
%     end
%     diary(diaryFile)
%     diary on
%     % text values
%     disp('======================================================================================================================')
%     disp('[1-S2] Text values for gamma/HbT/reflectance changes')
%     disp('======================================================================================================================')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     % cortical MUA/LFP
%     [~,index] = max(data.Contra.mean_CortMUA);
%     disp(['Contra stim Cort gamma MUA P/P (%): ' num2str(round(data.Contra.mean_CortMUA(index),1)) ' +/- ' num2str(round(data.Contra.std_CortMUA(index),1))]); disp(' ')
%     [~,index] = max(data.Ipsi.mean_CortMUA);
%     disp(['Ipsil stim Cort gamma MUA P/P (%): ' num2str(round(data.Ipsi.mean_CortMUA(index),1)) ' +/- ' num2str(round(data.Ipsi.std_CortMUA(index),1))]); disp(' ')
%     [~,index] = max(data.Auditory.mean_CortMUA);
%     disp(['Audit stim Cort gamma MUA P/P (%): ' num2str(round(data.Auditory.mean_CortMUA(index),1)) ' +/- ' num2str(round(data.Auditory.std_CortMUA(index),1))]); disp(' ')
%     % cortical LFP
%     disp(['Contra stim Cort gamma LFP P/P (%): ' num2str(round(data.Contra.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.Contra.std_CortS_Gam,1))]); disp(' ')
%     disp(['Ipsil stim Cort gamma LFP P/P (%): ' num2str(round(data.Ipsi.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.Ipsi.std_CortS_Gam,1))]); disp(' ')
%     disp(['Audit stim Cort gamma LFP P/P (%): ' num2str(round(data.Auditory.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.Auditory.std_CortS_Gam,1))]); disp(' ')
%     % hippocampal MUA
%     [~,index] = max(data.Contra.mean_HipMUA);
%     disp(['Contra stim Hip gamma MUA P/P (%): ' num2str(round(data.Contra.mean_HipMUA(index),1)) ' +/- ' num2str(round(data.Contra.std_HipMUA(index),1))]); disp(' ')
%     [~,index] = max(data.Ipsi.mean_HipMUA);
%     disp(['Ipsil stim Hip gamma MUA P/P (%): ' num2str(round(data.Ipsi.mean_HipMUA(index),1)) ' +/- ' num2str(round(data.Ipsi.std_HipMUA(index),1))]); disp(' ')
%     [~,index] = max(data.Auditory.mean_HipMUA);
%     disp(['Audit stim Hip gamma MUA P/P (%): ' num2str(round(data.Auditory.mean_HipMUA(index),1)) ' +/- ' num2str(round(data.Auditory.std_HipMUA(index),1))]); disp(' ')
%     % hippocampal LFP
%     disp(['Contra stim Hip gamma LFP P/P (%): ' num2str(round(data.Contra.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Contra.std_HipS_Gam,1))]); disp(' ')
%     disp(['Ipsil stim Hip gamma LFP P/P (%): ' num2str(round(data.Ipsi.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Ipsi.std_HipS_Gam,1))]); disp(' ')
%     disp(['Auditory stim Hip gamma LFP P/P (%): ' num2str(round(data.Auditory.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Auditory.std_HipS_Gam,1))]); disp(' ')
%     % HbT
%     [~,index] = max(data.Contra.mean_HbT);
%     disp(['Contra stim [HbT] (uM): ' num2str(round(data.Contra.mean_HbT(index),1)) ' +/- ' num2str(round(data.Contra.std_HbT(index),1))]); disp(' ')
%     [~,index] = max(data.Ipsi.mean_HbT);
%     disp(['Ipsil stim [HbT] (uM): ' num2str(round(data.Ipsi.mean_HbT(index),1)) ' +/- ' num2str(round(data.Ipsi.std_HbT(index),1))]); disp(' ')
%     [~,index] = max(data.Auditory.mean_HbT);
%     disp(['Audit stim [HbT] (uM): ' num2str(round(data.Auditory.mean_HbT(index),1)) ' +/- ' num2str(round(data.Auditory.std_HbT(index),1))]); disp(' ')
%     % R/R
%     [~,index] = min(data.Contra.mean_CBV);
%     disp(['Contra stim refl R/R (%): ' num2str(round(data.Contra.mean_CBV(index),1)) ' +/- ' num2str(round(data.Contra.std_CBV(index),1))]); disp(' ')
%     [~,index] = min(data.Ipsi.mean_CBV);
%     disp(['Ipsil stim refl R/R (%): ' num2str(round(data.Ipsi.mean_CBV(index),1)) ' +/- ' num2str(round(data.Ipsi.std_CBV(index),1))]); disp(' ')
%     [~,index] = min(data.Auditory.mean_CBV);
%     disp(['Audit stim refl R/R (%): ' num2str(round(data.Auditory.mean_CBV(index),1)) ' +/- ' num2str(round(data.Auditory.std_CBV(index),1))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     diary off
% end

end
