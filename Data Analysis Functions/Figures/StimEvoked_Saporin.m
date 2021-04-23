function [AnalysisResults] = StimEvoked_Saporin(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: 
%________________________________________________________________________________________________________________________

%% set-up
animalIDs = {'T141','T155','T156','T157','T142','T144','T159','T172','T150','T165','T166','T177','T179','T186','T187','T188','T189'};
C57BL6J_IDs = {'T141','T155','T156','T157','T186','T187','T188','T189'};
SSP_SAP_IDs = {'T142','T144','T159','T172'};
Blank_SAP_IDs = {'T150','T165','T166','T177','T179'};
solenoidNames = {'LPadSol','RPadSol','AudSol'};
compDataTypes = {'Ipsi','Contra','Auditory'};
hemispheres = {'adjLH','adjRH'};
treatments = {'C57BL6J','SSP_SAP','Blank_SAP'};
data = [];
cortVariables = {'HbT','CBV','cortMUA','cortGam','cortS','cortS_Gam','cortT','cortF','timeVector','count'};
hipVariables = {'hipMUA','hipGam','hipS','hipS_Gam','hipT','hipF','timeVector'};
%% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    % recognize treatment based on animal group
    if ismember(animalIDs{1,aa},C57BL6J_IDs) == true
        treatment = 'C57BL6J';
    elseif ismember(animalIDs{1,aa},SSP_SAP_IDs) == true
        treatment = 'SSP_SAP';
    elseif ismember(animalIDs{1,aa},Blank_SAP_IDs) == true
        treatment = 'Blank_SAP';
    end
    for bb = 1:length(solenoidNames)
        % left, right hemishpere hemo & neural data
        for cc = 1:length(hemispheres)
            % pre-allocate necessary variable fields
            data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).dummCheck = 1;
            for dd = 1:length(cortVariables)
                if isfield(data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}),cortVariables{1,dd}) == false
                    data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).(cortVariables{1,dd}) = [];
                end
            end
            data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).HbT = cat(1,data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).HbT,AnalysisResults.(animalIDs{1,aa}).EvokedAvgs.Stim.(hemispheres{1,cc}).(solenoidNames{1,bb}).CBV_HbT.HbT);
            data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).CBV = cat(1,data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).CBV,AnalysisResults.(animalIDs{1,aa}).EvokedAvgs.Stim.(hemispheres{1,cc}).(solenoidNames{1,bb}).CBV.CBV);
            data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).cortMUA = cat(1,data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).cortMUA,AnalysisResults.(animalIDs{1,aa}).EvokedAvgs.Stim.(hemispheres{1,cc}).(solenoidNames{1,bb}).MUA.corticalData);
            data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).cortGam = cat(1,data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).cortGam,AnalysisResults.(animalIDs{1,aa}).EvokedAvgs.Stim.(hemispheres{1,cc}).(solenoidNames{1,bb}).Gam.corticalData);
            data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).cortS = cat(3,data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).cortS,AnalysisResults.(animalIDs{1,aa}).EvokedAvgs.Stim.(hemispheres{1,cc}).(solenoidNames{1,bb}).LFP.corticalS);
            data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).cortS_Gam = cat(3,data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).cortS_Gam,AnalysisResults.(animalIDs{1,aa}).EvokedAvgs.Stim.(hemispheres{1,cc}).(solenoidNames{1,bb}).LFP.corticalS(49:end,20:23));
            data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).cortT = cat(1,data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).cortT,AnalysisResults.(animalIDs{1,aa}).EvokedAvgs.Stim.(hemispheres{1,cc}).(solenoidNames{1,bb}).LFP.T);
            data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).cortF = cat(1,data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).cortF,AnalysisResults.(animalIDs{1,aa}).EvokedAvgs.Stim.(hemispheres{1,cc}).(solenoidNames{1,bb}).LFP.F);
            data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).timeVector = cat(1,data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).timeVector,AnalysisResults.(animalIDs{1,aa}).EvokedAvgs.Stim.(hemispheres{1,cc}).(solenoidNames{1,bb}).timeVector);
            data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).count = cat(1,data.(treatment).(solenoidNames{1,bb}).(hemispheres{1,cc}).count,AnalysisResults.(animalIDs{1,aa}).EvokedAvgs.Stim.(hemispheres{1,cc}).(solenoidNames{1,bb}).count);
            % hippocampal neural data - preallocate necessary variable fields
            data.(treatment).(solenoidNames{1,bb}).Hip.dummyCheck = 1;
            for ee = 1:length(hipVariables)
                if isfield(data.(treatment).(solenoidNames{1,bb}).Hip,hipVariables{1,ee}) == false
                    data.(treatment).(solenoidNames{1,bb}).Hip.(hipVariables{1,ee}) = [];
                end
            end
            data.(treatment).(solenoidNames{1,bb}).Hip.hipMUA = cat(1,data.(treatment).(solenoidNames{1,bb}).Hip.hipMUA,AnalysisResults.(animalIDs{1,aa}).EvokedAvgs.Stim.adjLH.(solenoidNames{1,bb}).MUA.hippocampalData);
            data.(treatment).(solenoidNames{1,bb}).Hip.hipGam = cat(1,data.(treatment).(solenoidNames{1,bb}).Hip.hipGam,AnalysisResults.(animalIDs{1,aa}).EvokedAvgs.Stim.adjLH.(solenoidNames{1,bb}).Gam.hippocampalData);
            data.(treatment).(solenoidNames{1,bb}).Hip.hipS = cat(3,data.(treatment).(solenoidNames{1,bb}).Hip.hipS,AnalysisResults.(animalIDs{1,aa}).EvokedAvgs.Stim.adjLH.(solenoidNames{1,bb}).LFP.hippocampalS);
            data.(treatment).(solenoidNames{1,bb}).Hip.hipS_Gam = cat(3,data.(treatment).(solenoidNames{1,bb}).Hip.hipS_Gam,AnalysisResults.(animalIDs{1,aa}).EvokedAvgs.Stim.adjLH.(solenoidNames{1,bb}).LFP.hippocampalS(49:end,20:23));
            data.(treatment).(solenoidNames{1,bb}).Hip.hipT = cat(1,data.(treatment).(solenoidNames{1,bb}).Hip.hipT,AnalysisResults.(animalIDs{1,aa}).EvokedAvgs.Stim.adjLH.(solenoidNames{1,bb}).LFP.T);
            data.(treatment).(solenoidNames{1,bb}).Hip.hipF = cat(1,data.(treatment).(solenoidNames{1,bb}).Hip.hipF,AnalysisResults.(animalIDs{1,aa}).EvokedAvgs.Stim.adjLH.(solenoidNames{1,bb}).LFP.F);
            data.(treatment).(solenoidNames{1,bb}).Hip.timeVector = cat(1,data.(treatment).(solenoidNames{1,bb}).Hip.timeVector,AnalysisResults.(animalIDs{1,aa}).EvokedAvgs.Stim.adjLH.(solenoidNames{1,bb}).timeVector);
        end
    end
end
%% concatenate the data from the contra and ipsi data
% contra
for ff = 1:length(treatments)
    for gg = 1:length(cortVariables)
        data.(treatments{1,ff}).Contra.adjLH.(cortVariables{1,gg}) = data.(treatments{1,ff}).RPadSol.adjLH.(cortVariables{1,gg});
        data.(treatments{1,ff}).Contra.adjRH.(cortVariables{1,gg}) = data.(treatments{1,ff}).LPadSol.adjRH.(cortVariables{1,gg});
    end
    % hip
    for hh = 1:length(hipVariables)
        data.(treatments{1,ff}).Contra.Hip.(hipVariables{1,hh}) = data.(treatments{1,ff}).RPadSol.Hip.(hipVariables{1,hh});
    end
end
% Ipsi
for ff = 1:length(treatments)
    for gg = 1:length(cortVariables)
        data.(treatments{1,ff}).Ipsi.adjLH.(cortVariables{1,gg}) = data.(treatments{1,ff}).RPadSol.adjRH.(cortVariables{1,gg});
        data.(treatments{1,ff}).Ipsi.adjRH.(cortVariables{1,gg}) = data.(treatments{1,ff}).LPadSol.adjLH.(cortVariables{1,gg});
    end
    % hip
    for hh = 1:length(hipVariables)
        data.(treatments{1,ff}).Ipsi.Hip.(hipVariables{1,hh}) = data.(treatments{1,ff}).LPadSol.Hip.(hipVariables{1,hh});
    end
end
% auditory
for ff = 1:length(treatments)
    for gg = 1:length(cortVariables)
        data.(treatments{1,ff}).Auditory.adjLH.(cortVariables{1,gg}) = data.(treatments{1,ff}).AudSol.adjLH.(cortVariables{1,gg});
        data.(treatments{1,ff}).Auditory.adjRH.(cortVariables{1,gg}) = data.(treatments{1,ff}).AudSol.adjRH.(cortVariables{1,gg});
    end
    % hip
    for hh = 1:length(hipVariables)
        data.(treatments{1,ff}).Auditory.Hip.(hipVariables{1,hh}) = data.(treatments{1,ff}).AudSol.Hip.(hipVariables{1,hh});
    end
end
%% take the averages of each field through the proper dimension
for ee = 1:length(treatments)
    treatment = treatments{1,ee};
    for gg = 1:length(hemispheres)
        hemisphere = hemispheres{1,gg};
        for ff = 1:length(compDataTypes)
            compDataType = compDataTypes{1,ff};
            data.(treatment).(compDataType).(hemisphere).meanHbT = mean(data.(treatment).(compDataType).(hemisphere).HbT,1);
            data.(treatment).(compDataType).(hemisphere).stdHbT = std(data.(treatment).(compDataType).(hemisphere).HbT,0,1);
            data.(treatment).(compDataType).(hemisphere).meanCBV = mean(data.(treatment).(compDataType).(hemisphere).CBV,1);
            data.(treatment).(compDataType).(hemisphere).stdCBV = std(data.(treatment).(compDataType).(hemisphere).CBV,0,1);
            data.(treatment).(compDataType).(hemisphere).meanCortMUA = mean(data.(treatment).(compDataType).(hemisphere).cortMUA,1);
            data.(treatment).(compDataType).(hemisphere).stdCortMUA = std(data.(treatment).(compDataType).(hemisphere).cortMUA,0,1);
            data.(treatment).(compDataType).(hemisphere).meanCortGam = mean(data.(treatment).(compDataType).(hemisphere).cortGam,1);
            data.(treatment).(compDataType).(hemisphere).stdCortGam = std(data.(treatment).(compDataType).(hemisphere).cortGam,0,1);
            data.(treatment).(compDataType).(hemisphere).meanCortS = mean(data.(treatment).(compDataType).(hemisphere).cortS,3).*100;
            data.(treatment).(compDataType).(hemisphere).meanCortS_Gam = mean(mean(mean(data.(treatment).(compDataType).(hemisphere).cortS_Gam.*100,1),2),3);
            data.(treatment).(compDataType).(hemisphere).stdCortS_Gam = std(mean(mean(data.(treatment).(compDataType).(hemisphere).cortS_Gam.*100,1),2),0,3);
            data.(treatment).(compDataType).(hemisphere).meanCortT = mean(data.(treatment).(compDataType).(hemisphere).cortT,1);
            data.(treatment).(compDataType).(hemisphere).meanCortF = mean(data.(treatment).(compDataType).(hemisphere).cortF,1);
            data.(treatment).(compDataType).(hemisphere).meanTimeVector = mean(data.(treatment).(compDataType).(hemisphere).timeVector,1);
            data.(treatment).(compDataType).(hemisphere).meanCount = mean(data.(treatment).(compDataType).(hemisphere).count,1);
            data.(treatment).(compDataType).(hemisphere).stdCount = std(data.(treatment).(compDataType).(hemisphere).count,0,1);
            % hip
            data.(treatment).(compDataType).Hip.meanHipMUA = mean(data.(treatment).(compDataType).Hip.hipMUA,1);
            data.(treatment).(compDataType).Hip.stdHipMUA = std(data.(treatment).(compDataType).Hip.hipMUA,0,1);
            data.(treatment).(compDataType).Hip.meanHipGam = mean(data.(treatment).(compDataType).Hip.hipGam,1);
            data.(treatment).(compDataType).Hip.stdHipGam = std(data.(treatment).(compDataType).Hip.hipGam,0,1);
            data.(treatment).(compDataType).Hip.meanHipS = mean(data.(treatment).(compDataType).Hip.hipS,3).*100;
            data.(treatment).(compDataType).Hip.meanHipS_Gam = mean(mean(mean(data.(treatment).(compDataType).Hip.hipS_Gam.*100,1),2),3);
            data.(treatment).(compDataType).Hip.stdHipS_Gam = std(mean(mean(data.(treatment).(compDataType).Hip.hipS_Gam.*100,1),2),0,3);
            data.(treatment).(compDataType).Hip.meanT = mean(data.(treatment).(compDataType).Hip.hipT,1);
            data.(treatment).(compDataType).Hip.meanF = mean(data.(treatment).(compDataType).Hip.hipF,1);
            data.(treatment).(compDataType).Hip.meanTimeVector = mean(data.(treatment).(compDataType).Hip.timeVector,1); 
        end
    end
end
%% average stim-evoked figures
summaryFigure1 = figure;
sgtitle('Stimulus-evoked \DeltaHbT repsonses')
%% LH contra stim
ax1 = subplot(3,2,1);
% C57BL6Js
p1 = plot(data.C57BL6J.Contra.adjLH.meanTimeVector,data.C57BL6J.Contra.adjLH.meanHbT,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.C57BL6J.Contra.adjLH.meanTimeVector,data.C57BL6J.Contra.adjLH.meanHbT + data.C57BL6J.Contra.adjLH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.C57BL6J.Contra.adjLH.meanTimeVector,data.C57BL6J.Contra.adjLH.meanHbT - data.C57BL6J.Contra.adjLH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
p2 = plot(data.Blank_SAP.Contra.adjLH.meanTimeVector,data.Blank_SAP.Contra.adjLH.meanHbT,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.Contra.adjLH.meanTimeVector,data.Blank_SAP.Contra.adjLH.meanHbT + data.Blank_SAP.Contra.adjLH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.Contra.adjLH.meanTimeVector,data.Blank_SAP.Contra.adjLH.meanHbT - data.Blank_SAP.Contra.adjLH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
p3 = plot(data.SSP_SAP.Contra.adjLH.meanTimeVector,data.SSP_SAP.Contra.adjLH.meanHbT,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.Contra.adjLH.meanTimeVector,data.SSP_SAP.Contra.adjLH.meanHbT + data.SSP_SAP.Contra.adjLH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.Contra.adjLH.meanTimeVector,data.SSP_SAP.Contra.adjLH.meanHbT - data.SSP_SAP.Contra.adjLH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
title('LH (UnRx) \DeltaHbT - Contra Stim')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH contra stim
ax2 = subplot(3,2,2);
% C57BL6Js
plot(data.C57BL6J.Contra.adjRH.meanTimeVector,data.C57BL6J.Contra.adjRH.meanHbT,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.C57BL6J.Contra.adjRH.meanTimeVector,data.C57BL6J.Contra.adjRH.meanHbT + data.C57BL6J.Contra.adjRH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.C57BL6J.Contra.adjRH.meanTimeVector,data.C57BL6J.Contra.adjRH.meanHbT - data.C57BL6J.Contra.adjRH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
plot(data.Blank_SAP.Contra.adjRH.meanTimeVector,data.Blank_SAP.Contra.adjRH.meanHbT,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.Contra.adjRH.meanTimeVector,data.Blank_SAP.Contra.adjRH.meanHbT + data.Blank_SAP.Contra.adjRH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.Contra.adjRH.meanTimeVector,data.Blank_SAP.Contra.adjRH.meanHbT - data.Blank_SAP.Contra.adjRH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.Contra.adjRH.meanTimeVector,data.SSP_SAP.Contra.adjRH.meanHbT,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.Contra.adjRH.meanTimeVector,data.SSP_SAP.Contra.adjRH.meanHbT + data.SSP_SAP.Contra.adjRH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.Contra.adjRH.meanTimeVector,data.SSP_SAP.Contra.adjRH.meanHbT - data.SSP_SAP.Contra.adjRH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
title('RH (Rx) \DeltaHbT - Contra Stim')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% LH ipsi stim
ax3 = subplot(3,2,3);
% C57BL6Js
plot(data.C57BL6J.Ipsi.adjLH.meanTimeVector,data.C57BL6J.Ipsi.adjLH.meanHbT,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.C57BL6J.Ipsi.adjLH.meanTimeVector,data.C57BL6J.Ipsi.adjLH.meanHbT + data.C57BL6J.Ipsi.adjLH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.C57BL6J.Ipsi.adjLH.meanTimeVector,data.C57BL6J.Ipsi.adjLH.meanHbT - data.C57BL6J.Ipsi.adjLH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
plot(data.Blank_SAP.Ipsi.adjLH.meanTimeVector,data.Blank_SAP.Ipsi.adjLH.meanHbT,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.Ipsi.adjLH.meanTimeVector,data.Blank_SAP.Ipsi.adjLH.meanHbT + data.Blank_SAP.Ipsi.adjLH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.Ipsi.adjLH.meanTimeVector,data.Blank_SAP.Ipsi.adjLH.meanHbT - data.Blank_SAP.Ipsi.adjLH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.Ipsi.adjLH.meanTimeVector,data.SSP_SAP.Ipsi.adjLH.meanHbT,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.Ipsi.adjLH.meanTimeVector,data.SSP_SAP.Ipsi.adjLH.meanHbT + data.SSP_SAP.Ipsi.adjLH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.Ipsi.adjLH.meanTimeVector,data.SSP_SAP.Ipsi.adjLH.meanHbT - data.SSP_SAP.Ipsi.adjLH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
title('LH (UnRx) \DeltaHbT - Ipsi Stim')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH ipsi stim
ax4 = subplot(3,2,4);
% C57BL6Js
plot(data.C57BL6J.Ipsi.adjRH.meanTimeVector,data.C57BL6J.Ipsi.adjRH.meanHbT,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.C57BL6J.Ipsi.adjRH.meanTimeVector,data.C57BL6J.Ipsi.adjRH.meanHbT + data.C57BL6J.Ipsi.adjRH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.C57BL6J.Ipsi.adjRH.meanTimeVector,data.C57BL6J.Ipsi.adjRH.meanHbT - data.C57BL6J.Ipsi.adjRH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
plot(data.Blank_SAP.Ipsi.adjRH.meanTimeVector,data.Blank_SAP.Ipsi.adjRH.meanHbT,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.Ipsi.adjRH.meanTimeVector,data.Blank_SAP.Ipsi.adjRH.meanHbT + data.Blank_SAP.Ipsi.adjRH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.Ipsi.adjRH.meanTimeVector,data.Blank_SAP.Ipsi.adjRH.meanHbT - data.Blank_SAP.Ipsi.adjRH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.Ipsi.adjRH.meanTimeVector,data.SSP_SAP.Ipsi.adjRH.meanHbT,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.Ipsi.adjRH.meanTimeVector,data.SSP_SAP.Ipsi.adjRH.meanHbT + data.SSP_SAP.Ipsi.adjRH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.Ipsi.adjRH.meanTimeVector,data.SSP_SAP.Ipsi.adjRH.meanHbT - data.SSP_SAP.Ipsi.adjRH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
title('RH (Rx) \DeltaHbT - Ipsi Stim')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% LH auditory stim
ax5 = subplot(3,2,5);
% C57BL6Js
plot(data.C57BL6J.Auditory.adjLH.meanTimeVector,data.C57BL6J.Auditory.adjLH.meanHbT,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.C57BL6J.Auditory.adjLH.meanTimeVector,data.C57BL6J.Auditory.adjLH.meanHbT + data.C57BL6J.Auditory.adjLH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.C57BL6J.Auditory.adjLH.meanTimeVector,data.C57BL6J.Auditory.adjLH.meanHbT - data.C57BL6J.Auditory.adjLH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
plot(data.Blank_SAP.Auditory.adjLH.meanTimeVector,data.Blank_SAP.Auditory.adjLH.meanHbT,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.Auditory.adjLH.meanTimeVector,data.Blank_SAP.Auditory.adjLH.meanHbT + data.Blank_SAP.Auditory.adjLH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.Auditory.adjLH.meanTimeVector,data.Blank_SAP.Auditory.adjLH.meanHbT - data.Blank_SAP.Auditory.adjLH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.Auditory.adjLH.meanTimeVector,data.SSP_SAP.Auditory.adjLH.meanHbT,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.Auditory.adjLH.meanTimeVector,data.SSP_SAP.Auditory.adjLH.meanHbT + data.SSP_SAP.Auditory.adjLH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.Auditory.adjLH.meanTimeVector,data.SSP_SAP.Auditory.adjLH.meanHbT - data.SSP_SAP.Auditory.adjLH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
title('LH (UnRx) \DeltaHbT - Auditory Stim')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH auditory stim
ax6 = subplot(3,2,6);
% C57BL6Js
plot(data.C57BL6J.Auditory.adjRH.meanTimeVector,data.C57BL6J.Auditory.adjRH.meanHbT,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.C57BL6J.Auditory.adjRH.meanTimeVector,data.C57BL6J.Auditory.adjRH.meanHbT + data.C57BL6J.Auditory.adjRH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.C57BL6J.Auditory.adjRH.meanTimeVector,data.C57BL6J.Auditory.adjRH.meanHbT - data.C57BL6J.Auditory.adjRH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
plot(data.Blank_SAP.Auditory.adjRH.meanTimeVector,data.Blank_SAP.Auditory.adjRH.meanHbT,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.Auditory.adjRH.meanTimeVector,data.Blank_SAP.Auditory.adjRH.meanHbT + data.Blank_SAP.Auditory.adjRH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.Auditory.adjRH.meanTimeVector,data.Blank_SAP.Auditory.adjRH.meanHbT - data.Blank_SAP.Auditory.adjRH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.Auditory.adjRH.meanTimeVector,data.SSP_SAP.Auditory.adjRH.meanHbT,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.Auditory.adjRH.meanTimeVector,data.SSP_SAP.Auditory.adjRH.meanHbT + data.SSP_SAP.Auditory.adjRH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.Auditory.adjRH.meanTimeVector,data.SSP_SAP.Auditory.adjRH.meanHbT - data.SSP_SAP.Auditory.adjRH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
title('RH (Rx) \DeltaHbT - Auditory Stim')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% figure characteristics
linkaxes([ax1,ax2],'xy')
linkaxes([ax3,ax4],'xy')
linkaxes([ax5,ax6],'xy')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure1,[dirpath 'Stim_Evoked_HbT']);
    set(summaryFigure1,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Stim_Evoked_HbT'])
end
%% individual stim-evoked figures
summaryFigure2 = figure;
sgtitle('Stimulus-evoked \DeltaHbT repsonses - individual animals')
%% LH contra stim
ax1 = subplot(3,2,1);
% C57BL6Js
for aa = 1:size(data.C57BL6J.Contra.adjLH.HbT,1)
    plot(data.C57BL6J.Contra.adjLH.meanTimeVector,data.C57BL6J.Contra.adjLH.HbT(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Contra.adjLH.HbT,1)
    plot(data.Blank_SAP.Contra.adjLH.meanTimeVector,data.Blank_SAP.Contra.adjLH.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Contra.adjLH.HbT,1)
    plot(data.SSP_SAP.Contra.adjLH.meanTimeVector,data.SSP_SAP.Contra.adjLH.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('LH (UnRx) \DeltaHbT - Contra Stim')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH contra stim
ax2 = subplot(3,2,2);
% C57BL6Js
for aa = 1:size(data.C57BL6J.Contra.adjRH.HbT,1)
    plot(data.C57BL6J.Contra.adjRH.meanTimeVector,data.C57BL6J.Contra.adjRH.HbT(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Contra.adjRH.HbT,1)
    plot(data.Blank_SAP.Contra.adjRH.meanTimeVector,data.Blank_SAP.Contra.adjRH.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Contra.adjRH.HbT,1)
    plot(data.SSP_SAP.Contra.adjRH.meanTimeVector,data.SSP_SAP.Contra.adjRH.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('RH (Rx) \DeltaHbT - Contra Stim')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% LH ipsi stim
ax3 = subplot(3,2,3);
% C57BL6Js
for aa = 1:size(data.C57BL6J.Ipsi.adjLH.HbT,1)
    plot(data.C57BL6J.Ipsi.adjLH.meanTimeVector,data.C57BL6J.Ipsi.adjLH.HbT(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Ipsi.adjLH.HbT,1)
    plot(data.Blank_SAP.Ipsi.adjLH.meanTimeVector,data.Blank_SAP.Ipsi.adjLH.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Ipsi.adjLH.HbT,1)
    plot(data.SSP_SAP.Ipsi.adjLH.meanTimeVector,data.SSP_SAP.Ipsi.adjLH.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('LH (UnRx) \DeltaHbT - Ipsi Stim')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH ipsi stim
ax4 = subplot(3,2,4);
% C57BL6Js
for aa = 1:size(data.C57BL6J.Ipsi.adjRH.HbT,1)
    plot(data.C57BL6J.Ipsi.adjRH.meanTimeVector,data.C57BL6J.Ipsi.adjRH.HbT(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Ipsi.adjRH.HbT,1)
    plot(data.Blank_SAP.Ipsi.adjRH.meanTimeVector,data.Blank_SAP.Ipsi.adjRH.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Ipsi.adjRH.HbT,1)
    plot(data.SSP_SAP.Ipsi.adjRH.meanTimeVector,data.SSP_SAP.Ipsi.adjRH.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('RH (Rx) \DeltaHbT - Ipsi Stim')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% LH auditory stim
ax5 = subplot(3,2,5);
% C57BL6Js
for aa = 1:size(data.C57BL6J.Auditory.adjLH.HbT,1)
    plot(data.C57BL6J.Auditory.adjLH.meanTimeVector,data.C57BL6J.Auditory.adjLH.HbT(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Auditory.adjLH.HbT,1)
    plot(data.Blank_SAP.Auditory.adjLH.meanTimeVector,data.Blank_SAP.Auditory.adjLH.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Auditory.adjLH.HbT,1)
    plot(data.SSP_SAP.Auditory.adjLH.meanTimeVector,data.SSP_SAP.Auditory.adjLH.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('LH (UnRx) \DeltaHbT - Auditory Stim')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH auditory stim
ax6 = subplot(3,2,6);
% C57BL6Js
for aa = 1:size(data.C57BL6J.Auditory.adjRH.HbT,1)
    plot(data.C57BL6J.Auditory.adjRH.meanTimeVector,data.C57BL6J.Auditory.adjRH.HbT(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Auditory.adjRH.HbT,1)
    plot(data.Blank_SAP.Auditory.adjRH.meanTimeVector,data.Blank_SAP.Auditory.adjRH.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Auditory.adjRH.HbT,1)
    plot(data.SSP_SAP.Auditory.adjRH.meanTimeVector,data.SSP_SAP.Auditory.adjRH.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('RH (Rx) \DeltaHbT - Auditory Stim')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% figure characteristics
linkaxes([ax1,ax2],'xy')
linkaxes([ax3,ax4],'xy')
linkaxes([ax5,ax6],'xy')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure2,[dirpath 'indStim_Evoked_HbT']);
    set(summaryFigure2,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'indStim_Evoked_HbT'])
end
%% average whisk-evoked figures
summaryFigure3 = figure;
sgtitle('Stimulus-evoked cortical MUA [300-3000 Hz]  repsonses')
%% LH short whisks
ax1 = subplot(3,2,1);
% C57BL6Js
p1 = plot(data.C57BL6J.Contra.adjLH.meanTimeVector,data.C57BL6J.Contra.adjLH.meanCortMUA,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.C57BL6J.Contra.adjLH.meanTimeVector,data.C57BL6J.Contra.adjLH.meanCortMUA + data.C57BL6J.Contra.adjLH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.C57BL6J.Contra.adjLH.meanTimeVector,data.C57BL6J.Contra.adjLH.meanCortMUA - data.C57BL6J.Contra.adjLH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
p2 = plot(data.Blank_SAP.Contra.adjLH.meanTimeVector,data.Blank_SAP.Contra.adjLH.meanCortMUA,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.Contra.adjLH.meanTimeVector,data.Blank_SAP.Contra.adjLH.meanCortMUA + data.Blank_SAP.Contra.adjLH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.Contra.adjLH.meanTimeVector,data.Blank_SAP.Contra.adjLH.meanCortMUA - data.Blank_SAP.Contra.adjLH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
p3 = plot(data.SSP_SAP.Contra.adjLH.meanTimeVector,data.SSP_SAP.Contra.adjLH.meanCortMUA,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.Contra.adjLH.meanTimeVector,data.SSP_SAP.Contra.adjLH.meanCortMUA + data.SSP_SAP.Contra.adjLH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.Contra.adjLH.meanTimeVector,data.SSP_SAP.Contra.adjLH.meanCortMUA - data.SSP_SAP.Contra.adjLH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
title('LH (UnRx) MUA [300-3000 Hz] - Contra Stim')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH short whisks
ax2 = subplot(3,2,2);
% C57BL6Js
plot(data.C57BL6J.Contra.adjRH.meanTimeVector,data.C57BL6J.Contra.adjRH.meanCortMUA,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.C57BL6J.Contra.adjRH.meanTimeVector,data.C57BL6J.Contra.adjRH.meanCortMUA + data.C57BL6J.Contra.adjRH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.C57BL6J.Contra.adjRH.meanTimeVector,data.C57BL6J.Contra.adjRH.meanCortMUA - data.C57BL6J.Contra.adjRH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
plot(data.Blank_SAP.Contra.adjRH.meanTimeVector,data.Blank_SAP.Contra.adjRH.meanCortMUA,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.Contra.adjRH.meanTimeVector,data.Blank_SAP.Contra.adjRH.meanCortMUA + data.Blank_SAP.Contra.adjRH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.Contra.adjRH.meanTimeVector,data.Blank_SAP.Contra.adjRH.meanCortMUA - data.Blank_SAP.Contra.adjRH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.Contra.adjRH.meanTimeVector,data.SSP_SAP.Contra.adjRH.meanCortMUA,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.Contra.adjRH.meanTimeVector,data.SSP_SAP.Contra.adjRH.meanCortMUA + data.SSP_SAP.Contra.adjRH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.Contra.adjRH.meanTimeVector,data.SSP_SAP.Contra.adjRH.meanCortMUA - data.SSP_SAP.Contra.adjRH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
title('RH (Rx) MUA [300-3000 Hz] - Contra Stim')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% LH intermediate whisks
ax3 = subplot(3,2,3);
% C57BL6Js
plot(data.C57BL6J.Ipsi.adjLH.meanTimeVector,data.C57BL6J.Ipsi.adjLH.meanCortMUA,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.C57BL6J.Ipsi.adjLH.meanTimeVector,data.C57BL6J.Ipsi.adjLH.meanCortMUA + data.C57BL6J.Ipsi.adjLH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.C57BL6J.Ipsi.adjLH.meanTimeVector,data.C57BL6J.Ipsi.adjLH.meanCortMUA - data.C57BL6J.Ipsi.adjLH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
plot(data.Blank_SAP.Ipsi.adjLH.meanTimeVector,data.Blank_SAP.Ipsi.adjLH.meanCortMUA,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.Ipsi.adjLH.meanTimeVector,data.Blank_SAP.Ipsi.adjLH.meanCortMUA + data.Blank_SAP.Ipsi.adjLH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.Ipsi.adjLH.meanTimeVector,data.Blank_SAP.Ipsi.adjLH.meanCortMUA - data.Blank_SAP.Ipsi.adjLH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.Ipsi.adjLH.meanTimeVector,data.SSP_SAP.Ipsi.adjLH.meanCortMUA,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.Ipsi.adjLH.meanTimeVector,data.SSP_SAP.Ipsi.adjLH.meanCortMUA + data.SSP_SAP.Ipsi.adjLH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.Ipsi.adjLH.meanTimeVector,data.SSP_SAP.Ipsi.adjLH.meanCortMUA - data.SSP_SAP.Ipsi.adjLH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
title('LH (UnRx) MUA [300-3000 Hz] - Ipsi Stim')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH intermediate whisks
ax4 = subplot(3,2,4);
% C57BL6Js
plot(data.C57BL6J.Ipsi.adjRH.meanTimeVector,data.C57BL6J.Ipsi.adjRH.meanCortMUA,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.C57BL6J.Ipsi.adjRH.meanTimeVector,data.C57BL6J.Ipsi.adjRH.meanCortMUA + data.C57BL6J.Ipsi.adjRH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.C57BL6J.Ipsi.adjRH.meanTimeVector,data.C57BL6J.Ipsi.adjRH.meanCortMUA - data.C57BL6J.Ipsi.adjRH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
plot(data.Blank_SAP.Ipsi.adjRH.meanTimeVector,data.Blank_SAP.Ipsi.adjRH.meanCortMUA,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.Ipsi.adjRH.meanTimeVector,data.Blank_SAP.Ipsi.adjRH.meanCortMUA + data.Blank_SAP.Ipsi.adjRH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.Ipsi.adjRH.meanTimeVector,data.Blank_SAP.Ipsi.adjRH.meanCortMUA - data.Blank_SAP.Ipsi.adjRH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.Ipsi.adjRH.meanTimeVector,data.SSP_SAP.Ipsi.adjRH.meanCortMUA,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.Ipsi.adjRH.meanTimeVector,data.SSP_SAP.Ipsi.adjRH.meanCortMUA + data.SSP_SAP.Ipsi.adjRH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.Ipsi.adjRH.meanTimeVector,data.SSP_SAP.Ipsi.adjRH.meanCortMUA - data.SSP_SAP.Ipsi.adjRH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
title('RH (Rx) MUA [300-3000 Hz] - Ipsi Stim')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% LH long whisks
ax5 = subplot(3,2,5);
% C57BL6Js
plot(data.C57BL6J.Auditory.adjLH.meanTimeVector,data.C57BL6J.Auditory.adjLH.meanCortMUA,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.C57BL6J.Auditory.adjLH.meanTimeVector,data.C57BL6J.Auditory.adjLH.meanCortMUA + data.C57BL6J.Auditory.adjLH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.C57BL6J.Auditory.adjLH.meanTimeVector,data.C57BL6J.Auditory.adjLH.meanCortMUA - data.C57BL6J.Auditory.adjLH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
plot(data.Blank_SAP.Auditory.adjLH.meanTimeVector,data.Blank_SAP.Auditory.adjLH.meanCortMUA,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.Auditory.adjLH.meanTimeVector,data.Blank_SAP.Auditory.adjLH.meanCortMUA + data.Blank_SAP.Auditory.adjLH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.Auditory.adjLH.meanTimeVector,data.Blank_SAP.Auditory.adjLH.meanCortMUA - data.Blank_SAP.Auditory.adjLH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.Auditory.adjLH.meanTimeVector,data.SSP_SAP.Auditory.adjLH.meanCortMUA,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.Auditory.adjLH.meanTimeVector,data.SSP_SAP.Auditory.adjLH.meanCortMUA + data.SSP_SAP.Auditory.adjLH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.Auditory.adjLH.meanTimeVector,data.SSP_SAP.Auditory.adjLH.meanCortMUA - data.SSP_SAP.Auditory.adjLH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
title('LH (UnRx) MUA [300-3000 Hz] - Auditory Stim')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH long whisks
ax6 = subplot(3,2,6);
% C57BL6Js
plot(data.C57BL6J.Auditory.adjRH.meanTimeVector,data.C57BL6J.Auditory.adjRH.meanCortMUA,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.C57BL6J.Auditory.adjRH.meanTimeVector,data.C57BL6J.Auditory.adjRH.meanCortMUA + data.C57BL6J.Auditory.adjRH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.C57BL6J.Auditory.adjRH.meanTimeVector,data.C57BL6J.Auditory.adjRH.meanCortMUA - data.C57BL6J.Auditory.adjRH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
plot(data.Blank_SAP.Auditory.adjRH.meanTimeVector,data.Blank_SAP.Auditory.adjRH.meanCortMUA,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.Auditory.adjRH.meanTimeVector,data.Blank_SAP.Auditory.adjRH.meanCortMUA + data.Blank_SAP.Auditory.adjRH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.Auditory.adjRH.meanTimeVector,data.Blank_SAP.Auditory.adjRH.meanCortMUA - data.Blank_SAP.Auditory.adjRH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.Auditory.adjRH.meanTimeVector,data.SSP_SAP.Auditory.adjRH.meanCortMUA,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.Auditory.adjRH.meanTimeVector,data.SSP_SAP.Auditory.adjRH.meanCortMUA + data.SSP_SAP.Auditory.adjRH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.Auditory.adjRH.meanTimeVector,data.SSP_SAP.Auditory.adjRH.meanCortMUA - data.SSP_SAP.Auditory.adjRH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
title('RH (Rx) MUA [300-3000 Hz] - Auditory Stim')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% figure characteristics
linkaxes([ax1,ax2],'xy')
linkaxes([ax3,ax4],'xy')
linkaxes([ax5,ax6],'xy')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure3,[dirpath 'Stim_Evoked_MUA']);
    set(summaryFigure3,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Stim_Evoked_MUA'])
end
%% individual whisk-evoked figures
summaryFigure4 = figure;
sgtitle('Stimulus-evoked cortical MUA [300-3000 Hz] repsonses - individual animals')
%% LH short whisks
ax1 = subplot(3,2,1);
% C57BL6Js
for aa = 1:size(data.C57BL6J.Contra.adjLH.cortMUA,1)
    plot(data.C57BL6J.Contra.adjLH.meanTimeVector,data.C57BL6J.Contra.adjLH.cortMUA(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Contra.adjLH.cortMUA,1)
    plot(data.Blank_SAP.Contra.adjLH.meanTimeVector,data.Blank_SAP.Contra.adjLH.cortMUA(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Contra.adjLH.cortMUA,1)
    plot(data.SSP_SAP.Contra.adjLH.meanTimeVector,data.SSP_SAP.Contra.adjLH.cortMUA(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('LH (UnRx) MUA [300-3000 Hz] - Contra Stim')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
xlim([-2,10])
%% RH short whisks
ax2 = subplot(3,2,2);
% C57BL6Js
for aa = 1:size(data.C57BL6J.Contra.adjRH.cortMUA,1)
    plot(data.C57BL6J.Contra.adjRH.meanTimeVector,data.C57BL6J.Contra.adjRH.cortMUA(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Contra.adjRH.cortMUA,1)
    plot(data.Blank_SAP.Contra.adjRH.meanTimeVector,data.Blank_SAP.Contra.adjRH.cortMUA(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Contra.adjRH.cortMUA,1)
    plot(data.SSP_SAP.Contra.adjRH.meanTimeVector,data.SSP_SAP.Contra.adjRH.cortMUA(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('RH (Rx) MUA [300-3000 Hz] - Contra Stim')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% LH intermediate whisks
ax3 = subplot(3,2,3);
% C57BL6Js
for aa = 1:size(data.C57BL6J.Ipsi.adjLH.cortMUA,1)
    plot(data.C57BL6J.Ipsi.adjLH.meanTimeVector,data.C57BL6J.Ipsi.adjLH.cortMUA(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Ipsi.adjLH.cortMUA,1)
    plot(data.Blank_SAP.Ipsi.adjLH.meanTimeVector,data.Blank_SAP.Ipsi.adjLH.cortMUA(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Ipsi.adjLH.cortMUA,1)
    plot(data.SSP_SAP.Ipsi.adjLH.meanTimeVector,data.SSP_SAP.Ipsi.adjLH.cortMUA(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('LH (UnRx) MUA [300-3000 Hz] - Ipsi Stim')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH intermediate whisks
ax4 = subplot(3,2,4);
% C57BL6Js
for aa = 1:size(data.C57BL6J.Ipsi.adjRH.cortMUA,1)
    plot(data.C57BL6J.Ipsi.adjRH.meanTimeVector,data.C57BL6J.Ipsi.adjRH.cortMUA(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Ipsi.adjRH.cortMUA,1)
    plot(data.Blank_SAP.Ipsi.adjRH.meanTimeVector,data.Blank_SAP.Ipsi.adjRH.cortMUA(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Ipsi.adjRH.cortMUA,1)
    plot(data.SSP_SAP.Ipsi.adjRH.meanTimeVector,data.SSP_SAP.Ipsi.adjRH.cortMUA(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('RH (Rx) MUA [300-3000 Hz] - Ipsi Stim')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% LH long whisks
ax5 = subplot(3,2,5);
% C57BL6Js
for aa = 1:size(data.C57BL6J.Auditory.adjLH.cortMUA,1)
    plot(data.C57BL6J.Auditory.adjLH.meanTimeVector,data.C57BL6J.Auditory.adjLH.cortMUA(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Auditory.adjLH.cortMUA,1)
    plot(data.Blank_SAP.Auditory.adjLH.meanTimeVector,data.Blank_SAP.Auditory.adjLH.cortMUA(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Auditory.adjLH.cortMUA,1)
    plot(data.SSP_SAP.Auditory.adjLH.meanTimeVector,data.SSP_SAP.Auditory.adjLH.cortMUA(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('LH (UnRx) MUA [300-3000 Hz] - Auditory Stim')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH long whisks
ax6 = subplot(3,2,6);
% C57BL6Js
for aa = 1:size(data.C57BL6J.Auditory.adjRH.cortMUA,1)
    plot(data.C57BL6J.Auditory.adjRH.meanTimeVector,data.C57BL6J.Auditory.adjRH.cortMUA(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.Auditory.adjRH.cortMUA,1)
    plot(data.Blank_SAP.Auditory.adjRH.meanTimeVector,data.Blank_SAP.Auditory.adjRH.cortMUA(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.Auditory.adjRH.cortMUA,1)
    plot(data.SSP_SAP.Auditory.adjRH.meanTimeVector,data.SSP_SAP.Auditory.adjRH.cortMUA(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('RH (Rx) MUA [300-3000 Hz] - Auditory Stim')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'C57BL6J','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% figure characteristics
linkaxes([ax1,ax2],'xy')
linkaxes([ax3,ax4],'xy')
linkaxes([ax5,ax6],'xy')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure4,[dirpath 'indStim_Evoked_MUA']);
    set(summaryFigure4,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'indStim_Evoked_MUA'])
end
%% average neural responses
summaryFigure5 = figure;
sgtitle('Stimulus-evoked cortical neural (LFP) repsonses')
%% stim cortical LFP
subplot(3,6,1);
imagesc(data.C57BL6J.Contra.adjLH.meanCortT,data.C57BL6J.Contra.adjLH.meanCortF,data.C57BL6J.Contra.adjLH.meanCortS)
title('C57BL6J Contra Stim LH')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
caxis([-50,100])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% stim cortical LFP
subplot(3,6,2);
imagesc(data.C57BL6J.Contra.adjRH.meanCortT,data.C57BL6J.Contra.adjRH.meanCortF,data.C57BL6J.Contra.adjRH.meanCortS)
title('C57BL6J Contra Stim RH')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
caxis([-50,100])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% stim cortical LFP
subplot(3,6,3);
imagesc(data.SSP_SAP.Contra.adjLH.meanCortT,data.SSP_SAP.Contra.adjLH.meanCortF,data.SSP_SAP.Contra.adjLH.meanCortS)
title('SSP-SAP Contra Stim LH')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
caxis([-50,100])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% stim cortical LFP
subplot(3,6,4);
imagesc(data.SSP_SAP.Contra.adjRH.meanCortT,data.SSP_SAP.Contra.adjRH.meanCortF,data.SSP_SAP.Contra.adjRH.meanCortS)
title('SSP-SAP Contra Stim RH')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
caxis([-50,100])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% stim cortical LFP
subplot(3,6,5);
imagesc(data.Blank_SAP.Contra.adjLH.meanCortT,data.Blank_SAP.Contra.adjLH.meanCortF,data.Blank_SAP.Contra.adjLH.meanCortS)
title('Blank-SAP Contra Stim LH')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
caxis([-50,100])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% stim cortical LFP
subplot(3,6,6);
imagesc(data.Blank_SAP.Contra.adjRH.meanCortT,data.Blank_SAP.Contra.adjRH.meanCortF,data.Blank_SAP.Contra.adjRH.meanCortS)
title('Blank-SAP Contra Stim RH')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
caxis([-50,100])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% Ipsi cortical LFP
subplot(3,6,7);
imagesc(data.C57BL6J.Ipsi.adjLH.meanCortT,data.C57BL6J.Ipsi.adjLH.meanCortF,data.C57BL6J.Ipsi.adjLH.meanCortS)
title('C57BL6J Ipsi Stim LH')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
caxis([-50,100])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% Ipsi cortical LFP
subplot(3,6,8);
imagesc(data.C57BL6J.Ipsi.adjRH.meanCortT,data.C57BL6J.Ipsi.adjRH.meanCortF,data.C57BL6J.Ipsi.adjRH.meanCortS)
title('C57BL6J Ipsi Stim RH')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
caxis([-50,100])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% Ipsi cortical LFP
subplot(3,6,9);
imagesc(data.SSP_SAP.Ipsi.adjLH.meanCortT,data.SSP_SAP.Ipsi.adjLH.meanCortF,data.SSP_SAP.Ipsi.adjLH.meanCortS)
title('SSP-SAP Ipsi Stim LH')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
caxis([-50,100])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% Ipsi cortical LFP
subplot(3,6,10);
imagesc(data.SSP_SAP.Ipsi.adjRH.meanCortT,data.SSP_SAP.Ipsi.adjRH.meanCortF,data.SSP_SAP.Ipsi.adjRH.meanCortS)
title('SSP-SAP Ipsi Stim RH')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
caxis([-50,100])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% Ipsi cortical LFP
subplot(3,6,11);
imagesc(data.Blank_SAP.Ipsi.adjLH.meanCortT,data.Blank_SAP.Ipsi.adjLH.meanCortF,data.Blank_SAP.Ipsi.adjLH.meanCortS)
title('Blank-SAP Ipsi Stim LH')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
caxis([-50,100])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% Ipsi cortical LFP
subplot(3,6,12);
imagesc(data.Blank_SAP.Ipsi.adjRH.meanCortT,data.Blank_SAP.Ipsi.adjRH.meanCortF,data.Blank_SAP.Ipsi.adjRH.meanCortS)
title('Blank-SAP Ipsi Stim RH')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
caxis([-50,100])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% Auditory cortical LFP
subplot(3,6,13);
imagesc(data.C57BL6J.Auditory.adjLH.meanCortT,data.C57BL6J.Auditory.adjLH.meanCortF,data.C57BL6J.Auditory.adjLH.meanCortS)
title('C57BL6J Auditory Stim LH')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
caxis([-50,100])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% stim cortical LFP
subplot(3,6,14);
imagesc(data.C57BL6J.Auditory.adjRH.meanCortT,data.C57BL6J.Auditory.adjRH.meanCortF,data.C57BL6J.Auditory.adjRH.meanCortS)
title('C57BL6J Auditory Stim RH')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
caxis([-50,100])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% stim cortical LFP
subplot(3,6,15);
imagesc(data.SSP_SAP.Auditory.adjLH.meanCortT,data.SSP_SAP.Auditory.adjLH.meanCortF,data.SSP_SAP.Auditory.adjLH.meanCortS)
title('SSP-SAP Auditory Stim LH')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
caxis([-50,100])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% stim cortical LFP
subplot(3,6,16);
imagesc(data.SSP_SAP.Auditory.adjRH.meanCortT,data.SSP_SAP.Auditory.adjRH.meanCortF,data.SSP_SAP.Auditory.adjRH.meanCortS)
title('SSP-SAP Auditory Stim RH')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
caxis([-50,100])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% stim cortical LFP
subplot(3,6,17);
imagesc(data.Blank_SAP.Auditory.adjLH.meanCortT,data.Blank_SAP.Auditory.adjLH.meanCortF,data.Blank_SAP.Auditory.adjLH.meanCortS)
title('Blank-SAP Auditory Stim LH')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
caxis([-50,100])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% stim cortical LFP
subplot(3,6,18);
imagesc(data.Blank_SAP.Auditory.adjRH.meanCortT,data.Blank_SAP.Auditory.adjRH.meanCortF,data.Blank_SAP.Auditory.adjRH.meanCortS)
title('Blank-SAP Auditory Stim RH')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
caxis([-50,100])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure5,[dirpath 'Stim_Evoked_LFP']);
    set(summaryFigure5,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Stim_Evoked_LFP'])
end

end
