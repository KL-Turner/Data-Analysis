function [AnalysisResults] = StimEvoked_Saporin(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S2 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
%% set-up
animalIDs = {'T141','T155','T156','T157','T142','T144','T159','T172','T150','T165','T166'};
C57BL6J_IDs = {'T141','T155','T156','T157'};
SSP_SAP_IDs = {'T142','T144','T159','T172'};
Blank_SAP_IDs = {'T150','T165','T166'};
solenoidNames = {'LPadSol','RPadSol','AudSol'};
compDataTypes = {'Ipsi','Contra','Auditory'};
hemispheres = {'adjLH','adjRH'};
treatments = {'C57BL6J','Blank_SAP','SSP_SAP'};
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
            data.(treatment).(compDataType).(hemisphere).meanT = mean(data.(treatment).(compDataType).(hemisphere).cortT,1);
            data.(treatment).(compDataType).(hemisphere).meanF = mean(data.(treatment).(compDataType).(hemisphere).cortF,1);
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
%% Figures
summaryFigure = figure;
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

linkaxes([ax1,ax2],'xy')
linkaxes([ax3,ax4],'xy')
linkaxes([ax5,ax6],'xy')

end

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
% linkaxes([ax5,ax6,ax7,ax8],'xy')
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
