function [] = WhiskEvoked_Bilateral_IOS(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

resultsStruct = 'Results_Evoked';
load(resultsStruct);
%% set-up
expGroups = {'Naive','SSP-SAP','Blank-SAP'};
setName = 'IOS Set A';
animalIDs.all = {};
for aa = 1:length(expGroups)
    folderList = dir([expGroups{1,aa} delim setName]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs.all = horzcat(animalIDs.all,{folderList.name});
    animalIDs.(strrep(expGroups{1,aa},'-','_')) = {folderList.name};
end
treatments = {'Naive','SSP_SAP','Blank_SAP'};
whiskDataTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
hemispheres = {'adjLH','adjRH'};
data = [];
cortVariables = {'HbT','CBV','cortMUA','cortGam','cortS','cortS_Gam','cortT','cortF','timeVector'};
hipVariables = {'hipMUA','hipGam','hipS','hipS_Gam','hipT','hipF','timeVector'};
%% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs.all)
    animalID = animalIDs.all{1,aa};
    % recognize treatment based on animal group
    if ismember(animalID,animalIDs.Naive) == true
        treatment = 'Naive';
    elseif ismember(animalID,animalIDs.SSP_SAP) == true
        treatment = 'SSP_SAP';
    elseif ismember(animalID,animalIDs.Blank_SAP) == true
        treatment = 'Blank_SAP';
    end
    % short - intermediate - long whisks
    for bb = 1:length(whiskDataTypes)
        % left, right hemishpere hemo & neural data
        for cc = 1:length(hemispheres)
            % pre-allocate necessary variable fields
            data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}).dummCheck = 1;
            for dd = 1:length(cortVariables)
                if isfield(data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}),cortVariables{1,dd}) == false
                    data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}).(cortVariables{1,dd}) = [];
                end
            end
            data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}).HbT = cat(1,data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}).HbT,Results_Evoked.(animalIDs.all{1,aa}).Whisk.(hemispheres{1,cc}).(whiskDataTypes{1,bb}).CBV_HbT.HbT);
            data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}).CBV = cat(1,data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}).CBV,Results_Evoked.(animalIDs.all{1,aa}).Whisk.(hemispheres{1,cc}).(whiskDataTypes{1,bb}).CBV.CBV);
            data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}).cortMUA = cat(1,data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}).cortMUA,Results_Evoked.(animalIDs.all{1,aa}).Whisk.(hemispheres{1,cc}).(whiskDataTypes{1,bb}).MUA.corticalData);
            data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}).cortGam = cat(1,data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}).cortGam,Results_Evoked.(animalIDs.all{1,aa}).Whisk.(hemispheres{1,cc}).(whiskDataTypes{1,bb}).Gam.corticalData);
            data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}).cortS = cat(3,data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}).cortS,Results_Evoked.(animalIDs.all{1,aa}).Whisk.(hemispheres{1,cc}).(whiskDataTypes{1,bb}).LFP.corticalS);
            data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}).cortS_Gam = cat(3,data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}).cortS_Gam,Results_Evoked.(animalIDs.all{1,aa}).Whisk.(hemispheres{1,cc}).(whiskDataTypes{1,bb}).LFP.corticalS(49:end,20:23));
            data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}).cortT = cat(1,data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}).cortT,Results_Evoked.(animalIDs.all{1,aa}).Whisk.(hemispheres{1,cc}).(whiskDataTypes{1,bb}).LFP.T);
            data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}).cortF = cat(1,data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}).cortF,Results_Evoked.(animalIDs.all{1,aa}).Whisk.(hemispheres{1,cc}).(whiskDataTypes{1,bb}).LFP.F);
            data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}).timeVector = cat(1,data.(treatment).(whiskDataTypes{1,bb}).(hemispheres{1,cc}).timeVector,Results_Evoked.(animalIDs.all{1,aa}).Whisk.(hemispheres{1,cc}).(whiskDataTypes{1,bb}).timeVector);
        end
        % hippocampal neural data - preallocate necessary variable fields
        data.(treatment).(whiskDataTypes{1,bb}).Hip.dummyCheck = 1;
        for ee = 1:length(hipVariables)
            if isfield(data.(treatment).(whiskDataTypes{1,bb}).Hip,hipVariables{1,ee}) == false
                data.(treatment).(whiskDataTypes{1,bb}).Hip.(hipVariables{1,ee}) = [];
            end
        end
        data.(treatment).(whiskDataTypes{1,bb}).Hip.hipMUA = cat(1,data.(treatment).(whiskDataTypes{1,bb}).Hip.hipMUA,Results_Evoked.(animalIDs.all{1,aa}).Whisk.adjLH.(whiskDataTypes{1,bb}).MUA.hippocampalData);
        data.(treatment).(whiskDataTypes{1,bb}).Hip.hipGam = cat(1,data.(treatment).(whiskDataTypes{1,bb}).Hip.hipGam,Results_Evoked.(animalIDs.all{1,aa}).Whisk.adjLH.(whiskDataTypes{1,bb}).Gam.hippocampalData);
        data.(treatment).(whiskDataTypes{1,bb}).Hip.hipS = cat(3,data.(treatment).(whiskDataTypes{1,bb}).Hip.hipS,Results_Evoked.(animalIDs.all{1,aa}).Whisk.adjLH.(whiskDataTypes{1,bb}).LFP.hippocampalS);
        data.(treatment).(whiskDataTypes{1,bb}).Hip.hipS_Gam = cat(3,data.(treatment).(whiskDataTypes{1,bb}).Hip.hipS_Gam,Results_Evoked.(animalIDs.all{1,aa}).Whisk.adjLH.(whiskDataTypes{1,bb}).LFP.hippocampalS(49:end,20:23));
        data.(treatment).(whiskDataTypes{1,bb}).Hip.hipT = cat(1,data.(treatment).(whiskDataTypes{1,bb}).Hip.hipT,Results_Evoked.(animalIDs.all{1,aa}).Whisk.adjLH.(whiskDataTypes{1,bb}).LFP.T);
        data.(treatment).(whiskDataTypes{1,bb}).Hip.hipF = cat(1,data.(treatment).(whiskDataTypes{1,bb}).Hip.hipF,Results_Evoked.(animalIDs.all{1,aa}).Whisk.adjLH.(whiskDataTypes{1,bb}).LFP.F);
        data.(treatment).(whiskDataTypes{1,bb}).Hip.timeVector = cat(1,data.(treatment).(whiskDataTypes{1,bb}).Hip.timeVector,Results_Evoked.(animalIDs.all{1,aa}).Whisk.adjLH.(whiskDataTypes{1,bb}).timeVector);
    end
end
%% concatenate the data from the contra and ipsi data
for ff = 1:length(treatments)
    for gg = 1:length(whiskDataTypes)
        for hh = 1:length(hemispheres)
            data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).meanHbT = mean(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).HbT,1);
            data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).stdHbT = std(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).HbT,0,1);
            data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).meanCBV = mean(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).CBV,1);
            data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).stdCBV = std(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).CBV,0,1);
            data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).meanCortMUA = mean(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).cortMUA,1);
            data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).stdCortMUA = std(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).cortMUA,0,1);
            data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).meanCortGam = mean(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).cortGam,1);
            data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).stdCortGam = std(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).cortGam,0,1);
            data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).meanCortS = mean(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).cortS,3).*100;
            data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).mean_CortS_Gam = mean(mean(mean(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).cortS_Gam.*100,1),2),3);
            data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).std_CortS_Gam = std(mean(mean(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).cortS_Gam.*100,1),2),0,3);
            data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).meanCortT = mean(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).cortT,1);
            data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).meanCortF = mean(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).cortF,1);
            data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).meanTimeVector = mean(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).(hemispheres{1,hh}).timeVector,1);
        end
        data.(treatments{1,ff}).(whiskDataTypes{1,gg}).Hip.meanHipMUA = mean(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).Hip.hipMUA,1);
        data.(treatments{1,ff}).(whiskDataTypes{1,gg}).Hip.stdHipMUA = std(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).Hip.hipMUA,0,1);
        data.(treatments{1,ff}).(whiskDataTypes{1,gg}).Hip.meanHipGam = mean(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).Hip.hipGam,1);
        data.(treatments{1,ff}).(whiskDataTypes{1,gg}).Hip.stdHipGam = std(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).Hip.hipGam,0,1);
        data.(treatments{1,ff}).(whiskDataTypes{1,gg}).Hip.meanHipS = mean(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).Hip.hipS,3).*100;
        data.(treatments{1,ff}).(whiskDataTypes{1,gg}).Hip.mean_HipS_Gam = mean(mean(mean(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).Hip.hipS_Gam.*100,1),2),3);
        data.(treatments{1,ff}).(whiskDataTypes{1,gg}).Hip.std_HipS_Gam = std(mean(mean(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).Hip.hipS_Gam.*100,1),2),0,3);
        data.(treatments{1,ff}).(whiskDataTypes{1,gg}).Hip.meanHipT = mean(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).Hip.hipT,1);
        data.(treatments{1,ff}).(whiskDataTypes{1,gg}).Hip.meanHipF = mean(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).Hip.hipF,1);
        data.(treatments{1,ff}).(whiskDataTypes{1,gg}).Hip.meanTimeVector = mean(data.(treatments{1,ff}).(whiskDataTypes{1,gg}).Hip.timeVector,1);
    end
end
%% average whisk-evoked figures
summaryFigure1 = figure;
sgtitle('Whisking-evoked \DeltaHbT repsonses')
%% LH short whisks
ax1 = subplot(3,2,1);
% Naives
p1 = plot(data.Naive.ShortWhisks.adjLH.meanTimeVector,data.Naive.ShortWhisks.adjLH.meanHbT,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.Naive.ShortWhisks.adjLH.meanTimeVector,data.Naive.ShortWhisks.adjLH.meanHbT + data.Naive.ShortWhisks.adjLH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.Naive.ShortWhisks.adjLH.meanTimeVector,data.Naive.ShortWhisks.adjLH.meanHbT - data.Naive.ShortWhisks.adjLH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
p2 = plot(data.Blank_SAP.ShortWhisks.adjLH.meanTimeVector,data.Blank_SAP.ShortWhisks.adjLH.meanHbT,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.ShortWhisks.adjLH.meanTimeVector,data.Blank_SAP.ShortWhisks.adjLH.meanHbT + data.Blank_SAP.ShortWhisks.adjLH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.ShortWhisks.adjLH.meanTimeVector,data.Blank_SAP.ShortWhisks.adjLH.meanHbT - data.Blank_SAP.ShortWhisks.adjLH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
p3 = plot(data.SSP_SAP.ShortWhisks.adjLH.meanTimeVector,data.SSP_SAP.ShortWhisks.adjLH.meanHbT,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.ShortWhisks.adjLH.meanTimeVector,data.SSP_SAP.ShortWhisks.adjLH.meanHbT + data.SSP_SAP.ShortWhisks.adjLH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.ShortWhisks.adjLH.meanTimeVector,data.SSP_SAP.ShortWhisks.adjLH.meanHbT - data.SSP_SAP.ShortWhisks.adjLH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
title('LH (UnRx) \DeltaHbT - Short Whisks')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH short whisks
ax2 = subplot(3,2,2);
% Naives
plot(data.Naive.ShortWhisks.adjRH.meanTimeVector,data.Naive.ShortWhisks.adjRH.meanHbT,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.Naive.ShortWhisks.adjRH.meanTimeVector,data.Naive.ShortWhisks.adjRH.meanHbT + data.Naive.ShortWhisks.adjRH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.Naive.ShortWhisks.adjRH.meanTimeVector,data.Naive.ShortWhisks.adjRH.meanHbT - data.Naive.ShortWhisks.adjRH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
plot(data.Blank_SAP.ShortWhisks.adjRH.meanTimeVector,data.Blank_SAP.ShortWhisks.adjRH.meanHbT,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.ShortWhisks.adjRH.meanTimeVector,data.Blank_SAP.ShortWhisks.adjRH.meanHbT + data.Blank_SAP.ShortWhisks.adjRH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.ShortWhisks.adjRH.meanTimeVector,data.Blank_SAP.ShortWhisks.adjRH.meanHbT - data.Blank_SAP.ShortWhisks.adjRH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.ShortWhisks.adjRH.meanTimeVector,data.SSP_SAP.ShortWhisks.adjRH.meanHbT,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.ShortWhisks.adjRH.meanTimeVector,data.SSP_SAP.ShortWhisks.adjRH.meanHbT + data.SSP_SAP.ShortWhisks.adjRH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.ShortWhisks.adjRH.meanTimeVector,data.SSP_SAP.ShortWhisks.adjRH.meanHbT - data.SSP_SAP.ShortWhisks.adjRH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
title('RH (Rx) \DeltaHbT - Short Whisks')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% LH intermediate whisks
ax3 = subplot(3,2,3);
% Naives
plot(data.Naive.IntermediateWhisks.adjLH.meanTimeVector,data.Naive.IntermediateWhisks.adjLH.meanHbT,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.Naive.IntermediateWhisks.adjLH.meanTimeVector,data.Naive.IntermediateWhisks.adjLH.meanHbT + data.Naive.IntermediateWhisks.adjLH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.Naive.IntermediateWhisks.adjLH.meanTimeVector,data.Naive.IntermediateWhisks.adjLH.meanHbT - data.Naive.IntermediateWhisks.adjLH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
plot(data.Blank_SAP.IntermediateWhisks.adjLH.meanTimeVector,data.Blank_SAP.IntermediateWhisks.adjLH.meanHbT,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.IntermediateWhisks.adjLH.meanTimeVector,data.Blank_SAP.IntermediateWhisks.adjLH.meanHbT + data.Blank_SAP.IntermediateWhisks.adjLH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.IntermediateWhisks.adjLH.meanTimeVector,data.Blank_SAP.IntermediateWhisks.adjLH.meanHbT - data.Blank_SAP.IntermediateWhisks.adjLH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.IntermediateWhisks.adjLH.meanTimeVector,data.SSP_SAP.IntermediateWhisks.adjLH.meanHbT,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.IntermediateWhisks.adjLH.meanTimeVector,data.SSP_SAP.IntermediateWhisks.adjLH.meanHbT + data.SSP_SAP.IntermediateWhisks.adjLH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.IntermediateWhisks.adjLH.meanTimeVector,data.SSP_SAP.IntermediateWhisks.adjLH.meanHbT - data.SSP_SAP.IntermediateWhisks.adjLH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
title('LH (UnRx) \DeltaHbT - Intermediate Whisks')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH intermediate whisks
ax4 = subplot(3,2,4);
% Naives
plot(data.Naive.IntermediateWhisks.adjRH.meanTimeVector,data.Naive.IntermediateWhisks.adjRH.meanHbT,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.Naive.IntermediateWhisks.adjRH.meanTimeVector,data.Naive.IntermediateWhisks.adjRH.meanHbT + data.Naive.IntermediateWhisks.adjRH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.Naive.IntermediateWhisks.adjRH.meanTimeVector,data.Naive.IntermediateWhisks.adjRH.meanHbT - data.Naive.IntermediateWhisks.adjRH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
plot(data.Blank_SAP.IntermediateWhisks.adjRH.meanTimeVector,data.Blank_SAP.IntermediateWhisks.adjRH.meanHbT,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.IntermediateWhisks.adjRH.meanTimeVector,data.Blank_SAP.IntermediateWhisks.adjRH.meanHbT + data.Blank_SAP.IntermediateWhisks.adjRH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.IntermediateWhisks.adjRH.meanTimeVector,data.Blank_SAP.IntermediateWhisks.adjRH.meanHbT - data.Blank_SAP.IntermediateWhisks.adjRH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.IntermediateWhisks.adjRH.meanTimeVector,data.SSP_SAP.IntermediateWhisks.adjRH.meanHbT,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.IntermediateWhisks.adjRH.meanTimeVector,data.SSP_SAP.IntermediateWhisks.adjRH.meanHbT + data.SSP_SAP.IntermediateWhisks.adjRH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.IntermediateWhisks.adjRH.meanTimeVector,data.SSP_SAP.IntermediateWhisks.adjRH.meanHbT - data.SSP_SAP.IntermediateWhisks.adjRH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
title('RH (Rx) \DeltaHbT - Intermediate Whisks')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% LH long whisks
ax5 = subplot(3,2,5);
% Naives
plot(data.Naive.LongWhisks.adjLH.meanTimeVector,data.Naive.LongWhisks.adjLH.meanHbT,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.Naive.LongWhisks.adjLH.meanTimeVector,data.Naive.LongWhisks.adjLH.meanHbT + data.Naive.LongWhisks.adjLH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.Naive.LongWhisks.adjLH.meanTimeVector,data.Naive.LongWhisks.adjLH.meanHbT - data.Naive.LongWhisks.adjLH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
plot(data.Blank_SAP.LongWhisks.adjLH.meanTimeVector,data.Blank_SAP.LongWhisks.adjLH.meanHbT,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.LongWhisks.adjLH.meanTimeVector,data.Blank_SAP.LongWhisks.adjLH.meanHbT + data.Blank_SAP.LongWhisks.adjLH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.LongWhisks.adjLH.meanTimeVector,data.Blank_SAP.LongWhisks.adjLH.meanHbT - data.Blank_SAP.LongWhisks.adjLH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.LongWhisks.adjLH.meanTimeVector,data.SSP_SAP.LongWhisks.adjLH.meanHbT,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.LongWhisks.adjLH.meanTimeVector,data.SSP_SAP.LongWhisks.adjLH.meanHbT + data.SSP_SAP.LongWhisks.adjLH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.LongWhisks.adjLH.meanTimeVector,data.SSP_SAP.LongWhisks.adjLH.meanHbT - data.SSP_SAP.LongWhisks.adjLH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
title('LH (UnRx) \DeltaHbT - Long Whisks')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH long whisks
ax6 = subplot(3,2,6);
% Naives
plot(data.Naive.LongWhisks.adjRH.meanTimeVector,data.Naive.LongWhisks.adjRH.meanHbT,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.Naive.LongWhisks.adjRH.meanTimeVector,data.Naive.LongWhisks.adjRH.meanHbT + data.Naive.LongWhisks.adjRH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.Naive.LongWhisks.adjRH.meanTimeVector,data.Naive.LongWhisks.adjRH.meanHbT - data.Naive.LongWhisks.adjRH.stdHbT,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
plot(data.Blank_SAP.LongWhisks.adjRH.meanTimeVector,data.Blank_SAP.LongWhisks.adjRH.meanHbT,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.LongWhisks.adjRH.meanTimeVector,data.Blank_SAP.LongWhisks.adjRH.meanHbT + data.Blank_SAP.LongWhisks.adjRH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.LongWhisks.adjRH.meanTimeVector,data.Blank_SAP.LongWhisks.adjRH.meanHbT - data.Blank_SAP.LongWhisks.adjRH.stdHbT,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.LongWhisks.adjRH.meanTimeVector,data.SSP_SAP.LongWhisks.adjRH.meanHbT,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.LongWhisks.adjRH.meanTimeVector,data.SSP_SAP.LongWhisks.adjRH.meanHbT + data.SSP_SAP.LongWhisks.adjRH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.LongWhisks.adjRH.meanTimeVector,data.SSP_SAP.LongWhisks.adjRH.meanHbT - data.SSP_SAP.LongWhisks.adjRH.stdHbT,'color',colors('electric purple'),'LineWidth',0.5)
title('RH (Rx) \DeltaHbT - Long Whisks')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% figure characteristics
linkaxes([ax1,ax2],'xy')
linkaxes([ax3,ax4],'xy')
linkaxes([ax5,ax6],'xy')
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Whisking Evoked - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure1,[dirpath 'AverageWhiskEvoked_HbT']);
    set(summaryFigure1,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'AverageWhiskEvoked_HbT'])
end
%% individual whisk-evoked figures
summaryFigure2 = figure;
sgtitle('Whisking-evoked \DeltaHbT repsonses - individual animals')
%% LH short whisks
ax1 = subplot(3,2,1);
% Naives
for aa = 1:size(data.Naive.ShortWhisks.adjLH.HbT,1)
    plot(data.Naive.ShortWhisks.adjLH.meanTimeVector,data.Naive.ShortWhisks.adjLH.HbT(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.ShortWhisks.adjLH.HbT,1)
    plot(data.Blank_SAP.ShortWhisks.adjLH.meanTimeVector,data.Blank_SAP.ShortWhisks.adjLH.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.ShortWhisks.adjLH.HbT,1)
    plot(data.SSP_SAP.ShortWhisks.adjLH.meanTimeVector,data.SSP_SAP.ShortWhisks.adjLH.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('LH (UnRx) \DeltaHbT - Short Whisks')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
%% RH short whisks
ax2 = subplot(3,2,2);
% Naives
for aa = 1:size(data.Naive.ShortWhisks.adjRH.HbT,1)
    plot(data.Naive.ShortWhisks.adjRH.meanTimeVector,data.Naive.ShortWhisks.adjRH.HbT(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.ShortWhisks.adjRH.HbT,1)
    plot(data.Blank_SAP.ShortWhisks.adjRH.meanTimeVector,data.Blank_SAP.ShortWhisks.adjRH.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.ShortWhisks.adjRH.HbT,1)
    plot(data.SSP_SAP.ShortWhisks.adjRH.meanTimeVector,data.SSP_SAP.ShortWhisks.adjRH.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('RH (Rx) \DeltaHbT - Short Whisks')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% LH intermediate whisks
ax3 = subplot(3,2,3);
% Naives
for aa = 1:size(data.Naive.IntermediateWhisks.adjLH.HbT,1)
    plot(data.Naive.IntermediateWhisks.adjLH.meanTimeVector,data.Naive.IntermediateWhisks.adjLH.HbT(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.IntermediateWhisks.adjLH.HbT,1)
    plot(data.Blank_SAP.IntermediateWhisks.adjLH.meanTimeVector,data.Blank_SAP.IntermediateWhisks.adjLH.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.IntermediateWhisks.adjLH.HbT,1)
    plot(data.SSP_SAP.IntermediateWhisks.adjLH.meanTimeVector,data.SSP_SAP.IntermediateWhisks.adjLH.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('LH (UnRx) \DeltaHbT - Intermediate Whisks')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH intermediate whisks
ax4 = subplot(3,2,4);
% Naives
for aa = 1:size(data.Naive.IntermediateWhisks.adjRH.HbT,1)
    plot(data.Naive.IntermediateWhisks.adjRH.meanTimeVector,data.Naive.IntermediateWhisks.adjRH.HbT(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.IntermediateWhisks.adjRH.HbT,1)
    plot(data.Blank_SAP.IntermediateWhisks.adjRH.meanTimeVector,data.Blank_SAP.IntermediateWhisks.adjRH.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.IntermediateWhisks.adjRH.HbT,1)
    plot(data.SSP_SAP.IntermediateWhisks.adjRH.meanTimeVector,data.SSP_SAP.IntermediateWhisks.adjRH.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('RH (Rx) \DeltaHbT - Intermediate Whisks')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% LH long whisks
ax5 = subplot(3,2,5);
% Naives
for aa = 1:size(data.Naive.LongWhisks.adjLH.HbT,1)
    plot(data.Naive.LongWhisks.adjLH.meanTimeVector,data.Naive.LongWhisks.adjLH.HbT(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.LongWhisks.adjLH.HbT,1)
    plot(data.Blank_SAP.LongWhisks.adjLH.meanTimeVector,data.Blank_SAP.LongWhisks.adjLH.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.LongWhisks.adjLH.HbT,1)
    plot(data.SSP_SAP.LongWhisks.adjLH.meanTimeVector,data.SSP_SAP.LongWhisks.adjLH.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('LH (UnRx) \DeltaHbT - Long Whisks')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH long whisks
ax6 = subplot(3,2,6);
% Naives
for aa = 1:size(data.Naive.LongWhisks.adjRH.HbT,1)
    plot(data.Naive.LongWhisks.adjRH.meanTimeVector,data.Naive.LongWhisks.adjRH.HbT(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.LongWhisks.adjRH.HbT,1)
    plot(data.Blank_SAP.LongWhisks.adjRH.meanTimeVector,data.Blank_SAP.LongWhisks.adjRH.HbT(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.LongWhisks.adjRH.HbT,1)
    plot(data.SSP_SAP.LongWhisks.adjRH.meanTimeVector,data.SSP_SAP.LongWhisks.adjRH.HbT(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('RH (Rx) \DeltaHbT - Long Whisks')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% figure characteristics
linkaxes([ax1,ax2],'xy')
linkaxes([ax3,ax4],'xy')
linkaxes([ax5,ax6],'xy')
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Whisking Evoked - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure2,[dirpath 'IndividualWhiskEvoked_HbT']);
    set(summaryFigure2,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'IndividualWhiskEvoked_HbT'])
end
%% average whisk-evoked figures
summaryFigure3 = figure;
sgtitle('Whisking-evoked cortical MUA [300-3000 Hz]  repsonses')
%% LH short whisks
ax1 = subplot(3,2,1);
% Naives
p1 = plot(data.Naive.ShortWhisks.adjLH.meanTimeVector,data.Naive.ShortWhisks.adjLH.meanCortMUA,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.Naive.ShortWhisks.adjLH.meanTimeVector,data.Naive.ShortWhisks.adjLH.meanCortMUA + data.Naive.ShortWhisks.adjLH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.Naive.ShortWhisks.adjLH.meanTimeVector,data.Naive.ShortWhisks.adjLH.meanCortMUA - data.Naive.ShortWhisks.adjLH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
p2 = plot(data.Blank_SAP.ShortWhisks.adjLH.meanTimeVector,data.Blank_SAP.ShortWhisks.adjLH.meanCortMUA,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.ShortWhisks.adjLH.meanTimeVector,data.Blank_SAP.ShortWhisks.adjLH.meanCortMUA + data.Blank_SAP.ShortWhisks.adjLH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.ShortWhisks.adjLH.meanTimeVector,data.Blank_SAP.ShortWhisks.adjLH.meanCortMUA - data.Blank_SAP.ShortWhisks.adjLH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
p3 = plot(data.SSP_SAP.ShortWhisks.adjLH.meanTimeVector,data.SSP_SAP.ShortWhisks.adjLH.meanCortMUA,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.ShortWhisks.adjLH.meanTimeVector,data.SSP_SAP.ShortWhisks.adjLH.meanCortMUA + data.SSP_SAP.ShortWhisks.adjLH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.ShortWhisks.adjLH.meanTimeVector,data.SSP_SAP.ShortWhisks.adjLH.meanCortMUA - data.SSP_SAP.ShortWhisks.adjLH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
title('LH (UnRx) MUA [300-3000 Hz] - Short Whisks')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH short whisks
ax2 = subplot(3,2,2);
% Naives
plot(data.Naive.ShortWhisks.adjRH.meanTimeVector,data.Naive.ShortWhisks.adjRH.meanCortMUA,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.Naive.ShortWhisks.adjRH.meanTimeVector,data.Naive.ShortWhisks.adjRH.meanCortMUA + data.Naive.ShortWhisks.adjRH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.Naive.ShortWhisks.adjRH.meanTimeVector,data.Naive.ShortWhisks.adjRH.meanCortMUA - data.Naive.ShortWhisks.adjRH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
plot(data.Blank_SAP.ShortWhisks.adjRH.meanTimeVector,data.Blank_SAP.ShortWhisks.adjRH.meanCortMUA,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.ShortWhisks.adjRH.meanTimeVector,data.Blank_SAP.ShortWhisks.adjRH.meanCortMUA + data.Blank_SAP.ShortWhisks.adjRH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.ShortWhisks.adjRH.meanTimeVector,data.Blank_SAP.ShortWhisks.adjRH.meanCortMUA - data.Blank_SAP.ShortWhisks.adjRH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.ShortWhisks.adjRH.meanTimeVector,data.SSP_SAP.ShortWhisks.adjRH.meanCortMUA,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.ShortWhisks.adjRH.meanTimeVector,data.SSP_SAP.ShortWhisks.adjRH.meanCortMUA + data.SSP_SAP.ShortWhisks.adjRH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.ShortWhisks.adjRH.meanTimeVector,data.SSP_SAP.ShortWhisks.adjRH.meanCortMUA - data.SSP_SAP.ShortWhisks.adjRH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
title('RH (Rx) MUA [300-3000 Hz] - Short Whisks')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% LH intermediate whisks
ax3 = subplot(3,2,3);
% Naives
plot(data.Naive.IntermediateWhisks.adjLH.meanTimeVector,data.Naive.IntermediateWhisks.adjLH.meanCortMUA,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.Naive.IntermediateWhisks.adjLH.meanTimeVector,data.Naive.IntermediateWhisks.adjLH.meanCortMUA + data.Naive.IntermediateWhisks.adjLH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.Naive.IntermediateWhisks.adjLH.meanTimeVector,data.Naive.IntermediateWhisks.adjLH.meanCortMUA - data.Naive.IntermediateWhisks.adjLH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
plot(data.Blank_SAP.IntermediateWhisks.adjLH.meanTimeVector,data.Blank_SAP.IntermediateWhisks.adjLH.meanCortMUA,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.IntermediateWhisks.adjLH.meanTimeVector,data.Blank_SAP.IntermediateWhisks.adjLH.meanCortMUA + data.Blank_SAP.IntermediateWhisks.adjLH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.IntermediateWhisks.adjLH.meanTimeVector,data.Blank_SAP.IntermediateWhisks.adjLH.meanCortMUA - data.Blank_SAP.IntermediateWhisks.adjLH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.IntermediateWhisks.adjLH.meanTimeVector,data.SSP_SAP.IntermediateWhisks.adjLH.meanCortMUA,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.IntermediateWhisks.adjLH.meanTimeVector,data.SSP_SAP.IntermediateWhisks.adjLH.meanCortMUA + data.SSP_SAP.IntermediateWhisks.adjLH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.IntermediateWhisks.adjLH.meanTimeVector,data.SSP_SAP.IntermediateWhisks.adjLH.meanCortMUA - data.SSP_SAP.IntermediateWhisks.adjLH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
title('LH (UnRx) MUA [300-3000 Hz] - Intermediate Whisks')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH intermediate whisks
ax4 = subplot(3,2,4);
% Naives
plot(data.Naive.IntermediateWhisks.adjRH.meanTimeVector,data.Naive.IntermediateWhisks.adjRH.meanCortMUA,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.Naive.IntermediateWhisks.adjRH.meanTimeVector,data.Naive.IntermediateWhisks.adjRH.meanCortMUA + data.Naive.IntermediateWhisks.adjRH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.Naive.IntermediateWhisks.adjRH.meanTimeVector,data.Naive.IntermediateWhisks.adjRH.meanCortMUA - data.Naive.IntermediateWhisks.adjRH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
plot(data.Blank_SAP.IntermediateWhisks.adjRH.meanTimeVector,data.Blank_SAP.IntermediateWhisks.adjRH.meanCortMUA,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.IntermediateWhisks.adjRH.meanTimeVector,data.Blank_SAP.IntermediateWhisks.adjRH.meanCortMUA + data.Blank_SAP.IntermediateWhisks.adjRH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.IntermediateWhisks.adjRH.meanTimeVector,data.Blank_SAP.IntermediateWhisks.adjRH.meanCortMUA - data.Blank_SAP.IntermediateWhisks.adjRH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.IntermediateWhisks.adjRH.meanTimeVector,data.SSP_SAP.IntermediateWhisks.adjRH.meanCortMUA,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.IntermediateWhisks.adjRH.meanTimeVector,data.SSP_SAP.IntermediateWhisks.adjRH.meanCortMUA + data.SSP_SAP.IntermediateWhisks.adjRH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.IntermediateWhisks.adjRH.meanTimeVector,data.SSP_SAP.IntermediateWhisks.adjRH.meanCortMUA - data.SSP_SAP.IntermediateWhisks.adjRH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
title('RH (Rx) MUA [300-3000 Hz] - Intermediate Whisks')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% LH long whisks
ax5 = subplot(3,2,5);
% Naives
plot(data.Naive.LongWhisks.adjLH.meanTimeVector,data.Naive.LongWhisks.adjLH.meanCortMUA,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.Naive.LongWhisks.adjLH.meanTimeVector,data.Naive.LongWhisks.adjLH.meanCortMUA + data.Naive.LongWhisks.adjLH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.Naive.LongWhisks.adjLH.meanTimeVector,data.Naive.LongWhisks.adjLH.meanCortMUA - data.Naive.LongWhisks.adjLH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
plot(data.Blank_SAP.LongWhisks.adjLH.meanTimeVector,data.Blank_SAP.LongWhisks.adjLH.meanCortMUA,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.LongWhisks.adjLH.meanTimeVector,data.Blank_SAP.LongWhisks.adjLH.meanCortMUA + data.Blank_SAP.LongWhisks.adjLH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.LongWhisks.adjLH.meanTimeVector,data.Blank_SAP.LongWhisks.adjLH.meanCortMUA - data.Blank_SAP.LongWhisks.adjLH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.LongWhisks.adjLH.meanTimeVector,data.SSP_SAP.LongWhisks.adjLH.meanCortMUA,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.LongWhisks.adjLH.meanTimeVector,data.SSP_SAP.LongWhisks.adjLH.meanCortMUA + data.SSP_SAP.LongWhisks.adjLH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.LongWhisks.adjLH.meanTimeVector,data.SSP_SAP.LongWhisks.adjLH.meanCortMUA - data.SSP_SAP.LongWhisks.adjLH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
title('LH (UnRx) MUA [300-3000 Hz] - Long Whisks')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH long whisks
ax6 = subplot(3,2,6);
% Naives
plot(data.Naive.LongWhisks.adjRH.meanTimeVector,data.Naive.LongWhisks.adjRH.meanCortMUA,'color',colors('sapphire'),'LineWidth',2);
hold on
plot(data.Naive.LongWhisks.adjRH.meanTimeVector,data.Naive.LongWhisks.adjRH.meanCortMUA + data.Naive.LongWhisks.adjRH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
plot(data.Naive.LongWhisks.adjRH.meanTimeVector,data.Naive.LongWhisks.adjRH.meanCortMUA - data.Naive.LongWhisks.adjRH.stdCortMUA,'color',colors('sapphire'),'LineWidth',0.5)
% Blank-SAP
plot(data.Blank_SAP.LongWhisks.adjRH.meanTimeVector,data.Blank_SAP.LongWhisks.adjRH.meanCortMUA,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.LongWhisks.adjRH.meanTimeVector,data.Blank_SAP.LongWhisks.adjRH.meanCortMUA + data.Blank_SAP.LongWhisks.adjRH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.LongWhisks.adjRH.meanTimeVector,data.Blank_SAP.LongWhisks.adjRH.meanCortMUA - data.Blank_SAP.LongWhisks.adjRH.stdCortMUA,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
plot(data.SSP_SAP.LongWhisks.adjRH.meanTimeVector,data.SSP_SAP.LongWhisks.adjRH.meanCortMUA,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.LongWhisks.adjRH.meanTimeVector,data.SSP_SAP.LongWhisks.adjRH.meanCortMUA + data.SSP_SAP.LongWhisks.adjRH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.LongWhisks.adjRH.meanTimeVector,data.SSP_SAP.LongWhisks.adjRH.meanCortMUA - data.SSP_SAP.LongWhisks.adjRH.stdCortMUA,'color',colors('electric purple'),'LineWidth',0.5)
title('RH (Rx) MUA [300-3000 Hz] - Long Whisks')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% figure characteristics
linkaxes([ax1,ax2],'xy')
linkaxes([ax3,ax4],'xy')
linkaxes([ax5,ax6],'xy')
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Whisking Evoked - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure3,[dirpath 'AverageWhiskEvoked_MUA']);
    set(summaryFigure3,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'AverageWhiskEvoked_MUA'])
end
%% individual whisk-evoked figures
summaryFigure4 = figure;
sgtitle('Whisking-evoked cortical MUA [300-3000 Hz] repsonses - individual animals')
%% LH short whisks
ax1 = subplot(3,2,1);
% Naives
for aa = 1:size(data.Naive.ShortWhisks.adjLH.cortMUA,1)
    plot(data.Naive.ShortWhisks.adjLH.meanTimeVector,data.Naive.ShortWhisks.adjLH.cortMUA(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.ShortWhisks.adjLH.cortMUA,1)
    plot(data.Blank_SAP.ShortWhisks.adjLH.meanTimeVector,data.Blank_SAP.ShortWhisks.adjLH.cortMUA(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.ShortWhisks.adjLH.cortMUA,1)
    plot(data.SSP_SAP.ShortWhisks.adjLH.meanTimeVector,data.SSP_SAP.ShortWhisks.adjLH.cortMUA(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('LH (UnRx) MUA [300-3000 Hz] - Short Whisks')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
set(gca,'box','off')
xlim([-2,10])
%% RH short whisks
ax2 = subplot(3,2,2);
% Naives
for aa = 1:size(data.Naive.ShortWhisks.adjRH.cortMUA,1)
    plot(data.Naive.ShortWhisks.adjRH.meanTimeVector,data.Naive.ShortWhisks.adjRH.cortMUA(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.ShortWhisks.adjRH.cortMUA,1)
    plot(data.Blank_SAP.ShortWhisks.adjRH.meanTimeVector,data.Blank_SAP.ShortWhisks.adjRH.cortMUA(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.ShortWhisks.adjRH.cortMUA,1)
    plot(data.SSP_SAP.ShortWhisks.adjRH.meanTimeVector,data.SSP_SAP.ShortWhisks.adjRH.cortMUA(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('RH (Rx) MUA [300-3000 Hz] - Short Whisks')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% LH intermediate whisks
ax3 = subplot(3,2,3);
% Naives
for aa = 1:size(data.Naive.IntermediateWhisks.adjLH.cortMUA,1)
    plot(data.Naive.IntermediateWhisks.adjLH.meanTimeVector,data.Naive.IntermediateWhisks.adjLH.cortMUA(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.IntermediateWhisks.adjLH.cortMUA,1)
    plot(data.Blank_SAP.IntermediateWhisks.adjLH.meanTimeVector,data.Blank_SAP.IntermediateWhisks.adjLH.cortMUA(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.IntermediateWhisks.adjLH.cortMUA,1)
    plot(data.SSP_SAP.IntermediateWhisks.adjLH.meanTimeVector,data.SSP_SAP.IntermediateWhisks.adjLH.cortMUA(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('LH (UnRx) MUA [300-3000 Hz] - Intermediate Whisks')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH intermediate whisks
ax4 = subplot(3,2,4);
% Naives
for aa = 1:size(data.Naive.IntermediateWhisks.adjRH.cortMUA,1)
    plot(data.Naive.IntermediateWhisks.adjRH.meanTimeVector,data.Naive.IntermediateWhisks.adjRH.cortMUA(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.IntermediateWhisks.adjRH.cortMUA,1)
    plot(data.Blank_SAP.IntermediateWhisks.adjRH.meanTimeVector,data.Blank_SAP.IntermediateWhisks.adjRH.cortMUA(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.IntermediateWhisks.adjRH.cortMUA,1)
    plot(data.SSP_SAP.IntermediateWhisks.adjRH.meanTimeVector,data.SSP_SAP.IntermediateWhisks.adjRH.cortMUA(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('RH (Rx) MUA [300-3000 Hz] - Intermediate Whisks')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% LH long whisks
ax5 = subplot(3,2,5);
% Naives
for aa = 1:size(data.Naive.LongWhisks.adjLH.cortMUA,1)
    plot(data.Naive.LongWhisks.adjLH.meanTimeVector,data.Naive.LongWhisks.adjLH.cortMUA(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.LongWhisks.adjLH.cortMUA,1)
    plot(data.Blank_SAP.LongWhisks.adjLH.meanTimeVector,data.Blank_SAP.LongWhisks.adjLH.cortMUA(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.LongWhisks.adjLH.cortMUA,1)
    plot(data.SSP_SAP.LongWhisks.adjLH.meanTimeVector,data.SSP_SAP.LongWhisks.adjLH.cortMUA(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('LH (UnRx) MUA [300-3000 Hz] - Long Whisks')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% RH long whisks
ax6 = subplot(3,2,6);
% Naives
for aa = 1:size(data.Naive.LongWhisks.adjRH.cortMUA,1)
    plot(data.Naive.LongWhisks.adjRH.meanTimeVector,data.Naive.LongWhisks.adjRH.cortMUA(aa,:),'color',colors('sapphire'),'LineWidth',0.5);
    hold on
end
% Blank-SAP
for aa = 1:size(data.Blank_SAP.LongWhisks.adjRH.cortMUA,1)
    plot(data.Blank_SAP.LongWhisks.adjRH.meanTimeVector,data.Blank_SAP.LongWhisks.adjRH.cortMUA(aa,:),'color',colors('north texas green'),'LineWidth',0.5);
    hold on
end
% SSP-SAP
for aa = 1:size(data.SSP_SAP.LongWhisks.adjRH.cortMUA,1)
    plot(data.SSP_SAP.LongWhisks.adjRH.meanTimeVector,data.SSP_SAP.LongWhisks.adjRH.cortMUA(aa,:),'color',colors('electric purple'),'LineWidth',0.5);
    hold on
end
title('RH (Rx) MUA [300-3000 Hz] - Long Whisks')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
xlim([-2,10])
%% figure characteristics
linkaxes([ax1,ax2],'xy')
linkaxes([ax3,ax4],'xy')
linkaxes([ax5,ax6],'xy')
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Whisking Evoked - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure4,[dirpath 'IndividualWhiskEvoked_MUA']);
    set(summaryFigure4,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'IndividualWhiskEvoked_MUA'])
end
%% average neural responses
summaryFigure5 = figure;
sgtitle('Whisking-evoked cortical neural (LFP) repsonses')
%% brief whisks cortical LFP
ax1 = subplot(3,6,1);
imagesc(data.Naive.ShortWhisks.adjLH.meanCortT,data.Naive.ShortWhisks.adjLH.meanCortF,data.Naive.ShortWhisks.adjLH.meanCortS)
title('Naive short whisk LH')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% brief whisks cortical LFP
subplot(3,6,2);
imagesc(data.Naive.ShortWhisks.adjRH.meanCortT,data.Naive.ShortWhisks.adjRH.meanCortF,data.Naive.ShortWhisks.adjRH.meanCortS)
title('Naive short whisk RH')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% brief whisks cortical LFP
subplot(3,6,3);
imagesc(data.SSP_SAP.ShortWhisks.adjLH.meanCortT,data.SSP_SAP.ShortWhisks.adjLH.meanCortF,data.SSP_SAP.ShortWhisks.adjLH.meanCortS)
title('SSP-SAP short whisk LH')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% brief whisks cortical LFP
subplot(3,6,4);
imagesc(data.SSP_SAP.ShortWhisks.adjRH.meanCortT,data.SSP_SAP.ShortWhisks.adjRH.meanCortF,data.SSP_SAP.ShortWhisks.adjRH.meanCortS)
title('SSP-SAP short whisk RH')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% brief whisks cortical LFP
subplot(3,6,5);
imagesc(data.Blank_SAP.ShortWhisks.adjLH.meanCortT,data.Blank_SAP.ShortWhisks.adjLH.meanCortF,data.Blank_SAP.ShortWhisks.adjLH.meanCortS)
title('Blank-SAP short whisk LH')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% brief whisks cortical LFP
ax2 = subplot(3,6,6);
imagesc(data.Blank_SAP.ShortWhisks.adjRH.meanCortT,data.Blank_SAP.ShortWhisks.adjRH.meanCortF,data.Blank_SAP.ShortWhisks.adjRH.meanCortS)
title('Blank-SAP short whisk RH')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% intermediate whisks cortical LFP
subplot(3,6,7);
imagesc(data.Naive.IntermediateWhisks.adjLH.meanCortT,data.Naive.IntermediateWhisks.adjLH.meanCortF,data.Naive.IntermediateWhisks.adjLH.meanCortS)
title('Naive intermediate whisk LH')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% intermediate whisks cortical LFP
subplot(3,6,8);
imagesc(data.Naive.IntermediateWhisks.adjRH.meanCortT,data.Naive.IntermediateWhisks.adjRH.meanCortF,data.Naive.IntermediateWhisks.adjRH.meanCortS)
title('Naive intermediate whisk RH')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% intermediate whisks cortical LFP
subplot(3,6,9);
imagesc(data.SSP_SAP.IntermediateWhisks.adjLH.meanCortT,data.SSP_SAP.IntermediateWhisks.adjLH.meanCortF,data.SSP_SAP.IntermediateWhisks.adjLH.meanCortS)
title('SSP-SAP intermediate whisk LH')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% intermediate whisks cortical LFP
subplot(3,6,10);
imagesc(data.SSP_SAP.IntermediateWhisks.adjRH.meanCortT,data.SSP_SAP.IntermediateWhisks.adjRH.meanCortF,data.SSP_SAP.IntermediateWhisks.adjRH.meanCortS)
title('SSP-SAP intermediate whisk RH')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% intermediate whisks cortical LFP
subplot(3,6,11);
imagesc(data.Blank_SAP.IntermediateWhisks.adjLH.meanCortT,data.Blank_SAP.IntermediateWhisks.adjLH.meanCortF,data.Blank_SAP.IntermediateWhisks.adjLH.meanCortS)
title('Blank-SAP intermediate whisk LH')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% intermediate whisks cortical LFP
ax3 = subplot(3,6,12);
imagesc(data.Blank_SAP.IntermediateWhisks.adjRH.meanCortT,data.Blank_SAP.IntermediateWhisks.adjRH.meanCortF,data.Blank_SAP.IntermediateWhisks.adjRH.meanCortS)
title('Blank-SAP intermediate whisk RH')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
c2 = colorbar;
ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% long whisks cortical LFP
subplot(3,6,13);
imagesc(data.Naive.LongWhisks.adjLH.meanCortT,data.Naive.LongWhisks.adjLH.meanCortF,data.Naive.LongWhisks.adjLH.meanCortS)
title('Naive long whisk LH')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% long whisks cortical LFP
subplot(3,6,14);
imagesc(data.Naive.LongWhisks.adjRH.meanCortT,data.Naive.LongWhisks.adjRH.meanCortF,data.Naive.LongWhisks.adjRH.meanCortS)
title('Naive long whisk RH')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% long whisks cortical LFP
subplot(3,6,15);
imagesc(data.SSP_SAP.LongWhisks.adjLH.meanCortT,data.SSP_SAP.LongWhisks.adjLH.meanCortF,data.SSP_SAP.LongWhisks.adjLH.meanCortS)
title('SSP-SAP long whisk LH')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% long whisks cortical LFP
subplot(3,6,16);
imagesc(data.SSP_SAP.LongWhisks.adjRH.meanCortT,data.SSP_SAP.LongWhisks.adjRH.meanCortF,data.SSP_SAP.LongWhisks.adjRH.meanCortS)
title('SSP-SAP long whisk RH')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% long whisks cortical LFP
subplot(3,6,17);
imagesc(data.Blank_SAP.LongWhisks.adjLH.meanCortT,data.Blank_SAP.LongWhisks.adjLH.meanCortF,data.Blank_SAP.LongWhisks.adjLH.meanCortS)
title('Blank-SAP long whisk LH')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% long whisks cortical LFP
ax4 = subplot(3,6,18);
imagesc(data.Blank_SAP.LongWhisks.adjRH.meanCortT,data.Blank_SAP.LongWhisks.adjRH.meanCortF,data.Blank_SAP.LongWhisks.adjRH.meanCortS)
title('Blank-SAP long whisk RH')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis xy
set(gca,'box','off')
%% adjust colorbar positions
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax2Pos(3:4) = ax1Pos(3:4);
ax3Pos(3:4) = ax1Pos(3:4);
ax4Pos(3:4) = ax1Pos(3:4);
set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
set(ax4,'position',ax4Pos);
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Whisking Evoked - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure5,[dirpath 'AverageWhiskEvoked_LFP']);
    set(summaryFigure5,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'AverageWhiskEvoked_LFP'])
end

end
