function [] = FigS6_nNOS(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------

cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Transitions_Ephys';
load(resultsStruct);
cd(rootFolder)
groups = {'Blank_SAP','SSP_SAP'};
transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
%% IOS mean transitions between each arousal-state
% cd through each animal's directory and extract the appropriate analysis results
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Transitions_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(transitions)
            transition = transitions{1,cc};
            % pre-allocate necessary variable fields
            transitionData.(group).(transition).dummCheck = 1;
            if isfield(transitionData.(group).(transition),'HbT') == false
                transitionData.(group).(transition).HbT = [];
                transitionData.(group).(transition).difference = [];
                transitionData.(group).(transition).animalID = {};
                transitionData.(group).(transition).group = {};
            end
            transitionData.(group).(transition).HbT = cat(1,transitionData.(group).(transition).HbT,Results_Transitions_Ephys.(group).(animalID).RH.(transition).HbT);
            leadHbT = mean(Results_Transitions_Ephys.(group).(animalID).RH.(transition).HbT(1:30*20 + 1));
            lagHbT = mean(Results_Transitions_Ephys.(group).(animalID).RH.(transition).HbT(end - 30*20:end));
            transitionData.(group).(transition).difference = cat(1,transitionData.(group).(transition).difference,leadHbT - lagHbT);
            transitionData.(group).(transition).animalID = cat(1,transitionData.(group).(transition).animalID,animalID);
            transitionData.(group).(transition).group = cat(1,transitionData.(group).(transition).group,group);
        end
    end
end
% take average for each behavioral transition
for qq = 1:length(groups)
    group = groups{1,qq};
    for cc = 1:length(transitions)
        transition = transitions{1,cc};
        transitionData.(group).(transition).meanHbT = mean(transitionData.(group).(transition).HbT,1);
        transitionData.(group).(transition).stdHbT = std(transitionData.(group).(transition).HbT,0,1);
        transitionData.(group).(transition).meanDifference = mean(transitionData.(group).(transition).difference,1);
        transitionData.(group).(transition).stdDifference = std(transitionData.(group).(transition).difference,0,1);
    end
end

% GLME comparing peak correlation
for bb = 1:length(transitions)
    transition = transitions{1,bb};
    transitionStats.(transition).tableSize = cat(1,transitionData.Blank_SAP.(transition).difference,transitionData.SSP_SAP.(transition).difference);
    transitionStats.(transition).Table = table('Size',[size(transitionStats.(transition).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'AnimalID','Treatment','Difference'});
    transitionStats.(transition).Table.AnimalID = cat(1,transitionData.Blank_SAP.(transition).animalID,transitionData.SSP_SAP.(transition).animalID);
    transitionStats.(transition).Table.Treatment = cat(1,transitionData.Blank_SAP.(transition).group,transitionData.SSP_SAP.(transition).group);
    transitionStats.(transition).Table.Difference = cat(1,transitionData.Blank_SAP.(transition).difference,transitionData.SSP_SAP.(transition).difference);
    transitionStats.(transition).FitFormula = 'Difference ~ 1 + Treatment + (1|AnimalID)';
    transitionStats.(transition).Stats = fitglme(transitionStats.(transition).Table,transitionStats.(transition).FitFormula);
end
%% Arousal-state transitions
T1 = -30 + (1/30):(1/30):30;
figure;
for aa = 1:length(transitions)
    transition = transitions{1,aa};
    subplot(2,2,aa);
    p1 = plot(T1,transitionData.Blank_SAP.(transition).meanHbT,'-','color',colors('black'),'LineWidth',2);
    hold on
    plot(T1,transitionData.Blank_SAP.(transition).meanHbT + transitionData.Blank_SAP.(transition).stdHbT,'-','color',colors('black'),'LineWidth',0.5);
    plot(T1,transitionData.Blank_SAP.(transition).meanHbT - transitionData.Blank_SAP.(transition).stdHbT,'-','color',colors('black'),'LineWidth',0.5);
    p2 = plot(T1,transitionData.SSP_SAP.(transition).meanHbT,'-','color',colors('dark candy apple red'),'LineWidth',2);
    plot(T1,transitionData.SSP_SAP.(transition).meanHbT + transitionData.SSP_SAP.(transition).stdHbT,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
    plot(T1,transitionData.SSP_SAP.(transition).meanHbT - transitionData.SSP_SAP.(transition).stdHbT,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
    ylabel('\Delta[HbT] (\muM)')
    title(transition)
    legend([p1,p2],'Blank','SP')
    xlim([-30,30])
end

disp(transitionStats.AWAKEtoNREM.Stats)
disp(transitionStats.NREMtoAWAKE.Stats)
disp(transitionStats.NREMtoREM.Stats)
disp(transitionStats.REMtoAWAKE.Stats)

