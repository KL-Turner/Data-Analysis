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
            data.(group).(transition).dummCheck = 1;
            if isfield(data.(group).(transition),'HbT') == false
                data.(group).(transition).HbT = [];
                data.(group).(transition).difference = [];
            end
            data.(group).(transition).HbT = cat(1,data.(group).(transition).HbT,Results_Transitions_Ephys.(group).(animalID).RH.(transition).HbT);
            leadHbT = Results_Transitions_Ephys.(group).(animalID).RH.(transition).HbT(1:30*20 + 1);
            lagHbT = Results_Transitions_Ephys.(group).(animalID).RH.(transition).HbT(end - 30*20:end);
            data.(group).(transition).difference = cat(1,data.(group).(transition).difference,leadHbT - lagHbT);
        end
    end
end
% take average for each behavioral transition
for qq = 1:length(groups)
    group = groups{1,qq};
    for cc = 1:length(transitions)
        transition = transitions{1,cc};
        data.(group).(transition).meanHbT = mean(data.(group).(transition).HbT,1);
        data.(group).(transition).stdHbT = std(data.(group).(transition).HbT,0,1);
        data.(group).(transition).meanDifference = mean(data.(group).(transition).difference,1);
        data.(group).(transition).stdDifference = std(data.(group).(transition).difference,0,1);
    end
end
T1 = -30 + (1/30):(1/30):30;
%% Arousal-state transitions
figure;
for aa = 1:length(transitions)
    transition = transitions{1,aa};
    subplot(2,2,aa);
    p1 = plot(T1,data.Blank_SAP.(transition).meanHbT,'-','color',colors('black'),'LineWidth',2);
    hold on
    plot(T1,data.Blank_SAP.(transition).meanHbT + data.Blank_SAP.(transition).stdHbT,'-','color',colors('black'),'LineWidth',0.5);
    plot(T1,data.Blank_SAP.(transition).meanHbT - data.Blank_SAP.(transition).stdHbT,'-','color',colors('black'),'LineWidth',0.5);
    p2 = plot(T1,data.SSP_SAP.(transition).meanHbT,'-','color',colors('dark candy apple red'),'LineWidth',2);
    plot(T1,data.SSP_SAP.(transition).meanHbT + data.SSP_SAP.(transition).stdHbT,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
    plot(T1,data.SSP_SAP.(transition).meanHbT - data.SSP_SAP.(transition).stdHbT,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
    ylabel('\Delta[HbT] (\muM)')
    title(transition)
    legend([p1,p2],'Blank','SP')
    xlim([-30,30])
end

%% Arousal-state transitions differencefigure;
for aa = 1:length(transitions)
    transition = transitions{1,aa};
    subplot(2,2,aa);
    p1 = plot(T1,data.Blank_SAP.(transition).meanHbT,'-','color',colors('black'),'LineWidth',2);
    hold on
    plot(T1,data.Blank_SAP.(transition).meanHbT + data.Blank_SAP.(transition).stdHbT,'-','color',colors('black'),'LineWidth',0.5);
    plot(T1,data.Blank_SAP.(transition).meanHbT - data.Blank_SAP.(transition).stdHbT,'-','color',colors('black'),'LineWidth',0.5);
    p2 = plot(T1,data.SSP_SAP.(transition).meanHbT,'-','color',colors('dark candy apple red'),'LineWidth',2);
    plot(T1,data.SSP_SAP.(transition).meanHbT + data.SSP_SAP.(transition).stdHbT,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
    plot(T1,data.SSP_SAP.(transition).meanHbT - data.SSP_SAP.(transition).stdHbT,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
    ylabel('\Delta[HbT] (\muM)')
    title(transition)
    legend([p1,p2],'Blank','SP')
    xlim([-30,30])
end