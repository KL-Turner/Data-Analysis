function [] = StimEvoked_Pulse_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Evoked_Pulse';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Blank_SAP','SSP_SAP'};
solenoids = {'LPadSol','RPadSol','AudSol'};
dataTypes = {'HbT','timeVector'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Evoked_Pulse.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(solenoids)
            solenoid = solenoids{1,cc};
            data.(group).(solenoid).dummCheck = 1;
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                if isfield(data.(group).(solenoid),dataType) == false
                    data.(group).(solenoid).(dataType) = [];
                end
            end
            data.(group).(solenoid).HbT = cat(1,data.(group).(solenoid).HbT,Results_Evoked_Pulse.(group).(animalID).Stim.(solenoid).HbT);
            data.(group).(solenoid).timeVector = cat(1,data.(group).(solenoid).timeVector,Results_Evoked_Pulse.(group).(animalID).Stim.(solenoid).timeVector);
        end
    end
end
% pair stimulation types with hemisphere
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(solenoids)
        solenoid = solenoids{1,bb};
        [comparison] = FindSolenoidComparison('RH',solenoid);
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            data.(group).(comparison).(dataType) = data.(group).(solenoid).(dataType);
            data.(group).(comparison).(['mean_'  dataType]) = mean(data.(group).(solenoid).(dataType),1);
            data.(group).(comparison).(['stdErr_' dataType]) = std(data.(group).(solenoid).(dataType),1)./sqrt(size(data.(group).(solenoid).(dataType),1));
        end
    end
end
% figure
comparisons = {'ipsi','contra','aud'};
summaryFigure = figure;
sgtitle('Whisker stimulation [Pulse]')
for aa = 1:length(comparisons)
    comparison = comparisons{1,aa};
    p1 = plot(data.Blank_SAP.(comparison).mean_timeVector,data.Blank_SAP.(comparison).mean_HbT,'color',colors('north texas green'),'LineWidth',2);
    hold on;
    plot(data.Blank_SAP.(comparison).mean_timeVector,data.Blank_SAP.(comparison).mean_HbT + data.Blank_SAP.(comparison).stdErr_HbT,'color',colors('north texas green'),'LineWidth',0.25)
    plot(data.Blank_SAP.(comparison).mean_timeVector,data.Blank_SAP.(comparison).mean_HbT - data.Blank_SAP.(comparison).stdErr_HbT,'color',colors('north texas green'),'LineWidth',0.25)
    p2 = plot(data.SSP_SAP.(comparison).mean_timeVector,data.SSP_SAP.(comparison).mean_HbT,'color',colors('electric purple'),'LineWidth',2);
    plot(data.SSP_SAP.(comparison).mean_timeVector,data.SSP_SAP.(comparison).mean_HbT + data.SSP_SAP.(comparison).stdErr_HbT,'color',colors('electric purple'),'LineWidth',0.25)
    plot(data.SSP_SAP.(comparison).mean_timeVector,data.SSP_SAP.(comparison).mean_HbT - data.SSP_SAP.(comparison).stdErr_HbT,'color',colors('electric purple'),'LineWidth',0.25)
    title(comparison)
    ylabel('\DeltaHbT (\muM)')
    xlabel('Peri-stimulus time (s)')
    legend([p1,p2],'Blank-SAP','SSP-SAP')
    set(gca,'box','off')
    xlim([-2,10])
    axis square
    % save figure(s)
    if saveFigs == true
        dirpath = [rootFolder delim 'Summary Figures' delim 'Stimulus Evoked' delim];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(summaryFigure,[dirpath 'StimEvoked_Pulse_' comparison]);
    end
end