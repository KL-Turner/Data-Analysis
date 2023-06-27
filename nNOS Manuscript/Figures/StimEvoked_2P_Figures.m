function [] = StimEvoked_2P_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Evoked_2P';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'SSP_SAP','Blank_SAP'};
solenoids = {'LPadSol','RPadSol','AudSol'};
dataTypes = {'diameter','baseline','timeVector','count'};
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Evoked_2P.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        vIDs = fieldnames(Results_Evoked_2P.(group).(animalID));
        for cc = 1:length(vIDs)
            vID = vIDs{cc,1};
            for dd = 1:length(solenoids)
                solenoid = solenoids{1,dd};
                data.(group).(solenoid).dummCheck = 1;
                for ee = 1:length(dataTypes)
                    dataType = dataTypes{1,ee};
                    if isfield(data.(group).(solenoid),dataType) == false
                        data.(group).(solenoid).(dataType) = [];
                    end
                end
                data.(group).(solenoid).diameter = cat(1,data.(group).(solenoid).diameter,Results_Evoked_2P.(group).(animalID).(vID).Stim.(solenoid).diameter);
                data.(group).(solenoid).baseline = cat(1,data.(group).(solenoid).baseline,Results_Evoked_2P.(group).(animalID).(vID).Stim.(solenoid).baseline);
                data.(group).(solenoid).count = cat(1,data.(group).(solenoid).count,Results_Evoked_2P.(group).(animalID).(vID).Stim.(solenoid).count);
                data.(group).(solenoid).timeVector = cat(1,data.(group).(solenoid).timeVector,Results_Evoked_2P.(group).(animalID).(vID).Stim.(solenoid).timeVector);
            end
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
            data.(group).(comparison).(['stdErr_' dataType]) = std(data.(group).(solenoid).(dataType),0,1)./sqrt(size(data.(group).(solenoid).(dataType),1));
        end
    end
end
% figure
comparisons = {'ipsi','contra','aud'};
summaryFigure = figure;
sgtitle('Whisker stimulation [2P]')
for aa = 1:length(comparisons)
    comparison = comparisons{1,aa};
    subplot(1,3,aa)
    p1 = plot(data.Blank_SAP.(comparison).mean_timeVector,data.Blank_SAP.(comparison).mean_diameter,'color',colors('north texas green'),'LineWidth',2);
    hold on;
    plot(data.Blank_SAP.(comparison).mean_timeVector,data.Blank_SAP.(comparison).mean_diameter + data.Blank_SAP.(comparison).stdErr_diameter,'color',colors('north texas green'),'LineWidth',0.25)
    plot(data.Blank_SAP.(comparison).mean_timeVector,data.Blank_SAP.(comparison).mean_diameter - data.Blank_SAP.(comparison).stdErr_diameter,'color',colors('north texas green'),'LineWidth',0.25)
    p2 = plot(data.SSP_SAP.(comparison).mean_timeVector,data.SSP_SAP.(comparison).mean_diameter,'color',colors('electric purple'),'LineWidth',2);
    plot(data.SSP_SAP.(comparison).mean_timeVector,data.SSP_SAP.(comparison).mean_diameter + data.SSP_SAP.(comparison).stdErr_diameter,'color',colors('electric purple'),'LineWidth',0.25)
    plot(data.SSP_SAP.(comparison).mean_timeVector,data.SSP_SAP.(comparison).mean_diameter - data.SSP_SAP.(comparison).stdErr_diameter,'color',colors('electric purple'),'LineWidth',0.25)
    title(comparison)
    ylabel('\DeltaD/D (%)')
    xlabel('Peri-stimulus time (s)')
    legend([p1,p2],'Blank-SAP','SSP-SAP')
    set(gca,'box','off')
    xlim([-2,10])
    axis square
end
% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures' delim 'Stimulus Evoked' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'StimEvoked_2P']);
end