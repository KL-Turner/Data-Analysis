function [] = PupilStimEvoked_Ephys_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PupilEvoked_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Naive','Blank_SAP','SSP_SAP'};
dataTypes = {'mmArea','mmDiameter','zArea','zDiameter'};
solenoids = {'LPadSol','RPadSol','AudSol'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_PupilEvoked_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(solenoids)
                solenoid = solenoids{1,dd};
                data.(group).(dataType).(solenoid).dummCheck = 1;
                if isfield(data.(group).(dataType).(solenoid),'data') == false
                    data.(group).(dataType).(solenoid).data = [];
                    data.(group).(dataType).(solenoid).group = {};
                    data.(group).(dataType).(solenoid).animalID = {};
                end
                data.(group).(dataType).(solenoid).data = cat(1,data.(group).(dataType).(solenoid).data,Results_PupilEvoked_Ephys.(group).(animalID).(dataType).Stim.(solenoid).mean);
                data.(group).(dataType).(solenoid).group = cat(1,data.(group).(dataType).(solenoid).group,group);
                data.(group).(dataType).(solenoid).animalID = cat(1,data.(group).(dataType).(solenoid).animalID,animalID);
            end
        end
    end
end
% pair stimulation with hemisphere
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for cc = 1:length(solenoids)
            solenoid = solenoids{1,cc};
            [comparison] = FindSolenoidComparison('Both',solenoid);
            data.(group).(dataType).(comparison).meanData = mean(data.(group).(dataType).(solenoid).data,1);
            data.(group).(dataType).(comparison).stdErrData = std(data.(group).(dataType).(solenoid).data,0,1)./sqrt(size(data.(group).(dataType).(solenoid).data,1));
        end
    end
end
% figure
timeVector = (0:12*30)/30 - 2;
comparisons = {'contra','aud'};
labels = {'mm^2','mm','z-units','z-units'};
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    summaryFigure = figure;
    sgtitle(['Pupil whisker stimlation ' dataType ' [Ephys]'])
    for bb = 1:length(comparisons)
        subplot(1,2,bb);
        comparison = comparisons{1,bb};
        p1 = plot(timeVector,data.Naive.(dataType).(comparison).meanData,'color',colors('sapphire'),'LineWidth',2);
        hold on;
        plot(timeVector,data.Naive.(dataType).(comparison).meanData + data.Naive.(dataType).(comparison).stdErrData,'color',colors('sapphire'),'LineWidth',0.25)
        plot(timeVector,data.Naive.(dataType).(comparison).meanData - data.Naive.(dataType).(comparison).stdErrData,'color',colors('sapphire'),'LineWidth',0.25)
        p2 = plot(timeVector,data.Blank_SAP.(dataType).(comparison).meanData,'color',colors('north texas green'),'LineWidth',2);
        plot(timeVector,data.Blank_SAP.(dataType).(comparison).meanData + data.Blank_SAP.(dataType).(comparison).stdErrData,'color',colors('north texas green'),'LineWidth',0.25)
        plot(timeVector,data.Blank_SAP.(dataType).(comparison).meanData - data.Blank_SAP.(dataType).(comparison).stdErrData,'color',colors('north texas green'),'LineWidth',0.25)
        p3 = plot(timeVector,data.SSP_SAP.(dataType).(comparison).meanData,'color',colors('electric purple'),'LineWidth',2);
        plot(timeVector,data.SSP_SAP.(dataType).(comparison).meanData + data.SSP_SAP.(dataType).(comparison).stdErrData,'color',colors('electric purple'),'LineWidth',0.25)
        plot(timeVector,data.SSP_SAP.(dataType).(comparison).meanData - data.SSP_SAP.(dataType).(comparison).stdErrData,'color',colors('electric purple'),'LineWidth',0.25)
        title(comparison)
        ylabel(['\Delta' labels{1,aa}])
        xlabel('Peri-stimulus time (s)')
        if bb == 1
            legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
        end
        set(gca,'box','off')
        xlim([-2,10])
        axis square
    end
    linkaxes
    % save figure(s)
    if saveFigs == true
        dirpath = [rootFolder delim 'Summary Figures' delim 'Behavior' delim];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(summaryFigure,[dirpath 'PupilStimEvoked_Ephys_' dataType]);
    end
end