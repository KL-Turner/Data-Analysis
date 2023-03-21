function [] = StimEvoked_GCaMP_Figures(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Evoked_GCaMP';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'SSP_SAP','Blank_SAP'};
solenoids = {'LPadSol','RPadSol','AudSol'};
hemispheres = {'LH','RH','frontalLH','frontalRH'};
dataTypes = {'HbT','GCaMP','timeVector'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Evoked_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(solenoids)
                solenoid = solenoids{1,dd};
                data.(group).(hemisphere).(solenoid).dummCheck = 1;
                for ee = 1:length(dataTypes)
                    dataType = dataTypes{1,ee};
                    if isfield(data.(group).(hemisphere).(solenoid),dataType) == false
                        data.(group).(hemisphere).(solenoid).(dataType) = [];
                    end
                end
                data.(group).(hemisphere).(solenoid).HbT = cat(1,data.(group).(hemisphere).(solenoid).HbT,Results_Evoked_GCaMP.(group).(animalID).Stim.(hemisphere).(solenoid).HbT);
                data.(group).(hemisphere).(solenoid).GCaMP = cat(1,data.(group).(hemisphere).(solenoid).GCaMP,Results_Evoked_GCaMP.(group).(animalID).Stim.(hemisphere).(solenoid).GCaMP);
                data.(group).(hemisphere).(solenoid).timeVector = cat(1,data.(group).(hemisphere).(solenoid).timeVector,Results_Evoked_GCaMP.(group).(animalID).Stim.(hemisphere).(solenoid).timeVector);
            end
        end
    end
end
% pair stimulation with hemisphere
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(solenoids)
            solenoid = solenoids{1,cc};
            [comparison] = FindSolenoidComparison(hemisphere,solenoid);
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                data.(group).(hemisphere).(comparison).(dataType) = data.(group).(hemisphere).(solenoid).(dataType);
                data.(group).(hemisphere).(comparison).(['mean_' dataType]) = mean(data.(group).(hemisphere).(solenoid).(dataType),1);
                data.(group).(hemisphere).(comparison).(['stdErr_' dataType]) = std(data.(group).(hemisphere).(solenoid).(dataType),1)./sqrt(size(data.(group).(hemisphere).(solenoid).(dataType),1));
            end
        end
    end
end
% figure
figure;
sgtitle('GCaMP contralateral whisker stimlation')
for aa = 1:length(hemispheres)*2
    if aa <= length(hemispheres)
        ax(aa) = subplot(2,4,aa);
        hemisphere = hemispheres{1,aa};
        dataType = 'HbT';
        label = '\DeltaHbT (\muM)';
    else
        bx(aa) = subplot(2,4,aa);
        hemisphere = hemispheres{1,aa - length(hemispheres)};
        dataType = 'GCaMP';
        label = '\DeltaF/F (%)';
    end
    x0 = xline(0,'k');
    hold on
    x1 = xline(5,'r');
    % Blank-SAP
    p1 = plot(data.Blank_SAP.(hemisphere).contra.mean_timeVector,data.Blank_SAP.(hemisphere).contra.(['mean_' dataType]),'color',colors('north texas green'),'LineWidth',2);
    plot(data.Blank_SAP.(hemisphere).contra.mean_timeVector,data.Blank_SAP.(hemisphere).contra.(['mean_' dataType]) + data.Blank_SAP.(hemisphere).contra.(['stdErr_' dataType]),'color',colors('north texas green'),'LineWidth',0.5)
    plot(data.Blank_SAP.(hemisphere).contra.mean_timeVector,data.Blank_SAP.(hemisphere).contra.(['mean_' dataType]) - data.Blank_SAP.(hemisphere).contra.(['stdErr_' dataType]),'color',colors('north texas green'),'LineWidth',0.5)
    % SSP-SAP
    p2 = plot(data.SSP_SAP.(hemisphere).contra.mean_timeVector,data.SSP_SAP.(hemisphere).contra.(['mean_' dataType]),'color',colors('electric purple'),'LineWidth',2);
    plot(data.SSP_SAP.(hemisphere).contra.mean_timeVector,data.SSP_SAP.(hemisphere).contra.(['mean_' dataType]) + data.SSP_SAP.(hemisphere).contra.(['stdErr_' dataType]),'color',colors('electric purple'),'LineWidth',0.5)
    plot(data.SSP_SAP.(hemisphere).contra.mean_timeVector,data.SSP_SAP.(hemisphere).contra.(['mean_' dataType]) - data.SSP_SAP.(hemisphere).contra.(['stdErr_' dataType]),'color',colors('electric purple'),'LineWidth',0.5)
    title([hemisphere ' ' dataType])
    ylabel(label)
    xlabel('Peri-stimulus time (s)')
    if aa == 1
        legend([p1,p2,x0,x1],'Blank-SAP','SSP-SAP','stimOn','stimOff')
    end
    set(gca,'box','off')
    xlim([-2,10])
    axis square
end
linkaxes(ax)
linkaxes(bx)
