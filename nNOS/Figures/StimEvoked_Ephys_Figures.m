function [] = StimEvoked_Ephys_Figures(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: 
%________________________________________________________________________________________________________________________

cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Evoked_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Naive','SSP_SAP','Blank_SAP'};
solenoids = {'LPadSol','RPadSol','AudSol'};
hemispheres = {'LH','RH'};
dataTypes = {'HbT','cortMUA','cortGam','cortS','T','F','timeVector'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Evoked_Ephys.(group));
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
                data.(group).(hemisphere).(solenoid).HbT = cat(1,data.(group).(hemisphere).(solenoid).HbT,Results_Evoked_Ephys.(group).(animalID).Stim.(hemisphere).(solenoid).HbT);
                data.(group).(hemisphere).(solenoid).cortMUA = cat(1,data.(group).(hemisphere).(solenoid).cortMUA,Results_Evoked_Ephys.(group).(animalID).Stim.(hemisphere).(solenoid).cortMUA);
                data.(group).(hemisphere).(solenoid).cortGam = cat(1,data.(group).(hemisphere).(solenoid).cortGam,Results_Evoked_Ephys.(group).(animalID).Stim.(hemisphere).(solenoid).cortGam);
                data.(group).(hemisphere).(solenoid).cortS = cat(3,data.(group).(hemisphere).(solenoid).cortS,Results_Evoked_Ephys.(group).(animalID).Stim.(hemisphere).(solenoid).cortLFP);
                data.(group).(hemisphere).(solenoid).T = cat(1,data.(group).(hemisphere).(solenoid).T,Results_Evoked_Ephys.(group).(animalID).Stim.(hemisphere).(solenoid).T);
                data.(group).(hemisphere).(solenoid).F = cat(1,data.(group).(hemisphere).(solenoid).F,Results_Evoked_Ephys.(group).(animalID).Stim.(hemisphere).(solenoid).F);
                data.(group).(hemisphere).(solenoid).timeVector = cat(1,data.(group).(hemisphere).(solenoid).timeVector,Results_Evoked_Ephys.(group).(animalID).Stim.(hemisphere).(solenoid).timeVector);
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
                if strcmp(dataType,'cortS')
                    data.(group).(hemisphere).(comparison).(dataType) = data.(group).(hemisphere).(solenoid).(dataType);
                    data.(group).(hemisphere).(comparison).(['mean_' dataType]) = mean(data.(group).(hemisphere).(solenoid).(dataType),3);
                else
                    data.(group).(hemisphere).(comparison).(dataType) = data.(group).(hemisphere).(solenoid).(dataType);
                    data.(group).(hemisphere).(comparison).(['mean_' dataType]) = mean(data.(group).(hemisphere).(solenoid).(dataType),1);
                    data.(group).(hemisphere).(comparison).(['stdErr_' dataType]) = std(data.(group).(hemisphere).(solenoid).(dataType),1)./sqrt(size(data.(group).(hemisphere).(solenoid).(dataType),1));
                end
            end
        end
    end
end
% figure
figure;
sgtitle('Ephys contralateral whisker stimlation (HbT)')
for aa = 1:2
    ax(aa) = subplot(1,2,aa);
    hemisphere = hemispheres{1,aa};
    dataType = 'HbT';
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
    % Naive
    p3 = plot(data.Naive.(hemisphere).contra.mean_timeVector,data.Naive.(hemisphere).contra.(['mean_' dataType]),'color',colors('sapphire'),'LineWidth',2);
    plot(data.Naive.(hemisphere).contra.mean_timeVector,data.Naive.(hemisphere).contra.(['mean_' dataType]) + data.Naive.(hemisphere).contra.(['stdErr_' dataType]),'color',colors('sapphire'),'LineWidth',0.5)
    plot(data.Naive.(hemisphere).contra.mean_timeVector,data.Naive.(hemisphere).contra.(['mean_' dataType]) - data.Naive.(hemisphere).contra.(['stdErr_' dataType]),'color',colors('sapphire'),'LineWidth',0.5)
    title([hemisphere ' ' dataType])
    label = '\DeltaHbT (\muM)';
    ylabel(label)
    xlabel('Peri-stimulus time (s)')
    if aa == 1
        legend([p1,p2,p3,x0,x1],'Blank-SAP','SSP-SAP','Naive','stimOn','stimOff')
    end
    set(gca,'box','off')
    xlim([-2,10])
    axis square
end
linkaxes(ax)
% figure
figure;
sgtitle('Ephys contralateral whisker stimlation (gamma)')
for aa = 1:2
    bx(aa) = subplot(1,2,aa);
    hemisphere = hemispheres{1,aa};
    dataType = 'cortGam';
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
    % Naive
    p3 = plot(data.Naive.(hemisphere).contra.mean_timeVector,data.Naive.(hemisphere).contra.(['mean_' dataType]),'color',colors('sapphire'),'LineWidth',2);
    plot(data.Naive.(hemisphere).contra.mean_timeVector,data.Naive.(hemisphere).contra.(['mean_' dataType]) + data.Naive.(hemisphere).contra.(['stdErr_' dataType]),'color',colors('sapphire'),'LineWidth',0.5)
    plot(data.Naive.(hemisphere).contra.mean_timeVector,data.Naive.(hemisphere).contra.(['mean_' dataType]) - data.Naive.(hemisphere).contra.(['stdErr_' dataType]),'color',colors('sapphire'),'LineWidth',0.5)
    title([hemisphere ' ' dataType])
    label = '\DeltaP/P (%)';
    ylabel(label)
    xlabel('Peri-stimulus time (s)')
    if aa == 1
        legend([p1,p2,p3,x0,x1],'Blank-SAP','SSP-SAP','Naive','stimOn','stimOff')
    end
    set(gca,'box','off')
    xlim([-2,10])
    axis square
end
linkaxes(bx)
% figure
figure;
sgtitle('Ephys contralateral whisker stimlation (MUA)')
for aa = 1:2
    cx(aa) = subplot(1,2,aa);
    hemisphere = hemispheres{1,aa};
    dataType = 'cortMUA';
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
    % Naive
    p3 = plot(data.Naive.(hemisphere).contra.mean_timeVector,data.Naive.(hemisphere).contra.(['mean_' dataType]),'color',colors('sapphire'),'LineWidth',2);
    plot(data.Naive.(hemisphere).contra.mean_timeVector,data.Naive.(hemisphere).contra.(['mean_' dataType]) + data.Naive.(hemisphere).contra.(['stdErr_' dataType]),'color',colors('sapphire'),'LineWidth',0.5)
    plot(data.Naive.(hemisphere).contra.mean_timeVector,data.Naive.(hemisphere).contra.(['mean_' dataType]) - data.Naive.(hemisphere).contra.(['stdErr_' dataType]),'color',colors('sapphire'),'LineWidth',0.5)
    title([hemisphere ' ' dataType])
    label = '\DeltaP/P (%)';
    ylabel(label)
    xlabel('Peri-stimulus time (s)')
    if aa == 1
        legend([p1,p2,p3,x0,x1],'Blank-SAP','SSP-SAP','Naive','stimOn','stimOff')
    end
    set(gca,'box','off')
    xlim([-2,10])
    axis square
end
linkaxes(cx)
% figure
groupNames = {'Naive','Blank_SAP','SSP_SAP','Naive','Blank_SAP','SSP_SAP'};
hemNames = {'LH','LH','LH','RH','RH','RH'};
figure;
sgtitle('Ephys contralateral whisker stimlation (LFP)')
for aa = 1:6
    subplot(2,3,aa);
    hemisphere = hemNames{1,aa};
    groupName = groupNames{1,aa};
    imagesc(data.(groupName).(hemisphere).contra.mean_T,data.(groupName).(hemisphere).contra.mean_F,data.(groupName).(hemisphere).contra.mean_cortS)
    title([groupName ' ' hemisphere])
    ylabel('Freq (Hz)')
    xlabel('Peri-stimulus time (s)')
    colorbar
    caxis([-1,5])
    set(gca,'Ticklength',[0,0])
    axis xy
    set(gca,'box','off')
    axis square
end
