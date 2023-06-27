function [] = PowerSpectrum_Ephys_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PowerSpec_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Naive','Blank_SAP','SSP_SAP'};
dataTypes = {'HbT','gammaBandPower','deltaBandPower'};
hemispheres = {'LH','RH'};
behaviors = {'Rest','NREM','REM','Alert','Asleep','All'};
variables = {'S','normS','f'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_PowerSpec_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                data.(group).(hemisphere).(dataType).dummCheck = 1;
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    if isfield(data.(group).(hemisphere).(dataType),behavior) == false
                        data.(group).(hemisphere).(dataType).(behavior).S = [];
                        data.(group).(hemisphere).(dataType).(behavior).f = [];
                        data.(group).(hemisphere).(dataType).(behavior).group = {};
                        data.(group).(hemisphere).(dataType).(behavior).animalID = {};
                    end
                    if isempty(Results_PowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).S) == false
                        data.(group).(hemisphere).(dataType).(behavior).S = cat(1,data.(group).(hemisphere).(dataType).(behavior).S,Results_PowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).S');
                        data.(group).(hemisphere).(dataType).(behavior).f = cat(1,data.(group).(hemisphere).(dataType).(behavior).f,Results_PowerSpec_Ephys.(group).(animalID).(hemisphere).(dataType).(behavior).f);
                        data.(group).(hemisphere).(dataType).(behavior).group = cat(1,data.(group).(hemisphere).(dataType).(behavior).group,group);
                        data.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,data.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
                    end
                end
            end
        end
    end
end
% find the peak of the resting PSD
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:size(data.(group).(hemisphere).(dataType).Rest.S,2)
                data.(group).(hemisphere).(dataType).baseline(dd,1) = max(data.(group).(hemisphere).(dataType).Rest.S(:,dd));
            end
        end
    end
end
% DC-shift each animal/hemisphere/behavior PSD with respect to the resting peak
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                for ee = 1:size(data.(group).(hemisphere).(dataType).(behavior).S,1)
                    data.(group).(hemisphere).(dataType).(behavior).normS(ee,:) = (data.(group).(hemisphere).(dataType).(behavior).S(ee,:))*(1/(data.(group).(hemisphere).(dataType).baseline(ee,1)));
                end
            end
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(hemispheres)
        hemisphere = hemispheres{1,bb};
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                for ee = 1:length(variables)
                    variable = variables{1,ee};
                    data.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(data.(group).(hemisphere).(dataType).(behavior).(variable),1);
                    data.(group).(hemisphere).(dataType).(behavior).(['stdErr_' variable]) = std(data.(group).(hemisphere).(dataType).(behavior).(variable),0,1)./sqrt(size(data.(group).(hemisphere).(dataType).(behavior).(variable),1));
                end
            end
        end
    end
end
% power figures
xlimits = ({[1/10,0.5],[1/30,0.5],[1/60,0.5],[0.003,0.5],[0.003,0.5],[0.003,0.5]});
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        summaryFigure = figure;
        sgtitle([hemisphere ' ' dataType ' power spectrum [Ephys]'])
        for cc = 1:length(behaviors)
            behavior = behaviors{1,cc};
            subplot(2,3,cc);
            p1 = loglog(data.Naive.(hemisphere).(dataType).(behavior).mean_f,data.Naive.(hemisphere).(dataType).(behavior).mean_normS,'color',colors('sapphire'),'LineWidth',2);
            hold on;
            loglog(data.Naive.(hemisphere).(dataType).(behavior).mean_f,data.Naive.(hemisphere).(dataType).(behavior).mean_normS + data.Naive.(hemisphere).(dataType).(behavior).stdErr_normS,'color',colors('sapphire'),'LineWidth',0.25);
            loglog(data.Naive.(hemisphere).(dataType).(behavior).mean_f,data.Naive.(hemisphere).(dataType).(behavior).mean_normS - data.Naive.(hemisphere).(dataType).(behavior).stdErr_normS,'color',colors('sapphire'),'LineWidth',0.25);
            p2 = loglog(data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_f,data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_normS,'color',colors('north texas green'),'LineWidth',2);
            loglog(data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_f,data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_normS + data.Blank_SAP.(hemisphere).(dataType).(behavior).stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
            loglog(data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_f,data.Blank_SAP.(hemisphere).(dataType).(behavior).mean_normS - data.Blank_SAP.(hemisphere).(dataType).(behavior).stdErr_normS,'color',colors('north texas green'),'LineWidth',0.25);
            p3 = loglog(data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_f,data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_normS,'color',colors('electric purple'),'LineWidth',2);
            loglog(data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_f,data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_normS + data.SSP_SAP.(hemisphere).(dataType).(behavior).stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
            loglog(data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_f,data.SSP_SAP.(hemisphere).(dataType).(behavior).mean_normS - data.SSP_SAP.(hemisphere).(dataType).(behavior).stdErr_normS,'color',colors('electric purple'),'LineWidth',0.25);
            ylabel('Power (a.u.)')
            xlabel('Freq (Hz)')
            title(behavior)
            xlim(xlimits{1,cc})
            if cc == 1
                legend([p1,p2,p3],'Naive','Blank-SAP','SSP-SAP')
            end
            set(gca,'box','off')
            axis square
        end
        % save figure(s)
        if saveFigs == true
            dirpath = [rootFolder delim 'Summary Figures' delim 'Power Spectrum' delim];
            if ~exist(dirpath,'dir')
                mkdir(dirpath);
            end
            savefig(summaryFigure,[dirpath 'PowerSpectrum_Ephys_' hemisphere '_' dataType]);
        end
    end
end