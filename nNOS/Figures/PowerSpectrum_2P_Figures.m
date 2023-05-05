function [] = PowerSpectrum_2P_Figures(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PowerSpec_2P';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
variables = {'S','f'};
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_PowerSpec_2P.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        vIDs = fieldnames(Results_PowerSpec_2P.(group).(animalID));
        for cc = 1:length(vIDs)
            vID = vIDs{cc,1};
            data.(group).dummCheck = 1;
            for dd = 1:length(variables)
                if isfield(data.(group),(variables{1,dd})) == false
                    data.(group).(variables{1,dd}) = [];
                end
            end
            data.(group).S = cat(1,data.(group).S,Results_PowerSpec_2P.(group).(animalID).(vID).S');
            data.(group).f = cat(1,data.(group).f,Results_PowerSpec_2P.(group).(animalID).(vID).f);
        end
    end
end
% take the averages of each field through the proper dimension
for ee = 1:length(groups)
    group = groups{1,ee};
    data.(group).mean_S = mean(data.(group).S,1);
    data.(group).stdErr_S = std(data.(group).S,0,1)./sqrt(size(data.(group).S,1));
    data.(group).mean_f = mean(data.(group).f,1);
end
% figure
summaryFigure = figure;
sgtitle('2P arteriole diameter power spectra')
L1 = loglog(data.Blank_SAP.mean_f,data.Blank_SAP.mean_S,'color',colors('north texas green'),'LineWidth',2);
hold on;
loglog(data.Blank_SAP.mean_f,data.Blank_SAP.mean_S + data.Blank_SAP.stdErr_S,'color',colors('north texas green'),'LineWidth',0.25);
loglog(data.Blank_SAP.mean_f,data.Blank_SAP.mean_S - data.Blank_SAP.stdErr_S,'color',colors('north texas green'),'LineWidth',0.25);
L2 = loglog(data.SSP_SAP.mean_f,data.SSP_SAP.mean_S,'color',colors('electric purple'),'LineWidth',2);
loglog(data.SSP_SAP.mean_f,data.SSP_SAP.mean_S + data.SSP_SAP.stdErr_S,'color',colors('electric purple'),'LineWidth',0.25);
loglog(data.SSP_SAP.mean_f,data.SSP_SAP.mean_S - data.SSP_SAP.stdErr_S,'color',colors('electric purple'),'LineWidth',0.25);
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
xlim([0.003,0.5]);
legend([L1,L2],'Blank-SAP','SSP-SAP')
set(gca,'box','off')
axis square
    % save figure(s)
    if saveFigs == true
        dirpath = [rootFolder delim 'Summary Figures' delim 'Power Spectrum' delim];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(summaryFigure,[dirpath 'PowerSpectrum_ArterioleDiameter_2P']);
    end