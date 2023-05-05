function [] = RunningSpectroscopy_Figures(rootFolder,saveFigs,delim)

% data analysis
path = [rootFolder delim 'Results_Zhang'];
cd(path);
% get the demographic info of all animals   
experiment.Blank = {
'T192';
'T205';
'T206';
'T208';
'T209';
'T211';
'T225';
};
experiment.SSP = {
'T200';
'T212';
'T213';
'T215';
'T216';
'T217';
'T218';
'T219';
};
% set groups based on toxin injection
group.Blank.HR = [];
group.Blank.HbT = [];
group.Blank.HbD = [];
group.Blank.RestVar_HR = [];
group.Blank.RestVar_HbT = [];
group.Blank.RestVar_HbD = [];
group.SSP = group.Blank;
% fieldnames
fields = fieldnames(experiment);
for a0 = 1:numel(fields)
    expCondition = fields{a0};
    for a1 = 1:numel(experiment.(expCondition))
        animalID = experiment.(expCondition){a1};
        searchfolder = fullfile(path, 'Results', animalID);
        ind_targetFile = getfilenames(searchfolder, [expCondition,'-SAP', '*.mat']);
        if ~isempty(ind_targetFile)
            out.(animalID) = genFigure_individual_SAP(ind_targetFile);
        end
        group.(expCondition).HR = [group.(expCondition).HR;nanmean(out.(animalID).HR.LTA,1)]; %#ok<*NANMEAN> 
        group.(expCondition).HbT = [group.(expCondition).HbT;nanmean(out.(animalID).HbT.PC.LTA,1)];
        group.(expCondition).HbD = [group.(expCondition).HbD;nanmean(out.(animalID).HbD.PC.LTA,1)];
        
        group.(expCondition).RestVar_HR = [group.(expCondition).RestVar_HR; nanmean(out.(animalID).RestVar.PC.HR)];
        group.(expCondition).RestVar_HbT = [group.(expCondition).RestVar_HbT; nanmean(out.(animalID).RestVar.PC.HbT)];
        group.(expCondition).RestVar_HbD = [group.(expCondition).RestVar_HbD; nanmean(out.(animalID).RestVar.PC.HbD)];
    end
end
% HR
[summaryFigure] = figure_LTA(nanmean(group.Blank.HR,1),nanmean(group.SSP.HR,1),...
           nanstd(group.Blank.HR,[],1)/sqrt(size(group.Blank.HR,1)),...
           nanstd(group.SSP.HR,[],1)/sqrt(size(group.SSP.HR,1))); %#ok<*NANSTD> 
legend({'Blank-SAP','SSP-SAP'});
xlabel('Time (s)')
ylabel('Heart rate (Hz)');
% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures' delim 'Running Spectroscopy' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'RunningSpectroscopy_HeartRate']);
end
% oxygen       
[summaryFigure] = figure_LTA(nanmean(group.Blank.HbD,1),nanmean(group.SSP.HbD,1),...
           nanstd(group.Blank.HbD,[],1)/sqrt(size(group.Blank.HbD,1)),...
           nanstd(group.SSP.HbD,[],1)/sqrt(size(group.SSP.HbD,1)));
legend({'Blank-SAP','SSP-SAP'});
xlabel('Time (s)')
ylabel('Oxygen, PC (uM)')
% save figure(s)
if saveFigs == true
    savefig(summaryFigure,[dirpath 'RunningSpectroscopy_Oxygen']);
end
% HbT      
[summaryFigure] = figure_LTA(nanmean(group.Blank.HbT,1),nanmean(group.SSP.HbT,1),...
           nanstd(group.Blank.HbT,[],1)/sqrt(size(group.Blank.HbT,1)),...
           nanstd(group.SSP.HbT,[],1)/sqrt(size(group.SSP.HbT,1)));
legend({'Blank-SAP','SSP-SAP'});
xlabel('Time (s)')
ylabel('HbT, PC (uM)')
% save figure(s)
if saveFigs == true
    savefig(summaryFigure,[dirpath 'RunningSpectroscopy_HbT']);
end