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
runningGroup.Blank.HR = [];
runningGroup.Blank.HbT = [];
runningGroup.Blank.HbD = [];
runningGroup.Blank.RestVar_HR = [];
runningGroup.Blank.RestVar_HbT = [];
runningGroup.Blank.RestVar_HbD = [];
runningGroup.SSP = runningGroup.Blank;
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
        runningGroup.(expCondition).HR = [runningGroup.(expCondition).HR;nanmean(out.(animalID).HR.LTA,1)]; %#ok<*NANMEAN> 
        runningGroup.(expCondition).HbT = [runningGroup.(expCondition).HbT;nanmean(out.(animalID).HbT.PC.LTA,1)];
        runningGroup.(expCondition).HbD = [runningGroup.(expCondition).HbD;nanmean(out.(animalID).HbD.PC.LTA,1)];
        
        runningGroup.(expCondition).RestVar_HR = [runningGroup.(expCondition).RestVar_HR; nanmean(out.(animalID).RestVar.PC.HR)];
        runningGroup.(expCondition).RestVar_HbT = [runningGroup.(expCondition).RestVar_HbT; nanmean(out.(animalID).RestVar.PC.HbT)];
        runningGroup.(expCondition).RestVar_HbD = [runningGroup.(expCondition).RestVar_HbD; nanmean(out.(animalID).RestVar.PC.HbD)];
    end
end
% HR
[summaryFigure] = figure_LTA(nanmean(runningGroup.Blank.HR,1),nanmean(runningGroup.SSP.HR,1),...
           nanstd(runningGroup.Blank.HR,[],1)/sqrt(size(runningGroup.Blank.HR,1)),...
           nanstd(runningGroup.SSP.HR,[],1)/sqrt(size(runningGroup.SSP.HR,1))); %#ok<*NANSTD> 
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
[summaryFigure] = figure_LTA(nanmean(runningGroup.Blank.HbD,1),nanmean(runningGroup.SSP.HbD,1),...
           nanstd(runningGroup.Blank.HbD,[],1)/sqrt(size(runningGroup.Blank.HbD,1)),...
           nanstd(runningGroup.SSP.HbD,[],1)/sqrt(size(runningGroup.SSP.HbD,1)));
legend({'Blank-SAP','SSP-SAP'});
xlabel('Time (s)')
ylabel('Oxygen, PC (uM)')
% save figure(s)
if saveFigs == true
    savefig(summaryFigure,[dirpath 'RunningSpectroscopy_Oxygen']);
end
% HbT      
[summaryFigure] = figure_LTA(nanmean(runningGroup.Blank.HbT,1),nanmean(runningGroup.SSP.HbT,1),...
           nanstd(runningGroup.Blank.HbT,[],1)/sqrt(size(runningGroup.Blank.HbT,1)),...
           nanstd(runningGroup.SSP.HbT,[],1)/sqrt(size(runningGroup.SSP.HbT,1)));
legend({'Blank-SAP','SSP-SAP'});
xlabel('Time (s)')
ylabel('HbT, PC (uM)')
% save figure(s)
if saveFigs == true
    savefig(summaryFigure,[dirpath 'RunningSpectroscopy_HbT']);
end