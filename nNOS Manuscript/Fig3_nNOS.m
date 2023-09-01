function [] = Fig3_nNOS(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------

%% IOS ephys brief whisker stim
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Evoked_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Naive','Blank_SAP','SSP_SAP'};
solenoids = {'LPadSol','RPadSol','AudSol'};
hemispheres = {'LH','RH'};
dataTypes = {'HbT','cortMUA','cortGam','cortS','T','F','timeVector','group','animalID'};
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
                iosEphysData.(group).(hemisphere).(solenoid).dummCheck = 1;
                for ee = 1:length(dataTypes)
                    dataType = dataTypes{1,ee};
                    if isfield(iosEphysData.(group).(hemisphere).(solenoid),dataType) == false
                        if any(strcmp(dataType,{'group','animalID'})) == true
                            iosEphysData.(group).(hemisphere).(solenoid).(dataType) = {};
                        else
                            iosEphysData.(group).(hemisphere).(solenoid).(dataType) = [];
                        end
                    end
                end
                iosEphysData.(group).(hemisphere).(solenoid).HbT = cat(1,iosEphysData.(group).(hemisphere).(solenoid).HbT,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).HbT);
                iosEphysData.(group).(hemisphere).(solenoid).cortMUA = cat(1,iosEphysData.(group).(hemisphere).(solenoid).cortMUA,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).cortMUA);
                iosEphysData.(group).(hemisphere).(solenoid).cortGam = cat(1,iosEphysData.(group).(hemisphere).(solenoid).cortGam,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).cortGam);
                iosEphysData.(group).(hemisphere).(solenoid).cortS = cat(3,iosEphysData.(group).(hemisphere).(solenoid).cortS,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).cortLFP);
                iosEphysData.(group).(hemisphere).(solenoid).T = cat(1,iosEphysData.(group).(hemisphere).(solenoid).T,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).T);
                iosEphysData.(group).(hemisphere).(solenoid).F = cat(1,iosEphysData.(group).(hemisphere).(solenoid).F,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).F);
                iosEphysData.(group).(hemisphere).(solenoid).timeVector = cat(1,iosEphysData.(group).(hemisphere).(solenoid).timeVector,Results_Evoked_Ephys.(group).(animalID).(hemisphere).Stim.(solenoid).timeVector);
                iosEphysData.(group).(hemisphere).(solenoid).group = cat(1,iosEphysData.(group).(hemisphere).(solenoid).group,group);
                iosEphysData.(group).(hemisphere).(solenoid).animalID = cat(1,iosEphysData.(group).(hemisphere).(solenoid).animalID,animalID);
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
            iosEphysData.(group).(hemisphere).(comparison).group = {};
            iosEphysData.(group).(hemisphere).(comparison).animalID = {};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                if any(strcmp(dataType,{'group','animalID'})) == true
                    iosEphysData.(group).(hemisphere).(comparison).(dataType) = iosEphysData.(group).(hemisphere).(solenoid).(dataType);
                elseif strcmp(dataType,'cortS')
                    iosEphysData.(group).(hemisphere).(comparison).(dataType) = iosEphysData.(group).(hemisphere).(solenoid).(dataType);
                    iosEphysData.(group).(hemisphere).(comparison).(['mean_' dataType]) = mean(iosEphysData.(group).(hemisphere).(solenoid).(dataType),3).*100;
                else
                    iosEphysData.(group).(hemisphere).(comparison).(dataType) = iosEphysData.(group).(hemisphere).(solenoid).(dataType);
                    iosEphysData.(group).(hemisphere).(comparison).(['mean_' dataType]) = mean(iosEphysData.(group).(hemisphere).(solenoid).(dataType),1);
                    iosEphysData.(group).(hemisphere).(comparison).(['stdErr_' dataType]) = std(iosEphysData.(group).(hemisphere).(solenoid).(dataType),1)./sqrt(size(iosEphysData.(group).(hemisphere).(solenoid).(dataType),1));
                    if strcmp(dataType,'HbT') == true
                        for ee = 1:size(iosEphysData.(group).(hemisphere).(comparison).HbT,1)
                            startIdx = find(iosEphysData.(group).(hemisphere).(solenoid).timeVector(ee,:) == 2);
                            endIdx =  find(iosEphysData.(group).(hemisphere).(solenoid).timeVector(ee,:) == 4);
                            timeSnip = iosEphysData.(group).(hemisphere).(solenoid).timeVector(ee,:);
                            hbtSnip = iosEphysData.(group).(hemisphere).(comparison).HbT(ee,:);
                            % iosEphysData.(group).(hemisphere).(comparison).AUC(ee,1) = trapz(timeSnip(startIdx:endIdx),hbtSnip(startIdx:endIdx));
                            iosEphysData.(group).(hemisphere).(comparison).AUC(ee,1) = mean(hbtSnip(startIdx:endIdx));
                        end
                    end
                end
            end
        end
    end
end
% statistics - generalized linear mixed effects model
iosEphysStats.tableSize = cat(1,iosEphysData.Blank_SAP.RH.contra.AUC,iosEphysData.SSP_SAP.RH.contra.AUC);
iosEphysStats.Table = table('Size',[size(iosEphysStats.tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Group','AUC'});
iosEphysStats.Table.Mouse = cat(1,iosEphysData.Blank_SAP.RH.contra.animalID,iosEphysData.SSP_SAP.RH.contra.animalID);
iosEphysStats.Table.Group = cat(1,iosEphysData.Blank_SAP.RH.contra.group,iosEphysData.SSP_SAP.RH.contra.group);
iosEphysStats.Table.AUC = cat(1,iosEphysData.Blank_SAP.RH.contra.AUC,iosEphysData.SSP_SAP.RH.contra.AUC);
iosEphysStats.FitFormula = 'AUC ~ 1 + Group + (1|Mouse)';
iosEphysStats.Stats = fitglme(iosEphysStats.Table,iosEphysStats.FitFormula);

%% running spectroscopy
% data analysis
path = [rootFolder delim 'Results_Zhang'];
cd(path);
% get the demographic info of all animals
experiment.Blank = {'T192';'T205';'T206';'T208';'T209';'T211';'T225'};
runningBlankGroup(1:length(experiment.Blank)) = {'Blank_SAP'};
experiment.SSP = {'T200';'T212';'T213';'T216';'T217';'T218';'T219'};
runningSSPGroup(1:length(experiment.SSP)) = {'SSP_SAP'};
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
        searchfolder = fullfile(path,'Results',animalID);
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
time = -89/30:1/30:5;
blankRunningMean = mean(runningGroup.Blank.HbT,1);
sapRunningMean = mean(runningGroup.SSP.HbT,1);
blankRunningStdErr = std(runningGroup.Blank.HbT,[],1)/sqrt(size(runningGroup.Blank.HbT,1));
sapRunningStdErr = std(runningGroup.SSP.HbT,[],1)/sqrt(size(runningGroup.SSP.HbT,1));
qzGroups = {'Blank','SSP'};
for aa = 1:length(qzGroups)
    qzGroup = qzGroups{1,aa};
    for ee = 1:size(runningGroup.(qzGroup).HbT,1)
        startIdx = find(time == 1.5);
        endIdx =  find(time == 5);
        timeSnip = time;
        hbtSnip = runningGroup.(qzGroup).HbT(ee,:);
        % runningStats.(qzGroup).AUC(ee,1) = trapz(timeSnip(startIdx:endIdx),hbtSnip(startIdx:endIdx));
        runningStats.(qzGroup).AUC(ee,1) = mean(hbtSnip(startIdx:endIdx));
    end
end
% statistics - generalized linear mixed effects model
runningStats.tableSize = cat(1,runningStats.Blank.AUC,runningStats.SSP.AUC);
runningStats.Table = table('Size',[size(runningStats.tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Group','AUC'});
runningStats.Table.Mouse = cat(1,experiment.Blank,experiment.SSP);
runningStats.Table.Group = cat(1,runningBlankGroup',runningSSPGroup');
runningStats.Table.AUC = cat(1,runningStats.Blank.AUC,runningStats.SSP.AUC);
runningStats.FitFormula = 'AUC ~ 1 + Group + (1|Mouse)';
runningStats.Stats = fitglme(runningStats.Table,runningStats.FitFormula);

%% two photon long whisker stim
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Evoked_2P';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'SSP_SAP','Blank_SAP'};
solenoids = {'LPadSol','RPadSol','AudSol'};
dataTypes = {'diameter','baseline','timeVector','count','group','animalID'};
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
                twoPdata.(group).(solenoid).dummCheck = 1;
                for ee = 1:length(dataTypes)
                    dataType = dataTypes{1,ee};
                    if isfield(twoPdata.(group).(solenoid),dataType) == false
                        if any(strcmp(dataType,{'group','animalID'})) == true
                            twoPdata.(group).(solenoid).(dataType) = {};
                        else
                            twoPdata.(group).(solenoid).(dataType) = [];
                        end
                    end
                end
                twoPdata.(group).(solenoid).diameter = cat(1,twoPdata.(group).(solenoid).diameter,Results_Evoked_2P.(group).(animalID).(vID).Stim.(solenoid).diameter*100);
                twoPdata.(group).(solenoid).baseline = cat(1,twoPdata.(group).(solenoid).baseline,Results_Evoked_2P.(group).(animalID).(vID).Stim.(solenoid).baseline);
                twoPdata.(group).(solenoid).count = cat(1,twoPdata.(group).(solenoid).count,Results_Evoked_2P.(group).(animalID).(vID).Stim.(solenoid).count);
                twoPdata.(group).(solenoid).timeVector = cat(1,twoPdata.(group).(solenoid).timeVector,Results_Evoked_2P.(group).(animalID).(vID).Stim.(solenoid).timeVector);
                twoPdata.(group).(solenoid).group = cat(1,twoPdata.(group).(solenoid).group,group);
                twoPdata.(group).(solenoid).animalID = cat(1,twoPdata.(group).(solenoid).animalID,animalID);
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
            twoPdata.(group).(comparison).(dataType) = twoPdata.(group).(solenoid).(dataType);
            if any(strcmp(dataType,{'group','animalID'})) == false
                twoPdata.(group).(comparison).(['mean_'  dataType]) = mean(twoPdata.(group).(solenoid).(dataType),1);
                twoPdata.(group).(comparison).(['stdErr_' dataType]) = std(twoPdata.(group).(solenoid).(dataType),0,1)./sqrt(size(twoPdata.(group).(solenoid).(dataType),1));
            end
        end
    end
end
for aa = 1:length(groups)
    group = groups{1,aa};
    for ee = 1:size(twoPdata.(group).contra.diameter,1)
        startIdx = find(twoPdata.(group).contra.timeVector(ee,:) == 3);
        endIdx =  find(twoPdata.(group).contra.timeVector(ee,:) == 7);
        timeSnip = twoPdata.(group).contra.timeVector(ee,:);
        diameterSnip = twoPdata.(group).contra.diameter(ee,:);
        % twoPdata.(group).contra.AUC(ee,1) = trapz(timeSnip(startIdx:endIdx),diameterSnip(startIdx:endIdx));
        twoPdata.(group).contra.AUC(ee,1) = mean(diameterSnip(startIdx:endIdx));
    end
end
% statistics - generalized linear mixed effects model
diameterStats.tableSize = cat(1,twoPdata.Blank_SAP.contra.AUC,twoPdata.SSP_SAP.contra.AUC);
diameterStats.Table = table('Size',[size(diameterStats.tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Group','AUC'});
diameterStats.Table.Mouse = cat(1,twoPdata.Blank_SAP.contra.animalID,twoPdata.SSP_SAP.contra.animalID);
diameterStats.Table.Group = cat(1,twoPdata.Blank_SAP.contra.group,twoPdata.SSP_SAP.contra.group);
diameterStats.Table.AUC = cat(1,twoPdata.Blank_SAP.contra.AUC,twoPdata.SSP_SAP.contra.AUC);
diameterStats.FitFormula = 'AUC ~ 1 + Group + (1|Mouse)';
diameterStats.Stats = fitglme(diameterStats.Table,diameterStats.FitFormula);

%% IOS long whisker stim
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Evoked_Pulse';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Blank_SAP','SSP_SAP'};
solenoids = {'LPadSol','RPadSol','AudSol'};
dataTypes = {'HbT','timeVector','animalID','group'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Evoked_Pulse.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(solenoids)
            solenoid = solenoids{1,cc};
            pulseData.(group).(solenoid).dummCheck = 1;
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                if isfield(pulseData.(group).(solenoid),dataType) == false
                    if any(strcmp(dataType,{'group','animalID'})) == true
                        pulseData.(group).(solenoid).(dataType) = {};
                    else
                        pulseData.(group).(solenoid).(dataType) = [];
                    end
                end
            end
            pulseData.(group).(solenoid).HbT = cat(1,pulseData.(group).(solenoid).HbT,Results_Evoked_Pulse.(group).(animalID).Stim.(solenoid).HbT);
            pulseData.(group).(solenoid).timeVector = cat(1,pulseData.(group).(solenoid).timeVector,Results_Evoked_Pulse.(group).(animalID).Stim.(solenoid).timeVector);
            pulseData.(group).(solenoid).group = cat(1,pulseData.(group).(solenoid).group,group);
            pulseData.(group).(solenoid).animalID = cat(1,pulseData.(group).(solenoid).animalID,animalID);
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
            pulseData.(group).(comparison).(dataType) = pulseData.(group).(solenoid).(dataType);
            if any(strcmp(dataType,{'group','animalID'})) == false
                pulseData.(group).(comparison).(['mean_'  dataType]) = mean(pulseData.(group).(solenoid).(dataType),1);
                pulseData.(group).(comparison).(['stdErr_' dataType]) = std(pulseData.(group).(solenoid).(dataType),1)./sqrt(size(pulseData.(group).(solenoid).(dataType),1));
            end
        end
    end
end
for aa = 1:length(groups)
    group = groups{1,aa};
    for ee = 1:size(pulseData.(group).contra.HbT,1)
        startIdx = find(pulseData.(group).contra.timeVector(ee,:) == 1.5);
        endIdx =  find(pulseData.(group).contra.timeVector(ee,:) == 6.5);
        timeSnip = pulseData.(group).contra.timeVector(ee,:);
        pulseSnip = pulseData.(group).contra.HbT(ee,:);
        % pulseData.(group).contra.AUC(ee,1) = trapz(timeSnip(startIdx:endIdx),pulseSnip(startIdx:endIdx));
        pulseData.(group).contra.AUC(ee,1) = mean(pulseSnip(startIdx:endIdx));
    end
end
% statistics - generalized linear mixed effects model
pulseStats.tableSize = cat(1,pulseData.Blank_SAP.contra.AUC,pulseData.SSP_SAP.contra.AUC);
pulseStats.Table = table('Size',[size(pulseStats.tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Group','AUC'});
pulseStats.Table.Mouse = cat(1,pulseData.Blank_SAP.contra.animalID,pulseData.SSP_SAP.contra.animalID);
pulseStats.Table.Group = cat(1,pulseData.Blank_SAP.contra.group,pulseData.SSP_SAP.contra.group);
pulseStats.Table.AUC = cat(1,pulseData.Blank_SAP.contra.AUC,pulseData.SSP_SAP.contra.AUC);
pulseStats.FitFormula = 'AUC ~ 1 + Group + (1|Mouse)';
pulseStats.Stats = fitglme(pulseStats.Table,pulseStats.FitFormula);

%% IOS GCaMP
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Evoked_GCaMP';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'Blank_SAP','SSP_SAP'};
solenoids = {'LPadSol','RPadSol','AudSol'};
hemispheres = {'LH','RH'};
dataTypes = {'HbT','HbO','HbR','GCaMP','timeVector','group','animalID'};
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
                gcampdata.(group).(hemisphere).(solenoid).dummCheck = 1;
                for ee = 1:length(dataTypes)
                    dataType = dataTypes{1,ee};
                    if isfield(gcampdata.(group).(hemisphere).(solenoid),dataType) == false
                        if any(strcmp(dataType,{'group','animalID'})) == true
                            gcampdata.(group).(hemisphere).(solenoid).(dataType) = {};
                        else
                            gcampdata.(group).(hemisphere).(solenoid).(dataType) = [];
                        end
                    end
                end
                gcampdata.(group).(hemisphere).(solenoid).HbT = cat(1,gcampdata.(group).(hemisphere).(solenoid).HbT,Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).HbT.mean);
                gcampdata.(group).(hemisphere).(solenoid).HbO = cat(1,gcampdata.(group).(hemisphere).(solenoid).HbO,Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).HbO.mean);
                gcampdata.(group).(hemisphere).(solenoid).HbR = cat(1,gcampdata.(group).(hemisphere).(solenoid).HbR,Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).HbR.mean);
                gcampdata.(group).(hemisphere).(solenoid).GCaMP = cat(1,gcampdata.(group).(hemisphere).(solenoid).GCaMP,Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).GCaMP.mean*100);
                gcampdata.(group).(hemisphere).(solenoid).timeVector = cat(1,gcampdata.(group).(hemisphere).(solenoid).timeVector,Results_Evoked_GCaMP.(group).(animalID).(hemisphere).Stim.(solenoid).HbT.timeVector);
                gcampdata.(group).(hemisphere).(solenoid).group = cat(1,gcampdata.(group).(hemisphere).(solenoid).group,group);
                gcampdata.(group).(hemisphere).(solenoid).animalID = cat(1,gcampdata.(group).(hemisphere).(solenoid).animalID,animalID);
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
            gcampdata.(group).(hemisphere).(comparison).group = {};
            gcampdata.(group).(hemisphere).(comparison).animalID = {};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                if any(strcmp(dataType,{'group','animalID'})) == true
                    gcampdata.(group).(hemisphere).(comparison).(dataType) = gcampdata.(group).(hemisphere).(solenoid).(dataType);
                else
                    gcampdata.(group).(hemisphere).(comparison).(dataType) = gcampdata.(group).(hemisphere).(solenoid).(dataType);
                    gcampdata.(group).(hemisphere).(comparison).(['mean_' dataType]) = mean(gcampdata.(group).(hemisphere).(solenoid).(dataType),1);
                    gcampdata.(group).(hemisphere).(comparison).(['stdErr_' dataType]) = std(gcampdata.(group).(hemisphere).(solenoid).(dataType),1)./sqrt(size(gcampdata.(group).(hemisphere).(solenoid).(dataType),1));
                    if strcmp(dataType,'timeVector') == false
                        for ee = 1:size(gcampdata.(group).(hemisphere).(comparison).(dataType),1)
                            if strcmp(dataType,'GCaMP') == true
                                startIdx = find(gcampdata.(group).(hemisphere).(solenoid).timeVector(ee,:) == 2);
                                endIdx =  find(gcampdata.(group).(hemisphere).(solenoid).timeVector(ee,:) == 5);
                            else
                                startIdx = find(gcampdata.(group).(hemisphere).(solenoid).timeVector(ee,:) == 1.5);
                                endIdx =  find(gcampdata.(group).(hemisphere).(solenoid).timeVector(ee,:) == 6.5);
                            end
                            timeSnip = gcampdata.(group).(hemisphere).(solenoid).timeVector(ee,:);
                            dataSnip = gcampdata.(group).(hemisphere).(comparison).(dataType)(ee,:);
                            % gcampdata.(group).(hemisphere).(comparison).(['AUC_' dataType])(ee,1) = trapz(timeSnip(startIdx:endIdx),dataSnip(startIdx:endIdx));
                            gcampdata.(group).(hemisphere).(comparison).(['AUC_' dataType])(ee,1) = mean(dataSnip(startIdx:endIdx));
                        end
                    end
                end
            end
        end
    end
end
% statistics - generalized linear mixed effects model
specDT = {'HbT','HbO','HbR','GCaMP'};
for aa = 1:length(specDT)
    gcampStats.(specDT{1,aa}).tableSize = cat(1,gcampdata.Blank_SAP.RH.contra.(['AUC_' specDT{1,aa}]),gcampdata.SSP_SAP.RH.contra.(['AUC_' specDT{1,aa}]));
    gcampStats.(specDT{1,aa}).Table = table('Size',[size(gcampStats.(specDT{1,aa}).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Group','AUC'});
    gcampStats.(specDT{1,aa}).Table.Mouse = cat(1,gcampdata.Blank_SAP.RH.contra.animalID,gcampdata.SSP_SAP.RH.contra.animalID);
    gcampStats.(specDT{1,aa}).Table.Group = cat(1,gcampdata.Blank_SAP.RH.contra.group,gcampdata.SSP_SAP.RH.contra.group);
    gcampStats.(specDT{1,aa}).Table.AUC = cat(1,gcampdata.Blank_SAP.RH.contra.(['AUC_' specDT{1,aa}]),gcampdata.SSP_SAP.RH.contra.(['AUC_' specDT{1,aa}]));
    gcampStats.(specDT{1,aa}).FitFormula = 'AUC ~ 1 + Group + (1|Mouse)';
    gcampStats.(specDT{1,aa}).Stats = fitglme(gcampStats.(specDT{1,aa}).Table,gcampStats.(specDT{1,aa}).FitFormula);
end

%% IOS variance signals
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_IntSig_Ephys';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Naive','Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
behaviors = {'Rest','Whisk','Stim','NREM','REM','Iso'};
dataTypes = {'HbT'};
variables = {'avg','p2p','vari'};
fs = 30;
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_IntSig_Ephys.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                iosSigData.(group).(hemisphere).(dataType).dummCheck = 1;
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    if isfield(iosSigData.(group).(hemisphere).(dataType),behavior) == false
                        iosSigData.(group).(hemisphere).(dataType).(behavior).avg = [];
                        iosSigData.(group).(hemisphere).(dataType).(behavior).p2p = [];
                        iosSigData.(group).(hemisphere).(dataType).(behavior).vari = [];
                        iosSigData.(group).(hemisphere).(dataType).(behavior).group = {};
                        iosSigData.(group).(hemisphere).(dataType).(behavior).animalID = {};
                    end
                    animalVar = [];
                    animalP2P = [];
                    for ff = 1:length(Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).indHbT)
                        if strcmp(behavior,'Rest') == true
                            dataArray = Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).indHbT{ff,1}(2*fs:end);
                        elseif strcmp(behavior,'Stim') == true
                            dataArray = Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).indHbT{ff,1} - mean(cell2mat(Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).indHbT),1);
                        else
                            dataArray = Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).indHbT{ff,1};
                        end
                        animalVar(ff,1) = var(dataArray);
                        animalP2P(ff,1) = max(dataArray) - min(dataArray);
                    end
                    iosSigData.(group).(hemisphere).(dataType).(behavior).avg = cat(1,iosSigData.(group).(hemisphere).(dataType).(behavior).avg,mean(Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).HbT));
                    iosSigData.(group).(hemisphere).(dataType).(behavior).p2p = cat(1,iosSigData.(group).(hemisphere).(dataType).(behavior).p2p,mean(animalP2P));
                    iosSigData.(group).(hemisphere).(dataType).(behavior).vari = cat(1,iosSigData.(group).(hemisphere).(dataType).(behavior).vari,mean(animalVar));
                    iosSigData.(group).(hemisphere).(dataType).(behavior).group = cat(1,iosSigData.(group).(hemisphere).(dataType).(behavior).group,group);
                    iosSigData.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,iosSigData.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
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
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                for ee = 1:length(variables)
                    variable = variables{1,ee};
                    iosSigData.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(iosSigData.(group).(hemisphere).(dataType).(behavior).(variable),1);
                    iosSigData.(group).(hemisphere).(dataType).(behavior).(['std_' variable]) = std(iosSigData.(group).(hemisphere).(dataType).(behavior).(variable),1);
                end
            end
        end
    end
end

%% IOS pulse variance
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_IntSig_Pulse';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
variables = {'avg','p2p','vari'};
fs = 30;
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_IntSig_Pulse.(group));
    pulseSigData.(group).dummCheck = 1;
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        if isfield(pulseSigData.(group),'avg') == false
            pulseSigData.(group).avg = [];
            pulseSigData.(group).p2p = [];
            pulseSigData.(group).vari = [];
            pulseSigData.(group).group = {};
            pulseSigData.(group).animalID = {};
        end
        animalVar = [];
        animalP2P = [];
        if isfield(Results_IntSig_Pulse.(group).(animalID),'Rest') == true
            for ff = 1:length(Results_IntSig_Pulse.(group).(animalID).Rest.indHbT)
                dataArray = Results_IntSig_Pulse.(group).(animalID).Rest.indHbT{ff,1}(2*fs:end);
                animalVar(ff,1) = var(dataArray);
                animalP2P(ff,1) = max(dataArray) - min(dataArray);
            end
            pulseSigData.(group).avg = cat(1,pulseSigData.(group).avg,mean(Results_IntSig_Pulse.(group).(animalID).Rest.HbT,'omitnan'));
            pulseSigData.(group).p2p = cat(1,pulseSigData.(group).p2p,mean(animalP2P,'omitnan'));
            pulseSigData.(group).vari = cat(1,pulseSigData.(group).vari,mean(animalVar,'omitnan'));
            pulseSigData.(group).group = cat(1,pulseSigData.(group).group,group);
            pulseSigData.(group).animalID = cat(1,pulseSigData.(group).animalID,animalID);
        end
    end
end
% mean/std
for aa = 1:length(groups)
    group = groups{1,aa};
    for ee = 1:length(variables)
        variable = variables{1,ee};
        pulseSigData.(group).(['mean_' variable]) = mean(pulseSigData.(group).(variable),1);
        pulseSigData.(group).(['std_' variable]) = std(pulseSigData.(group).(variable),1);
    end
end

%% GCaMP variance signals
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_IntSig_GCaMP';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
behaviors = {'Rest','Whisk','Stim','NREM','REM'};
dataTypes = {'HbT','HbO','HbR','GCaMP'};
variables = {'avg','p2p','vari'};
fs = 10;
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_IntSig_GCaMP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                gcampSigdata.(group).(hemisphere).(dataType).dummCheck = 1;
                for ee = 1:length(behaviors)
                    behavior = behaviors{1,ee};
                    if isfield(gcampSigdata.(group).(hemisphere).(dataType),behavior) == false
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).avg = [];
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).p2p = [];
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).vari = [];
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).group = {};
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).animalID = {};
                    end
                    animalVar = [];
                    animalP2P = [];
                    for ff = 1:length(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).indData)
                        if strcmp(behavior,'Rest') == true
                            dataArray = Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).indData{ff,1}(2*fs:end);
                        elseif strcmp(behavior,'Stim') == true
                            dataArray = Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).indData{ff,1} - mean(cell2mat(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).indData),1);
                        else
                            dataArray = Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).indData{ff,1};
                        end
                        animalVar(ff,1) = var(dataArray);
                        animalP2P(ff,1) = max(dataArray) - min(dataArray);
                    end
                    try
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).avg = cat(1,gcampSigdata.(group).(hemisphere).(dataType).(behavior).avg,mean(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).mean));
                    catch
                        gcampSigdata.(group).(hemisphere).(dataType).(behavior).avg = cat(1,gcampSigdata.(group).(hemisphere).(dataType).(behavior).avg,mean(cell2mat(Results_IntSig_GCaMP.(group).(animalID).(hemisphere).(dataType).(behavior).mean)));
                    end
                    gcampSigdata.(group).(hemisphere).(dataType).(behavior).p2p = cat(1,gcampSigdata.(group).(hemisphere).(dataType).(behavior).p2p,mean(animalP2P));
                    gcampSigdata.(group).(hemisphere).(dataType).(behavior).vari = cat(1,gcampSigdata.(group).(hemisphere).(dataType).(behavior).vari,mean(animalVar));
                    gcampSigdata.(group).(hemisphere).(dataType).(behavior).group = cat(1,gcampSigdata.(group).(hemisphere).(dataType).(behavior).group,group);
                    gcampSigdata.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,gcampSigdata.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
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
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                for ee = 1:length(variables)
                    variable = variables{1,ee};
                    gcampSigdata.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(gcampSigdata.(group).(hemisphere).(dataType).(behavior).(variable),1);
                    gcampSigdata.(group).(hemisphere).(dataType).(behavior).(['std_' variable]) = std(gcampSigdata.(group).(hemisphere).(dataType).(behavior).(variable),1);
                end
            end
        end
    end
end

%% LFP power
cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_PowerSpec_LFP';
load(resultsStruct);
cd(rootFolder)
% loop parameters
groups = {'Naive','Blank_SAP','SSP_SAP'};
hemispheres = {'LH','RH'};
behaviors = {'Alert','Asleep','All'};
variables = {'S','f','deltaS'};
dimensions = [2,1,1];
% extract analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_PowerSpec_LFP.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(hemispheres)
            hemisphere = hemispheres{1,cc};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                data.(group).(hemisphere).(behavior).dummCheck = 1;
                for ee = 1:length(variables)
                    variable = variables{1,ee};
                    dimension = dimensions(ee);
                    if isfield(data.(group).(hemisphere).(behavior),(variable)) == false
                        data.(group).(hemisphere).(behavior).(variable) = [];
                        data.(group).(hemisphere).(behavior).group = {};
                        data.(group).(hemisphere).(behavior).animalID = {};
                        data.(group).(hemisphere).(behavior).hemisphere = {};
                        data.(group).(hemisphere).(behavior).behavior = {};
                    end
                    % pull data if field isn't empty
                    if isempty(Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).S) == false
                        if strcmp(variable,'deltaS') == true
                            index = find(round(Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).f,2) == 4);
                            deltaIndex = index(end);
                            data.(group).(hemisphere).(behavior).(variable) = cat(dimension,data.(group).(hemisphere).(behavior).(variable),mean(Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).S(1:deltaIndex)));
                            % for stats
                            data.(group).(hemisphere).(behavior).group = cat(1,data.(group).(hemisphere).(behavior).group,group);
                            data.(group).(hemisphere).(behavior).animalID = cat(1,data.(group).(hemisphere).(behavior).animalID,animalID);
                            data.(group).(hemisphere).(behavior).hemisphere = cat(1,data.(group).(hemisphere).(behavior).hemisphere,hemisphere);
                            data.(group).(hemisphere).(behavior).behavior = cat(1,data.(group).(hemisphere).(behavior).behavior,behavior);
                        else
                            data.(group).(hemisphere).(behavior).(variable) = cat(dimension,data.(group).(hemisphere).(behavior).(variable),Results_PowerSpec_LFP.(group).(animalID).(hemisphere).(behavior).(variable));
                        end
                    end
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
            for dd = 1:length(variables)
                variable = variables{1,dd};
                dimension = dimensions(dd);
                data.(group).(hemisphere).(behavior).(['mean_' variable]) = mean(data.(group).(hemisphere).(behavior).(variable),dimension);
                data.(group).(hemisphere).(behavior).(['stdErr_' variable]) = std(data.(group).(hemisphere).(behavior).(variable),0,dimension)./sqrt(size(data.(group).(hemisphere).(behavior).(variable),dimension));
            end
        end
    end
end
% statistics - generalized linear mixed effects model
for aa = 1:length(hemispheres)
    hemisphere = hemispheres{1,aa};
    for bb = 1:length(behaviors)
        behavior = behaviors{1,bb};
        lfpStats.(hemisphere).(behavior).tableSize = cat(1,data.Blank_SAP.(hemisphere).(behavior).deltaS,data.SSP_SAP.(hemisphere).(behavior).deltaS,data.Naive.(hemisphere).(behavior).deltaS);
        lfpStats.(hemisphere).(behavior).Table = table('Size',[size(lfpStats.(hemisphere).(behavior).tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'group','animalID','behavior','deltaS'});
        lfpStats.(hemisphere).(behavior).Table.group = cat(1,data.Blank_SAP.(hemisphere).(behavior).group,data.SSP_SAP.(hemisphere).(behavior).group,data.Naive.(hemisphere).(behavior).group);
        lfpStats.(hemisphere).(behavior).Table.animalID = cat(1,data.Blank_SAP.(hemisphere).(behavior).animalID,data.SSP_SAP.(hemisphere).(behavior).animalID,data.Naive.(hemisphere).(behavior).animalID);
        lfpStats.(hemisphere).(behavior).Table.behavior = cat(1,data.Blank_SAP.(hemisphere).(behavior).behavior,data.SSP_SAP.(hemisphere).(behavior).behavior,data.Naive.(hemisphere).(behavior).behavior);
        lfpStats.(hemisphere).(behavior).Table.deltaS = cat(1,data.Blank_SAP.(hemisphere).(behavior).deltaS,data.SSP_SAP.(hemisphere).(behavior).deltaS,data.Naive.(hemisphere).(behavior).deltaS);
        lfpStats.(hemisphere).(behavior).FitFormula = 'deltaS ~ 1 + group + behavior + (1|animalID)';
        lfpStats.(hemisphere).(behavior).Stats = fitglme(lfpStats.(hemisphere).(behavior).Table,lfpStats.(hemisphere).(behavior).FitFormula);
    end
end

%% HbT variance statistics
blankVarData = cat(1,iosSigData.Blank_SAP.RH.HbT.Rest.vari,gcampSigdata.Blank_SAP.RH.HbT.Rest.vari,pulseSigData.Blank_SAP.vari);
sspVarData = cat(1,iosSigData.SSP_SAP.RH.HbT.Rest.vari,gcampSigdata.SSP_SAP.RH.HbT.Rest.vari,pulseSigData.SSP_SAP.vari);
restVarStats.tableSize = cat(1,iosSigData.Blank_SAP.RH.HbT.Rest.vari,iosSigData.SSP_SAP.RH.HbT.Rest.vari,gcampSigdata.Blank_SAP.RH.HbT.Rest.vari,gcampSigdata.SSP_SAP.RH.HbT.Rest.vari,pulseSigData.Blank_SAP.vari,pulseSigData.SSP_SAP.vari);
restVarStats.Table = table('Size',[size(restVarStats.tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Group','Variance'});
restVarStats.Table.Mouse = cat(1,iosSigData.Blank_SAP.RH.HbT.Rest.animalID,iosSigData.SSP_SAP.RH.HbT.Rest.animalID,gcampSigdata.Blank_SAP.RH.HbT.Rest.animalID,gcampSigdata.SSP_SAP.RH.HbT.Rest.animalID,pulseSigData.Blank_SAP.animalID,pulseSigData.SSP_SAP.animalID);
restVarStats.Table.Group = cat(1,iosSigData.Blank_SAP.RH.HbT.Rest.group,iosSigData.SSP_SAP.RH.HbT.Rest.group,gcampSigdata.Blank_SAP.RH.HbT.Rest.group,gcampSigdata.SSP_SAP.RH.HbT.Rest.group,pulseSigData.Blank_SAP.group,pulseSigData.SSP_SAP.group);
restVarStats.Table.Variance = cat(1,iosSigData.Blank_SAP.RH.HbT.Rest.vari,iosSigData.SSP_SAP.RH.HbT.Rest.vari,gcampSigdata.Blank_SAP.RH.HbT.Rest.vari,gcampSigdata.SSP_SAP.RH.HbT.Rest.vari,pulseSigData.Blank_SAP.vari,pulseSigData.SSP_SAP.vari);
restVarStats.FitFormula = 'Variance ~ 1 + Group + (1|Mouse)';
restVarStats.Stats = fitglme(restVarStats.Table,restVarStats.FitFormula);

%% figure
Fig3 = figure('Name','Figure 3','units','normalized','outerposition',[0 0 1 1]);

% ephys stimulation
subplot(3,3,1)
p1 = plot(iosEphysData.Blank_SAP.RH.contra.mean_timeVector,iosEphysData.Blank_SAP.RH.contra.mean_HbT,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(iosEphysData.Blank_SAP.RH.contra.mean_timeVector,iosEphysData.Blank_SAP.RH.contra.mean_HbT + iosEphysData.Blank_SAP.RH.contra.stdErr_HbT,'color',colors('north texas green'),'LineWidth',0.25)
plot(iosEphysData.Blank_SAP.RH.contra.mean_timeVector,iosEphysData.Blank_SAP.RH.contra.mean_HbT - iosEphysData.Blank_SAP.RH.contra.stdErr_HbT,'color',colors('north texas green'),'LineWidth',0.25)
p2 = plot(iosEphysData.SSP_SAP.RH.contra.mean_timeVector,iosEphysData.SSP_SAP.RH.contra.mean_HbT,'color',colors('electric purple'),'LineWidth',2);
plot(iosEphysData.SSP_SAP.RH.contra.mean_timeVector,iosEphysData.SSP_SAP.RH.contra.mean_HbT + iosEphysData.SSP_SAP.RH.contra.stdErr_HbT,'color',colors('electric purple'),'LineWidth',0.25)
plot(iosEphysData.SSP_SAP.RH.contra.mean_timeVector,iosEphysData.SSP_SAP.RH.contra.mean_HbT - iosEphysData.SSP_SAP.RH.contra.stdErr_HbT,'color',colors('electric purple'),'LineWidth',0.25)
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2],'Blank-SAP','SSP-SAP')
set(gca,'box','off')
axis square
axis tight
xlim([-2,5]);

subplot(3,3,2)
loglog(data.Blank_SAP.RH.All.mean_f,data.Blank_SAP.RH.All.mean_S,'color',colors('north texas green'),'LineWidth',2);
hold on
loglog(data.Blank_SAP.RH.All.mean_f,data.Blank_SAP.RH.All.mean_S + data.Blank_SAP.RH.All.stdErr_S,'color',colors('north texas green'),'LineWidth',0.25);
loglog(data.Blank_SAP.RH.All.mean_f,data.Blank_SAP.RH.All.mean_S - data.Blank_SAP.RH.All.stdErr_S,'color',colors('north texas green'),'LineWidth',0.25);
loglog(data.SSP_SAP.RH.All.mean_f,data.SSP_SAP.RH.All.mean_S,'color',colors('electric purple'),'LineWidth',2);
loglog(data.SSP_SAP.RH.All.mean_f,data.SSP_SAP.RH.All.mean_S + data.SSP_SAP.RH.All.stdErr_S,'color',colors('electric purple'),'LineWidth',0.25);
loglog(data.SSP_SAP.RH.All.mean_f,data.SSP_SAP.RH.All.mean_S - data.SSP_SAP.RH.All.stdErr_S,'color',colors('electric purple'),'LineWidth',0.25);
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
set(gca,'box','off')
axis square
axis tight
xlim([1,100]);

% ephys rest variance
subplot(3,3,3)
xInds = ones(1,length(blankVarData));
scatter(xInds*1,blankVarData,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e1 = errorbar(1,mean(blankVarData,'omitnan'),std(blankVarData,0,1,'omitnan'),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(sspVarData));
scatter(xInds*2,sspVarData,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(2,mean(sspVarData,'omitnan'),std(sspVarData,0,1,'omitnan'),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
ylabel('\Delta[HbT]^2 (\muM)')
set(gca,'box','off')
set(gca,'xtick',[])
axis square
axis tight
xlim([0,3]);

% GCaMP HbT
subplot(3,3,4);
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbT,'color',colors('north texas green'),'LineWidth',2);
hold on;
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbT + gcampdata.Blank_SAP.RH.contra.stdErr_HbT,'color',colors('north texas green'),'LineWidth',0.25)
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbT - gcampdata.Blank_SAP.RH.contra.stdErr_HbT,'color',colors('north texas green'),'LineWidth',0.25)
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbT,'color',colors('electric purple'),'LineWidth',2);
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbT + gcampdata.SSP_SAP.RH.contra.stdErr_HbT,'color',colors('electric purple'),'LineWidth',0.25)
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbT - gcampdata.SSP_SAP.RH.contra.stdErr_HbT,'color',colors('electric purple'),'LineWidth',0.25)
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
axis square
axis tight
xlim([-2,10]);

subplot(3,3,5)
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_GCaMP,'color',colors('north texas green'),'LineWidth',2);
hold on;
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_GCaMP + gcampdata.Blank_SAP.RH.contra.stdErr_GCaMP,'color',colors('north texas green'),'LineWidth',0.25)
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_GCaMP - gcampdata.Blank_SAP.RH.contra.stdErr_GCaMP,'color',colors('north texas green'),'LineWidth',0.25)
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_GCaMP,'color',colors('electric purple'),'LineWidth',2);
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_GCaMP + gcampdata.SSP_SAP.RH.contra.stdErr_GCaMP,'color',colors('electric purple'),'LineWidth',0.25)
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_GCaMP - gcampdata.SSP_SAP.RH.contra.stdErr_GCaMP,'color',colors('electric purple'),'LineWidth',0.25)
ylabel('\DeltaF/F (%)')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
axis square
axis tight
xlim([-2,10]);

% GCaMP HbO
subplot(3,3,6);
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbO,'color',colors('north texas green'),'LineWidth',2);
hold on;
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbO + gcampdata.Blank_SAP.RH.contra.stdErr_HbO,'color',colors('north texas green'),'LineWidth',0.25)
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbO - gcampdata.Blank_SAP.RH.contra.stdErr_HbO,'color',colors('north texas green'),'LineWidth',0.25)
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbO,'color',colors('electric purple'),'LineWidth',2);
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbO + gcampdata.SSP_SAP.RH.contra.stdErr_HbO,'color',colors('electric purple'),'LineWidth',0.25)
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbO - gcampdata.SSP_SAP.RH.contra.stdErr_HbO,'color',colors('electric purple'),'LineWidth',0.25)
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbR,'color',colors('north texas green'),'LineWidth',2);
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbR + gcampdata.Blank_SAP.RH.contra.stdErr_HbR,'color',colors('north texas green'),'LineWidth',0.25)
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbR - gcampdata.Blank_SAP.RH.contra.stdErr_HbR,'color',colors('north texas green'),'LineWidth',0.25)
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbR,'color',colors('electric purple'),'LineWidth',2);
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbR + gcampdata.SSP_SAP.RH.contra.stdErr_HbR,'color',colors('electric purple'),'LineWidth',0.25)
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbR - gcampdata.SSP_SAP.RH.contra.stdErr_HbR,'color',colors('electric purple'),'LineWidth',0.25)
ylabel('\DeltaHbO-R (\muM)')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
axis square
axis tight
xlim([-2,10]);

%% pulse stimulation
subplot(3,3,7)
plot(pulseData.Blank_SAP.contra.mean_timeVector,pulseData.Blank_SAP.contra.mean_HbT,'color',colors('north texas green'),'LineWidth',2);
hold on;
plot(pulseData.Blank_SAP.contra.mean_timeVector,pulseData.Blank_SAP.contra.mean_HbT + pulseData.Blank_SAP.contra.stdErr_HbT,'color',colors('north texas green'),'LineWidth',0.25)
plot(pulseData.Blank_SAP.contra.mean_timeVector,pulseData.Blank_SAP.contra.mean_HbT - pulseData.Blank_SAP.contra.stdErr_HbT,'color',colors('north texas green'),'LineWidth',0.25)
plot(pulseData.SSP_SAP.contra.mean_timeVector,pulseData.SSP_SAP.contra.mean_HbT,'color',colors('electric purple'),'LineWidth',2);
plot(pulseData.SSP_SAP.contra.mean_timeVector,pulseData.SSP_SAP.contra.mean_HbT + pulseData.SSP_SAP.contra.stdErr_HbT,'color',colors('electric purple'),'LineWidth',0.25)
plot(pulseData.SSP_SAP.contra.mean_timeVector,pulseData.SSP_SAP.contra.mean_HbT - pulseData.SSP_SAP.contra.stdErr_HbT,'color',colors('electric purple'),'LineWidth',0.25)
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
axis square
axis tight
xlim([-2,10]);
% 2p stimulation
subplot(3,3,8)
plot(twoPdata.Blank_SAP.contra.mean_timeVector,twoPdata.Blank_SAP.contra.mean_diameter,'color',colors('north texas green'),'LineWidth',2);
hold on;
plot(twoPdata.Blank_SAP.contra.mean_timeVector,twoPdata.Blank_SAP.contra.mean_diameter + twoPdata.Blank_SAP.contra.stdErr_diameter,'color',colors('north texas green'),'LineWidth',0.25)
plot(twoPdata.Blank_SAP.contra.mean_timeVector,twoPdata.Blank_SAP.contra.mean_diameter - twoPdata.Blank_SAP.contra.stdErr_diameter,'color',colors('north texas green'),'LineWidth',0.25)
plot(twoPdata.SSP_SAP.contra.mean_timeVector,twoPdata.SSP_SAP.contra.mean_diameter,'color',colors('electric purple'),'LineWidth',2);
plot(twoPdata.SSP_SAP.contra.mean_timeVector,twoPdata.SSP_SAP.contra.mean_diameter + twoPdata.SSP_SAP.contra.stdErr_diameter,'color',colors('electric purple'),'LineWidth',0.25)
plot(twoPdata.SSP_SAP.contra.mean_timeVector,twoPdata.SSP_SAP.contra.mean_diameter - twoPdata.SSP_SAP.contra.stdErr_diameter,'color',colors('electric purple'),'LineWidth',0.25)
ylabel('\DeltaD/D (%)')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
axis square
axis tight
xlim([-2,10]);
% running spectroscopy
subplot(3,3,9)
plot(time,blankRunningMean,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(time,blankRunningMean + blankRunningStdErr,'color',colors('north texas green'),'LineWidth',0.25);
plot(time,blankRunningMean - blankRunningStdErr,'color',colors('north texas green'),'LineWidth',0.25);
plot(time,sapRunningMean,'color',colors('electric purple'),'LineWidth',2);
plot(time,sapRunningMean + sapRunningStdErr,'color',colors('electric purple'),'LineWidth',0.25);
plot(time,sapRunningMean - sapRunningStdErr,'color',colors('electric purple'),'LineWidth',0.25);
xlabel('Time (s)')
ylabel('\DeltaHbT (\muM)')
set(gca,'box','off')
axis square
axis tight
xlim([-2,5]);

%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'MATLAB Figure Panels' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(Fig3,[dirpath 'Fig3']);
    set(Fig3,'PaperPositionMode','auto');
    print('-vector','-dpdf','-fillpage',[dirpath 'Fig3'])
    diaryFile = [dirpath 'Fig3_Readout.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    % statistical diary
    diary(diaryFile)
    diary on
   
    % IOS brief stim AUC
    disp('AUC 0.1 s stim, n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(iosEphysData.Blank_SAP.RH.contra.AUC)) ' +/- ' num2str(std(iosEphysData.Blank_SAP.RH.contra.AUC,0,1)./sqrt(size(iosEphysData.Blank_SAP.RH.contra.AUC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(iosEphysData.SSP_SAP.RH.contra.AUC)) ' +/- ' num2str(std(iosEphysData.SSP_SAP.RH.contra.AUC,0,1)./sqrt(size(iosEphysData.SSP_SAP.RH.contra.AUC,1)))]); disp(' ')
    disp('======================================================================================================================')
    disp('GLME statistics for IOS AUC (t = 2:4 sec)')
    disp('======================================================================================================================')
    disp(iosEphysStats.Stats)
   
    % LFP delta power
    disp('LFP delta power, n = 9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(data.Blank_SAP.RH.All.deltaS)) ' +/- ' num2str(std(data.Blank_SAP.RH.All.deltaS,0,1)./sqrt(size(data.Blank_SAP.RH.All.deltaS,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(data.SSP_SAP.RH.All.deltaS)) ' +/- ' num2str(std(data.SSP_SAP.RH.All.deltaS,0,1)./sqrt(size(data.SSP_SAP.RH.All.deltaS,1)))]); disp(' ')
    disp('======================================================================================================================')
    disp('GLME statistics for LFP delta (1:4 Hz)')
    disp('======================================================================================================================')
    disp(lfpStats.RH.All.Stats)
  
    % IOS resting variance
    disp('IOS resting variance, n = 24 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(blankVarData)) ' +/- ' num2str(std(blankVarData,0,1)./sqrt(size(blankVarData,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(sspVarData)) ' +/- ' num2str(std(sspVarData,0,1)./sqrt(size(sspVarData,1)))]); disp(' ')
    disp('======================================================================================================================')
    disp('GLME statistics for resting HbT variance')
    disp('======================================================================================================================')
    disp(restVarStats.Stats)
   
    % IOS long stim HbT
    disp('IOS 5s stim HbT, n = 7-8 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(gcampdata.Blank_SAP.RH.contra.AUC_HbT)) ' +/- ' num2str(std(gcampdata.Blank_SAP.RH.contra.AUC_HbT,0,1)./sqrt(size(gcampdata.Blank_SAP.RH.contra.AUC_HbT,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(gcampdata.SSP_SAP.RH.contra.AUC_HbT)) ' +/- ' num2str(std(gcampdata.SSP_SAP.RH.contra.AUC_HbT,0,1)./sqrt(size(gcampdata.SSP_SAP.RH.contra.AUC_HbT,1)))]); disp(' ')
    disp('======================================================================================================================')
    disp('GLME statistics for GCaMP (HbT) (t = 1.5:6.5 sec)')
    disp('======================================================================================================================')
    disp(gcampStats.HbT.Stats)
   
    % IOS long stim GCaMP
    disp('IOS 5s stim GCaMP, n = 7-8 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(gcampdata.Blank_SAP.RH.contra.AUC_GCaMP)) ' +/- ' num2str(std(gcampdata.Blank_SAP.RH.contra.AUC_GCaMP,0,1)./sqrt(size(gcampdata.Blank_SAP.RH.contra.AUC_GCaMP,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(gcampdata.SSP_SAP.RH.contra.AUC_GCaMP)) ' +/- ' num2str(std(gcampdata.SSP_SAP.RH.contra.AUC_GCaMP,0,1)./sqrt(size(gcampdata.SSP_SAP.RH.contra.AUC_GCaMP,1)))]); disp(' ')
    disp('======================================================================================================================')
    disp('GLME statistics for GCaMP (GCaMP) (t = 2:5 sec)')
    disp('======================================================================================================================')
    disp(gcampStats.GCaMP.Stats)
  
    % IOS long stim HbO
    disp('IOS 5s stim HbO, n = 7-8 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(gcampdata.Blank_SAP.RH.contra.AUC_HbO)) ' +/- ' num2str(std(gcampdata.Blank_SAP.RH.contra.AUC_HbO,0,1)./sqrt(size(gcampdata.Blank_SAP.RH.contra.AUC_HbO,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(gcampdata.SSP_SAP.RH.contra.AUC_HbO)) ' +/- ' num2str(std(gcampdata.SSP_SAP.RH.contra.AUC_HbO,0,1)./sqrt(size(gcampdata.SSP_SAP.RH.contra.AUC_HbO,1)))]); disp(' ')
    disp('======================================================================================================================')
    disp('GLME statistics for GCaMP (HbO) (t = 1.5:6.5 sec)')
    disp('======================================================================================================================')
    disp(gcampStats.HbO.Stats)
    
    % IOS long stim HbR
    disp('IOS 5s stim HbR, n = 7-8 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(gcampdata.Blank_SAP.RH.contra.AUC_HbR)) ' +/- ' num2str(std(gcampdata.Blank_SAP.RH.contra.AUC_HbR,0,1)./sqrt(size(gcampdata.Blank_SAP.RH.contra.AUC_HbR,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(gcampdata.SSP_SAP.RH.contra.AUC_HbR)) ' +/- ' num2str(std(gcampdata.SSP_SAP.RH.contra.AUC_HbR,0,1)./sqrt(size(gcampdata.SSP_SAP.RH.contra.AUC_HbR,1)))]); disp(' ')
    disp('======================================================================================================================')
    disp('GLME statistics for GCaMP (HbR) (t = 1.5:6.5 sec)')
    disp('======================================================================================================================')
    disp(gcampStats.HbR.Stats)
    
    % IOS long stim
    disp('IOS 5s stim AUC, n = 7-8 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(pulseData.Blank_SAP.contra.AUC)) ' +/- ' num2str(std(pulseData.Blank_SAP.contra.AUC,0,1)./sqrt(size(pulseData.Blank_SAP.contra.AUC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(pulseData.SSP_SAP.contra.AUC)) ' +/- ' num2str(std(pulseData.SSP_SAP.contra.AUC,0,1)./sqrt(size(pulseData.SSP_SAP.contra.AUC,1)))]); disp(' ')
    disp('======================================================================================================================')
    disp('GLME statistics for IOS pulse AUC (t = 1.5:6.5 sec)')
    disp('======================================================================================================================')
    disp(pulseStats.Stats)
   
    % IOS running
    disp('Running AUC, n = 7 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(runningStats.Blank.AUC)) ' +/- ' num2str(std(runningStats.Blank.AUC,0,1)./sqrt(size(runningStats.Blank.AUC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(runningStats.SSP.AUC)) ' +/- ' num2str(std(runningStats.SSP.AUC,0,1)./sqrt(size(runningStats.SSP.AUC,1)))]); disp(' ')
    disp('======================================================================================================================')
    disp('GLME statistics for running AUC (t = 1.5:2.5 sec)')
    disp('======================================================================================================================')
    disp(runningStats.Stats)
    
    % 2P long stim
    disp('2P 5s stim AUC, n = 7-9 mice per group, mean +/- SEM'); disp(' ')
    disp(['Blank-SAP ' num2str(mean(twoPdata.Blank_SAP.contra.AUC)) ' +/- ' num2str(std(twoPdata.Blank_SAP.contra.AUC,0,1)./sqrt(size(twoPdata.Blank_SAP.contra.AUC,1)))]); disp(' ')
    disp(['SSP-SAP ' num2str(mean(twoPdata.SSP_SAP.contra.AUC)) ' +/- ' num2str(std(twoPdata.SSP_SAP.contra.AUC,0,1)./sqrt(size(twoPdata.SSP_SAP.contra.AUC,1)))]); disp(' ')
    disp('======================================================================================================================')
    disp('GLME statistics for two photon AUC (t = 3:7 sec)')
    disp('======================================================================================================================')
    disp(diameterStats.Stats)
    
    diary off
end