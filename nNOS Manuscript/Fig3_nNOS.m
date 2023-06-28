function [] = Fig3_nNOS(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
%% IOS ephys
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
        if any(strcmp(animalID,{'1T142','1T172'})) == false
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
                            iosEphysData.(group).(hemisphere).(comparison).AUC(ee,1) = trapz(timeSnip(startIdx:endIdx),hbtSnip(startIdx:endIdx));
                        end
                    end
                end
            end
        end
    end
end
[iosEphysStats.h,iosEphysStats.p] = ttest2(iosEphysData.Blank_SAP.RH.contra.AUC,iosEphysData.SSP_SAP.RH.contra.AUC);
%% running spectroscopy
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
time = -89/30:1/30:5;
blankRunningMean = mean(runningGroup.Blank.HbT,1);
sapRunningMean = mean(runningGroup.SSP.HbT,1);
blankRunningStdErr = std(runningGroup.Blank.HbT,[],1)/sqrt(size(runningGroup.Blank.HbT,1));
sapRunningStdErr = std(runningGroup.SSP.HbT,[],1)/sqrt(size(runningGroup.SSP.HbT,1));
qzGroups = {'Blank','SSP'};
for aa = 1:length(qzGroups)
    qzGroup = qzGroups{1,aa};
    for ee = 1:size(runningGroup.(qzGroup).HbT,1)
        startIdx = find(time == 1);
        endIdx =  find(time == 5);
        timeSnip = time;
        hbtSnip = runningGroup.(qzGroup).HbT(ee,:);
        runningStats.(qzGroup).AUC(ee,1) = trapz(timeSnip(startIdx:endIdx),hbtSnip(startIdx:endIdx));
    end
end
[runningStats.h,runningStats.p] = ttest2(runningStats.Blank.AUC,runningStats.SSP.AUC);
%% two photon
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
        if strcmp(animalID,'T212') == false
            vIDs = fieldnames(Results_Evoked_2P.(group).(animalID));
            for cc = 1:length(vIDs)
                vID = vIDs{cc,1};
                for dd = 1:length(solenoids)
                    solenoid = solenoids{1,dd};
                    twoPdata.(group).(solenoid).dummCheck = 1;
                    for ee = 1:length(dataTypes)
                        dataType = dataTypes{1,ee};
                        if isfield(twoPdata.(group).(solenoid),dataType) == false
                            twoPdata.(group).(solenoid).(dataType) = [];
                        end
                    end
                    twoPdata.(group).(solenoid).diameter = cat(1,twoPdata.(group).(solenoid).diameter,Results_Evoked_2P.(group).(animalID).(vID).Stim.(solenoid).diameter*100);
                    twoPdata.(group).(solenoid).baseline = cat(1,twoPdata.(group).(solenoid).baseline,Results_Evoked_2P.(group).(animalID).(vID).Stim.(solenoid).baseline);
                    twoPdata.(group).(solenoid).count = cat(1,twoPdata.(group).(solenoid).count,Results_Evoked_2P.(group).(animalID).(vID).Stim.(solenoid).count);
                    twoPdata.(group).(solenoid).timeVector = cat(1,twoPdata.(group).(solenoid).timeVector,Results_Evoked_2P.(group).(animalID).(vID).Stim.(solenoid).timeVector);
                end
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
            twoPdata.(group).(comparison).(['mean_'  dataType]) = mean(twoPdata.(group).(solenoid).(dataType),1);
            twoPdata.(group).(comparison).(['stdErr_' dataType]) = std(twoPdata.(group).(solenoid).(dataType),0,1)./sqrt(size(twoPdata.(group).(solenoid).(dataType),1));
        end
    end
end
for aa = 1:length(groups)
    group = groups{1,aa};
    for ee = 1:size(twoPdata.(group).contra.diameter,1)
        startIdx = find(twoPdata.(group).contra.timeVector(ee,:) == 2);
        endIdx =  find(twoPdata.(group).contra.timeVector(ee,:) == 7);
        timeSnip = twoPdata.(group).contra.timeVector(ee,:);
        diameterSnip = twoPdata.(group).contra.diameter(ee,:);
        twoPdata.(group).contra.AUC(ee,1) = trapz(timeSnip(startIdx:endIdx),diameterSnip(startIdx:endIdx));
    end
end
[twoPStats.h,twoPStats.p] = ttest2(twoPdata.Blank_SAP.contra.AUC,twoPdata.SSP_SAP.contra.AUC);
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
                            startIdx = find(gcampdata.(group).(hemisphere).(solenoid).timeVector(ee,:) == 2);
                            endIdx =  find(gcampdata.(group).(hemisphere).(solenoid).timeVector(ee,:) == 7);
                            timeSnip = gcampdata.(group).(hemisphere).(solenoid).timeVector(ee,:);
                            dataSnip = gcampdata.(group).(hemisphere).(comparison).(dataType)(ee,:);
                            gcampdata.(group).(hemisphere).(comparison).(['AUC_' dataType])(ee,1) = trapz(timeSnip(startIdx:endIdx),dataSnip(startIdx:endIdx));
                        end
                    end
                end
            end
        end
    end
end
[gcampStats.HbT.h,gcampStats.HbT.p] = ttest2(gcampdata.Blank_SAP.RH.contra.AUC_HbT,gcampdata.SSP_SAP.RH.contra.AUC_HbT);
[gcampStats.HbO.h,gcampStats.HbO.p] = ttest2(gcampdata.Blank_SAP.RH.contra.AUC_HbO,gcampdata.SSP_SAP.RH.contra.AUC_HbO);
[gcampStats.HbR.h,gcampStats.HbR.p] = ttest2(gcampdata.Blank_SAP.RH.contra.AUC_HbR,gcampdata.SSP_SAP.RH.contra.AUC_HbR);
%% IOS hbt signals
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
        if any(strcmp(animalID,{'T142','T172'})) == false
            for cc = 1:length(hemispheres)
                hemisphere = hemispheres{1,cc};
                for dd = 1:length(dataTypes)
                    dataType = dataTypes{1,dd};
                    iossigdata.(group).(hemisphere).(dataType).dummCheck = 1;
                    for ee = 1:length(behaviors)
                        behavior = behaviors{1,ee};
                        if isfield(iossigdata.(group).(hemisphere).(dataType),behavior) == false
                            iossigdata.(group).(hemisphere).(dataType).(behavior).avg = [];
                            iossigdata.(group).(hemisphere).(dataType).(behavior).p2p = [];
                            iossigdata.(group).(hemisphere).(dataType).(behavior).vari = [];
                            iossigdata.(group).(hemisphere).(dataType).(behavior).group = {};
                            iossigdata.(group).(hemisphere).(dataType).(behavior).animalID = {};
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
                        iossigdata.(group).(hemisphere).(dataType).(behavior).avg = cat(1,iossigdata.(group).(hemisphere).(dataType).(behavior).avg,mean(Results_IntSig_Ephys.(group).(animalID).(hemisphere).(behavior).HbT));
                        iossigdata.(group).(hemisphere).(dataType).(behavior).p2p = cat(1,iossigdata.(group).(hemisphere).(dataType).(behavior).p2p,mean(animalP2P));
                        iossigdata.(group).(hemisphere).(dataType).(behavior).vari = cat(1,iossigdata.(group).(hemisphere).(dataType).(behavior).vari,mean(animalVar));
                        iossigdata.(group).(hemisphere).(dataType).(behavior).group = cat(1,iossigdata.(group).(hemisphere).(dataType).(behavior).group,group);
                        iossigdata.(group).(hemisphere).(dataType).(behavior).animalID = cat(1,iossigdata.(group).(hemisphere).(dataType).(behavior).animalID,animalID);
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
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for dd = 1:length(behaviors)
                behavior = behaviors{1,dd};
                for ee = 1:length(variables)
                    variable = variables{1,ee};
                    iossigdata.(group).(hemisphere).(dataType).(behavior).(['mean_' variable]) = mean(iossigdata.(group).(hemisphere).(dataType).(behavior).(variable),1);
                    iossigdata.(group).(hemisphere).(dataType).(behavior).(['std_' variable]) = std(iossigdata.(group).(hemisphere).(dataType).(behavior).(variable),1);
                end
            end
        end
    end
end
[iosSigStats.Rest.h,iosSigStats.Rest.p] = ttest2(iossigdata.Blank_SAP.RH.HbT.Rest.vari,iossigdata.SSP_SAP.RH.HbT.Rest.vari);
[iosSigStats.Stim.h,iosSigStats.Stim.p] = ttest2(iossigdata.Blank_SAP.RH.HbT.Stim.vari,iossigdata.SSP_SAP.RH.HbT.Stim.vari);
%% GCaMP signals
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
[gcampSigStats.HbT.Rest.h,gcampSigStats.HbT.Rest.p] = ttest2(gcampSigdata.Blank_SAP.RH.HbT.Rest.vari,gcampSigdata.SSP_SAP.RH.HbT.Rest.vari);
[gcampSigStats.HbO.Rest.h,gcampSigStats.HbO.Rest.p] = ttest2(gcampSigdata.Blank_SAP.RH.HbO.Rest.vari,gcampSigdata.SSP_SAP.RH.HbO.Rest.vari);
[gcampSigStats.HbR.Rest.h,gcampSigStats.HbR.Rest.p] = ttest2(gcampSigdata.Blank_SAP.RH.HbR.Rest.vari,gcampSigdata.SSP_SAP.RH.HbR.Rest.vari);
% [gcampSigStats.HbT.Stim.h,gcampSigStats.HbT.Stim.p] = ttest2(gcampSigdata.Blank_SAP.RH.HbT.Stim.vari,gcampSigdata.SSP_SAP.RH.HbT.Stim.vari);
% [gcampSigStats.HbO.Stim.h,gcampSigStats.HbO.Stim.p] = ttest2(gcampSigdata.Blank_SAP.RH.HbO.Stim.vari,gcampSigdata.SSP_SAP.RH.HbO.Stim.vari);
% [gcampSigStats.HbR.Stim.h,gcampSigStats.HbR.Stim.p] = ttest2(gcampSigdata.Blank_SAP.RH.HbR.Stim.vari,gcampSigdata.SSP_SAP.RH.HbR.Stim.vari);
%% figure
Fig3 = figure('Name','Figure 3');
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
xlim([-2,5])
ylim([-10,20])
axis square
subplot(3,3,2)
% running spectroscopy
plot(time,blankRunningMean,'color',colors('north texas green'),'LineWidth',2);
hold on
plot(time,blankRunningMean + blankRunningStdErr,'color',colors('north texas green'),'LineWidth',0.25);
plot(time,blankRunningMean - blankRunningStdErr,'color',colors('north texas green'),'LineWidth',0.25);
plot(time,sapRunningMean,'color',colors('electric purple'),'LineWidth',2);
plot(time,sapRunningMean + sapRunningStdErr,'color',colors('electric purple'),'LineWidth',0.25);
plot(time,sapRunningMean - sapRunningStdErr,'color',colors('electric purple'),'LineWidth',0.25);
xlim([-2,5]);
ylim([-5,25])
xlabel('Time (s)')
ylabel('\DeltaHbT (\muM)')
set(gca,'box','off')
axis square
% 2p stimulation
subplot(3,3,3)
plot(twoPdata.Blank_SAP.contra.mean_timeVector,twoPdata.Blank_SAP.contra.mean_diameter,'color',colors('north texas green'),'LineWidth',2);
hold on;
plot(twoPdata.Blank_SAP.contra.mean_timeVector,twoPdata.Blank_SAP.contra.mean_diameter + twoPdata.Blank_SAP.contra.stdErr_diameter,'color',colors('north texas green'),'LineWidth',0.25)
plot(twoPdata.Blank_SAP.contra.mean_timeVector,twoPdata.Blank_SAP.contra.mean_diameter - twoPdata.Blank_SAP.contra.stdErr_diameter,'color',colors('north texas green'),'LineWidth',0.25)
plot(twoPdata.SSP_SAP.contra.mean_timeVector,twoPdata.SSP_SAP.contra.mean_diameter,'color',colors('electric purple'),'LineWidth',2);
plot(twoPdata.SSP_SAP.contra.mean_timeVector,twoPdata.SSP_SAP.contra.mean_diameter + twoPdata.SSP_SAP.contra.stdErr_diameter,'color',colors('electric purple'),'LineWidth',0.25)
plot(twoPdata.SSP_SAP.contra.mean_timeVector,twoPdata.SSP_SAP.contra.mean_diameter - twoPdata.SSP_SAP.contra.stdErr_diameter,'color',colors('electric purple'),'LineWidth',0.25)
ylabel('\DeltaD/D (%)')
xlabel('Peri-stimulus time (s)')
xlim([-2,10])
ylim([-3,18])
set(gca,'box','off')
axis square
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
xlim([-2,10])
ylim([-5,30])
axis square
% GCaMP HbO
subplot(3,3,5);
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbO,'color',colors('north texas green'),'LineWidth',2);
hold on;
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbO + gcampdata.Blank_SAP.RH.contra.stdErr_HbO,'color',colors('north texas green'),'LineWidth',0.25)
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbO - gcampdata.Blank_SAP.RH.contra.stdErr_HbO,'color',colors('north texas green'),'LineWidth',0.25)
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbO,'color',colors('electric purple'),'LineWidth',2);
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbO + gcampdata.SSP_SAP.RH.contra.stdErr_HbO,'color',colors('electric purple'),'LineWidth',0.25)
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbO - gcampdata.SSP_SAP.RH.contra.stdErr_HbO,'color',colors('electric purple'),'LineWidth',0.25)
ylabel('\DeltaHbO (\muM)')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
xlim([-2,10])
ylim([-5,40])
axis square
% GCaMP HbR
subplot(3,3,6);
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbR,'color',colors('north texas green'),'LineWidth',2);
hold on;
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbR + gcampdata.Blank_SAP.RH.contra.stdErr_HbR,'color',colors('north texas green'),'LineWidth',0.25)
plot(gcampdata.Blank_SAP.RH.contra.mean_timeVector,gcampdata.Blank_SAP.RH.contra.mean_HbR - gcampdata.Blank_SAP.RH.contra.stdErr_HbR,'color',colors('north texas green'),'LineWidth',0.25)
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbR,'color',colors('electric purple'),'LineWidth',2);
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbR + gcampdata.SSP_SAP.RH.contra.stdErr_HbR,'color',colors('electric purple'),'LineWidth',0.25)
plot(gcampdata.SSP_SAP.RH.contra.mean_timeVector,gcampdata.SSP_SAP.RH.contra.mean_HbR - gcampdata.SSP_SAP.RH.contra.stdErr_HbR,'color',colors('electric purple'),'LineWidth',0.25)
ylabel('\DeltaHbR (\muM)')
xlabel('Peri-stimulus time (s)')
set(gca,'box','off')
xlim([-2,10])
ylim([-10,2])
axis square
% ephys rest variance
subplot(3,5,11)
xInds = ones(1,length(iossigdata.Blank_SAP.RH.HbT.Rest.vari));
scatter(xInds*1,iossigdata.Blank_SAP.RH.HbT.Rest.vari,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e1 = errorbar(1,iossigdata.Blank_SAP.RH.HbT.Rest.mean_vari,iossigdata.Blank_SAP.RH.HbT.Rest.std_vari,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(iossigdata.SSP_SAP.RH.HbT.Rest.vari));
scatter(xInds*2,iossigdata.SSP_SAP.RH.HbT.Rest.vari,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(2,iossigdata.SSP_SAP.RH.HbT.Rest.mean_vari,iossigdata.SSP_SAP.RH.HbT.Rest.std_vari,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
ylabel('\Delta[HbT]^2 (\muM)')
xlim([0,3])
set(gca,'box','off')
set(gca,'xtick',[])
axis square
% ephys stim residuals variance
subplot(3,5,12)
xInds = ones(1,length(iossigdata.Blank_SAP.RH.HbT.Stim.vari));
scatter(xInds*1,iossigdata.Blank_SAP.RH.HbT.Stim.vari,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e1 = errorbar(1,iossigdata.Blank_SAP.RH.HbT.Stim.mean_vari,iossigdata.Blank_SAP.RH.HbT.Stim.std_vari,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(iossigdata.SSP_SAP.RH.HbT.Stim.vari));
scatter(xInds*2,iossigdata.SSP_SAP.RH.HbT.Stim.vari,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(2,iossigdata.SSP_SAP.RH.HbT.Stim.mean_vari,iossigdata.SSP_SAP.RH.HbT.Stim.std_vari,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
ylabel('\Delta[HbT]^2 (\muM)')
xlim([0,3])
set(gca,'box','off')
set(gca,'xtick',[])
axis square
% gcamp rest variance (HbT)
subplot(3,5,13)
xInds = ones(1,length(gcampSigdata.Blank_SAP.RH.HbT.Rest.vari));
scatter(xInds*1,gcampSigdata.Blank_SAP.RH.HbT.Rest.vari,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e1 = errorbar(1,gcampSigdata.Blank_SAP.RH.HbT.Rest.mean_vari,gcampSigdata.Blank_SAP.RH.HbT.Rest.std_vari,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(gcampSigdata.SSP_SAP.RH.HbT.Rest.vari));
scatter(xInds*2,gcampSigdata.SSP_SAP.RH.HbT.Rest.vari,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(2,gcampSigdata.SSP_SAP.RH.HbT.Rest.mean_vari,gcampSigdata.SSP_SAP.RH.HbT.Rest.std_vari,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
ylabel('\Delta[HbT]^2 (\muM)')
xlim([0,3])
set(gca,'box','off')
set(gca,'xtick',[])
axis square
%  GCaMP HbO variance residuals
subplot(3,5,14)
xInds = ones(1,length(gcampSigdata.Blank_SAP.RH.HbO.Rest.vari));
scatter(xInds*1,gcampSigdata.Blank_SAP.RH.HbO.Rest.vari,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e1 = errorbar(1,gcampSigdata.Blank_SAP.RH.HbO.Rest.mean_vari,gcampSigdata.Blank_SAP.RH.HbO.Rest.std_vari,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(gcampSigdata.SSP_SAP.RH.HbO.Rest.vari));
scatter(xInds*2,gcampSigdata.SSP_SAP.RH.HbO.Rest.vari,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(2,gcampSigdata.SSP_SAP.RH.HbO.Rest.mean_vari,gcampSigdata.SSP_SAP.RH.HbO.Rest.std_vari,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
ylabel('\Delta[HbO]^2 (\muM)')
xlim([0,3])
set(gca,'box','off')
set(gca,'xtick',[])
axis square
% GCaMP HbR variance residuals
subplot(3,5,15)
xInds = ones(1,length(gcampSigdata.Blank_SAP.RH.HbR.Rest.vari));
scatter(xInds*1,gcampSigdata.Blank_SAP.RH.HbR.Rest.vari,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e1 = errorbar(1,gcampSigdata.Blank_SAP.RH.HbR.Rest.mean_vari,gcampSigdata.Blank_SAP.RH.HbR.Rest.std_vari,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(gcampSigdata.SSP_SAP.RH.HbR.Rest.vari));
scatter(xInds*2,gcampSigdata.SSP_SAP.RH.HbR.Rest.vari,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(2,gcampSigdata.SSP_SAP.RH.HbR.Rest.mean_vari,gcampSigdata.SSP_SAP.RH.HbR.Rest.std_vari,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
ylabel('\Delta[HbR]^2 (\muM)')
xlim([0,3])
set(gca,'box','off')
set(gca,'xtick',[])
axis square
% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Figure Panels' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(Fig3,[dirpath 'Fig3']);
end