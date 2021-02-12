function [AnalysisResults] = Fig7_Test(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 7 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

% colorBlack = [(0/256),(0/256),(0/256)];
% colorGrey = [(209/256),(211/256),(212/256)];
% colorRfcAwake = [(0/256),(64/256),(64/256)];
% colorRfcNREM = [(0/256),(174/256),(239/256)];
% colorRfcREM = [(190/256),(30/256),(45/256)];
colorRest = [(0/256),(166/256),(81/256)];
colorWhisk = [(31/256),(120/256),(179/256)];
% colorStim = [(255/256),(28/256),(206/256)];
colorNREM = [(191/256),(0/256),(255/256)];
colorREM = [(254/256),(139/256),(0/256)];
colorAlert = [(255/256),(191/256),(0/256)];
colorAsleep = [(0/256),(128/256),(255/256)];
colorAll = [(183/256),(115/256),(51/256)];
% colorIso = [(0/256),(256/256),(256/256)];
%% set-up and process data
IOS_animalIDs = {'T135','T141','T142','T144','T151','T155','T156','T157','T159'};
treatments = {'C57BL6J','SSP_SAP'};
C57BL6J_IDs = {'T141','T155','T156','T157'};
SSP_SAP_IDs = {'T135','T142','T144','T151','T159'};
behavFields = {'Rest','NREM','REM','Awake','Sleep','All'};
behavFields2 = {'Rest','NREM','REM','Awake','All'};
behavFields3 = {'Rest','Whisk','NREM','REM','Awake','Sleep','All'};
dataTypes = {'CBV_HbT','gammaBandPower'};
%% average coherence during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.C57BL6J.Coherr = [];
data.SSP_SAP.Coherr = [];
for a = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,a};
    if ismember(animalID,C57BL6J_IDs) == true
        treatment = 'C57BL6J';
        for b = 1:length(behavFields)
            behavField = behavFields{1,b};
            % create the behavior folder for the first iteration of the loop
            if isfield(data.(treatment).Coherr,behavField) == false
                data.(treatment).Coherr.(behavField) = [];
            end
            for c = 1:length(dataTypes)
                dataType = dataTypes{1,c};
                % don't concatenate empty arrays where there was no data for this behavior
                if isempty(AnalysisResults.(animalID).Coherence.(behavField).(dataType).C) == false
                    % create the data type folder for the first iteration of the loop
                    if isfield(data.(treatment).Coherr.(behavField),dataType) == false
                        data.(treatment).Coherr.(behavField).(dataType).C = [];
                        data.(treatment).Coherr.(behavField).(dataType).f = [];
                        data.(treatment).Coherr.(behavField).(dataType).confC = [];
                    end
                    % concatenate C/f for existing data - exclude any empty sets
                    data.(treatment).Coherr.(behavField).(dataType).C = cat(2,data.(treatment).Coherr.(behavField).(dataType).C,AnalysisResults.(animalID).Coherence.(behavField).(dataType).C);
                    data.(treatment).Coherr.(behavField).(dataType).f = cat(1,data.(treatment).Coherr.(behavField).(dataType).f,AnalysisResults.(animalID).Coherence.(behavField).(dataType).f);
                    data.(treatment).Coherr.(behavField).(dataType).confC = cat(1,data.(treatment).Coherr.(behavField).(dataType).confC,AnalysisResults.(animalID).Coherence.(behavField).(dataType).confC);
                end
            end
        end
    elseif ismember(animalID,SSP_SAP_IDs) == true
        treatment = 'SSP_SAP';
        for b = 1:length(behavFields)
            behavField = behavFields{1,b};
            % create the behavior folder for the first iteration of the loop
            if isfield(data.(treatment).Coherr,behavField) == false
                data.(treatment).Coherr.(behavField) = [];
            end
            for c = 1:length(dataTypes)
                dataType = dataTypes{1,c};
                % don't concatenate empty arrays where there was no data for this behavior
                if isempty(AnalysisResults.(animalID).Coherence.(behavField).(dataType).C) == false
                    % create the data type folder for the first iteration of the loop
                    if isfield(data.(treatment).Coherr.(behavField),dataType) == false
                        data.(treatment).Coherr.(behavField).(dataType).C = [];
                        data.(treatment).Coherr.(behavField).(dataType).f = [];
                        data.(treatment).Coherr.(behavField).(dataType).confC = [];
                    end
                    % concatenate C/f for existing data - exclude any empty sets
                    data.(treatment).Coherr.(behavField).(dataType).C = cat(2,data.(treatment).Coherr.(behavField).(dataType).C,AnalysisResults.(animalID).Coherence.(behavField).(dataType).C);
                    data.(treatment).Coherr.(behavField).(dataType).f = cat(1,data.(treatment).Coherr.(behavField).(dataType).f,AnalysisResults.(animalID).Coherence.(behavField).(dataType).f);
                    data.(treatment).Coherr.(behavField).(dataType).confC = cat(1,data.(treatment).Coherr.(behavField).(dataType).confC,AnalysisResults.(animalID).Coherence.(behavField).(dataType).confC);
                end
            end
        end
    end
end
% take mean/StD of C/f and determine confC line
for qq = 1:length(treatments)
    treatment = treatments{1,qq};
    for ee = 1:length(behavFields)
        behavField = behavFields{1,ee};
        for ff = 1:length(dataTypes)
            dataType = dataTypes{1,ff};
            data.(treatment).Coherr.(behavField).(dataType).meanC = mean(data.(treatment).Coherr.(behavField).(dataType).C,2);
            data.(treatment).Coherr.(behavField).(dataType).stdC = std(data.(treatment).Coherr.(behavField).(dataType).C,0,2);
            data.(treatment).Coherr.(behavField).(dataType).meanf = mean(data.(treatment).Coherr.(behavField).(dataType).f,1);
            data.(treatment).Coherr.(behavField).(dataType).maxConfC = geomean(data.(treatment).Coherr.(behavField).(dataType).confC);
            data.(treatment).Coherr.(behavField).(dataType).maxConfC_Y = ones(length(data.(treatment).Coherr.(behavField).(dataType).meanf),1)*data.(treatment).Coherr.(behavField).(dataType).maxConfC;
        end
    end
end
% %% power spectra during different behaviors
% % cd through each animal's directory and extract the appropriate analysis results
% for a = 1:length(IOS_animalIDs)
%     animalID = IOS_animalIDs{1,a};
%     for b = 1:length(behavFields)
%         behavField = behavFields{1,b};
%         for c = 1:length(dataTypes)
%             dataType = dataTypes{1,c};
%             data.PowerSpec.(behavField).(dataType).adjLH.S{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjLH.S;
%             data.PowerSpec.(behavField).(dataType).adjLH.f{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjLH.f;
%             data.PowerSpec.(behavField).(dataType).adjRH.S{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjRH.S;
%             data.PowerSpec.(behavField).(dataType).adjRH.f{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjRH.f;
%         end
%     end
% end
% % find the peak of the resting PSD for each animal/hemisphere
% for a = 1:length(IOS_animalIDs)
%     for c = 1:length(dataTypes)
%         dataType = dataTypes{1,c};
%         data.PowerSpec.baseline.(dataType).LH{a,1} = max(data.PowerSpec.Rest.(dataType).adjLH.S{a,1});
%         data.PowerSpec.baseline.(dataType).RH{a,1} = max(data.PowerSpec.Rest.(dataType).adjRH.S{a,1});
%     end
% end
% % DC-shift each animal/hemisphere/behavior PSD with respect to the resting peak
% for a = 1:length(IOS_animalIDs)
%     for dd = 1:length(behavFields)
%         behavField = behavFields{1,dd};
%         for j = 1:length(dataTypes)
%             dataType = dataTypes{1,j};
%             for ee = 1:size(data.PowerSpec.(behavField).(dataType).adjLH.S,2)
%                 data.PowerSpec.(behavField).(dataType).normLH{a,1} = (data.PowerSpec.(behavField).(dataType).adjLH.S{a,1})*(1/(data.PowerSpec.baseline.(dataType).LH{a,1}));
%                 data.PowerSpec.(behavField).(dataType).normRH{a,1} = (data.PowerSpec.(behavField).(dataType).adjRH.S{a,1})*(1/(data.PowerSpec.baseline.(dataType).RH{a,1}));
%             end
%         end
%     end
% end
% % concatenate the data from the left and right hemispheres - removes any empty data
% for ee = 1:length(behavFields)
%     behavField = behavFields{1,ee};
%     for ff = 1:length(dataTypes)
%         dataType = dataTypes{1,ff};
%         data.PowerSpec.(behavField).(dataType).cat_S = [];
%         data.PowerSpec.(behavField).(dataType).cat_f = [];
%         for z = 1:length(data.PowerSpec.(behavField).(dataType).normLH)
%             data.PowerSpec.(behavField).(dataType).cat_S = cat(2,data.PowerSpec.(behavField).(dataType).cat_S,data.PowerSpec.(behavField).(dataType).normLH{z,1},data.PowerSpec.(behavField).(dataType).normRH{z,1});
%             data.PowerSpec.(behavField).(dataType).cat_f = cat(1,data.PowerSpec.(behavField).(dataType).cat_f,data.PowerSpec.(behavField).(dataType).adjLH.f{z,1},data.PowerSpec.(behavField).(dataType).adjRH.f{z,1});
%         end
%     end
% end
% % take mean/StD of S/f
% for h = 1:length(behavFields)
%     behavField = behavFields{1,h};
%     for j = 1:length(dataTypes)
%         dataType = dataTypes{1,j};
%         data.PowerSpec.(behavField).(dataType).meanCortS = mean(data.PowerSpec.(behavField).(dataType).cat_S,2);
%         data.PowerSpec.(behavField).(dataType).stdCortS = std(data.PowerSpec.(behavField).(dataType).cat_S,0,2);
%         data.PowerSpec.(behavField).(dataType).meanCortf = mean(data.PowerSpec.(behavField).(dataType).cat_f,1);
%     end
% end
% %% vessel power spectra of different behaviors
% % cd through each animal's directory and extract the appropriate analysis results
% data.VesselPowerSpec = [];
% for aa = 1:length(TwoP_animalIDs)
%     animalID = TwoP_animalIDs{1,aa};
%     for bb = 1:length(behavFields2)
%         behavField = behavFields2{1,bb};
%         if isfield(AnalysisResults.(animalID).PowerSpectra,behavField) == true
%             vesselIDs = fieldnames(AnalysisResults.(animalID).PowerSpectra.(behavField));
%             for cc = 1:length(vesselIDs)
%                 vesselID = vesselIDs{cc,1};
%                 data.VesselPowerSpec.(animalID).(behavField).(vesselID).S = AnalysisResults.(animalID).PowerSpectra.(behavField).(vesselID).S;
%                 data.VesselPowerSpec.(animalID).(behavField).(vesselID).f = AnalysisResults.(animalID).PowerSpectra.(behavField).(vesselID).f;
%             end
%         end
%     end
% end
% % find the peak of the resting PSD for each arteriole
% for aa = 1:length(TwoP_animalIDs)
%     animalID = TwoP_animalIDs{1,aa};
%     vesselIDs = fieldnames(data.VesselPowerSpec.(animalID).Rest);
%     for cc = 1:length(vesselIDs)
%         vesselID = vesselIDs{cc,1};
%         data.VesselPowerSpec.(animalID).baseline.(vesselID) = max(data.VesselPowerSpec.(animalID).Rest.(vesselID).S);
%     end
% end
% % DC-shift each arteriole's PSD with respect to the resting peak
% for aa = 1:length(TwoP_animalIDs)
%     animalID = TwoP_animalIDs{1,aa};
%     for bb = 1:length(behavFields2)
%         behavField = behavFields2{1,bb};
%         if isfield(data.VesselPowerSpec.(animalID),behavField) == true
%             vesselIDs = fieldnames(data.VesselPowerSpec.(animalID).(behavField));
%             for cc = 1:length(vesselIDs)
%                 vesselID = vesselIDs{cc,1};
%                 data.VesselPowerSpec.(animalID).(behavField).(vesselID).normS = data.VesselPowerSpec.(animalID).(behavField).(vesselID).S*(1/data.VesselPowerSpec.(animalID).baseline.(vesselID));
%             end
%         end
%     end
% end
% % concatenate the data from all arterioles - removes any empty data
% for aa = 1:length(TwoP_animalIDs)
%     animalID = TwoP_animalIDs{1,aa};
%     for bb = 1:length(behavFields2)
%         behavField = behavFields2{1,bb};
%         if isfield(data.VesselPowerSpec,behavField) == false
%             data.VesselPowerSpec.(behavField).S = [];
%             data.VesselPowerSpec.(behavField).f = [];
%         end
%         if isfield(data.VesselPowerSpec.(animalID),behavField) == true
%             vesselIDs = fieldnames(data.VesselPowerSpec.(animalID).(behavField));
%             for cc = 1:length(vesselIDs)
%                 vesselID = vesselIDs{cc,1};
%                 data.VesselPowerSpec.(behavField).S = cat(2,data.VesselPowerSpec.(behavField).S,data.VesselPowerSpec.(animalID).(behavField).(vesselID).normS);
%                 data.VesselPowerSpec.(behavField).f = cat(1,data.VesselPowerSpec.(behavField).f,data.VesselPowerSpec.(animalID).(behavField).(vesselID).f);
%             end
%         end
%     end
% end
% % take mean/StD of S/f
% for dd = 1:length(behavFields2)
%     behavField = behavFields2{1,dd};
%     data.VesselPowerSpec.(behavField).meanS = mean(data.VesselPowerSpec.(behavField).S,2);
%     data.VesselPowerSpec.(behavField).StDS = std(data.VesselPowerSpec.(behavField).S,0,2);
%     data.VesselPowerSpec.(behavField).meanf = mean(data.VesselPowerSpec.(behavField).f,1);
% end
%% Pearson's correlations during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.C57BL6J.CorrCoef = [];
data.SSP_SAP.CorrCoef = [];
for a = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,a};
    if ismember(animalID,C57BL6J_IDs) == true
        treatment = 'C57BL6J';
        for b = 1:length(behavFields3)
            behavField = behavFields3{1,b};
            % create the behavior folder for the first iteration of the loop
            if isfield(data.(treatment).CorrCoef,behavField) == false
                data.(treatment).CorrCoef.(behavField) = [];
            end
            for c = 1:length(dataTypes)
                dataType = dataTypes{1,c};
                % don't concatenate empty arrays where there was no data for this behavior
                if isempty(AnalysisResults.(animalID).CorrCoeff.(behavField).(dataType).meanR) == false
                    % create the data type folder for the first iteration of the loop
                    if isfield(data.(treatment).CorrCoef.(behavField),dataType) == false
                        data.(treatment).CorrCoef.(behavField).(dataType).meanRs = [];
                        data.(treatment).CorrCoef.(behavField).(dataType).animalID = {};
                        data.(treatment).CorrCoef.(behavField).(dataType).behavior = {};
                    end
                    % concatenate mean R and the animalID/behavior for statistics table
                    data.(treatment).CorrCoef.(behavField).(dataType).meanRs = cat(1,data.(treatment).CorrCoef.(behavField).(dataType).meanRs,AnalysisResults.(animalID).CorrCoeff.(behavField).(dataType).meanR);
                    data.(treatment).CorrCoef.(behavField).(dataType).animalID = cat(1,data.(treatment).CorrCoef.(behavField).(dataType).animalID,animalID);
                    data.(treatment).CorrCoef.(behavField).(dataType).behavior = cat(1,data.(treatment).CorrCoef.(behavField).(dataType).behavior,behavField);
                end
            end
        end
    elseif ismember(animalID,SSP_SAP_IDs) == true
        treatment = 'SSP_SAP';
        for b = 1:length(behavFields3)
            behavField = behavFields3{1,b};
            % create the behavior folder for the first iteration of the loop
            if isfield(data.(treatment).CorrCoef,behavField) == false
                data.(treatment).CorrCoef.(behavField) = [];
            end
            for c = 1:length(dataTypes)
                dataType = dataTypes{1,c};
                % don't concatenate empty arrays where there was no data for this behavior
                if isempty(AnalysisResults.(animalID).CorrCoeff.(behavField).(dataType).meanR) == false
                    % create the data type folder for the first iteration of the loop
                    if isfield(data.(treatment).CorrCoef.(behavField),dataType) == false
                        data.(treatment).CorrCoef.(behavField).(dataType).meanRs = [];
                        data.(treatment).CorrCoef.(behavField).(dataType).animalID = {};
                        data.(treatment).CorrCoef.(behavField).(dataType).behavior = {};
                    end
                    % concatenate mean R and the animalID/behavior for statistics table
                    data.(treatment).CorrCoef.(behavField).(dataType).meanRs = cat(1,data.(treatment).CorrCoef.(behavField).(dataType).meanRs,AnalysisResults.(animalID).CorrCoeff.(behavField).(dataType).meanR);
                    data.(treatment).CorrCoef.(behavField).(dataType).animalID = cat(1,data.(treatment).CorrCoef.(behavField).(dataType).animalID,animalID);
                    data.(treatment).CorrCoef.(behavField).(dataType).behavior = cat(1,data.(treatment).CorrCoef.(behavField).(dataType).behavior,behavField);
                end
            end
        end
    end
end
% take mean/STD of R
for qq = 1:length(treatments)
    treatment = treatments{1,qq};
    for ee = 1:length(behavFields3)
        behavField = behavFields3{1,ee};
        for ff = 1:length(dataTypes)
            dataType = dataTypes{1,ff};
            data.(treatment).CorrCoef.(behavField).(dataType).meanR = mean(data.(treatment).CorrCoef.(behavField).(dataType).meanRs,1);
            data.(treatment).CorrCoef.(behavField).(dataType).stdR = std(data.(treatment).CorrCoef.(behavField).(dataType).meanRs,0,1);
        end
    end
end
% %% statistics - generalized linear mixed-effects model
% % HbT
% HbTtableSize = cat(1,data.CorrCoef.Rest.CBV_HbT.meanRs,data.CorrCoef.Whisk.CBV_HbT.meanRs,data.CorrCoef.NREM.CBV_HbT.meanRs,data.CorrCoef.REM.CBV_HbT.meanRs,...
%     data.CorrCoef.Awake.CBV_HbT.meanRs,data.CorrCoef.Sleep.CBV_HbT.meanRs,data.CorrCoef.All.CBV_HbT.meanRs);
% HbTTable = table('Size',[size(HbTtableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','CorrCoef','Behavior'});
% HbTTable.Mouse = cat(1,data.CorrCoef.Rest.CBV_HbT.animalID,data.CorrCoef.Whisk.CBV_HbT.animalID,data.CorrCoef.NREM.CBV_HbT.animalID,data.CorrCoef.REM.CBV_HbT.animalID,...
%     data.CorrCoef.Awake.CBV_HbT.animalID,data.CorrCoef.Sleep.CBV_HbT.animalID,data.CorrCoef.All.CBV_HbT.animalID);
% HbTTable.CorrCoef = cat(1,data.CorrCoef.Rest.CBV_HbT.meanRs,data.CorrCoef.Whisk.CBV_HbT.meanRs,data.CorrCoef.NREM.CBV_HbT.meanRs,data.CorrCoef.REM.CBV_HbT.meanRs,...
%     data.CorrCoef.Awake.CBV_HbT.meanRs,data.CorrCoef.Sleep.CBV_HbT.meanRs,data.CorrCoef.All.CBV_HbT.meanRs);
% HbTTable.Behavior = cat(1,data.CorrCoef.Rest.CBV_HbT.behavior,data.CorrCoef.Whisk.CBV_HbT.behavior,data.CorrCoef.NREM.CBV_HbT.behavior,data.CorrCoef.REM.CBV_HbT.behavior,...
%     data.CorrCoef.Awake.CBV_HbT.behavior,data.CorrCoef.Sleep.CBV_HbT.behavior,data.CorrCoef.All.CBV_HbT.behavior);
% HbTFitFormula = 'CorrCoef ~ 1 + Behavior + (1|Mouse)';
% HbTStats = fitglme(HbTTable,HbTFitFormula);
% % gamma-band power
% gammatableSize = cat(1,data.CorrCoef.Rest.gammaBandPower.meanRs,data.CorrCoef.Whisk.gammaBandPower.meanRs,data.CorrCoef.NREM.gammaBandPower.meanRs,data.CorrCoef.REM.gammaBandPower.meanRs,...
%     data.CorrCoef.Awake.gammaBandPower.meanRs,data.CorrCoef.Sleep.gammaBandPower.meanRs,data.CorrCoef.All.gammaBandPower.meanRs);
% gammaTable = table('Size',[size(gammatableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','CorrCoef','Behavior'});
% gammaTable.Mouse = cat(1,data.CorrCoef.Rest.gammaBandPower.animalID,data.CorrCoef.Whisk.gammaBandPower.animalID,data.CorrCoef.NREM.gammaBandPower.animalID,data.CorrCoef.REM.gammaBandPower.animalID,...
%     data.CorrCoef.Awake.gammaBandPower.animalID,data.CorrCoef.Sleep.gammaBandPower.animalID,data.CorrCoef.All.gammaBandPower.animalID);
% gammaTable.CorrCoef = cat(1,data.CorrCoef.Rest.gammaBandPower.meanRs,data.CorrCoef.Whisk.gammaBandPower.meanRs,data.CorrCoef.NREM.gammaBandPower.meanRs,data.CorrCoef.REM.gammaBandPower.meanRs,...
%     data.CorrCoef.Awake.gammaBandPower.meanRs,data.CorrCoef.Sleep.gammaBandPower.meanRs,data.CorrCoef.All.gammaBandPower.meanRs);
% gammaTable.Behavior = cat(1,data.CorrCoef.Rest.gammaBandPower.behavior,data.CorrCoef.Whisk.gammaBandPower.behavior,data.CorrCoef.NREM.gammaBandPower.behavior,data.CorrCoef.REM.gammaBandPower.behavior,...
%     data.CorrCoef.Awake.gammaBandPower.behavior,data.CorrCoef.Sleep.gammaBandPower.behavior,data.CorrCoef.All.gammaBandPower.behavior);
% gammaFitFormula = 'CorrCoef ~ 1 + Behavior + (1|Mouse)';
% gammaStats = fitglme(gammaTable,gammaFitFormula);
%% Fig. 7
summaryFigure = figure('Name','Fig7 (a-g)');
sgtitle('Figure 7 - Turner et al. 2020')
CC_xIndsA = ones(1,length(data.C57BL6J.CorrCoef.Rest.CBV_HbT.meanRs));
CC_xInds = ones(1,length(data.SSP_SAP.CorrCoef.Rest.CBV_HbT.meanRs));
CC_xInds2A = ones(1,length(data.C57BL6J.CorrCoef.Awake.CBV_HbT.animalID));
CC_xInds2 = ones(1,length(data.SSP_SAP.CorrCoef.Awake.CBV_HbT.animalID));
CC_xInds3A = ones(1,length(data.C57BL6J.CorrCoef.Sleep.CBV_HbT.animalID));
CC_xInds3 = ones(1,length(data.SSP_SAP.CorrCoef.Sleep.CBV_HbT.animalID));

%% [7a] power spectra of gamma-band power during different arousal-states
% ax1 = subplot(3,3,1);
% L1 = loglog(data.PowerSpec.Rest.gammaBandPower.meanCortf,data.PowerSpec.Rest.gammaBandPower.meanCortS,'color',colorRest,'LineWidth',2);
% hold on
% rectangle('Position',[0.005,0.15,0.1 - 0.005,95],'FaceColor','w','EdgeColor','w')
% L2 = loglog(data.PowerSpec.NREM.gammaBandPower.meanCortf,data.PowerSpec.NREM.gammaBandPower.meanCortS,'color',colorNREM,'LineWidth',2);
% rectangle('Position',[0.005,0.15,1/30 - 0.005,95],'FaceColor','w','EdgeColor','w')
% L3 = loglog(data.PowerSpec.REM.gammaBandPower.meanCortf,data.PowerSpec.REM.gammaBandPower.meanCortS,'color',colorREM,'LineWidth',2);
% rectangle('Position',[0.005,0.15,1/60 - 0.005,95],'FaceColor','w','EdgeColor','w')
% L4 = loglog(data.PowerSpec.Awake.gammaBandPower.meanCortf,data.PowerSpec.Awake.gammaBandPower.meanCortS,'color',colorAlert,'LineWidth',2);
% L5 = loglog(data.PowerSpec.Sleep.gammaBandPower.meanCortf,data.PowerSpec.Sleep.gammaBandPower.meanCortS,'color',colorAsleep,'LineWidth',2);
% L6 = loglog(data.PowerSpec.All.gammaBandPower.meanCortf,data.PowerSpec.All.gammaBandPower.meanCortS,'color',colorAll,'LineWidth',2);
% xline(1/10,'color','k');
% xline(1/30,'color','k');
% xline(1/60,'color','k');
% title({'[7a] Cortical power','Gamma-band [30-100 Hz]'})
% ylabel('Power (a.u.)')
% xlabel('Freq (Hz)')
% legend([L1,L2,L3,L4,L5,L6],'Rest','NREM','REM','Alert','Asleep','All','Location','SouthEast')
% axis square
% xlim([0.003,0.5])
% ylim([0.1,100])
% set(gca,'box','off')
% ax1.TickLength = [0.03,0.03];
%% [7b] coherence^2 between bilateral gamma-band power during different arousal-states
% ax2 = subplot(3,3,2);
% semilogx(data.Coherr.Rest.gammaBandPower.meanf,data.Coherr.Rest.gammaBandPower.meanC.^2,'color',colorRest,'LineWidth',2);
% hold on
% rectangle('Position',[0.005,0.05,0.1 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
% semilogx(data.Coherr.NREM.gammaBandPower.meanf,data.Coherr.NREM.gammaBandPower.meanC.^2,'color',colorNREM,'LineWidth',2);
% rectangle('Position',[0.005,0.05,1/30 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
% semilogx(data.Coherr.REM.gammaBandPower.meanf,data.Coherr.REM.gammaBandPower.meanC.^2,'color',colorREM,'LineWidth',2);
% rectangle('Position',[0.005,0.05,1/60 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
% semilogx(data.Coherr.Awake.gammaBandPower.meanf,data.Coherr.Awake.gammaBandPower.meanC.^2,'color',colorAlert,'LineWidth',2);
% semilogx(data.Coherr.Sleep.gammaBandPower.meanf,data.Coherr.Sleep.gammaBandPower.meanC.^2,'color',colorAsleep,'LineWidth',2);
% semilogx(data.Coherr.All.gammaBandPower.meanf,data.Coherr.All.gammaBandPower.meanC.^2,'color',colorAll,'LineWidth',2);
% % confidence lines
% % semilogx(data.Coherr.Rest.gammaBandPower.meanf,data.Coherr.Rest.gammaBandPower.maxConfC_Y.^2,'-','color',colorRest,'LineWidth',1);
% % semilogx(data.Coherr.NREM.gammaBandPower.meanf,data.Coherr.NREM.gammaBandPower.maxConfC_Y.^2,'-','color',colorNREM,'LineWidth',1);
% % semilogx(data.Coherr.REM.gammaBandPower.meanf,data.Coherr.REM.gammaBandPower.maxConfC_Y.^2,'-','color',colorREM,'LineWidth',1);
% % semilogx(data.Coherr.Awake.gammaBandPower.meanf,data.Coherr.Awake.gammaBandPower.maxConfC_Y.^2,'-','color',colorAlert,'LineWidth',1);
% % semilogx(data.Coherr.Sleep.gammaBandPower.meanf,data.Coherr.Sleep.gammaBandPower.maxConfC_Y.^2,'-','color',colorAsleep,'LineWidth',1);
% % semilogx(data.Coherr.All.gammaBandPower.meanf,data.Coherr.All.gammaBandPower.maxConfC_Y.^2,'-','color',colorAll,'LineWidth',1);
% xline(1/10,'color','k');
% xline(1/30,'color','k');
% xline(1/60,'color','k');
% ylabel('Coherence^2')
% xlabel('Freq (Hz)')
% title({'[7b] Bilateral coherence^2','Gamma-band [30-100 Hz]'})
% axis square
% xlim([0.003,0.5])
% ylim([0,1])
% set(gca,'box','off')
% ax2.TickLength = [0.03,0.03];
%% [7c] Pearson's correlations between bilateral gamma-band power during different arousal-states
% ax3 = subplot(3,3,3);
% s1 = scatter(CC_xInds*1,data.CorrCoef.Rest.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
% hold on
% e1 = errorbar(1,data.CorrCoef.Rest.gammaBandPower.meanR,data.CorrCoef.Rest.gammaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e1.Color = 'black';
% e1.MarkerSize = 10;
% e1.CapSize = 10;
% s2 = scatter(CC_xInds*2,data.CorrCoef.Whisk.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on', 'jitterAmount',0.25);
% e2 = errorbar(2,data.CorrCoef.Whisk.gammaBandPower.meanR,data.CorrCoef.Whisk.gammaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e2.Color = 'black';
% e2.MarkerSize = 10;
% e2.CapSize = 10;
% s3 = scatter(CC_xInds*3,data.CorrCoef.NREM.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
% e3 = errorbar(3,data.CorrCoef.NREM.gammaBandPower.meanR,data.CorrCoef.NREM.gammaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e3.Color = 'black';
% e3.MarkerSize = 10;
% e3.CapSize = 10;
% s4 = scatter(CC_xInds*4,data.CorrCoef.REM.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
% e4 = errorbar(4,data.CorrCoef.REM.gammaBandPower.meanR,data.CorrCoef.REM.gammaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e4.Color = 'black';
% e4.MarkerSize = 10;
% e4.CapSize = 10;
% s5 = scatter(CC_xInds2*5,data.CorrCoef.Awake.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
% e5 = errorbar(5,data.CorrCoef.Awake.gammaBandPower.meanR,data.CorrCoef.Awake.gammaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e5.Color = 'black';
% e5.MarkerSize = 10;
% e5.CapSize = 10;
% s6 = scatter(CC_xInds3*6,data.CorrCoef.Sleep.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
% e6 = errorbar(6,data.CorrCoef.Sleep.gammaBandPower.meanR,data.CorrCoef.Sleep.gammaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e6.Color = 'black';
% e6.MarkerSize = 10;
% e6.CapSize = 10;
% s7 = scatter(CC_xInds*7,data.CorrCoef.All.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
% e7 = errorbar(7,data.CorrCoef.All.gammaBandPower.meanR,data.CorrCoef.All.gammaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e7.Color = 'black';
% e7.MarkerSize = 10;
% e7.CapSize = 10;
% title({'[7c] Cortical Pearson''s corr. coef','Gamma-band [30-100 Hz]'})
% ylabel('Corr. coefficient')
% legend([s1,s2,s3,s4,s5,s6,s7],'Rest','Whisk','NREM','REM','Alert','Asleep','All')
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
% axis square
% xlim([0,length(behavFields3) + 1])
% ylim([-0.1,1])
% set(gca,'box','off')
% ax3.TickLength = [0.03,0.03];
%% [7d] power spectra of HbT during different arousal-states
% ax4 = subplot(1,2,1);
% loglog(data.PowerSpec.Rest.CBV_HbT.meanCortf,data.PowerSpec.Rest.CBV_HbT.meanCortS,'color',colorRest,'LineWidth',2);
% hold on
% rectangle('Position',[0.005,0.015,0.1 - 0.005,950],'FaceColor','w','EdgeColor','w')
% loglog(data.PowerSpec.NREM.CBV_HbT.meanCortf,data.PowerSpec.NREM.CBV_HbT.meanCortS,'color',colorNREM,'LineWidth',2);
% rectangle('Position',[0.005,0.015,1/30 - 0.005,950],'FaceColor','w','EdgeColor','w')
% loglog(data.PowerSpec.REM.CBV_HbT.meanCortf,data.PowerSpec.REM.CBV_HbT.meanCortS,'color',colorREM,'LineWidth',2);
% rectangle('Position',[0.005,0.015,1/60 - 0.005,950],'FaceColor','w','EdgeColor','w')
% loglog(data.PowerSpec.Awake.CBV_HbT.meanCortf,data.PowerSpec.Awake.CBV_HbT.meanCortS,'color',colorAlert,'LineWidth',2);
% loglog(data.PowerSpec.Sleep.CBV_HbT.meanCortf,data.PowerSpec.Sleep.CBV_HbT.meanCortS,'color',colorAsleep,'LineWidth',2);
% loglog(data.PowerSpec.All.CBV_HbT.meanCortf,data.PowerSpec.All.CBV_HbT.meanCortS,'color',colorAll,'LineWidth',2);
% xline(1/10,'color','k');
% xline(1/30,'color','k');
% xline(1/60,'color','k');
% title({'[7d] Cortical power','\Delta[HbT] (\muM)'})
% ylabel('Power (a.u.)')
% xlabel('Freq (Hz)')
% axis square
% xlim([0.003,0.5]);
% ylim([0.01,1000])
% set(gca,'box','off')
% ax4.TickLength = [0.03,0.03];
%% [7b] coherence^2 between bilateral HbT during different arousal-states
ax5 = subplot(1,2,1);
semilogx(data.C57BL6J.Coherr.Rest.CBV_HbT.meanf,data.C57BL6J.Coherr.Rest.CBV_HbT.meanC.^2,'color',colorRest,'LineWidth',2);
hold on
semilogx(data.SSP_SAP.Coherr.Rest.CBV_HbT.meanf,data.SSP_SAP.Coherr.Rest.CBV_HbT.meanC.^2,'color',colorRest,'LineWidth',1);
rectangle('Position',[0.005,0.05,0.1 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.C57BL6J.Coherr.NREM.CBV_HbT.meanf,data.C57BL6J.Coherr.NREM.CBV_HbT.meanC.^2,'color',colorNREM,'LineWidth',2);
semilogx(data.SSP_SAP.Coherr.NREM.CBV_HbT.meanf,data.SSP_SAP.Coherr.NREM.CBV_HbT.meanC.^2,'color',colorNREM,'LineWidth',1);
rectangle('Position',[0.005,0.05,1/30 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.C57BL6J.Coherr.REM.CBV_HbT.meanf,data.C57BL6J.Coherr.REM.CBV_HbT.meanC.^2,'color',colorREM,'LineWidth',2);
semilogx(data.SSP_SAP.Coherr.REM.CBV_HbT.meanf,data.SSP_SAP.Coherr.REM.CBV_HbT.meanC.^2,'color',colorREM,'LineWidth',1);
rectangle('Position',[0.005,0.05,1/60 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.C57BL6J.Coherr.Awake.CBV_HbT.meanf,data.C57BL6J.Coherr.Awake.CBV_HbT.meanC.^2,'color',colorAlert,'LineWidth',2);
semilogx(data.SSP_SAP.Coherr.Awake.CBV_HbT.meanf,data.SSP_SAP.Coherr.Awake.CBV_HbT.meanC.^2,'color',colorAlert,'LineWidth',1);
semilogx(data.C57BL6J.Coherr.Sleep.CBV_HbT.meanf,data.C57BL6J.Coherr.Sleep.CBV_HbT.meanC.^2,'color',colorAsleep,'LineWidth',2);
semilogx(data.SSP_SAP.Coherr.Sleep.CBV_HbT.meanf,data.SSP_SAP.Coherr.Sleep.CBV_HbT.meanC.^2,'color',colorAsleep,'LineWidth',1);
semilogx(data.C57BL6J.Coherr.All.CBV_HbT.meanf,data.C57BL6J.Coherr.All.CBV_HbT.meanC.^2,'color',colorAll,'LineWidth',2);
semilogx(data.SSP_SAP.Coherr.All.CBV_HbT.meanf,data.SSP_SAP.Coherr.All.CBV_HbT.meanC.^2,'color',colorAll,'LineWidth',1);
% confidence lines
% semilogx(data.Coherr.Rest.CBV_HbT.meanf,data.Coherr.Rest.CBV_HbT.maxConfC_Y.^2,'-','color',colorRest,'LineWidth',1);
% semilogx(data.Coherr.NREM.CBV_HbT.meanf,data.Coherr.NREM.CBV_HbT.maxConfC_Y.^2,'-','color',colorNREM,'LineWidth',1);
% semilogx(data.Coherr.REM.CBV_HbT.meanf,data.Coherr.REM.CBV_HbT.maxConfC_Y.^2,'-','color',colorREM,'LineWidth',1);
% semilogx(data.Coherr.Awake.CBV_HbT.meanf,data.Coherr.Awake.CBV_HbT.maxConfC_Y.^2,'-','color',colorAlert,'LineWidth',1);
% semilogx(data.Coherr.Sleep.CBV_HbT.meanf,data.Coherr.Sleep.CBV_HbT.maxConfC_Y.^2,'-','color',colorAsleep,'LineWidth',1);
% semilogx(data.Coherr.All.CBV_HbT.meanf,data.Coherr.All.CBV_HbT.maxConfC_Y.^2,'-','color',colorAll,'LineWidth',1);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[7e] Bilateral coherence^2','\Delta[HbT] \muM (%)'})
axis square
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [7e] Pearson's correlations between bilateral HbT during different arousal-states
ax6 = subplot(1,2,2);
scatter(CC_xIndsA*1,data.C57BL6J.CorrCoef.Rest.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.C57BL6J.CorrCoef.Rest.CBV_HbT.meanR,data.C57BL6J.CorrCoef.Rest.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(CC_xInds*2,data.SSP_SAP.CorrCoef.Rest.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e2 = errorbar(2,data.SSP_SAP.CorrCoef.Rest.CBV_HbT.meanR,data.SSP_SAP.CorrCoef.Rest.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;

scatter(CC_xIndsA*3,data.C57BL6J.CorrCoef.Whisk.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.C57BL6J.CorrCoef.Whisk.CBV_HbT.meanR,data.C57BL6J.CorrCoef.Whisk.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(CC_xInds*4,data.SSP_SAP.CorrCoef.Whisk.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.SSP_SAP.CorrCoef.Whisk.CBV_HbT.meanR,data.SSP_SAP.CorrCoef.Whisk.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;

scatter(CC_xIndsA*5,data.C57BL6J.CorrCoef.NREM.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.C57BL6J.CorrCoef.NREM.CBV_HbT.meanR,data.C57BL6J.CorrCoef.NREM.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(CC_xInds*6,data.SSP_SAP.CorrCoef.NREM.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.SSP_SAP.CorrCoef.NREM.CBV_HbT.meanR,data.SSP_SAP.CorrCoef.NREM.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;

scatter(CC_xIndsA*7,data.C57BL6J.CorrCoef.REM.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e7 = errorbar(7,data.C57BL6J.CorrCoef.REM.CBV_HbT.meanR,data.C57BL6J.CorrCoef.REM.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
scatter(CC_xInds*8,data.SSP_SAP.CorrCoef.REM.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e8 = errorbar(8,data.SSP_SAP.CorrCoef.REM.CBV_HbT.meanR,data.SSP_SAP.CorrCoef.REM.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
e8.MarkerSize = 10;
e8.CapSize = 10;

scatter(CC_xInds2A*9,data.C57BL6J.CorrCoef.Awake.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e9 = errorbar(9,data.C57BL6J.CorrCoef.Awake.CBV_HbT.meanR,data.C57BL6J.CorrCoef.Awake.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e9.Color = 'black';
e9.MarkerSize = 10;
e9.CapSize = 10;
scatter(CC_xInds2*10,data.SSP_SAP.CorrCoef.Awake.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e10 = errorbar(10,data.SSP_SAP.CorrCoef.Awake.CBV_HbT.meanR,data.SSP_SAP.CorrCoef.Awake.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e10.Color = 'black';
e10.MarkerSize = 10;
e10.CapSize = 10;

scatter(CC_xInds3A*11,data.C57BL6J.CorrCoef.Sleep.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e11 = errorbar(11,data.C57BL6J.CorrCoef.Sleep.CBV_HbT.meanR,data.C57BL6J.CorrCoef.Sleep.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e11.Color = 'black';
e11.MarkerSize = 10;
e11.CapSize = 10;
scatter(CC_xInds3*12,data.SSP_SAP.CorrCoef.Sleep.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e12 = errorbar(12,data.SSP_SAP.CorrCoef.Sleep.CBV_HbT.meanR,data.SSP_SAP.CorrCoef.Sleep.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e12.Color = 'black';
e12.MarkerSize = 10;
e12.CapSize = 10;

scatter(CC_xIndsA*13,data.C57BL6J.CorrCoef.All.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e13 = errorbar(13,data.C57BL6J.CorrCoef.All.CBV_HbT.meanR,data.C57BL6J.CorrCoef.All.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e13.Color = 'black';
e13.MarkerSize = 10;
e13.CapSize = 10;
scatter(CC_xInds*14,data.SSP_SAP.CorrCoef.All.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e14 = errorbar(14,data.SSP_SAP.CorrCoef.All.CBV_HbT.meanR,data.SSP_SAP.CorrCoef.All.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e14.Color = 'black';
e14.MarkerSize = 10;
e14.CapSize = 10;

title({'[7f] Cortical Pearson''s corr. coef','\DeltaHbT \muM (%)'})
ylabel('Corr. coefficient')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,15])
ylim([0,1])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
% %% [7g] arteriole power spectra of HbT during different arousal-states
% ax7 = subplot(3,3,7);
% loglog(data.VesselPowerSpec.Rest.meanf,data.VesselPowerSpec.Rest.meanS,'color',colorRest,'LineWidth',2);
% hold on
% rectangle('Position',[0.005,0.015,0.1 - 0.005,95],'FaceColor','w','EdgeColor','w')
% loglog(data.VesselPowerSpec.NREM.meanf,data.VesselPowerSpec.NREM.meanS,'color',colorNREM,'LineWidth',2);
% rectangle('Position',[0.005,0.015,1/30 - 0.005,95],'FaceColor','w','EdgeColor','w')
% loglog(data.VesselPowerSpec.REM.meanf,data.VesselPowerSpec.REM.meanS,'color',colorREM,'LineWidth',2);
% rectangle('Position',[0.005,0.015,1/60 - 0.005,95],'FaceColor','w','EdgeColor','w')
% loglog(data.VesselPowerSpec.Awake.meanf,data.VesselPowerSpec.Awake.meanS,'color',colorAlert,'LineWidth',2);
% loglog(data.VesselPowerSpec.All.meanf,data.VesselPowerSpec.All.meanS,'color',colorAll,'LineWidth',2);
% xline(1/10,'color','k');
% xline(1/30,'color','k');
% xline(1/60,'color','k');
% title({'[7g] Arteriole power','\DeltaD/D (%)'})
% ylabel('Power (a.u.)')
% xlabel('Freq (Hz)')
% axis square
% xlim([0.003,0.5]);
% ylim([0.01,100])
% set(gca,'box','off')
% ax7.TickLength = [0.03,0.03];
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig7']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig7'])
end
%     %% statistical diary
%     diaryFile = [dirpath 'Fig7_Statistics.txt'];
%     if exist(diaryFile,'file') == 2
%         delete(diaryFile)
%     end
%     diary(diaryFile)
%     diary on
%     % HbT statistical diary
%     disp('======================================================================================================================')
%     disp('[7c] Generalized linear mixed-effects model statistics for mean HbT corr. coef during Rest, Whisk, NREM, and REM')
%     disp('======================================================================================================================')
%     disp(gammaStats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp(['Rest  Gamma P/P R: ' num2str(round(data.CorrCoef.Rest.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Rest.gammaBandPower.stdR,2))]); disp(' ')
%     disp(['Whisk Gamma P/P R: ' num2str(round(data.CorrCoef.Whisk.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Whisk.gammaBandPower.stdR,2))]); disp(' ')
%     disp(['NREM  Gamma P/P R: ' num2str(round(data.CorrCoef.NREM.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.NREM.gammaBandPower.stdR,2))]); disp(' ')
%     disp(['REM   Gamma P/P R: ' num2str(round(data.CorrCoef.REM.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.REM.gammaBandPower.stdR,2))]); disp(' ')
%     disp(['Awake Gamma P/P R: ' num2str(round(data.CorrCoef.Awake.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Awake.gammaBandPower.stdR,2))]); disp(' ')
%     disp(['Sleep Gamma P/P R: ' num2str(round(data.CorrCoef.Sleep.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Sleep.gammaBandPower.stdR,2))]); disp(' ')
%     disp(['All   Gamma P/P R: ' num2str(round(data.CorrCoef.All.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.All.gammaBandPower.stdR,2))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     % gamma statistical diary
%     disp('======================================================================================================================')
%     disp('[7f] Generalized linear mixed-effects model statistics for mean gamma-band corr. coef during Rest, Whisk, NREM, and REM')
%     disp('======================================================================================================================')
%     disp(HbTStats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp(['Rest  [HbT] (uM) R: ' num2str(round(data.CorrCoef.Rest.CBV_HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Rest.CBV_HbT.stdR,2))]); disp(' ')
%     disp(['Whisk [HbT] (uM) R: ' num2str(round(data.CorrCoef.Whisk.CBV_HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Whisk.CBV_HbT.stdR,2))]); disp(' ')
%     disp(['NREM  [HbT] (uM) R: ' num2str(round(data.CorrCoef.NREM.CBV_HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.NREM.CBV_HbT.stdR,2))]); disp(' ')
%     disp(['REM   [HbT] (uM) R: ' num2str(round(data.CorrCoef.REM.CBV_HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.REM.CBV_HbT.stdR,2))]); disp(' ')
%     disp(['Awake [HbT] (uM) R: ' num2str(round(data.CorrCoef.Awake.CBV_HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Awake.CBV_HbT.stdR,2))]); disp(' ')
%     disp(['Sleep [HbT] (uM) R: ' num2str(round(data.CorrCoef.Sleep.CBV_HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Sleep.CBV_HbT.stdR,2))]); disp(' ')
%     disp(['All   [HbT] (uM) R: ' num2str(round(data.CorrCoef.All.CBV_HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.All.CBV_HbT.stdR,2))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     diary off
%     %% organized for supplemental table
%     % variable names
%     ColumnNames_R = {'Rest','Whisk','NREM','REM','Awake','Sleep','All'};
%     % gamma-band R
%     for aa = 1:length(ColumnNames_R)
%         Gamma_R_MeanStD{1,aa} = [num2str(round(data.CorrCoef.(ColumnNames_R{1,aa}).gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.(ColumnNames_R{1,aa}).gammaBandPower.stdR,2))]; %#ok<*AGROW>
%     end
%     % gamma-band R p-values
%     for aa = 1:length(ColumnNames_R)
%         if strcmp(ColumnNames_R{1,aa},'Rest') == true
%             Gamma_R_pVal{1,aa} = {' '};
%         else
%             Gamma_R_pVal{1,aa} = ['p < ' num2str(gammaStats.Coefficients.pValue(aa,1))];
%         end
%     end
%     % HbT R
%     for aa = 1:length(ColumnNames_R)
%         HbT_MeanStD{1,aa} = [num2str(round(data.CorrCoef.(ColumnNames_R{1,aa}).CBV_HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.(ColumnNames_R{1,aa}).CBV_HbT.stdR,2))];
%     end
%     % HbT R p-values
%     for aa = 1:length(ColumnNames_R)
%         if strcmp(ColumnNames_R{1,aa},'Rest') == true
%             HbT_R_pVal{1,aa} = {' '};
%         else
%             HbT_R_pVal{1,aa} = ['p < ' num2str(HbTStats.Coefficients.pValue(aa,1))];
%         end
%     end
%     %% save table data
%     if isfield(AnalysisResults,'CorrCoef') == false
%         AnalysisResults.CorrCoef = [];
%     end
%     if isfield(AnalysisResults.CorrCoef,'gammaBandPower') == false
%         AnalysisResults.CorrCoef.columnNames = ColumnNames_R;
%         AnalysisResults.CorrCoef.gammaBandPower.meanStD = Gamma_R_MeanStD;
%         AnalysisResults.CorrCoef.gammaBandPower.p = Gamma_R_pVal;
%         AnalysisResults.CorrCoef.CBV_HbT.meanStD = HbT_MeanStD;
%         AnalysisResults.CorrCoef.CBV_HbT.p = HbT_R_pVal;
%         cd(rootFolder)
%         save('AnalysisResults.mat','AnalysisResults')
%     end
% end

end
