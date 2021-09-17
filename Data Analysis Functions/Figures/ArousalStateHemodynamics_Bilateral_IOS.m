function [] = ArousalStateHemodynamics_Bilateral_IOS(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: 
%________________________________________________________________________________________________________________________

% colorBlack = [(0/256),(0/256),(0/256)];
% colorGrey = [(209/256),(211/256),(212/256)];
% colorRfcAwake = [(0/256),(64/256),(64/256)];
% colorRfcNREM = [(0/256),(174/256),(239/256)];
% colorRfcREM = [(190/256),(30/256),(45/256)];
colorWhisk = [(31/256),(120/256),(179/256)];
colorStim = [(255/256),(28/256),(206/256)];
colorNREM = [(191/256),(0/256),(255/256)];
colorREM = [(254/256),(139/256),(0/256)];
% colorAlert = [(255/256),(191/256),(0/256)];
% colorAsleep = [(0/256),(128/256),(255/256)];
% colorAll = [(183/256),(115/256),(51/256)];
colorIso = [(0/256),(256/256),(256/256)];
%% set-up and process data
resultsStruct = 'Results_BehavHbT';
load(resultsStruct);
expGroups = {'Naive','SSP-SAP','Blank-SAP'};
setName = 'IOS Set A';
animalIDs.all = {};
for aa = 1:length(expGroups)
    folderList = dir([expGroups{1,aa} delim setName]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs.all = horzcat(animalIDs.all,{folderList.name});
    animalIDs.(strrep(expGroups{1,aa},'-','_')) = {folderList.name};
end
treatments = {'Naive','SSP_SAP','Blank_SAP'};
behavFields = {'Rest','Whisk','Stim','NREM','REM'};
%% mean HbT comparison between behaviors
% pre-allocate the date for each day
for aa = 1:length(animalIDs.all)
    animalID = animalIDs.all{1,aa};
    whiskFileIDs = unique(Results_BehavHbT.(animalID).Whisk.FileIDs);
    whiskFileDates = [];
    % identify the unique days present for each animal using the whisking field.
    for bb = 1:length(whiskFileIDs)
        whiskFileDates{bb,1} = ConvertDate_IOS(whiskFileIDs{bb,1}); %#ok<*AGROW>
    end
    uniqueWhiskFileDates = unique(whiskFileDates);
    % put pre-allocate each date
    for dd = 1:length(behavFields)
        behavField = behavFields{1,dd};
        for ee = 1:length(uniqueWhiskFileDates)
            fileDate = uniqueWhiskFileDates{ee,1};
            data.(animalID).(behavField).(fileDate).MeanLH = [];
            data.(animalID).(behavField).(fileDate).MeanRH = [];
            data.(animalID).(behavField).(fileDate).IndLH = {};
            data.(animalID).(behavField).(fileDate).IndRH = {};
        end
        procData.(behavField).animalID{aa,1} = animalID;
        procData.(behavField).behavior{aa,1} = behavField;
        procData.(behavField).LH{aa,1} = 'LH';
        procData.(behavField).RH{aa,1} = 'RH';
    end
end
% put data into cell for each unique date
for ff = 1:length(animalIDs.all)
    animalID = animalIDs.all{1,ff};
    for gg = 1:length(behavFields)
        behavField = behavFields{1,gg};
        % data is structured slightly differently depending on class
        if strcmp(behavField,'Rest') == true || strcmp(behavField,'Whisk') == true
            fileIDs = Results_BehavHbT.(animalID).(behavField).FileIDs;
            for hh = 1:length(fileIDs)
                fileDate = ConvertDate_IOS(fileIDs{hh,1});
                data.(animalID).(behavField).(fileDate).MeanLH = cat(1,data.(animalID).(behavField).(fileDate).MeanLH,Results_BehavHbT.(animalID).(behavField).MeanAdjLH(hh,1));
                data.(animalID).(behavField).(fileDate).MeanRH = cat(1,data.(animalID).(behavField).(fileDate).MeanRH,Results_BehavHbT.(animalID).(behavField).MeanAdjRH(hh,1));
                data.(animalID).(behavField).(fileDate).IndLH = cat(1,data.(animalID).(behavField).(fileDate).IndLH,Results_BehavHbT.(animalID).(behavField).IndAdjLH{hh,1});
                data.(animalID).(behavField).(fileDate).IndRH = cat(1,data.(animalID).(behavField).(fileDate).IndRH,Results_BehavHbT.(animalID).(behavField).IndAdjRH{hh,1});
            end
        elseif strcmp(behavField,'Stim') == true
            % left hem stims
            fileIDs = Results_BehavHbT.(animalID).(behavField).LH_FileIDs;
            for hh = 1:length(fileIDs)
                fileDate = ConvertDate_IOS(fileIDs{hh,1});
                data.(animalID).(behavField).(fileDate).MeanLH = cat(1,data.(animalID).(behavField).(fileDate).MeanLH,Results_BehavHbT.(animalID).(behavField).MeanAdjLH(hh,1));
                data.(animalID).(behavField).(fileDate).IndLH = cat(1,data.(animalID).(behavField).(fileDate).IndLH,Results_BehavHbT.(animalID).(behavField).IndAdjLH{hh,1});
            end
            % right hem stims
            fileIDs = Results_BehavHbT.(animalID).(behavField).RH_FileIDs;
            for hh = 1:length(fileIDs)
                fileDate = ConvertDate_IOS(fileIDs{hh,1});
                data.(animalID).(behavField).(fileDate).MeanRH = cat(1,data.(animalID).(behavField).(fileDate).MeanRH,Results_BehavHbT.(animalID).(behavField).MeanAdjRH(hh,1));
                data.(animalID).(behavField).(fileDate).IndRH = cat(1,data.(animalID).(behavField).(fileDate).IndRH,Results_BehavHbT.(animalID).(behavField).IndAdjRH{hh,1});
            end
        else
            fileIDs = Results_BehavHbT.(animalID).(behavField).FileIDs;
            for ii = 1:length(fileIDs)
                fileDate = ConvertDate_IOS(fileIDs{ii,1});
                data.(animalID).(behavField).(fileDate).MeanLH = cat(1,data.(animalID).(behavField).(fileDate).MeanLH,Results_BehavHbT.(animalID).(behavField).MeanAdjLH(ii,1));
                data.(animalID).(behavField).(fileDate).MeanRH = cat(1,data.(animalID).(behavField).(fileDate).MeanRH,Results_BehavHbT.(animalID).(behavField).MeanAdjRH(ii,1));
                data.(animalID).(behavField).(fileDate).IndLH = cat(1,data.(animalID).(behavField).(fileDate).IndLH,Results_BehavHbT.(animalID).(behavField).IndAdjLH{ii,1});
                data.(animalID).(behavField).(fileDate).IndRH = cat(1,data.(animalID).(behavField).(fileDate).IndRH,Results_BehavHbT.(animalID).(behavField).IndAdjRH{ii,1});
            end
        end
    end
end
% find the mean of the 10-second resting periods from each day to determine a baseline
for jj = 1:length(animalIDs.all)
    animalID = animalIDs.all{1,jj};
    whiskFileIDs = unique(Results_BehavHbT.(animalID).Whisk.FileIDs);
    whiskFileDates = [];
    for kk = 1:length(whiskFileIDs)
        whiskFileDates{kk,1} = ConvertDate_IOS(whiskFileIDs{kk,1}); %#ok<*AGROW>
    end
    uniqueWhiskFileDates = unique(whiskFileDates);
    % take mean from each day. Days with no data will show up as NaN and be excluded
    for ll = 1:length(uniqueWhiskFileDates)
        fileDate = uniqueWhiskFileDates{ll,1};
        data.(animalID).Rest.(fileDate).baselineLH = mean(data.(animalID).Rest.(fileDate).MeanLH);
        data.(animalID).Rest.(fileDate).baselineRH = mean(data.(animalID).Rest.(fileDate).MeanRH);
    end
end
% subtract the 10-second resting baseline for each day from the other data types. If the day doesn't have resting data,
% exclude it from analysis
for mm = 1:length(animalIDs.all)
    animalID = animalIDs.all{1,mm};
    for nn = 1:length(behavFields)
        behavField = behavFields{1,nn};
        % subtract each day's 10-second baseline from each behavior field
        fileDates = fieldnames(data.(animalID).(behavField));
        for oo = 1:length(fileDates)
            fileDate = fileDates{oo,1};
            if strcmp(behavField,'Stim') == true || strcmp(behavField,'Whisk') == true
                data.(animalID).(behavField).(fileDate).CorrMeanLH = data.(animalID).(behavField).(fileDate).MeanLH;% - data.(animalID).Rest.(fileDate).baselineLH;
                data.(animalID).(behavField).(fileDate).CorrMeanRH = data.(animalID).(behavField).(fileDate).MeanRH;%; - data.(animalID).Rest.(fileDate).baselineRH;
                for pp = 1:length(data.(animalID).(behavField).(fileDate).IndLH)
                    data.(animalID).(behavField).(fileDate).CorrIndLH{pp,1} = data.(animalID).(behavField).(fileDate).IndLH{pp,1};% - data.(animalID).Rest.(fileDate).baselineLH;
                end
                for pp = 1:length(data.(animalID).(behavField).(fileDate).IndRH)
                    data.(animalID).(behavField).(fileDate).CorrIndRH{pp,1} = data.(animalID).(behavField).(fileDate).IndRH{pp,1};% - data.(animalID).Rest.(fileDate).baselineRH;
                end
            else
                data.(animalID).(behavField).(fileDate).CorrMeanLH = data.(animalID).(behavField).(fileDate).MeanLH - data.(animalID).Rest.(fileDate).baselineLH;
                data.(animalID).(behavField).(fileDate).CorrMeanRH = data.(animalID).(behavField).(fileDate).MeanRH - data.(animalID).Rest.(fileDate).baselineRH;
                for pp = 1:length(data.(animalID).(behavField).(fileDate).IndLH)
                    data.(animalID).(behavField).(fileDate).CorrIndLH{pp,1} = data.(animalID).(behavField).(fileDate).IndLH{pp,1} - data.(animalID).Rest.(fileDate).baselineLH;
                end
                for pp = 1:length(data.(animalID).(behavField).(fileDate).IndRH)
                    data.(animalID).(behavField).(fileDate).CorrIndRH{pp,1} = data.(animalID).(behavField).(fileDate).IndRH{pp,1} - data.(animalID).Rest.(fileDate).baselineRH;
                end
            end
        end
    end
end
% take the mean of the corrected data from each unique day
for qq = 1:length(animalIDs.all)
    animalID = animalIDs.all{1,qq};
    for rr = 1:length(behavFields)
        behavField = behavFields{1,rr};
        fileDates = fieldnames(data.(animalID).(behavField));
        for ss = 1:length(fileDates)
            fileDate = fileDates{ss,1};
            data.(animalID).(behavField).(fileDate).DayMeanLH = mean(data.(animalID).(behavField).(fileDate).CorrMeanLH);
            data.(animalID).(behavField).(fileDate).DayMeanRH = mean(data.(animalID).(behavField).(fileDate).CorrMeanRH);
            data.(animalID).(behavField).(fileDate).DayAllMeanLH = [];
            data.(animalID).(behavField).(fileDate).DayAllMeanRH = [];
            data.(animalID).(behavField).(fileDate).DayIndLH = [];
            data.(animalID).(behavField).(fileDate).DayIndRH = [];
            % concatenate individual trials into a single array for each unique day
            if isfield(data.(animalID).(behavField).(fileDate),'CorrIndLH') == true
                % left means - diff loop is necessary as STIM field has diff number of events
                for tt = 1:length(data.(animalID).(behavField).(fileDate).CorrMeanLH)
                    data.(animalID).(behavField).(fileDate).DayAllMeanLH = cat(2,data.(animalID).(behavField).(fileDate).DayIndLH,data.(animalID).(behavField).(fileDate).CorrIndLH{tt,1});
                end
                % right means
                for tt = 1:length(data.(animalID).(behavField).(fileDate).CorrMeanRH)
                    data.(animalID).(behavField).(fileDate).DayAllMeanRH = cat(2,data.(animalID).(behavField).(fileDate).DayIndRH,data.(animalID).(behavField).(fileDate).CorrIndRH{tt,1});
                end
                % left individual data pts
                for tt = 1:length(data.(animalID).(behavField).(fileDate).CorrIndLH)
                    data.(animalID).(behavField).(fileDate).DayIndLH = cat(2,data.(animalID).(behavField).(fileDate).DayIndLH,data.(animalID).(behavField).(fileDate).CorrIndLH{tt,1});
                end
                % right individual data pts
                for tt = 1:length(data.(animalID).(behavField).(fileDate).CorrIndRH)
                    data.(animalID).(behavField).(fileDate).DayIndRH = cat(2,data.(animalID).(behavField).(fileDate).DayIndRH,data.(animalID).(behavField).(fileDate).CorrIndRH{tt,1});
                end
            end
        end
    end
end
% put all the corrected means from each unique day into a single vector
nans = 1;
for uu = 1:length(animalIDs.all)
    animalID = animalIDs.all{1,uu};
    for vv = 1:length(behavFields)
        behavField = behavFields{1,vv};
        fileDates = fieldnames(data.(animalID).(behavField));
        procData.(animalID).(behavField).DayMeansLH = [];
        procData.(animalID).(behavField).DayMeansRH = [];
        procData.(animalID).(behavField).CatIndLH = [];
        procData.(animalID).(behavField).CatIndRH = [];
        for ww = 1:length(fileDates)
            fileDate = fileDates{ww,1};
            if isnan(data.(animalID).(behavField).(fileDate).DayMeanLH) == false
                procData.(animalID).(behavField).DayMeansLH = cat(1,procData.(animalID).(behavField).DayMeansLH,data.(animalID).(behavField).(fileDate).DayMeanLH);
                procData.(animalID).(behavField).DayMeansRH = cat(1,procData.(animalID).(behavField).DayMeansRH,data.(animalID).(behavField).(fileDate).DayMeanRH);
                procData.(animalID).(behavField).CatIndLH = cat(2,procData.(animalID).(behavField).CatIndLH,data.(animalID).(behavField).(fileDate).DayIndLH);
                procData.(animalID).(behavField).CatIndRH = cat(2,procData.(animalID).(behavField).CatIndRH,data.(animalID).(behavField).(fileDate).DayIndRH);
            else
                nans = nans + 1;
            end
        end
    end
end
behavFields2 = {'Rest','Whisk','Stim','NREM','REM','Iso'};
for zz = 1:length(animalIDs.all)
    animalID = animalIDs.all{1,zz};
    % recognize treatment based on animal group
    if ismember(animalIDs.all{1,zz},animalIDs.Naive) == true
        treatment = 'Naive';
    elseif ismember(animalIDs.all{1,zz},animalIDs.SSP_SAP) == true
        treatment = 'SSP_SAP';
    elseif ismember(animalIDs.all{1,zz},animalIDs.Blank_SAP) == true
        treatment = 'Blank_SAP';
    end
    for yy = 1:length(behavFields2)
        behavField = behavFields2{1,yy};
        testData.(treatment).(behavField).dummCheck = 1;
        if isfield(testData.(treatment).(behavField),'meanLH') == false
            testData.(treatment).(behavField).meanLH = [];
            testData.(treatment).(behavField).meanRH = [];
            testData.(treatment).(behavField).animalIDs = [];
            testData.(treatment).(behavField).treatment = [];
        end
        if strcmp(behavField,'Iso') == false
            testData.(treatment).(behavField).meanLH = cat(1,testData.(treatment).(behavField).meanLH,mean(procData.(animalID).(behavField).DayMeansLH));
            testData.(treatment).(behavField).meanRH = cat(1,testData.(treatment).(behavField).meanRH,mean(procData.(animalID).(behavField).DayMeansRH));
            testData.(treatment).(behavField).animalIDs = cat(1,testData.(treatment).(behavField).animalIDs,animalID);
            testData.(treatment).(behavField).treatment = cat(1,testData.(treatment).(behavField).treatment,{treatment});
        else
            testData.(treatment).(behavField).meanLH = cat(1,testData.(treatment).(behavField).meanLH,mean(Results_BehavHbT.(animalID).(behavField).adjLH));
            testData.(treatment).(behavField).meanRH = cat(1,testData.(treatment).(behavField).meanRH,mean(Results_BehavHbT.(animalID).(behavField).adjRH));
            testData.(treatment).(behavField).animalIDs = cat(1,testData.(treatment).(behavField).animalIDs,animalID);
            testData.(treatment).(behavField).treatment = cat(1,testData.(treatment).(behavField).treatment,{treatment});
        end
    end
end
% take the mean and stdev across animals
for qqq = 1:length(treatments)
    treatment = treatments{1,qqq};
    for aaa = 1:length(behavFields2)
        behavField = behavFields2{1,aaa};
        testData.(treatment).(behavField).LH_MeanCBV = nanmean(testData.(treatment).(behavField).meanLH,1);
        testData.(treatment).(behavField).LH_StdMeanCBV = nanstd(testData.(treatment).(behavField).meanLH,0,1);
        testData.(treatment).(behavField).RH_MeanCBV = nanmean(testData.(treatment).(behavField).meanRH,1);
        testData.(treatment).(behavField).RH_StdMeanCBV = nanstd(testData.(treatment).(behavField).meanRH,0,1);
    end
end
%% statistics - generalized linear mixed effects model
behavFields3 = {'Whisk','Stim','NREM','REM','Iso'};
for bb = 1:length(behavFields3)
    behavField = behavFields3{1,bb};
    % statistics - generalized linear mixed effects model
    Stats.(behavField).tableSize = cat(1,testData.Blank_SAP.(behavField).meanRH,testData.SSP_SAP.(behavField).meanRH);
    Stats.(behavField).Table = table('Size',[size(Stats.(behavField).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Treatment','HbT'});
    Stats.(behavField).Table.Mouse = cat(1,testData.Blank_SAP.(behavField).animalIDs,testData.SSP_SAP.(behavField).animalIDs);
    Stats.(behavField).Table.Treatment = cat(1,testData.Blank_SAP.(behavField).treatment,testData.SSP_SAP.(behavField).treatment);
    Stats.(behavField).Table.HbT = cat(1,testData.Blank_SAP.(behavField).meanRH,testData.SSP_SAP.(behavField).meanRH);
    Stats.(behavField).FitFormula = 'HbT ~ 1 + Treatment + (1|Mouse)';
    Stats.(behavField).Stats = fitglme(Stats.(behavField).Table,Stats.(behavField).FitFormula);
end
%% mean HbT during different behaviors
summaryFigure = figure;
C57_xInds = ones(1,length(testData.Naive.Rest.meanLH));
SSP_xInds = ones(1,length(testData.SSP_SAP.Rest.meanLH));
Blank_xInds = ones(1,length(testData.Blank_SAP.Rest.meanLH));
%% Whisk
% Naive
b1 = scatter(1,testData.Naive.Whisk.LH_MeanCBV,'d','MarkerFaceColor',colors('sapphire'),'MarkerEdgeColor','k','LineWidth',1.5);
hold on
scatter(2,testData.Naive.Whisk.RH_MeanCBV,'d','MarkerFaceColor',colors('sapphire'),'MarkerEdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.Naive.Whisk.meanLH)
    x = [C57_xInds(1,aa)*1,C57_xInds(1,aa)*2];
    y = [testData.Naive.Whisk.meanLH(aa,1),testData.Naive.Whisk.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorWhisk,'jitter','off', 'jitterAmount',0.25)
end
% SSP-SAP
b2 = scatter(3,testData.SSP_SAP.Whisk.LH_MeanCBV,'d','MarkerFaceColor',colors('electric purple'),'MarkerEdgeColor','k','LineWidth',1.5);
scatter(4,testData.SSP_SAP.Whisk.RH_MeanCBV,'d','MarkerFaceColor',colors('electric purple'),'MarkerEdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.SSP_SAP.Whisk.meanLH)
    x = [SSP_xInds(1,aa)*3,SSP_xInds(1,aa)*4];
    y = [testData.SSP_SAP.Whisk.meanLH(aa,1),testData.SSP_SAP.Whisk.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorWhisk,'jitter','off', 'jitterAmount',0.25)
end
% Blank-SAP
b3 = scatter(5,testData.Blank_SAP.Whisk.LH_MeanCBV,'d','MarkerFaceColor',colors('north texas green'),'MarkerEdgeColor','k','LineWidth',1.5);
scatter(6,testData.Blank_SAP.Whisk.RH_MeanCBV,'d','MarkerFaceColor',colors('north texas green'),'MarkerEdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.Blank_SAP.Whisk.meanLH)
    x = [Blank_xInds(1,aa)*5,Blank_xInds(1,aa)*6];
    y = [testData.Blank_SAP.Whisk.meanLH(aa,1),testData.Blank_SAP.Whisk.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorWhisk,'jitter','off', 'jitterAmount',0.25)
end
%% Stim
% Naive
scatter(8,testData.Naive.Stim.LH_MeanCBV,'d','MarkerFaceColor',colors('sapphire'),'MarkerEdgeColor','k','LineWidth',1.5)
hold on
scatter(9,testData.Naive.Stim.RH_MeanCBV,'d','MarkerFaceColor',colors('sapphire'),'MarkerEdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.Naive.Stim.meanLH)
    x = [C57_xInds(1,aa)*8,C57_xInds(1,aa)*9];
    y = [testData.Naive.Stim.meanLH(aa,1),testData.Naive.Stim.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorStim,'jitter','off', 'jitterAmount',0.25)
end
% SSP-SAP
scatter(10,testData.SSP_SAP.Stim.LH_MeanCBV,'d','MarkerFaceColor',colors('electric purple'),'MarkerEdgeColor','k','LineWidth',1.5)
scatter(11,testData.SSP_SAP.Stim.RH_MeanCBV,'d','MarkerFaceColor',colors('electric purple'),'MarkerEdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.SSP_SAP.Stim.meanLH)
    x = [SSP_xInds(1,aa)*10,SSP_xInds(1,aa)*11];
    y = [testData.SSP_SAP.Stim.meanLH(aa,1),testData.SSP_SAP.Stim.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorStim,'jitter','off', 'jitterAmount',0.25)
end
% Blank-SAP
scatter(12,testData.Blank_SAP.Stim.LH_MeanCBV,'d','MarkerFaceColor',colors('north texas green'),'MarkerEdgeColor','k','LineWidth',1.5)
scatter(13,testData.Blank_SAP.Stim.RH_MeanCBV,'d','MarkerFaceColor',colors('north texas green'),'MarkerEdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.Blank_SAP.Stim.meanLH)
    x = [Blank_xInds(1,aa)*12,Blank_xInds(1,aa)*13];
    y = [testData.Blank_SAP.Stim.meanLH(aa,1),testData.Blank_SAP.Stim.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorStim,'jitter','off', 'jitterAmount',0.25)
end
%% NREM
% Naive
scatter(15,testData.Naive.NREM.LH_MeanCBV,'d','MarkerFaceColor',colors('sapphire'),'MarkerEdgeColor','k','LineWidth',1.5)
hold on
scatter(16,testData.Naive.NREM.RH_MeanCBV,'d','MarkerFaceColor',colors('sapphire'),'MarkerEdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.Naive.NREM.meanLH)
    x = [C57_xInds(1,aa)*15,C57_xInds(1,aa)*16];
    y = [testData.Naive.NREM.meanLH(aa,1),testData.Naive.NREM.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorNREM,'jitter','off', 'jitterAmount',0.25)
end
% SSP-SAP
scatter(17,testData.SSP_SAP.NREM.LH_MeanCBV,'d','MarkerFaceColor',colors('electric purple'),'MarkerEdgeColor','k','LineWidth',1.5)
scatter(18,testData.SSP_SAP.NREM.RH_MeanCBV,'d','MarkerFaceColor',colors('electric purple'),'MarkerEdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.SSP_SAP.NREM.meanLH)
    x = [SSP_xInds(1,aa)*17,SSP_xInds(1,aa)*18];
    y = [testData.SSP_SAP.NREM.meanLH(aa,1),testData.SSP_SAP.NREM.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorNREM,'jitter','off', 'jitterAmount',0.25)
end
% Blank-SAP
scatter(19,testData.Blank_SAP.NREM.LH_MeanCBV,'d','MarkerFaceColor',colors('north texas green'),'MarkerEdgeColor','k','LineWidth',1.5)
scatter(20,testData.Blank_SAP.NREM.RH_MeanCBV,'d','MarkerFaceColor',colors('north texas green'),'MarkerEdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.Blank_SAP.NREM.meanLH)
    x = [Blank_xInds(1,aa)*19,Blank_xInds(1,aa)*20];
    y = [testData.Blank_SAP.NREM.meanLH(aa,1),testData.Blank_SAP.NREM.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorNREM,'jitter','off', 'jitterAmount',0.25)
end
%% REM
% Naive
scatter(22,testData.Naive.REM.LH_MeanCBV,'d','MarkerFaceColor',colors('sapphire'),'MarkerEdgeColor','k','LineWidth',1.5)
hold on
scatter(23,testData.Naive.REM.RH_MeanCBV,'d','MarkerFaceColor',colors('sapphire'),'MarkerEdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.Naive.REM.meanLH)
    x = [C57_xInds(1,aa)*22,C57_xInds(1,aa)*23];
    y = [testData.Naive.REM.meanLH(aa,1),testData.Naive.REM.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorREM,'jitter','off', 'jitterAmount',0.25)
end
% SSP-SAP
scatter(24,testData.SSP_SAP.REM.LH_MeanCBV,'d','MarkerFaceColor',colors('electric purple'),'MarkerEdgeColor','k','LineWidth',1.5)
scatter(25,testData.SSP_SAP.REM.RH_MeanCBV,'d','MarkerFaceColor',colors('electric purple'),'MarkerEdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.SSP_SAP.REM.meanLH)
    x = [SSP_xInds(1,aa)*24,SSP_xInds(1,aa)*25];
    y = [testData.SSP_SAP.REM.meanLH(aa,1),testData.SSP_SAP.REM.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorREM,'jitter','off', 'jitterAmount',0.25)
end
% Blank-SAP
scatter(26,testData.Blank_SAP.REM.LH_MeanCBV,'d','MarkerFaceColor',colors('north texas green'),'MarkerEdgeColor','k','LineWidth',1.5)
scatter(27,testData.Blank_SAP.REM.RH_MeanCBV,'d','MarkerFaceColor',colors('north texas green'),'MarkerEdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.Blank_SAP.REM.meanLH)
    x = [Blank_xInds(1,aa)*26,Blank_xInds(1,aa)*27];
    y = [testData.Blank_SAP.REM.meanLH(aa,1),testData.Blank_SAP.REM.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorREM,'jitter','off', 'jitterAmount',0.25)
end
%% Iso
% Naive
scatter(29,testData.Naive.Iso.LH_MeanCBV,'d','MarkerFaceColor',colors('sapphire'),'MarkerEdgeColor','k','LineWidth',1.5)
hold on
scatter(30,testData.Naive.Iso.RH_MeanCBV,'d','MarkerFaceColor',colors('sapphire'),'MarkerEdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.Naive.Iso.meanLH)
    x = [C57_xInds(1,aa)*29,C57_xInds(1,aa)*30];
    y = [testData.Naive.Iso.meanLH(aa,1),testData.Naive.Iso.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorIso,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorIso,'jitter','off', 'jitterAmount',0.25)
end
% SSP-SAP
scatter(31,testData.SSP_SAP.Iso.LH_MeanCBV,'d','MarkerFaceColor',colors('electric purple'),'MarkerEdgeColor','k','LineWidth',1.5)
scatter(32,testData.SSP_SAP.Iso.RH_MeanCBV,'d','MarkerFaceColor',colors('electric purple'),'MarkerEdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.SSP_SAP.Iso.meanLH)
    x = [SSP_xInds(1,aa)*31,SSP_xInds(1,aa)*32];
    y = [testData.SSP_SAP.Iso.meanLH(aa,1),testData.SSP_SAP.Iso.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorIso,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorIso,'jitter','off', 'jitterAmount',0.25)
end
% Blank-SAP
scatter(33,testData.Blank_SAP.Iso.LH_MeanCBV,'d','MarkerFaceColor',colors('north texas green'),'MarkerEdgeColor','k','LineWidth',1.5)
scatter(34,testData.Blank_SAP.Iso.RH_MeanCBV,'d','MarkerFaceColor',colors('north texas green'),'MarkerEdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.Blank_SAP.Iso.meanLH)
    x = [Blank_xInds(1,aa)*33,Blank_xInds(1,aa)*34];
    y = [testData.Blank_SAP.Iso.meanLH(aa,1),testData.Blank_SAP.Iso.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorIso,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorIso,'jitter','off', 'jitterAmount',0.25)
end
%% figure characteristics
title({'Mean \DeltaHbT (\muM)','during arousal-states'})
ylabel('\DeltaHbT (\muM)')
legend([b1,b2,b3],'Naive (UnRx black, RH red)','SSP-SAP','Blank-SAP','Location','NorthWest')
set(gca,'xtick',[3.5,10.5,17.5,24.5,31.5])
set(gca,'xticklabel',{'Whisk','Stim','NREM','REM','Isoflurane'})
xtickangle(45)
axis square
xlim([0,35])
ylim([-5,250])
set(gca,'box','off')
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Arousal State Hemodynamics - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Arousal_HbT']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Arousal_HbT'])
end

end
