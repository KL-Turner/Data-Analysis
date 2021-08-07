function [AnalysisResults] = MeanHbT_Saporin(rootFolder,saveFigs,delim,AnalysisResults)
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
expGroups = {'C57BL6J','SSP-SAP','Blank-SAP'};
setName = 'IOS Set A';
animalIDs.all = {};
for aa = 1:length(expGroups)
    folderList = dir([expGroups{1,aa} delim setName]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs.all = horzcat(animalIDs.all,{folderList.name});
    animalIDs.(strrep(expGroups{1,aa},'-','_')) = {folderList.name};
end
treatments = {'C57BL6J','SSP_SAP','Blank_SAP'};
behavFields = {'Rest','Whisk','Stim','NREM','REM'};
%% mean HbT comparison between behaviors
% pre-allocate the date for each day
for aa = 1:length(animalIDs.all)
    animalID = animalIDs.all{1,aa};
    whiskFileIDs = unique(AnalysisResults.(animalID).MeanCBV.Whisk.CBV_HbT.FileIDs);
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
            data.HbT.(animalID).(behavField).(fileDate).MeanLH = [];
            data.HbT.(animalID).(behavField).(fileDate).MeanRH = [];
            data.HbT.(animalID).(behavField).(fileDate).IndLH = {};
            data.HbT.(animalID).(behavField).(fileDate).IndRH = {};
        end
        procData.HbT.(behavField).animalID{aa,1} = animalID;
        procData.HbT.(behavField).behavior{aa,1} = behavField;
        procData.HbT.(behavField).LH{aa,1} = 'LH';
        procData.HbT.(behavField).RH{aa,1} = 'RH';
    end
end
% put data into cell for each unique date
for ff = 1:length(animalIDs.all)
    animalID = animalIDs.all{1,ff};
    for gg = 1:length(behavFields)
        behavField = behavFields{1,gg};
        % data is structured slightly differently depending on class
        if strcmp(behavField,'Rest') == true || strcmp(behavField,'Whisk') == true
            fileIDs = AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.FileIDs;
            for hh = 1:length(fileIDs)
                fileDate = ConvertDate_IOS(fileIDs{hh,1});
                data.HbT.(animalID).(behavField).(fileDate).MeanLH = cat(1,data.HbT.(animalID).(behavField).(fileDate).MeanLH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.MeanAdjLH(hh,1));
                data.HbT.(animalID).(behavField).(fileDate).MeanRH = cat(1,data.HbT.(animalID).(behavField).(fileDate).MeanRH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.MeanAdjRH(hh,1));
                data.HbT.(animalID).(behavField).(fileDate).IndLH = cat(1,data.HbT.(animalID).(behavField).(fileDate).IndLH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjLH{hh,1});
                data.HbT.(animalID).(behavField).(fileDate).IndRH = cat(1,data.HbT.(animalID).(behavField).(fileDate).IndRH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjRH{hh,1});
            end
        elseif strcmp(behavField,'Stim') == true
            % left hem stims
            fileIDs = AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.LH_FileIDs;
            for hh = 1:length(fileIDs)
                fileDate = ConvertDate_IOS(fileIDs{hh,1});
                data.HbT.(animalID).(behavField).(fileDate).MeanLH = cat(1,data.HbT.(animalID).(behavField).(fileDate).MeanLH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.MeanAdjLH(hh,1));
                data.HbT.(animalID).(behavField).(fileDate).IndLH = cat(1,data.HbT.(animalID).(behavField).(fileDate).IndLH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjLH{hh,1});
            end
            % right hem stims
            fileIDs = AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.RH_FileIDs;
            for hh = 1:length(fileIDs)
                fileDate = ConvertDate_IOS(fileIDs{hh,1});
                data.HbT.(animalID).(behavField).(fileDate).MeanRH = cat(1,data.HbT.(animalID).(behavField).(fileDate).MeanRH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.MeanAdjRH(hh,1));
                data.HbT.(animalID).(behavField).(fileDate).IndRH = cat(1,data.HbT.(animalID).(behavField).(fileDate).IndRH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjRH{hh,1});
            end
        else
            fileIDs = AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.FileIDs;
            for ii = 1:length(fileIDs)
                fileDate = ConvertDate_IOS(fileIDs{ii,1});
                data.HbT.(animalID).(behavField).(fileDate).MeanLH = cat(1,data.HbT.(animalID).(behavField).(fileDate).MeanLH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.MeanAdjLH(ii,1));
                data.HbT.(animalID).(behavField).(fileDate).MeanRH = cat(1,data.HbT.(animalID).(behavField).(fileDate).MeanRH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.MeanAdjRH(ii,1));
                data.HbT.(animalID).(behavField).(fileDate).IndLH = cat(1,data.HbT.(animalID).(behavField).(fileDate).IndLH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjLH{ii,1});
                data.HbT.(animalID).(behavField).(fileDate).IndRH = cat(1,data.HbT.(animalID).(behavField).(fileDate).IndRH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjRH{ii,1});
            end
        end
    end
end
% find the mean of the 10-second resting periods from each day to determine a baseline
for jj = 1:length(animalIDs.all)
    animalID = animalIDs.all{1,jj};
    whiskFileIDs = unique(AnalysisResults.(animalID).MeanCBV.Whisk.CBV_HbT.FileIDs);
    whiskFileDates = [];
    for kk = 1:length(whiskFileIDs)
        whiskFileDates{kk,1} = ConvertDate_IOS(whiskFileIDs{kk,1}); %#ok<*AGROW>
    end
    uniqueWhiskFileDates = unique(whiskFileDates);
    % take mean from each day. Days with no data will show up as NaN and be excluded
    for ll = 1:length(uniqueWhiskFileDates)
        fileDate = uniqueWhiskFileDates{ll,1};
        data.HbT.(animalID).Rest.(fileDate).baselineLH = mean(data.HbT.(animalID).Rest.(fileDate).MeanLH);
        data.HbT.(animalID).Rest.(fileDate).baselineRH = mean(data.HbT.(animalID).Rest.(fileDate).MeanRH);
    end
end
% subtract the 10-second resting baseline for each day from the other data types. If the day doesn't have resting data,
% exclude it from analysis
for mm = 1:length(animalIDs.all)
    animalID = animalIDs.all{1,mm};
    for nn = 1:length(behavFields)
        behavField = behavFields{1,nn};
        % subtract each day's 10-second baseline from each behavior field
        fileDates = fieldnames(data.HbT.(animalID).(behavField));
        for oo = 1:length(fileDates)
            fileDate = fileDates{oo,1};
            if strcmp(behavField,'Stim') == true || strcmp(behavField,'Whisk') == true
                data.HbT.(animalID).(behavField).(fileDate).CorrMeanLH = data.HbT.(animalID).(behavField).(fileDate).MeanLH;% - data.HbT.(animalID).Rest.(fileDate).baselineLH;
                data.HbT.(animalID).(behavField).(fileDate).CorrMeanRH = data.HbT.(animalID).(behavField).(fileDate).MeanRH;%; - data.HbT.(animalID).Rest.(fileDate).baselineRH;
                for pp = 1:length(data.HbT.(animalID).(behavField).(fileDate).IndLH)
                    data.HbT.(animalID).(behavField).(fileDate).CorrIndLH{pp,1} = data.HbT.(animalID).(behavField).(fileDate).IndLH{pp,1};% - data.HbT.(animalID).Rest.(fileDate).baselineLH;
                end
                for pp = 1:length(data.HbT.(animalID).(behavField).(fileDate).IndRH)
                    data.HbT.(animalID).(behavField).(fileDate).CorrIndRH{pp,1} = data.HbT.(animalID).(behavField).(fileDate).IndRH{pp,1};% - data.HbT.(animalID).Rest.(fileDate).baselineRH;
                end
            else
                data.HbT.(animalID).(behavField).(fileDate).CorrMeanLH = data.HbT.(animalID).(behavField).(fileDate).MeanLH - data.HbT.(animalID).Rest.(fileDate).baselineLH;
                data.HbT.(animalID).(behavField).(fileDate).CorrMeanRH = data.HbT.(animalID).(behavField).(fileDate).MeanRH - data.HbT.(animalID).Rest.(fileDate).baselineRH;
                for pp = 1:length(data.HbT.(animalID).(behavField).(fileDate).IndLH)
                    data.HbT.(animalID).(behavField).(fileDate).CorrIndLH{pp,1} = data.HbT.(animalID).(behavField).(fileDate).IndLH{pp,1} - data.HbT.(animalID).Rest.(fileDate).baselineLH;
                end
                for pp = 1:length(data.HbT.(animalID).(behavField).(fileDate).IndRH)
                    data.HbT.(animalID).(behavField).(fileDate).CorrIndRH{pp,1} = data.HbT.(animalID).(behavField).(fileDate).IndRH{pp,1} - data.HbT.(animalID).Rest.(fileDate).baselineRH;
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
        fileDates = fieldnames(data.HbT.(animalID).(behavField));
        for ss = 1:length(fileDates)
            fileDate = fileDates{ss,1};
            data.HbT.(animalID).(behavField).(fileDate).DayMeanLH = mean(data.HbT.(animalID).(behavField).(fileDate).CorrMeanLH);
            data.HbT.(animalID).(behavField).(fileDate).DayMeanRH = mean(data.HbT.(animalID).(behavField).(fileDate).CorrMeanRH);
            data.HbT.(animalID).(behavField).(fileDate).DayAllMeanLH = [];
            data.HbT.(animalID).(behavField).(fileDate).DayAllMeanRH = [];
            data.HbT.(animalID).(behavField).(fileDate).DayIndLH = [];
            data.HbT.(animalID).(behavField).(fileDate).DayIndRH = [];
            % concatenate individual trials into a single array for each unique day
            if isfield(data.HbT.(animalID).(behavField).(fileDate),'CorrIndLH') == true
                % left means - diff loop is necessary as STIM field has diff number of events
                for tt = 1:length(data.HbT.(animalID).(behavField).(fileDate).CorrMeanLH)
                    data.HbT.(animalID).(behavField).(fileDate).DayAllMeanLH = cat(2,data.HbT.(animalID).(behavField).(fileDate).DayIndLH,data.HbT.(animalID).(behavField).(fileDate).CorrIndLH{tt,1});
                end
                % right means
                for tt = 1:length(data.HbT.(animalID).(behavField).(fileDate).CorrMeanRH)
                    data.HbT.(animalID).(behavField).(fileDate).DayAllMeanRH = cat(2,data.HbT.(animalID).(behavField).(fileDate).DayIndRH,data.HbT.(animalID).(behavField).(fileDate).CorrIndRH{tt,1});
                end
                % left individual data pts
                for tt = 1:length(data.HbT.(animalID).(behavField).(fileDate).CorrIndLH)
                    data.HbT.(animalID).(behavField).(fileDate).DayIndLH = cat(2,data.HbT.(animalID).(behavField).(fileDate).DayIndLH,data.HbT.(animalID).(behavField).(fileDate).CorrIndLH{tt,1});
                end
                % right individual data pts
                for tt = 1:length(data.HbT.(animalID).(behavField).(fileDate).CorrIndRH)
                    data.HbT.(animalID).(behavField).(fileDate).DayIndRH = cat(2,data.HbT.(animalID).(behavField).(fileDate).DayIndRH,data.HbT.(animalID).(behavField).(fileDate).CorrIndRH{tt,1});
                end
            end
        end
    end
end
% P=put all the corrected means from each unique day into a single vector
nans = 1;
for uu = 1:length(animalIDs.all)
    animalID = animalIDs.all{1,uu};
    for vv = 1:length(behavFields)
        behavField = behavFields{1,vv};
        fileDates = fieldnames(data.HbT.(animalID).(behavField));
        procData.HbT.(animalID).(behavField).DayMeansLH = [];
        procData.HbT.(animalID).(behavField).DayMeansRH = [];
        procData.HbT.(animalID).(behavField).CatIndLH = [];
        procData.HbT.(animalID).(behavField).CatIndRH = [];
        for ww = 1:length(fileDates)
            fileDate = fileDates{ww,1};
            if isnan(data.HbT.(animalID).(behavField).(fileDate).DayMeanLH) == false
                procData.HbT.(animalID).(behavField).DayMeansLH = cat(1,procData.HbT.(animalID).(behavField).DayMeansLH,data.HbT.(animalID).(behavField).(fileDate).DayMeanLH);
                procData.HbT.(animalID).(behavField).DayMeansRH = cat(1,procData.HbT.(animalID).(behavField).DayMeansRH,data.HbT.(animalID).(behavField).(fileDate).DayMeanRH);
                procData.HbT.(animalID).(behavField).CatIndLH = cat(2,procData.HbT.(animalID).(behavField).CatIndLH,data.HbT.(animalID).(behavField).(fileDate).DayIndLH);
                procData.HbT.(animalID).(behavField).CatIndRH = cat(2,procData.HbT.(animalID).(behavField).CatIndRH,data.HbT.(animalID).(behavField).(fileDate).DayIndRH);
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
    if ismember(animalIDs.all{1,zz},animalIDs.C57BL6J) == true
        treatment = 'C57BL6J';
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
        end
        if strcmp(behavField,'Iso') == false
            testData.(treatment).(behavField).meanLH = cat(1,testData.(treatment).(behavField).meanLH,mean(procData.HbT.(animalID).(behavField).DayMeansLH));
            testData.(treatment).(behavField).meanRH = cat(1,testData.(treatment).(behavField).meanRH,mean(procData.HbT.(animalID).(behavField).DayMeansRH));
        else
            testData.(treatment).(behavField).meanLH = cat(1,testData.(treatment).(behavField).meanLH,mean(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.adjLH));
            testData.(treatment).(behavField).meanRH = cat(1,testData.(treatment).(behavField).meanRH,mean(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.adjRH));
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
%% mean HbT during different behaviors
summaryFigure = figure('Name','Fig5 (a-f)');
C57_xInds = ones(1,length(testData.C57BL6J.Rest.meanLH));
SSP_xInds = ones(1,length(testData.SSP_SAP.Rest.meanLH));
Blank_xInds = ones(1,length(testData.Blank_SAP.Rest.meanLH));
%% Whisk
% C57BL6J
b1 = bar(1,testData.C57BL6J.Whisk.LH_MeanCBV,'FaceColor',colors('sapphire'),'EdgeColor','k','LineWidth',1.5);
hold on
bar(2,testData.C57BL6J.Whisk.RH_MeanCBV,'FaceColor',colors('sapphire'),'EdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.C57BL6J.Whisk.meanLH)
    x = [C57_xInds(1,aa)*1,C57_xInds(1,aa)*2];
    y = [testData.C57BL6J.Whisk.meanLH(aa,1),testData.C57BL6J.Whisk.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorWhisk,'jitter','off', 'jitterAmount',0.25)
end
% SSP-SAP
b2 = bar(3,testData.SSP_SAP.Whisk.LH_MeanCBV,'FaceColor',colors('electric purple'),'EdgeColor','k','LineWidth',1.5);
bar(4,testData.SSP_SAP.Whisk.RH_MeanCBV,'FaceColor',colors('electric purple'),'EdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.SSP_SAP.Whisk.meanLH)
    x = [SSP_xInds(1,aa)*3,SSP_xInds(1,aa)*4];
    y = [testData.SSP_SAP.Whisk.meanLH(aa,1),testData.SSP_SAP.Whisk.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorWhisk,'jitter','off', 'jitterAmount',0.25)
end
% Blank-SAP
b3 = bar(5,testData.Blank_SAP.Whisk.LH_MeanCBV,'FaceColor',colors('north texas green'),'EdgeColor','k','LineWidth',1.5);
bar(6,testData.Blank_SAP.Whisk.RH_MeanCBV,'FaceColor',colors('north texas green'),'EdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.Blank_SAP.Whisk.meanLH)
    x = [Blank_xInds(1,aa)*5,Blank_xInds(1,aa)*6];
    y = [testData.Blank_SAP.Whisk.meanLH(aa,1),testData.Blank_SAP.Whisk.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorWhisk,'jitter','off', 'jitterAmount',0.25)
end
%% Stim
% C57BL6J
bar(8,testData.C57BL6J.Stim.LH_MeanCBV,'FaceColor',colors('sapphire'),'EdgeColor','k','LineWidth',1.5)
hold on
bar(9,testData.C57BL6J.Stim.RH_MeanCBV,'FaceColor',colors('sapphire'),'EdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.C57BL6J.Stim.meanLH)
    x = [C57_xInds(1,aa)*8,C57_xInds(1,aa)*9];
    y = [testData.C57BL6J.Stim.meanLH(aa,1),testData.C57BL6J.Stim.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorStim,'jitter','off', 'jitterAmount',0.25)
end
% SSP-SAP
bar(10,testData.SSP_SAP.Stim.LH_MeanCBV,'FaceColor',colors('electric purple'),'EdgeColor','k','LineWidth',1.5)
bar(11,testData.SSP_SAP.Stim.RH_MeanCBV,'FaceColor',colors('electric purple'),'EdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.SSP_SAP.Stim.meanLH)
    x = [SSP_xInds(1,aa)*10,SSP_xInds(1,aa)*11];
    y = [testData.SSP_SAP.Stim.meanLH(aa,1),testData.SSP_SAP.Stim.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorStim,'jitter','off', 'jitterAmount',0.25)
end
% Blank-SAP
bar(12,testData.Blank_SAP.Stim.LH_MeanCBV,'FaceColor',colors('north texas green'),'EdgeColor','k','LineWidth',1.5)
bar(13,testData.Blank_SAP.Stim.RH_MeanCBV,'FaceColor',colors('north texas green'),'EdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.Blank_SAP.Stim.meanLH)
    x = [Blank_xInds(1,aa)*12,Blank_xInds(1,aa)*13];
    y = [testData.Blank_SAP.Stim.meanLH(aa,1),testData.Blank_SAP.Stim.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorStim,'jitter','off', 'jitterAmount',0.25)
end
%% NREM
% C57BL6J
bar(15,testData.C57BL6J.NREM.LH_MeanCBV,'FaceColor',colors('sapphire'),'EdgeColor','k','LineWidth',1.5)
hold on
bar(16,testData.C57BL6J.NREM.RH_MeanCBV,'FaceColor',colors('sapphire'),'EdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.C57BL6J.NREM.meanLH)
    x = [C57_xInds(1,aa)*15,C57_xInds(1,aa)*16];
    y = [testData.C57BL6J.NREM.meanLH(aa,1),testData.C57BL6J.NREM.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorNREM,'jitter','off', 'jitterAmount',0.25)
end
% SSP-SAP
bar(17,testData.SSP_SAP.NREM.LH_MeanCBV,'FaceColor',colors('electric purple'),'EdgeColor','k','LineWidth',1.5)
bar(18,testData.SSP_SAP.NREM.RH_MeanCBV,'FaceColor',colors('electric purple'),'EdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.SSP_SAP.NREM.meanLH)
    x = [SSP_xInds(1,aa)*17,SSP_xInds(1,aa)*18];
    y = [testData.SSP_SAP.NREM.meanLH(aa,1),testData.SSP_SAP.NREM.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorNREM,'jitter','off', 'jitterAmount',0.25)
end
% Blank-SAP
bar(19,testData.Blank_SAP.NREM.LH_MeanCBV,'FaceColor',colors('north texas green'),'EdgeColor','k','LineWidth',1.5)
bar(20,testData.Blank_SAP.NREM.RH_MeanCBV,'FaceColor',colors('north texas green'),'EdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.Blank_SAP.NREM.meanLH)
    x = [Blank_xInds(1,aa)*19,Blank_xInds(1,aa)*20];
    y = [testData.Blank_SAP.NREM.meanLH(aa,1),testData.Blank_SAP.NREM.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorNREM,'jitter','off', 'jitterAmount',0.25)
end
%% REM
% C57BL6J
bar(22,testData.C57BL6J.REM.LH_MeanCBV,'FaceColor',colors('sapphire'),'EdgeColor','k','LineWidth',1.5)
hold on
bar(23,testData.C57BL6J.REM.RH_MeanCBV,'FaceColor',colors('sapphire'),'EdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.C57BL6J.REM.meanLH)
    x = [C57_xInds(1,aa)*22,C57_xInds(1,aa)*23];
    y = [testData.C57BL6J.REM.meanLH(aa,1),testData.C57BL6J.REM.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorREM,'jitter','off', 'jitterAmount',0.25)
end
% SSP-SAP
bar(24,testData.SSP_SAP.REM.LH_MeanCBV,'FaceColor',colors('electric purple'),'EdgeColor','k','LineWidth',1.5)
bar(25,testData.SSP_SAP.REM.RH_MeanCBV,'FaceColor',colors('electric purple'),'EdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.SSP_SAP.REM.meanLH)
    x = [SSP_xInds(1,aa)*24,SSP_xInds(1,aa)*25];
    y = [testData.SSP_SAP.REM.meanLH(aa,1),testData.SSP_SAP.REM.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorREM,'jitter','off', 'jitterAmount',0.25)
end
% Blank-SAP
bar(26,testData.Blank_SAP.REM.LH_MeanCBV,'FaceColor',colors('north texas green'),'EdgeColor','k','LineWidth',1.5)
bar(27,testData.Blank_SAP.REM.RH_MeanCBV,'FaceColor',colors('north texas green'),'EdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.Blank_SAP.REM.meanLH)
    x = [Blank_xInds(1,aa)*26,Blank_xInds(1,aa)*27];
    y = [testData.Blank_SAP.REM.meanLH(aa,1),testData.Blank_SAP.REM.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorREM,'jitter','off', 'jitterAmount',0.25)
end
%% Iso
% C57BL6J
bar(29,testData.C57BL6J.Iso.LH_MeanCBV,'FaceColor',colors('sapphire'),'EdgeColor','k','LineWidth',1.5)
hold on
bar(30,testData.C57BL6J.Iso.RH_MeanCBV,'FaceColor',colors('sapphire'),'EdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.C57BL6J.Iso.meanLH)
    x = [C57_xInds(1,aa)*29,C57_xInds(1,aa)*30];
    y = [testData.C57BL6J.Iso.meanLH(aa,1),testData.C57BL6J.Iso.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorIso,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorIso,'jitter','off', 'jitterAmount',0.25)
end
% SSP-SAP
bar(31,testData.SSP_SAP.Iso.LH_MeanCBV,'FaceColor',colors('electric purple'),'EdgeColor','k','LineWidth',1.5)
bar(32,testData.SSP_SAP.Iso.RH_MeanCBV,'FaceColor',colors('electric purple'),'EdgeColor','r','LineWidth',1.5)
for aa = 1:length(testData.SSP_SAP.Iso.meanLH)
    x = [SSP_xInds(1,aa)*31,SSP_xInds(1,aa)*32];
    y = [testData.SSP_SAP.Iso.meanLH(aa,1),testData.SSP_SAP.Iso.meanRH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colorIso,'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','r','MarkerFaceColor',colorIso,'jitter','off', 'jitterAmount',0.25)
end
% Blank-SAP
bar(33,testData.Blank_SAP.Iso.LH_MeanCBV,'FaceColor',colors('north texas green'),'EdgeColor','k','LineWidth',1.5)
bar(34,testData.Blank_SAP.Iso.RH_MeanCBV,'FaceColor',colors('north texas green'),'EdgeColor','r','LineWidth',1.5)
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
legend([b1,b2,b3],'C57BL6J (UnRx black, RH red)','SSP-SAP','Blank-SAP','Location','NorthWest')
set(gca,'xtick',[3.5,10.5,17.5,24.5,31.5])
set(gca,'xticklabel',{'Whisk','Stim','NREM','REM','Isoflurane'})
xtickangle(45)
axis square
xlim([0,35])
ylim([-5,250])
set(gca,'box','off')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Arousal_HbT']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Arousal_HbT'])
end

end
