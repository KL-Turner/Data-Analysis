function [AnalysisResults] = MeanHbT_Saporin(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 5 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

% colorBlack = [(0/256),(0/256),(0/256)];
% colorGrey = [(209/256),(211/256),(212/256)];
% colorRfcAwake = [(0/256),(64/256),(64/256)];
% colorRfcNREM = [(0/256),(174/256),(239/256)];
% colorRfcREM = [(190/256),(30/256),(45/256)];
colorRest = [(0/256),(166/256),(81/256)];
colorWhisk = [(31/256),(120/256),(179/256)];
colorStim = [(255/256),(28/256),(206/256)];
colorNREM = [(191/256),(0/256),(255/256)];
colorREM = [(254/256),(139/256),(0/256)];
% colorAlert = [(255/256),(191/256),(0/256)];
% colorAsleep = [(0/256),(128/256),(255/256)];
% colorAll = [(183/256),(115/256),(51/256)];
% colorIso = [(0/256),(256/256),(256/256)];
%% set-up and process data
animalIDs = {'T135','T141','T142','T144','T151','T155','T156','T157','T159'};
C57BL6J_IDs = {'T141','T155','T156','T157'};
SSP_SAP_IDs = {'T135','T142','T144','T151','T159'};
treatments = {'C57BL6J','SSP_SAP'};
%% mean HbT comparison between behaviors
% pre-allocate the date for each day
behavFields = {'Rest','Whisk','Stim','NREM','REM','Iso'};
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
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
for ff = 1:length(animalIDs)
    animalID = animalIDs{1,ff};
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
for jj = 1:length(animalIDs)
    animalID = animalIDs{1,jj};
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
for mm = 1:length(animalIDs)
    animalID = animalIDs{1,mm};
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
for qq = 1:length(animalIDs)
    animalID = animalIDs{1,qq};
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
for uu = 1:length(animalIDs)
    animalID = animalIDs{1,uu};
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
xxx = 1;
zzz = 1;
for zz = 1:length(animalIDs)
    animalID = animalIDs{1,zz};
    if ismember(animalID,C57BL6J_IDs) == true
        treatment = 'C57BL6J';
        for yy = 1:length(behavFields)
            behavField = behavFields{1,yy};
            testData.(treatment).(behavField).meanLH(xxx,1) = mean(procData.HbT.(animalID).(behavField).DayMeansLH);
            testData.(treatment).(behavField).meanRH(xxx,1) = mean(procData.HbT.(animalID).(behavField).DayMeansRH);
        end
        xxx = xxx + 1;
    elseif ismember(animalID,SSP_SAP_IDs) == true
        treatment = 'SSP_SAP';
        for yy = 1:length(behavFields)
            behavField = behavFields{1,yy};
            testData.(treatment).(behavField).meanLH(zzz,1) = mean(procData.HbT.(animalID).(behavField).DayMeansLH);
            testData.(treatment).(behavField).meanRH(zzz,1) = mean(procData.HbT.(animalID).(behavField).DayMeansRH);
        end
        zzz = zzz + 1;
    end
end
% take the mean and stdev across animals
for qqq = 1:length(treatments)
    treatment = treatments{1,qqq};
    for aaa = 1:length(behavFields)
        behavField = behavFields{1,aaa};
        testData.(treatment).(behavField).LH_MeanCBV = mean(testData.(treatment).(behavField).meanLH,1);
        testData.(treatment).(behavField).LH_StdMeanCBV = std(testData.(treatment).(behavField).meanLH,0,1);
        testData.(treatment).(behavField).RH_MeanCBV = mean(testData.(treatment).(behavField).meanRH,1);
        testData.(treatment).(behavField).RH_StdMeanCBV = std(testData.(treatment).(behavField).meanRH,0,1);
    end
end
%% mean HbT during different behaviors
summaryFigure = figure('Name','Fig5 (a-f)');
C57_xInds = ones(1,length(testData.C57BL6J.Rest.meanLH));
SSP_xInds = ones(1,length(testData.SSP_SAP.Rest.meanLH));
Blank_xInds = ones(1,length(testData.Blank_SAP.Rest.meanLH));
s1 = scatter(C57*1,testData.C57BL6J.Rest.meanLH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,testData.C57BL6J.Rest.LH_MeanCBV,testData.C57BL6J.Rest.LH_StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(C57*2,testData.C57BL6J.Rest.meanRH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,testData.C57BL6J.Rest.RH_MeanCBV,testData.C57BL6J.Rest.RH_StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(SAP*3,testData.SSP_SAP.Rest.meanLH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,testData.SSP_SAP.Rest.LH_MeanCBV,testData.SSP_SAP.Rest.LH_StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
s4 = scatter(SAP*4,testData.SSP_SAP.Rest.meanRH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,testData.SSP_SAP.Rest.RH_MeanCBV,testData.SSP_SAP.Rest.RH_StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
s5 = scatter(C57*5,testData.C57BL6J.Whisk.meanLH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,testData.C57BL6J.Whisk.LH_MeanCBV,testData.C57BL6J.Whisk.LH_StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
s6 = scatter(C57*6,testData.C57BL6J.Whisk.meanRH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e6 = errorbar(6,testData.C57BL6J.Whisk.RH_MeanCBV,testData.C57BL6J.Whisk.RH_StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
s7 = scatter(SAP*7,testData.SSP_SAP.Whisk.meanLH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e7 = errorbar(7,testData.SSP_SAP.Whisk.LH_MeanCBV,testData.SSP_SAP.Whisk.LH_StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
s8 = scatter(SAP*8,testData.SSP_SAP.Whisk.meanRH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e8 = errorbar(8,testData.SSP_SAP.Whisk.RH_MeanCBV,testData.SSP_SAP.Whisk.RH_StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
e8.MarkerSize = 10;
e8.CapSize = 10;
s9 = scatter(C57*9,testData.C57BL6J.Stim.meanLH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e9 = errorbar(9,testData.C57BL6J.Stim.LH_MeanCBV,testData.C57BL6J.Stim.LH_StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e9.Color = 'black';
e9.MarkerSize = 10;
e9.CapSize = 10;
s10 = scatter(C57*10,testData.C57BL6J.Stim.meanRH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e10 = errorbar(10,testData.C57BL6J.Stim.RH_MeanCBV,testData.C57BL6J.Stim.RH_StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e10.Color = 'black';
e10.MarkerSize = 10;
e10.CapSize = 10;
s11 = scatter(SAP*11,testData.SSP_SAP.Stim.meanLH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e11 = errorbar(11,testData.SSP_SAP.Stim.LH_MeanCBV,testData.SSP_SAP.Stim.LH_StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e11.Color = 'black';
e11.MarkerSize = 10;
e11.CapSize = 10;
s12 = scatter(SAP*12,testData.SSP_SAP.Stim.meanRH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e12 = errorbar(12,testData.SSP_SAP.Stim.RH_MeanCBV,testData.SSP_SAP.Stim.RH_StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e12.Color = 'black';
e12.MarkerSize = 10;
e12.CapSize = 10;
s13 = scatter(C57*13,testData.C57BL6J.NREM.meanLH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e13 = errorbar(13,testData.C57BL6J.NREM.LH_MeanCBV,testData.C57BL6J.NREM.LH_StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e13.Color = 'black';
e13.MarkerSize = 10;
e13.CapSize = 10;
s14 = scatter(C57*14,testData.C57BL6J.NREM.meanRH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e14 = errorbar(14,testData.C57BL6J.NREM.RH_MeanCBV,testData.C57BL6J.NREM.RH_StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e14.Color = 'black';
e14.MarkerSize = 10;
e14.CapSize = 10;
s15 = scatter(SAP*15,testData.SSP_SAP.NREM.meanLH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e15 = errorbar(15,testData.SSP_SAP.NREM.LH_MeanCBV,testData.SSP_SAP.NREM.LH_StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e15.Color = 'black';
e15.MarkerSize = 10;
e15.CapSize = 10;
s16 = scatter(SAP*16,testData.SSP_SAP.NREM.meanRH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e16 = errorbar(16,testData.SSP_SAP.NREM.RH_MeanCBV,testData.SSP_SAP.NREM.RH_StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e16.Color = 'black';
e16.MarkerSize = 10;
e16.CapSize = 10;
s17 = scatter(C57*17,testData.C57BL6J.REM.meanLH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e17 = errorbar(17,testData.C57BL6J.REM.LH_MeanCBV,testData.C57BL6J.REM.LH_StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e17.Color = 'black';
e17.MarkerSize = 10;
e17.CapSize = 10;
s18 = scatter(C57*18,testData.C57BL6J.REM.meanRH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e18 = errorbar(18,testData.C57BL6J.REM.RH_MeanCBV,testData.C57BL6J.REM.RH_StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e18.Color = 'black';
e18.MarkerSize = 10;
e18.CapSize = 10;
s19 = scatter(SAP*19,testData.SSP_SAP.REM.meanLH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e19 = errorbar(19,testData.SSP_SAP.REM.LH_MeanCBV,testData.SSP_SAP.REM.LH_StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e19.Color = 'black';
e19.MarkerSize = 10;
e19.CapSize = 10;
s20 = scatter(SAP*20,testData.SSP_SAP.REM.meanRH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e20 = errorbar(20,testData.SSP_SAP.REM.RH_MeanCBV,testData.SSP_SAP.REM.RH_StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e20.Color = 'black';
e20.MarkerSize = 10;
e20.CapSize = 10;
%% figure characteristics
title({'[5a] Mean \Delta[C57BL6J] (\muM)','during arousal-states'})
ylabel('\DeltaC57BL6J (\muM)')
legend([s1,s5,s9,s13,s17],'Rest','Whisk','Stim','NREM','REM','Location','NorthWest')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,21])
% ylim([-10,100])
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
