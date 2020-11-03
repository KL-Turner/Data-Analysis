function [AnalysisResults] = Fig5_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 5 for Turner_Gheres_Proctor_Drew_eLife2020
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
IOS_animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
TwoP_animalIDs = {'T115','T116','T117','T118','T125','T126'};
LDF_animalIDs = {'T108','T109','T110','T111','T119','T120','T121','T122'};
%% mean HbT comparison between behaviors
% pre-allocate the date for each day
IOS_behavFields = {'Rest','Whisk','Stim','NREM','REM'};
for aa = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,aa};
    whiskFileIDs = unique(AnalysisResults.(animalID).MeanCBV.Whisk.CBV_HbT.FileIDs);
    whiskFileDates = [];
    % identify the unique days present for each animal using the whisking field.
    for bb = 1:length(whiskFileIDs)
        whiskFileDates{bb,1} = ConvertDate_IOS_eLife2020(whiskFileIDs{bb,1}); %#ok<*AGROW>
    end
    uniqueWhiskFileDates = unique(whiskFileDates);
    % put pre-allocate each date
    for dd = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,dd};
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
for ff = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,ff};
    for gg = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,gg};
        % data is structured slightly differently depending on class
        if strcmp(behavField,'Rest') == true || strcmp(behavField,'Whisk') == true
            fileIDs = AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.FileIDs;
            for hh = 1:length(fileIDs)
                fileDate = ConvertDate_IOS_eLife2020(fileIDs{hh,1});
                data.HbT.(animalID).(behavField).(fileDate).MeanLH = cat(1,data.HbT.(animalID).(behavField).(fileDate).MeanLH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.MeanAdjLH(hh,1));
                data.HbT.(animalID).(behavField).(fileDate).MeanRH = cat(1,data.HbT.(animalID).(behavField).(fileDate).MeanRH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.MeanAdjRH(hh,1));
                data.HbT.(animalID).(behavField).(fileDate).IndLH = cat(1,data.HbT.(animalID).(behavField).(fileDate).IndLH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjLH{hh,1});
                data.HbT.(animalID).(behavField).(fileDate).IndRH = cat(1,data.HbT.(animalID).(behavField).(fileDate).IndRH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjRH{hh,1});
            end
        elseif strcmp(behavField,'Stim') == true
            % left hem stims
            fileIDs = AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.LH_FileIDs;
            for hh = 1:length(fileIDs)
                fileDate = ConvertDate_IOS_eLife2020(fileIDs{hh,1});
                data.HbT.(animalID).(behavField).(fileDate).MeanLH = cat(1,data.HbT.(animalID).(behavField).(fileDate).MeanLH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.MeanAdjLH(hh,1));
                data.HbT.(animalID).(behavField).(fileDate).IndLH = cat(1,data.HbT.(animalID).(behavField).(fileDate).IndLH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjLH{hh,1});
            end
            % right hem stims
            fileIDs = AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.RH_FileIDs;
            for hh = 1:length(fileIDs)
                fileDate = ConvertDate_IOS_eLife2020(fileIDs{hh,1});
                data.HbT.(animalID).(behavField).(fileDate).MeanRH = cat(1,data.HbT.(animalID).(behavField).(fileDate).MeanRH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.MeanAdjRH(hh,1));
                data.HbT.(animalID).(behavField).(fileDate).IndRH = cat(1,data.HbT.(animalID).(behavField).(fileDate).IndRH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjRH{hh,1});
            end
        else
            fileIDs = AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.FileIDs;
            for ii = 1:length(fileIDs)
                fileDate = ConvertDate_IOS_eLife2020(fileIDs{ii,1});
                data.HbT.(animalID).(behavField).(fileDate).MeanLH = cat(1,data.HbT.(animalID).(behavField).(fileDate).MeanLH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.MeanAdjLH(ii,1));
                data.HbT.(animalID).(behavField).(fileDate).MeanRH = cat(1,data.HbT.(animalID).(behavField).(fileDate).MeanRH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.MeanAdjRH(ii,1));
                data.HbT.(animalID).(behavField).(fileDate).IndLH = cat(1,data.HbT.(animalID).(behavField).(fileDate).IndLH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjLH{ii,1});
                data.HbT.(animalID).(behavField).(fileDate).IndRH = cat(1,data.HbT.(animalID).(behavField).(fileDate).IndRH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjRH{ii,1});
            end
        end
    end
end
% find the mean of the 10-second resting periods from each day to determine a baseline
for jj = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,jj};
    whiskFileIDs = unique(AnalysisResults.(animalID).MeanCBV.Whisk.CBV_HbT.FileIDs);
    whiskFileDates = [];
    for kk = 1:length(whiskFileIDs)
        whiskFileDates{kk,1} = ConvertDate_IOS_eLife2020(whiskFileIDs{kk,1}); %#ok<*AGROW>
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
for mm = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,mm};
    for nn = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,nn};
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
for qq = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,qq};
    for rr = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,rr};
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
for uu = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,uu};
    for vv = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,vv};
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
% put all the means (of the corrected means) from each unique day into a single vector
for yy = 1:length(IOS_behavFields)
    behavField = IOS_behavFields{1,yy};
    procData.HbT.(behavField).IndMeanCBV = [];
    procData.HbT.(behavField).CatCBV = [];
    procData.HbT.(behavField).meanLH = [];
    procData.HbT.(behavField).meanRH = [];
    for zz = 1:length(IOS_animalIDs)
        animalID = IOS_animalIDs{1,zz};
        procData.HbT.(behavField).IndMeanCBV = cat(1,procData.HbT.(behavField).IndMeanCBV,mean(procData.HbT.(animalID).(behavField).DayMeansLH),mean(procData.HbT.(animalID).(behavField).DayMeansRH));
        procData.HbT.(behavField).meanLH = cat(1,procData.HbT.(behavField).meanLH,mean(procData.HbT.(animalID).(behavField).DayMeansLH));
        procData.HbT.(behavField).meanRH = cat(1,procData.HbT.(behavField).meanRH,mean(procData.HbT.(animalID).(behavField).DayMeansRH));
        procData.HbT.(behavField).CatCBV = cat(2,procData.HbT.(behavField).CatCBV,procData.HbT.(animalID).(behavField).CatIndLH,procData.HbT.(animalID).(behavField).CatIndRH);
    end
end
% take the mean and stdev across animals
for aaa = 1:length(IOS_behavFields)
    behavField = IOS_behavFields{1,aaa};
    procData.HbT.(behavField).MeanCBV = mean(procData.HbT.(behavField).IndMeanCBV,1);
    procData.HbT.(behavField).StdMeanCBV = std(procData.HbT.(behavField).IndMeanCBV,0,1);
end
%% statistics - generalized linear mixed effects model
HbTtableSize = cat(1,procData.HbT.Rest.meanLH,procData.HbT.Rest.meanRH,procData.HbT.Whisk.meanLH,procData.HbT.Whisk.meanRH,...
    procData.HbT.Stim.meanLH,procData.HbT.Stim.meanRH,procData.HbT.NREM.meanLH,procData.HbT.NREM.meanRH,procData.HbT.REM.meanLH,procData.HbT.REM.meanRH);
HbTTable = table('Size',[size(HbTtableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','HbT','Behavior','Hemisphere'});
HbTTable.Mouse = cat(1,procData.HbT.Rest.animalID,procData.HbT.Rest.animalID,procData.HbT.Whisk.animalID,procData.HbT.Whisk.animalID,...
    procData.HbT.Stim.animalID,procData.HbT.Stim.animalID,procData.HbT.NREM.animalID,procData.HbT.NREM.animalID,procData.HbT.REM.animalID,procData.HbT.REM.animalID);
HbTTable.HbT = cat(1,procData.HbT.Rest.meanLH,procData.HbT.Rest.meanRH,procData.HbT.Whisk.meanLH,procData.HbT.Whisk.meanRH,...
    procData.HbT.Stim.meanLH,procData.HbT.Stim.meanRH,procData.HbT.NREM.meanLH,procData.HbT.NREM.meanRH,procData.HbT.REM.meanLH,procData.HbT.REM.meanRH);
HbTTable.Behavior = cat(1,procData.HbT.Rest.behavior,procData.HbT.Rest.behavior,procData.HbT.Whisk.behavior,procData.HbT.Whisk.behavior,...
    procData.HbT.Stim.behavior,procData.HbT.Stim.behavior,procData.HbT.NREM.behavior,procData.HbT.NREM.behavior,procData.HbT.REM.behavior,procData.HbT.REM.behavior);
HbTTable.Hemisphere = cat(1,procData.HbT.Rest.LH,procData.HbT.Rest.RH,procData.HbT.Whisk.LH,procData.HbT.Whisk.RH,...
    procData.HbT.Stim.LH,procData.HbT.Stim.RH,procData.HbT.NREM.LH,procData.HbT.NREM.RH,procData.HbT.REM.LH,procData.HbT.REM.RH);
HbTFitFormula = 'HbT ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
HbTStats = fitglme(HbTTable,HbTFitFormula); %#ok<*NASGU>
%% mean vessel diameter comparison between behaviors
% pre-allocate the date for each day
TwoP_behavFields = {'Rest','Whisk','NREM','REM'};
for aa = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,aa};
    whiskFileIDs = unique(AnalysisResults.(animalID).MeanVesselDiameter.Whisk.allFileIDs);
    whiskFileDates = [];
    % identify the unique days present for each animal using the whisking field.
    for bb = 1:length(whiskFileIDs)
        whiskFileDates{bb,1} = ConvertDate_2P_eLife2020(whiskFileIDs{bb,1}); %#ok<*AGROW>
    end
    uniqueWhiskFileDates = unique(whiskFileDates);
    % put pre-allocate each date
    for dd = 1:length(TwoP_behavFields)
        behavField = TwoP_behavFields{1,dd};
        if isfield(AnalysisResults.(animalID).MeanVesselDiameter,behavField) == true
            vIDs = fieldnames(AnalysisResults.(animalID).MeanVesselDiameter.(behavField));
            for qq = 1:length(vIDs)
                vID = vIDs{qq,1};
                if strcmp(vID,'allFileIDs') == false && strcmp(vID,'V1') == false
                    for ee = 1:length(uniqueWhiskFileDates)
                        fileDate = uniqueWhiskFileDates{ee,1};
                        data.TwoP.(animalID).(behavField).(vID).(fileDate).mean = [];
                        data.TwoP.(animalID).(behavField).(vID).(fileDate).max = [];
                        data.TwoP.(animalID).(behavField).(vID).(fileDate).indData = {};
                    end
                end
            end
        end
    end
end
% put data into cell for each unique date
for ff = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,ff};
    for gg = 1:length(TwoP_behavFields)
        behavField = TwoP_behavFields{1,gg};
        % data is structured slightly differently depending on class
        if isfield(data.TwoP.(animalID),behavField) == true
            vIDs = fieldnames(data.TwoP.(animalID).(behavField));
            for qq = 1:length(vIDs)
                vID = vIDs{qq,1};
                fileIDs = AnalysisResults.(animalID).MeanVesselDiameter.(behavField).(vID).fileIDs;
                for hh = 1:length(fileIDs)
                    fileDate = ConvertDate_2P_eLife2020(fileIDs{hh,1});
                    data.TwoP.(animalID).(behavField).(vID).(fileDate).mean = cat(1,data.TwoP.(animalID).(behavField).(vID).(fileDate).mean,AnalysisResults.(animalID).MeanVesselDiameter.(behavField).(vID).mean(hh,1));
                    data.TwoP.(animalID).(behavField).(vID).(fileDate).max = cat(1,data.TwoP.(animalID).(behavField).(vID).(fileDate).max,AnalysisResults.(animalID).MeanVesselDiameter.(behavField).(vID).max(hh,1));
                    data.TwoP.(animalID).(behavField).(vID).(fileDate).indData = cat(1,data.TwoP.(animalID).(behavField).(vID).(fileDate).indData,AnalysisResults.(animalID).MeanVesselDiameter.(behavField).(vID).indEvents{hh,1});
                end
            end
        end
    end
end
% find the mean of the 10-second resting periods from each day to determine a baseline
for jj = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,jj};
    whiskFileIDs = unique(AnalysisResults.(animalID).MeanVesselDiameter.Whisk.allFileIDs);
    whiskFileDates = [];
    for kk = 1:length(whiskFileIDs)
        whiskFileDates{kk,1} = ConvertDate_IOS_eLife2020(whiskFileIDs{kk,1}); %#ok<*AGROW>
    end
    uniqueWhiskFileDates = unique(whiskFileDates);
    % take mean from each day. Days with no data will show up as NaN and be excluded
    vIDs = fieldnames(data.TwoP.(animalID).Rest);
    for qq = 1:length(vIDs)
        vID = vIDs{qq,1};
        for ll = 1:length(uniqueWhiskFileDates)
            fileDate = uniqueWhiskFileDates{ll,1};
            data.TwoP.(animalID).Rest.(vID).(fileDate).baseline = mean(data.TwoP.(animalID).Rest.(vID).(fileDate).mean);
        end
    end
end
% subtract the 10-second resting baseline for each day from the other data types. If the day doesn't have resting data,
% exclude it from analysis
for mm = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,mm};
    for nn = 1:length(TwoP_behavFields)
        behavField = TwoP_behavFields{1,nn};
        % subtract each day's 10-second baseline from each behavior field
        if isfield(data.TwoP.(animalID),behavField) == true
            vIDs = fieldnames(data.TwoP.(animalID).(behavField));
            for qq = 1:length(vIDs)
                vID = vIDs{qq,1};
                fileDates = fieldnames(data.TwoP.(animalID).(behavField).(vID));
                for oo = 1:length(fileDates)
                    fileDate = fileDates{oo,1};
                    if strcmp(behavField,'Whisk') == true
                        data.TwoP.(animalID).(behavField).(vID).(fileDate).CorrMean = data.TwoP.(animalID).(behavField).(vID).(fileDate).mean;
                        data.TwoP.(animalID).(behavField).(vID).(fileDate).CorrMax = data.TwoP.(animalID).(behavField).(vID).(fileDate).max;
                        for pp = 1:length(data.TwoP.(animalID).(behavField).(vID).(fileDate).indData)
                            data.TwoP.(animalID).(behavField).(vID).(fileDate).CorrInd{pp,1} = data.TwoP.(animalID).(behavField).(vID).(fileDate).indData{pp,1};
                        end
                    else
                        data.TwoP.(animalID).(behavField).(vID).(fileDate).CorrMean = data.TwoP.(animalID).(behavField).(vID).(fileDate).mean - data.TwoP.(animalID).Rest.(vID).(fileDate).baseline;
                        data.TwoP.(animalID).(behavField).(vID).(fileDate).CorrMax = data.TwoP.(animalID).(behavField).(vID).(fileDate).max - data.TwoP.(animalID).Rest.(vID).(fileDate).baseline;
                        for pp = 1:length(data.TwoP.(animalID).(behavField).(vID).(fileDate).indData)
                            data.TwoP.(animalID).(behavField).(vID).(fileDate).CorrInd{pp,1} = data.TwoP.(animalID).(behavField).(vID).(fileDate).indData{pp,1} - data.TwoP.(animalID).Rest.(vID).(fileDate).baseline;
                        end
                    end
                end
            end
        end
    end
end
% take the mean of the corrected data from each unique day
for qqq = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,qqq};
    for rr = 1:length(TwoP_behavFields)
        behavField = TwoP_behavFields{1,rr};
        if isfield(data.TwoP.(animalID),behavField) == true
            vIDs = fieldnames(data.TwoP.(animalID).(behavField));
            for qq = 1:length(vIDs)
                vID = vIDs{qq,1};
                fileDates = fieldnames(data.TwoP.(animalID).(behavField).(vID));
                for ss = 1:length(fileDates)
                    fileDate = fileDates{ss,1};
                    data.TwoP.(animalID).(behavField).(vID).(fileDate).DayMean = mean(data.TwoP.(animalID).(behavField).(vID).(fileDate).CorrMean);
                    data.TwoP.(animalID).(behavField).(vID).(fileDate).DayMax = mean(data.TwoP.(animalID).(behavField).(vID).(fileDate).CorrMax);
                    data.TwoP.(animalID).(behavField).(vID).(fileDate).DayInd = [];
                    % concatenate individual trials into a single array for each unique day
                    if isfield(data.TwoP.(animalID).(behavField).(vID).(fileDate),'CorrInd') == true
                        for tt = 1:length(data.TwoP.(animalID).(behavField).(vID).(fileDate).CorrInd)
                            data.TwoP.(animalID).(behavField).(vID).(fileDate).DayInd = cat(2,data.TwoP.(animalID).(behavField).(vID).(fileDate).DayInd,data.TwoP.(animalID).(behavField).(vID).(fileDate).CorrInd{tt,1});
                        end
                    end
                end
            end
        end
    end
end
% put all the corrected means from each unique day into a single vector
for uu = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,uu};
    for vv = 1:length(TwoP_behavFields)
        behavField = TwoP_behavFields{1,vv};
        if isfield(data.TwoP.(animalID),behavField) == true
            vIDs = fieldnames(data.TwoP.(animalID).(behavField));
            for qq = 1:length(vIDs)
                vID = vIDs{qq,1};
                fileDates = fieldnames(data.TwoP.(animalID).(behavField).(vID));
                procData.TwoP.(animalID).(behavField).(vID).DayMeans = [];
                procData.TwoP.(animalID).(behavField).(vID).DayMaxs = [];
                procData.TwoP.(animalID).(behavField).(vID).CatInd = [];
                for ww = 1:length(fileDates)
                    fileDate = fileDates{ww,1};
                    if isnan(data.TwoP.(animalID).(behavField).(vID).(fileDate).DayMean) == false
                        procData.TwoP.(animalID).(behavField).(vID).DayMeans = cat(1,procData.TwoP.(animalID).(behavField).(vID).DayMeans,data.TwoP.(animalID).(behavField).(vID).(fileDate).DayMean);
                        procData.TwoP.(animalID).(behavField).(vID).DayMaxs = cat(1,procData.TwoP.(animalID).(behavField).(vID).DayMaxs,data.TwoP.(animalID).(behavField).(vID).(fileDate).DayMax);
                        procData.TwoP.(animalID).(behavField).(vID).CatInd = cat(2,procData.TwoP.(animalID).(behavField).(vID).CatInd,data.TwoP.(animalID).(behavField).(vID).(fileDate).DayInd);
                    end
                end
            end
        end
    end
end
% put all the means (of the corrected means) from each unique day into a single vector
for yy = 1:length(TwoP_behavFields)
    behavField = TwoP_behavFields{1,yy};
    procData.TwoP.(behavField).IndMeanDiam = [];
    procData.TwoP.(behavField).IndMaxDiam = [];
    procData.TwoP.(behavField).CatIndDiam = [];
    idx = 1;
    for zz = 1:length(TwoP_animalIDs)
        animalID = TwoP_animalIDs{1,zz};
        if isfield(procData.TwoP.(animalID),behavField) == true
            vIDs = fieldnames(procData.TwoP.(animalID).(behavField));
            for qq = 1:length(vIDs)
                vID = vIDs{qq,1};
                procData.TwoP.(behavField).IndMeanDiam = cat(1,procData.TwoP.(behavField).IndMeanDiam,nanmean(procData.TwoP.(animalID).(behavField).(vID).DayMeans));
                procData.TwoP.(behavField).IndMaxDiam = cat(1,procData.TwoP.(behavField).IndMaxDiam,nanmean(procData.TwoP.(animalID).(behavField).(vID).DayMaxs));
                procData.TwoP.(behavField).CatIndDiam = cat(2,procData.TwoP.(behavField).CatIndDiam,procData.TwoP.(animalID).(behavField).(vID).CatInd);
                procData.TwoP.(behavField).animalID{idx,1} = animalID;
                procData.TwoP.(behavField).behavior{idx,1} = behavField;
                procData.TwoP.(behavField).vID{idx,1} = vID;
                idx = idx + 1;
            end
        end
    end
end
% take the mean and stdev across animals
for aaa = 1:length(TwoP_behavFields)
    behavField = TwoP_behavFields{1,aaa};
    procData.TwoP.(behavField).MeanDiam = nanmean(procData.TwoP.(behavField).IndMeanDiam,1);
    procData.TwoP.(behavField).StdMeanDiam = nanstd(procData.TwoP.(behavField).IndMeanDiam,0,1);
    procData.TwoP.(behavField).MaxDiam = nanmean(procData.TwoP.(behavField).IndMaxDiam,1);
    procData.TwoP.(behavField).StdMaxDiam = nanstd(procData.TwoP.(behavField).IndMaxDiam,0,1);
end
% statistics - linear mixed effects model
tableSize = cat(1,procData.TwoP.Rest.animalID,procData.TwoP.Whisk.animalID,procData.TwoP.NREM.animalID,procData.TwoP.REM.animalID);
vesselDiameterTable = table('Size',[size(tableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','Diameter','Behavior','Vessel'});
vesselDiameterTable.Mouse = cat(1,procData.TwoP.Rest.animalID,procData.TwoP.Whisk.animalID,procData.TwoP.NREM.animalID,procData.TwoP.REM.animalID);
vesselDiameterTable.Diameter = cat(1,procData.TwoP.Rest.IndMaxDiam,procData.TwoP.Whisk.IndMaxDiam,procData.TwoP.NREM.IndMaxDiam,procData.TwoP.REM.IndMaxDiam);
vesselDiameterTable.Behavior = cat(1,procData.TwoP.Rest.behavior,procData.TwoP.Whisk.behavior,procData.TwoP.NREM.behavior,procData.TwoP.REM.behavior);
vesselDiameterTable.Vessel = cat(1,procData.TwoP.Rest.vID,procData.TwoP.Whisk.vID,procData.TwoP.NREM.vID,procData.TwoP.REM.vID);
vesselFitFormula = 'Diameter ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Vessel)';
vesselStats = fitglme(vesselDiameterTable,vesselFitFormula);
%% LDf comparison between behaviors
% pre-allocate the date for each day
LDF_behavFields = {'Rest','Whisk','NREM','REM'};
for aa = 1:length(LDF_animalIDs)
    animalID = LDF_animalIDs{1,aa};
    for dd = 1:length(LDF_behavFields)
        behavField = LDF_behavFields{1,dd};
        data.LDF.(animalID).(behavField).mean = AnalysisResults.(animalID).LDFlow.(behavField).mean;
        data.LDF.(animalID).(behavField).indData = AnalysisResults.(animalID).LDFlow.(behavField).indData;
        procData.LDF.(behavField).animalID{aa,1} = animalID;
        procData.LDF.(behavField).behavior{aa,1} = behavField;
    end
end
% find the mean of the 10-second resting periods from each day to determine a baseline
for jj = 1:length(LDF_animalIDs)
    animalID = LDF_animalIDs{1,jj};
    data.LDF.(animalID).Rest.baseline = mean(data.LDF.(animalID).Rest.mean);
end
% subtract the 10-second resting baseline for each day from the other data types. If the day doesn't have resting data,
% exclude it from analysis
for mm = 1:length(LDF_animalIDs)
    animalID = LDF_animalIDs{1,mm};
    for nn = 1:length(LDF_behavFields)
        behavField = LDF_behavFields{1,nn};
        % subtract each day's 10-second baseline from each behavior field
        if strcmp(behavField,'Whisk') == true
            data.LDF.(animalID).(behavField).CorrMean = data.LDF.(animalID).(behavField).mean;
            for pp = 1:length(data.LDF.(animalID).(behavField).indData)
                data.LDF.(animalID).(behavField).CorrInd{pp,1} = data.LDF.(animalID).(behavField).indData{pp,1};
            end
        else
            data.LDF.(animalID).(behavField).CorrMean = data.LDF.(animalID).(behavField).mean - data.LDF.(animalID).Rest.baseline;
            for pp = 1:length(data.LDF.(animalID).(behavField).indData)
                data.LDF.(animalID).(behavField).CorrInd{pp,1} = data.LDF.(animalID).(behavField).indData{pp,1} - data.LDF.(animalID).Rest.baseline;
            end
        end
    end
end
% take the mean of the corrected data from each unique day
for qq = 1:length(LDF_animalIDs)
    animalID = LDF_animalIDs{1,qq};
    for rr = 1:length(LDF_behavFields)
        behavField = LDF_behavFields{1,rr};
        procData.LDF.(animalID).(behavField).DayMean = mean(data.LDF.(animalID).(behavField).CorrMean);
        procData.LDF.(animalID).(behavField).CatInd = [];
        if isnan(procData.LDF.(animalID).(behavField).DayMean) == false
            for zz = 1:length(data.LDF.(animalID).(behavField).CorrInd)
                procData.LDF.(animalID).(behavField).CatInd = cat(2,procData.LDF.(animalID).(behavField).CatInd,data.LDF.(animalID).(behavField).CorrInd{zz,1});
            end
        end
    end
end
% put all the means (of the corrected means) from each unique day into a single vector
for yy = 1:length(LDF_behavFields)
    behavField = LDF_behavFields{1,yy};
    procData.LDF.(behavField).IndMeanLDF = [];
    procData.LDF.(behavField).CatLDF = [];
    for zz = 1:length(LDF_animalIDs)
        animalID = LDF_animalIDs{1,zz};
        procData.LDF.(behavField).IndMeanLDF = cat(1,procData.LDF.(behavField).IndMeanLDF,procData.LDF.(animalID).(behavField).DayMean);
        procData.LDF.(behavField).CatLDF = cat(2,procData.LDF.(behavField).CatLDF,procData.LDF.(animalID).(behavField).CatInd);
    end
end
% take the mean and stdev across animals
for aaa = 1:length(LDF_behavFields)
    behavField = LDF_behavFields{1,aaa};
    procData.LDF.(behavField).MeanLDF = mean(procData.LDF.(behavField).IndMeanLDF,1);
    procData.LDF.(behavField).StdLDF = std(procData.LDF.(behavField).IndMeanLDF,0,1);
end
% statistics - linear mixed effects model
tableSize = cat(1,procData.LDF.Rest.IndMeanLDF,procData.LDF.Whisk.IndMeanLDF,procData.LDF.NREM.IndMeanLDF,procData.LDF.REM.IndMeanLDF);
flowTable = table('Size',[size(tableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','Flow','Behavior'});
flowTable.Mouse = cat(1,procData.LDF.Rest.animalID,procData.LDF.Whisk.animalID,procData.LDF.NREM.animalID,procData.LDF.REM.animalID);
flowTable.Flow = cat(1,procData.LDF.Rest.IndMeanLDF,procData.LDF.Whisk.IndMeanLDF,procData.LDF.NREM.IndMeanLDF,procData.LDF.REM.IndMeanLDF);
flowTable.Behavior = cat(1,procData.LDF.Rest.behavior,procData.LDF.Whisk.behavior,procData.LDF.NREM.behavior,procData.LDF.REM.behavior);
flowFitFormula = 'Flow ~ 1 + Behavior + (1|Mouse)';
flowStats = fitglme(flowTable,flowFitFormula);
%% Fig. 5
summaryFigure = figure('Name','Fig5 (a-f)');
sgtitle('Figure 5 - Turner et al. 2020')
%% [5a] mean HbT during different behaviors
ax1 = subplot(2,3,1);
HbT_xInds = ones(1,length(IOS_animalIDs)*2);
s1 = scatter(HbT_xInds*1,procData.HbT.Rest.IndMeanCBV,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,procData.HbT.Rest.MeanCBV,procData.HbT.Rest.StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(HbT_xInds*2,procData.HbT.Whisk.IndMeanCBV,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,procData.HbT.Whisk.MeanCBV,procData.HbT.Whisk.StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(HbT_xInds*3,procData.HbT.Stim.IndMeanCBV,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,procData.HbT.Stim.MeanCBV,procData.HbT.Stim.StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
s4 = scatter(HbT_xInds*4,procData.HbT.NREM.IndMeanCBV,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,procData.HbT.NREM.MeanCBV,procData.HbT.NREM.StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
s5 = scatter(HbT_xInds*5,procData.HbT.REM.IndMeanCBV,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,procData.HbT.REM.MeanCBV,procData.HbT.REM.StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
title({'[5a] Mean \Delta[HbT] (\muM)','during arousal-states'})
ylabel('\DeltaHbT (\muM)')
legend([s1,s2,s3,s4,s5],'Rest','Whisk','Stim','NREM','REM','Location','NorthWest')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(IOS_behavFields) + 1])
ylim([-10,100])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [5b] mean vessel diameter during different behaviors
ax2 = subplot(2,3,2);
TwoP_xIndsRest = ones(1,length(procData.TwoP.Rest.IndMeanDiam));
TwoP_xIndsWhisk = ones(1,length(procData.TwoP.Whisk.IndMeanDiam));
TwoP_xIndsNREM = ones(1,length(procData.TwoP.NREM.IndMeanDiam));
TwoP_xIndsREM = ones(1,length(procData.TwoP.REM.IndMeanDiam));
scatter(TwoP_xIndsRest*1,procData.TwoP.Rest.IndMeanDiam,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,procData.TwoP.Rest.MeanDiam,procData.TwoP.Rest.StdMeanDiam,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(TwoP_xIndsWhisk*2,procData.TwoP.Whisk.IndMeanDiam,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,procData.TwoP.Whisk.MeanDiam,procData.TwoP.Whisk.StdMeanDiam,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(TwoP_xIndsNREM*3,procData.TwoP.NREM.IndMeanDiam,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,procData.TwoP.NREM.MeanDiam,procData.TwoP.NREM.StdMeanDiam,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(TwoP_xIndsREM*4,procData.TwoP.REM.IndMeanDiam,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,procData.TwoP.REM.MeanDiam,procData.TwoP.REM.StdMeanDiam,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
title({'[5b] Mean \DeltaD/D (%)','during arousal-states'})
ylabel('\DeltaD/D (%)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(TwoP_behavFields) + 1])
ylim([-10,60])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [5c] mean vessel diameter during different behaviors
ax3 = subplot(2,3,3);
LDF_xInds = ones(1,length(LDF_animalIDs));
scatter(LDF_xInds*1,procData.LDF.Rest.IndMeanLDF,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,procData.LDF.Rest.MeanLDF,procData.LDF.Rest.StdLDF,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(LDF_xInds*2,procData.LDF.Whisk.IndMeanLDF,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,procData.LDF.Whisk.MeanLDF,procData.LDF.Whisk.StdLDF,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(LDF_xInds*3,procData.LDF.NREM.IndMeanLDF,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,procData.LDF.NREM.MeanLDF,procData.LDF.NREM.StdLDF,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(LDF_xInds*4,procData.LDF.REM.IndMeanLDF,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,procData.LDF.REM.MeanLDF,procData.LDF.REM.StdLDF,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
title({'[5c] Mean \DeltaQ/Q (%)','during arousal-states'})
ylabel('\DeltaQ/Q (%)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(LDF_behavFields) + 1])
ylim([-10,60])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [5a bottom] mean HbT distribution during different behaviors
ax4 = subplot(2,3,4);
[xCurve1,yCurve1] = SmoothHistogramBinsFit_eLife2020(procData.HbT.Rest.CatCBV,11,'normal');
[xCurve2,yCurve2] = SmoothHistogramBinsFit_eLife2020(procData.HbT.Whisk.CatCBV,11,'normal');
[xCurve3,yCurve3] = SmoothHistogramBinsFit_eLife2020(procData.HbT.Stim.CatCBV,11,'normal');
[xCurve4,yCurve4] = SmoothHistogramBinsFit_eLife2020(procData.HbT.NREM.CatCBV,11,'normal');
[xCurve5,yCurve5] = SmoothHistogramBinsFit_eLife2020(procData.HbT.REM.CatCBV,11,'normal');
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
hold on
plot(xCurve2,yCurve2,'color',colorWhisk,'LineWidth',2)
plot(xCurve3,yCurve3,'color',colorStim,'LineWidth',2)
plot(xCurve4,yCurve4,'color',colorNREM,'LineWidth',2)
plot(xCurve5,yCurve5,'color',colorREM,'LineWidth',2)
title({'\Delta[HbT] (\muM)','arousal-state distribution'})
xlabel('\Delta[HbT] (\muM)')
ylabel('Probability')
axis square
set(gca,'box','off')
xlim([-50,150])
ylim([0,0.6])
ax4.TickLength = [0.03,0.03];
%% [5b bottom] vessel diameter distribution during different behaviors
ax5 = subplot(2,3,5);
[xCurve1,yCurve1] = SmoothHistogramBinsFit_eLife2020(procData.TwoP.Rest.CatIndDiam,10,'normal');
[xCurve2,yCurve2] = SmoothHistogramBinsFit_eLife2020(procData.TwoP.Whisk.CatIndDiam,10,'normal');
[xCurve3,yCurve3] = SmoothHistogramBinsFit_eLife2020(procData.TwoP.NREM.CatIndDiam,10,'normal');
[xCurve4,yCurve4] = SmoothHistogramBinsFit_eLife2020(procData.TwoP.REM.CatIndDiam,10,'normal');
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
hold on
plot(xCurve2,yCurve2,'color',colorWhisk,'LineWidth',2)
plot(xCurve3,yCurve3,'color',colorNREM,'LineWidth',2)
plot(xCurve4,yCurve4,'color',colorREM,'LineWidth',2)
title({'\DeltaD/D (%)','arousal-state distribution'})
xlabel('\DeltaD/D (%)')
ylabel('Probability')
axis square
set(gca,'box','off')
xlim([-30,80])
ylim([0,0.6])
ax5.TickLength = [0.03,0.03];
%% [5c bottom] LDF arousal-state vessel distribution
ax6 = subplot(2,3,6);
[xCurve1,yCurve1] = SmoothHistogramBinsFit_eLife2020(procData.LDF.Rest.CatLDF,4,'normal');
[xCurve2,yCurve2] = SmoothHistogramBinsFit_eLife2020(procData.LDF.Whisk.CatLDF,11,'normal');
[xCurve3,yCurve3] = SmoothHistogramBinsFit_eLife2020(procData.LDF.NREM.CatLDF,11,'normal');
[xCurve4,yCurve4] = SmoothHistogramBinsFit_eLife2020(procData.LDF.REM.CatLDF,11,'normal');
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
hold on
plot(xCurve2,yCurve2,'color',colorWhisk,'LineWidth',2)
plot(xCurve3,yCurve3,'color',colorNREM,'LineWidth',2)
plot(xCurve4,yCurve4,'color',colorREM,'LineWidth',2)
title({'\DeltaQ/Q (%)','arousal-state distribution'})
xlabel('\DeltaQ/Q (%)')
ylabel('Probability')
axis square
xlim([-40,80])
ylim([0,0.6])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    set(summaryFigure,'PaperPositionMode','auto');
    savefig(summaryFigure,[dirpath 'Fig5']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Fig5'])
    %% statistical diary
    diaryFile = [dirpath 'Fig5_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    % HbT statistical diary
    disp('======================================================================================================================')
    disp('[5a] Generalized linear mixed-effects model statistics for mean HbT during Rest, Whisk, Stim, NREM, and REM')
    disp('======================================================================================================================')
    disp(HbTStats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  [HbT] (uM): Set = 0']); disp(' ') %#ok<NBRAK>
    disp(['Whisk [HbT] (uM): ' num2str(round(procData.HbT.Whisk.MeanCBV,1)) ' +/- ' num2str(round(procData.HbT.Whisk.StdMeanCBV,1))]); disp(' ')
    disp(['Stim  [HbT] (uM): ' num2str(round(procData.HbT.Stim.MeanCBV,1)) ' +/- ' num2str(round(procData.HbT.Stim.StdMeanCBV,1))]); disp(' ')
    disp(['NREM  [HbT] (uM): ' num2str(round(procData.HbT.NREM.MeanCBV,1)) ' +/- ' num2str(round(procData.HbT.NREM.StdMeanCBV,1))]); disp(' ')
    disp(['REM   [HbT] (uM): ' num2str(round(procData.HbT.REM.MeanCBV,1)) ' +/- ' num2str(round(procData.HbT.REM.StdMeanCBV,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % Peak vessel diameter statistical diary
    disp('======================================================================================================================')
    disp('[5b] Generalized linear mixed-effects model statistics for mean vessel diameter during Rest, Whisk, NREM, and REM')
    disp('======================================================================================================================')
    disp(vesselStats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  D/D (%): Set = 0']); disp(' ') %#ok<NBRAK>
    disp(['Whisk D/D (%): ' num2str(round(procData.TwoP.Whisk.MeanDiam,1)) ' +/- ' num2str(round(procData.TwoP.Whisk.StdMeanDiam,1))]); disp(' ')
    disp(['NREM  D/D (%): ' num2str(round(procData.TwoP.NREM.MeanDiam,1)) ' +/- ' num2str(round(procData.TwoP.NREM.StdMeanDiam,1))]); disp(' ')
    disp(['REM   D/D (%): ' num2str(round(procData.TwoP.REM.MeanDiam,1)) ' +/- ' num2str(round(procData.TwoP.REM.StdMeanDiam,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % LDF flow statistical diary
    disp('======================================================================================================================')
    disp('[5c] Generalized linear mixed-effects model statistics for mean doppler flow during Rest, Whisk, NREM, and REM')
    disp('======================================================================================================================')
    disp(flowStats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  Q/Q (%): Set = 0']); disp(' ') %#ok<NBRAK>
    disp(['Whisk Q/Q (%): ' num2str(round(procData.LDF.Whisk.MeanLDF,1)) ' +/- ' num2str(round(procData.LDF.Whisk.StdLDF,1))]); disp(' ')
    disp(['NREM  Q/Q (%): ' num2str(round(procData.LDF.NREM.MeanLDF,1)) ' +/- ' num2str(round(procData.LDF.NREM.StdLDF,1))]); disp(' ')
    disp(['REM   Q/Q (%): ' num2str(round(procData.LDF.REM.MeanLDF,1)) ' +/- ' num2str(round(procData.LDF.REM.StdLDF,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    diary off
end

end
