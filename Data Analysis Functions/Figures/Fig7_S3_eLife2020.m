function [AnalysisResults] = Fig7_S3_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 7-S3 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

% colorBlack = [(0/256),(0/256),(0/256)];
% colorGrey = [(209/256),(211/256),(212/256)];
% colorRfcAwake = [(0/256),(64/256),(64/256)];
% colorRfcNREM = [(0/256),(174/256),(239/256)];
% colorRfcREM = [(190/256),(30/256),(45/256)];
colorRest = [(0/256),(166/256),(81/256)];
% colorWhisk = [(31/256),(120/256),(179/256)];
% colorStim = [(255/256),(28/256),(206/256)];
colorNREM = [(191/256),(0/256),(255/256)];
colorREM = [(254/256),(139/256),(0/256)];
colorAlert = [(255/256),(191/256),(0/256)];
colorAsleep = [(0/256),(128/256),(255/256)];
colorAll = [(183/256),(115/256),(51/256)];
% colorIso = [(0/256),(256/256),(256/256)];
%% set-up and process data
IOS_animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
behavFields = {'Rest','NREM','REM','Awake','Sleep','All'};
dataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower'};
%% average coherence during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.Coherr = [];
for aa = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        % create the behavior folder for the first iteration of the loop
        if isfield(data.Coherr,behavField) == false
            data.Coherr.(behavField) = [];
        end
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(AnalysisResults.(animalID).Coherence.(behavField).(dataType).C) == false
                % create the data type folder for the first iteration of the loop
                if isfield(data.Coherr.(behavField),dataType) == false
                    data.Coherr.(behavField).(dataType).C = [];
                    data.Coherr.(behavField).(dataType).f = [];
                    data.Coherr.(behavField).(dataType).confC = [];
                    data.Coherr.(behavField).(dataType).animalID = {};
                    data.Coherr.(behavField).(dataType).behavior = {};
                end
                % concatenate C/f for existing data - exclude any empty sets
                data.Coherr.(behavField).(dataType).C = cat(2,data.Coherr.(behavField).(dataType).C,AnalysisResults.(animalID).Coherence.(behavField).(dataType).C);
                data.Coherr.(behavField).(dataType).f = cat(1,data.Coherr.(behavField).(dataType).f,AnalysisResults.(animalID).Coherence.(behavField).(dataType).f);
                data.Coherr.(behavField).(dataType).confC = cat(1,data.Coherr.(behavField).(dataType).confC,AnalysisResults.(animalID).Coherence.(behavField).(dataType).confC);
                if isempty(AnalysisResults.(animalID).Coherence.(behavField).(dataType).C) == false
                    data.Coherr.(behavField).(dataType).animalID = cat(1,data.Coherr.(behavField).(dataType).animalID,animalID);
                    data.Coherr.(behavField).(dataType).behavior = cat(1,data.Coherr.(behavField).(dataType).behavior,behavField);
                end
            end
        end
    end
end
% find 0.1/0.01 Hz peaks in coherence
for dd = 1:length(behavFields)
    behavField = behavFields{1,dd};
    for ee = 1:length(dataTypes)
        dataType = dataTypes{1,ee};
        for ff = 1:size(data.Coherr.(behavField).(dataType).C,2)
            if strcmp(behavField,'Rest') == true
                f_short = data.Coherr.(behavField).(dataType).f(ff,:);
                C = data.Coherr.(behavField).(dataType).C(:,ff);
                f_long = 0:0.01:0.5;
                C_long = interp1(f_short,C,f_long);
                index01 = find(f_long == 0.1);
                data.Coherr.(behavField).(dataType).C01(ff,1) = C_long(index01).^2; %#ok<FNDSB>
            elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
                F = round(data.Coherr.(behavField).(dataType).f(ff,:),2);
                C = data.Coherr.(behavField).(dataType).C(:,ff);
                index01 = find(F == 0.1);
                data.Coherr.(behavField).(dataType).C01(ff,1) = C(index01(1)).^2;
            else
                F = round(data.Coherr.(behavField).(dataType).f(ff,:),3);
                C = data.Coherr.(behavField).(dataType).C(:,ff);
                index01 = find(F == 0.1);
                index001 = find(F == 0.01);
                data.Coherr.(behavField).(dataType).C01(ff,1) = C(index01(1)).^2;
                data.Coherr.(behavField).(dataType).C001(ff,1) = C(index001(1)).^2;
            end
        end
    end
end
% take mean/StD of peak C
for gg = 1:length(behavFields)
    behavField = behavFields{1,gg};
    for hh = 1:length(dataTypes)
        dataType = dataTypes{1,hh};
        if strcmp(behavField,'Rest') == true || strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
            data.Coherr.(behavField).(dataType).meanC01 = mean(data.Coherr.(behavField).(dataType).C01,1);
            data.Coherr.(behavField).(dataType).stdC01 = std(data.Coherr.(behavField).(dataType).C01,0,1);
        else
            data.Coherr.(behavField).(dataType).meanC01 = mean(data.Coherr.(behavField).(dataType).C01,1);
            data.Coherr.(behavField).(dataType).stdC01 = std(data.Coherr.(behavField).(dataType).C01,0,1);
            data.Coherr.(behavField).(dataType).meanC001 = mean(data.Coherr.(behavField).(dataType).C001,1);
            data.Coherr.(behavField).(dataType).stdC001 = std(data.Coherr.(behavField).(dataType).C001,0,1);
        end
    end
end
%% statistics - generalized linear mixed effects model
% delta PSD @ 0.1 Hz
Delta_Coh01_tableSize = cat(1,data.Coherr.Rest.deltaBandPower.C01,data.Coherr.NREM.deltaBandPower.C01,data.Coherr.REM.deltaBandPower.C01,...
    data.Coherr.Awake.deltaBandPower.C01,data.Coherr.Sleep.deltaBandPower.C01,data.Coherr.All.deltaBandPower.C01);
Delta_Coh01_Table = table('Size',[size(Delta_Coh01_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','Coh01'});
Delta_Coh01_Table.Mouse = cat(1,data.Coherr.Rest.deltaBandPower.animalID,data.Coherr.NREM.deltaBandPower.animalID,data.Coherr.REM.deltaBandPower.animalID,...
    data.Coherr.Awake.deltaBandPower.animalID,data.Coherr.Sleep.deltaBandPower.animalID,data.Coherr.All.deltaBandPower.animalID);
Delta_Coh01_Table.Behavior = cat(1,data.Coherr.Rest.deltaBandPower.behavior,data.Coherr.NREM.deltaBandPower.behavior,data.Coherr.REM.deltaBandPower.behavior,...
    data.Coherr.Awake.deltaBandPower.behavior,data.Coherr.Sleep.deltaBandPower.behavior,data.Coherr.All.deltaBandPower.behavior);
Delta_Coh01_Table.Coh01 = cat(1,data.Coherr.Rest.deltaBandPower.C01,data.Coherr.NREM.deltaBandPower.C01,data.Coherr.REM.deltaBandPower.C01,...
    data.Coherr.Awake.deltaBandPower.C01,data.Coherr.Sleep.deltaBandPower.C01,data.Coherr.All.deltaBandPower.C01);
Delta_Coh01_FitFormula = 'Coh01 ~ 1 + Behavior + (1|Mouse)';
Delta_Coh01_Stats = fitglme(Delta_Coh01_Table,Delta_Coh01_FitFormula);
% delta PSD @ 0.01 Hz
Delta_Coh001_tableSize = cat(1,data.Coherr.Awake.deltaBandPower.C001,data.Coherr.Sleep.deltaBandPower.C001,data.Coherr.All.deltaBandPower.C001);
Delta_Coh001_Table = table('Size',[size(Delta_Coh001_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','Coh001'});
Delta_Coh001_Table.Mouse = cat(1,data.Coherr.Awake.deltaBandPower.animalID,data.Coherr.Sleep.deltaBandPower.animalID,data.Coherr.All.deltaBandPower.animalID);
Delta_Coh001_Table.Behavior = cat(1,data.Coherr.Awake.deltaBandPower.behavior,data.Coherr.Sleep.deltaBandPower.behavior,data.Coherr.All.deltaBandPower.behavior);
Delta_Coh001_Table.Coh001 = cat(1,data.Coherr.Awake.deltaBandPower.C001,data.Coherr.Sleep.deltaBandPower.C001,data.Coherr.All.deltaBandPower.C001);
Delta_Coh001_FitFormula = 'Coh001 ~ 1 + Behavior + (1|Mouse)';
Delta_Coh001_Stats = fitglme(Delta_Coh001_Table,Delta_Coh001_FitFormula);
% theta PSD @ 0.1 Hz
Theta_Coh01_tableSize = cat(1,data.Coherr.Rest.thetaBandPower.C01,data.Coherr.NREM.thetaBandPower.C01,data.Coherr.REM.thetaBandPower.C01,...
    data.Coherr.Awake.thetaBandPower.C01,data.Coherr.Sleep.thetaBandPower.C01,data.Coherr.All.thetaBandPower.C01);
Theta_Coh01_Table = table('Size',[size(Theta_Coh01_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','Coh01'});
Theta_Coh01_Table.Mouse = cat(1,data.Coherr.Rest.thetaBandPower.animalID,data.Coherr.NREM.thetaBandPower.animalID,data.Coherr.REM.thetaBandPower.animalID,...
    data.Coherr.Awake.thetaBandPower.animalID,data.Coherr.Sleep.thetaBandPower.animalID,data.Coherr.All.thetaBandPower.animalID);
Theta_Coh01_Table.Behavior = cat(1,data.Coherr.Rest.thetaBandPower.behavior,data.Coherr.NREM.thetaBandPower.behavior,data.Coherr.REM.thetaBandPower.behavior,...
    data.Coherr.Awake.thetaBandPower.behavior,data.Coherr.Sleep.thetaBandPower.behavior,data.Coherr.All.thetaBandPower.behavior);
Theta_Coh01_Table.Coh01 = cat(1,data.Coherr.Rest.thetaBandPower.C01,data.Coherr.NREM.thetaBandPower.C01,data.Coherr.REM.thetaBandPower.C01,...
    data.Coherr.Awake.thetaBandPower.C01,data.Coherr.Sleep.thetaBandPower.C01,data.Coherr.All.thetaBandPower.C01);
Theta_Coh01_FitFormula = 'Coh01 ~ 1 + Behavior + (1|Mouse)';
Theta_Coh01_Stats = fitglme(Theta_Coh01_Table,Theta_Coh01_FitFormula);
% theta PSD @ 0.01 Hz
Theta_Coh001_tableSize = cat(1,data.Coherr.Awake.thetaBandPower.C001,data.Coherr.Sleep.thetaBandPower.C001,data.Coherr.All.thetaBandPower.C001);
Theta_Coh001_Table = table('Size',[size(Theta_Coh001_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','Coh001'});
Theta_Coh001_Table.Mouse = cat(1,data.Coherr.Awake.thetaBandPower.animalID,data.Coherr.Sleep.thetaBandPower.animalID,data.Coherr.All.thetaBandPower.animalID);
Theta_Coh001_Table.Behavior = cat(1,data.Coherr.Awake.thetaBandPower.behavior,data.Coherr.Sleep.thetaBandPower.behavior,data.Coherr.All.thetaBandPower.behavior);
Theta_Coh001_Table.Coh001 = cat(1,data.Coherr.Awake.thetaBandPower.C001,data.Coherr.Sleep.thetaBandPower.C001,data.Coherr.All.thetaBandPower.C001);
Theta_Coh001_FitFormula = 'Coh001 ~ 1 + Behavior + (1|Mouse)';
Theta_Coh001_Stats = fitglme(Theta_Coh001_Table,Theta_Coh001_FitFormula);
% alpha PSD @ 0.1 Hz
Alpha_Coh01_tableSize = cat(1,data.Coherr.Rest.alphaBandPower.C01,data.Coherr.NREM.alphaBandPower.C01,data.Coherr.REM.alphaBandPower.C01,...
    data.Coherr.Awake.alphaBandPower.C01,data.Coherr.Sleep.alphaBandPower.C01,data.Coherr.All.alphaBandPower.C01);
Alpha_Coh01_Table = table('Size',[size(Alpha_Coh01_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','Coh01'});
Alpha_Coh01_Table.Mouse = cat(1,data.Coherr.Rest.alphaBandPower.animalID,data.Coherr.NREM.alphaBandPower.animalID,data.Coherr.REM.alphaBandPower.animalID,...
    data.Coherr.Awake.alphaBandPower.animalID,data.Coherr.Sleep.alphaBandPower.animalID,data.Coherr.All.alphaBandPower.animalID);
Alpha_Coh01_Table.Behavior = cat(1,data.Coherr.Rest.alphaBandPower.behavior,data.Coherr.NREM.alphaBandPower.behavior,data.Coherr.REM.alphaBandPower.behavior,...
    data.Coherr.Awake.alphaBandPower.behavior,data.Coherr.Sleep.alphaBandPower.behavior,data.Coherr.All.alphaBandPower.behavior);
Alpha_Coh01_Table.Coh01 = cat(1,data.Coherr.Rest.alphaBandPower.C01,data.Coherr.NREM.alphaBandPower.C01,data.Coherr.REM.alphaBandPower.C01,...
    data.Coherr.Awake.alphaBandPower.C01,data.Coherr.Sleep.alphaBandPower.C01,data.Coherr.All.alphaBandPower.C01);
Alpha_Coh01_FitFormula = 'Coh01 ~ 1 + Behavior + (1|Mouse)';
Alpha_Coh01_Stats = fitglme(Alpha_Coh01_Table,Alpha_Coh01_FitFormula);
% alpha PSD @ 0.01 Hz
Alpha_Coh001_tableSize = cat(1,data.Coherr.Awake.alphaBandPower.C001,data.Coherr.Sleep.alphaBandPower.C001,data.Coherr.All.alphaBandPower.C001);
Alpha_Coh001_Table = table('Size',[size(Alpha_Coh001_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','Coh001'});
Alpha_Coh001_Table.Mouse = cat(1,data.Coherr.Awake.alphaBandPower.animalID,data.Coherr.Sleep.alphaBandPower.animalID,data.Coherr.All.alphaBandPower.animalID);
Alpha_Coh001_Table.Behavior = cat(1,data.Coherr.Awake.alphaBandPower.behavior,data.Coherr.Sleep.alphaBandPower.behavior,data.Coherr.All.alphaBandPower.behavior);
Alpha_Coh001_Table.Coh001 = cat(1,data.Coherr.Awake.alphaBandPower.C001,data.Coherr.Sleep.alphaBandPower.C001,data.Coherr.All.alphaBandPower.C001);
Alpha_Coh001_FitFormula = 'Coh001 ~ 1 + Behavior + (1|Mouse)';
Alpha_Coh001_Stats = fitglme(Alpha_Coh001_Table,Alpha_Coh001_FitFormula);
% beta PSD @ 0.1 Hz
Beta_Coh01_tableSize = cat(1,data.Coherr.Rest.betaBandPower.C01,data.Coherr.NREM.betaBandPower.C01,data.Coherr.REM.betaBandPower.C01,...
    data.Coherr.Awake.betaBandPower.C01,data.Coherr.Sleep.betaBandPower.C01,data.Coherr.All.betaBandPower.C01);
Beta_Coh01_Table = table('Size',[size(Beta_Coh01_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','Coh01'});
Beta_Coh01_Table.Mouse = cat(1,data.Coherr.Rest.betaBandPower.animalID,data.Coherr.NREM.betaBandPower.animalID,data.Coherr.REM.betaBandPower.animalID,...
    data.Coherr.Awake.betaBandPower.animalID,data.Coherr.Sleep.betaBandPower.animalID,data.Coherr.All.betaBandPower.animalID);
Beta_Coh01_Table.Behavior = cat(1,data.Coherr.Rest.betaBandPower.behavior,data.Coherr.NREM.betaBandPower.behavior,data.Coherr.REM.betaBandPower.behavior,...
    data.Coherr.Awake.betaBandPower.behavior,data.Coherr.Sleep.betaBandPower.behavior,data.Coherr.All.betaBandPower.behavior);
Beta_Coh01_Table.Coh01 = cat(1,data.Coherr.Rest.betaBandPower.C01,data.Coherr.NREM.betaBandPower.C01,data.Coherr.REM.betaBandPower.C01,...
    data.Coherr.Awake.betaBandPower.C01,data.Coherr.Sleep.betaBandPower.C01,data.Coherr.All.betaBandPower.C01);
Beta_Coh01_FitFormula = 'Coh01 ~ 1 + Behavior + (1|Mouse)';
Beta_Coh01_Stats = fitglme(Beta_Coh01_Table,Beta_Coh01_FitFormula);
% beta PSD @ 0.01 Hz
Beta_Coh001_tableSize = cat(1,data.Coherr.Awake.betaBandPower.C001,data.Coherr.Sleep.betaBandPower.C001,data.Coherr.All.betaBandPower.C001);
Beta_Coh001_Table = table('Size',[size(Beta_Coh001_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','Coh001'});
Beta_Coh001_Table.Mouse = cat(1,data.Coherr.Awake.betaBandPower.animalID,data.Coherr.Sleep.betaBandPower.animalID,data.Coherr.All.betaBandPower.animalID);
Beta_Coh001_Table.Behavior = cat(1,data.Coherr.Awake.betaBandPower.behavior,data.Coherr.Sleep.betaBandPower.behavior,data.Coherr.All.betaBandPower.behavior);
Beta_Coh001_Table.Coh001 = cat(1,data.Coherr.Awake.betaBandPower.C001,data.Coherr.Sleep.betaBandPower.C001,data.Coherr.All.betaBandPower.C001);
Beta_Coh001_FitFormula = 'Coh001 ~ 1 + Behavior + (1|Mouse)';
Beta_Coh001_Stats = fitglme(Beta_Coh001_Table,Beta_Coh001_FitFormula);
%% Power spectra during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            data.PowerSpec.(behavField).(dataType).adjLH.S{aa,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjLH.S;
            data.PowerSpec.(behavField).(dataType).adjLH.f{aa,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjLH.f;
            data.PowerSpec.(behavField).(dataType).adjRH.S{aa,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjRH.S;
            data.PowerSpec.(behavField).(dataType).adjRH.f{aa,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjRH.f;
        end
    end
end
% find the peak of the resting PSD for each animal/hemisphere
for aa = 1:length(IOS_animalIDs)
    for cc = 1:length(dataTypes)
        dataType = dataTypes{1,cc};
        data.PowerSpec.baseline.(dataType).LH{aa,1} = max(data.PowerSpec.Rest.(dataType).adjLH.S{aa,1});
        data.PowerSpec.baseline.(dataType).RH{aa,1} = max(data.PowerSpec.Rest.(dataType).adjRH.S{aa,1});
    end
end
% DC-shift each animal/hemisphere/behavior PSD with respect to the resting peak
for aa = 1:length(IOS_animalIDs)
    for dd = 1:length(behavFields)
        behavField = behavFields{1,dd};
        for j = 1:length(dataTypes)
            dataType = dataTypes{1,j};
            for ee = 1:size(data.PowerSpec.(behavField).(dataType).adjLH.S,2)
                data.PowerSpec.(behavField).(dataType).normLH{aa,1} = (data.PowerSpec.(behavField).(dataType).adjLH.S{aa,1})*(1/(data.PowerSpec.baseline.(dataType).LH{aa,1}));
                data.PowerSpec.(behavField).(dataType).normRH{aa,1} = (data.PowerSpec.(behavField).(dataType).adjRH.S{aa,1})*(1/(data.PowerSpec.baseline.(dataType).RH{aa,1}));
            end
        end
    end
end
% concatenate the data from the left and right hemispheres - removes any empty data
for dd = 1:length(behavFields)
    behavField = behavFields{1,dd};
    for ee = 1:length(dataTypes)
        dataType = dataTypes{1,ee};
        data.PowerSpec.(behavField).(dataType).cat_S = [];
        data.PowerSpec.(behavField).(dataType).cat_f = [];
        data.PowerSpec.(behavField).(dataType).animalID = {};
        data.PowerSpec.(behavField).(dataType).behavior = {};
        data.PowerSpec.(behavField).(dataType).hemisphere = {};
        for zz = 1:length(data.PowerSpec.(behavField).(dataType).normLH)
            data.PowerSpec.(behavField).(dataType).cat_S = cat(2,data.PowerSpec.(behavField).(dataType).cat_S,data.PowerSpec.(behavField).(dataType).normLH{zz,1},data.PowerSpec.(behavField).(dataType).normRH{zz,1});
            data.PowerSpec.(behavField).(dataType).cat_f = cat(1,data.PowerSpec.(behavField).(dataType).cat_f,data.PowerSpec.(behavField).(dataType).adjLH.f{zz,1},data.PowerSpec.(behavField).(dataType).adjRH.f{zz,1});
            if isempty(data.PowerSpec.(behavField).(dataType).normLH{zz,1}) == false
                data.PowerSpec.(behavField).(dataType).animalID = cat(1,data.PowerSpec.(behavField).(dataType).animalID,animalID,animalID);
                data.PowerSpec.(behavField).(dataType).behavior = cat(1,data.PowerSpec.(behavField).(dataType).behavior,behavField,behavField);
                data.PowerSpec.(behavField).(dataType).hemisphere = cat(1,data.PowerSpec.(behavField).(dataType).hemisphere,'LH','RH');
            end
        end
    end
end
% find 0.1/0.01 Hz peaks in PSD
for dd = 1:length(behavFields)
    behavField = behavFields{1,dd};
    for ee = 1:length(dataTypes)
        dataType = dataTypes{1,ee};
        for ff = 1:size(data.PowerSpec.(behavField).(dataType).cat_S,2)
            if strcmp(behavField,'Rest') == true
                f_short = data.PowerSpec.(behavField).(dataType).cat_f(ff,:);
                S = data.PowerSpec.(behavField).(dataType).cat_S(:,ff);
                f_long = 0:0.01:0.5;
                S_long = interp1(f_short,S,f_long);
                index01 = find(f_long == 0.1);
                data.PowerSpec.(behavField).(dataType).S01(ff,1) = S_long(index01); %#ok<FNDSB>
            elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
                F = round(data.PowerSpec.(behavField).(dataType).cat_f(ff,:),2);
                S = data.PowerSpec.(behavField).(dataType).cat_S(:,ff);
                index01 = find(F == 0.1);
                data.PowerSpec.(behavField).(dataType).S01(ff,1) = S(index01(1));
            else
                F = round(data.PowerSpec.(behavField).(dataType).cat_f(ff,:),3);
                S = data.PowerSpec.(behavField).(dataType).cat_S(:,ff);
                index01 = find(F == 0.1);
                index001 = find(F == 0.01);
                data.PowerSpec.(behavField).(dataType).S01(ff,1) = S(index01(1));
                data.PowerSpec.(behavField).(dataType).S001(ff,1) = S(index001(1));
            end
        end
    end
end
% take mean/StD of peak S
for dd = 1:length(behavFields)
    behavField = behavFields{1,dd};
    for ee = 1:length(dataTypes)
        dataType = dataTypes{1,ee};
        if strcmp(behavField,'Rest') == true || strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
            data.PowerSpec.(behavField).(dataType).meanS01 = mean(data.PowerSpec.(behavField).(dataType).S01,1);
            data.PowerSpec.(behavField).(dataType).stdS01 = std(data.PowerSpec.(behavField).(dataType).S01,0,1);
        else
            data.PowerSpec.(behavField).(dataType).meanS01 = mean(data.PowerSpec.(behavField).(dataType).S01,1);
            data.PowerSpec.(behavField).(dataType).stdS01 = std(data.PowerSpec.(behavField).(dataType).S01,0,1);
            data.PowerSpec.(behavField).(dataType).meanS001 = mean(data.PowerSpec.(behavField).(dataType).S001,1);
            data.PowerSpec.(behavField).(dataType).stdS001 = std(data.PowerSpec.(behavField).(dataType).S001,0,1);
        end
    end
end
%% statistics - generalized linear mixed effects model
% delta PSD @ 0.1 Hz
Delta_PSD01_tableSize = cat(1,data.PowerSpec.Rest.deltaBandPower.S01,data.PowerSpec.NREM.deltaBandPower.S01,data.PowerSpec.REM.deltaBandPower.S01,...
    data.PowerSpec.Awake.deltaBandPower.S01,data.PowerSpec.Sleep.deltaBandPower.S01,data.PowerSpec.All.deltaBandPower.S01);
Delta_PSD01_Table = table('Size',[size(Delta_PSD01_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','PSD01'});
Delta_PSD01_Table.Mouse = cat(1,data.PowerSpec.Rest.deltaBandPower.animalID,data.PowerSpec.NREM.deltaBandPower.animalID,data.PowerSpec.REM.deltaBandPower.animalID,...
    data.PowerSpec.Awake.deltaBandPower.animalID,data.PowerSpec.Sleep.deltaBandPower.animalID,data.PowerSpec.All.deltaBandPower.animalID);
Delta_PSD01_Table.Hemisphere = cat(1,data.PowerSpec.Rest.deltaBandPower.hemisphere,data.PowerSpec.NREM.deltaBandPower.hemisphere,data.PowerSpec.REM.deltaBandPower.hemisphere,...
    data.PowerSpec.Awake.deltaBandPower.hemisphere,data.PowerSpec.Sleep.deltaBandPower.hemisphere,data.PowerSpec.All.deltaBandPower.hemisphere);
Delta_PSD01_Table.Behavior = cat(1,data.PowerSpec.Rest.deltaBandPower.behavior,data.PowerSpec.NREM.deltaBandPower.behavior,data.PowerSpec.REM.deltaBandPower.behavior,...
    data.PowerSpec.Awake.deltaBandPower.behavior,data.PowerSpec.Sleep.deltaBandPower.behavior,data.PowerSpec.All.deltaBandPower.behavior);
Delta_PSD01_Table.PSD01 = cat(1,data.PowerSpec.Rest.deltaBandPower.S01,data.PowerSpec.NREM.deltaBandPower.S01,data.PowerSpec.REM.deltaBandPower.S01,...
    data.PowerSpec.Awake.deltaBandPower.S01,data.PowerSpec.Sleep.deltaBandPower.S01,data.PowerSpec.All.deltaBandPower.S01);
Delta_PSD01_FitFormula = 'PSD01 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
Delta_PSD01_Stats = fitglme(Delta_PSD01_Table,Delta_PSD01_FitFormula);
% delta PSD @ 0.01 Hz
Delta_PSD001_tableSize = cat(1,data.PowerSpec.Awake.deltaBandPower.S001,data.PowerSpec.Sleep.deltaBandPower.S001,data.PowerSpec.All.deltaBandPower.S001);
Delta_PSD001_Table = table('Size',[size(Delta_PSD001_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','PSD001'});
Delta_PSD001_Table.Mouse = cat(1,data.PowerSpec.Awake.deltaBandPower.animalID,data.PowerSpec.Sleep.deltaBandPower.animalID,data.PowerSpec.All.deltaBandPower.animalID);
Delta_PSD001_Table.Hemisphere = cat(1,data.PowerSpec.Awake.deltaBandPower.hemisphere,data.PowerSpec.Sleep.deltaBandPower.hemisphere,data.PowerSpec.All.deltaBandPower.hemisphere);
Delta_PSD001_Table.Behavior = cat(1,data.PowerSpec.Awake.deltaBandPower.behavior,data.PowerSpec.Sleep.deltaBandPower.behavior,data.PowerSpec.All.deltaBandPower.behavior);
Delta_PSD001_Table.PSD001 = cat(1,data.PowerSpec.Awake.deltaBandPower.S001,data.PowerSpec.Sleep.deltaBandPower.S001,data.PowerSpec.All.deltaBandPower.S001);
Delta_PSD001_FitFormula = 'PSD001 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
Delta_PSD001_Stats = fitglme(Delta_PSD001_Table,Delta_PSD001_FitFormula);
% theta PSD @ 0.1 Hz
Theta_PSD01_tableSize = cat(1,data.PowerSpec.Rest.thetaBandPower.S01,data.PowerSpec.NREM.thetaBandPower.S01,data.PowerSpec.REM.thetaBandPower.S01,...
    data.PowerSpec.Awake.thetaBandPower.S01,data.PowerSpec.Sleep.thetaBandPower.S01,data.PowerSpec.All.thetaBandPower.S01);
Theta_PSD01_Table = table('Size',[size(Theta_PSD01_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','PSD01'});
Theta_PSD01_Table.Mouse = cat(1,data.PowerSpec.Rest.thetaBandPower.animalID,data.PowerSpec.NREM.thetaBandPower.animalID,data.PowerSpec.REM.thetaBandPower.animalID,...
    data.PowerSpec.Awake.thetaBandPower.animalID,data.PowerSpec.Sleep.thetaBandPower.animalID,data.PowerSpec.All.thetaBandPower.animalID);
Theta_PSD01_Table.Hemisphere = cat(1,data.PowerSpec.Rest.thetaBandPower.hemisphere,data.PowerSpec.NREM.thetaBandPower.hemisphere,data.PowerSpec.REM.thetaBandPower.hemisphere,...
    data.PowerSpec.Awake.thetaBandPower.hemisphere,data.PowerSpec.Sleep.thetaBandPower.hemisphere,data.PowerSpec.All.thetaBandPower.hemisphere);
Theta_PSD01_Table.Behavior = cat(1,data.PowerSpec.Rest.thetaBandPower.behavior,data.PowerSpec.NREM.thetaBandPower.behavior,data.PowerSpec.REM.thetaBandPower.behavior,...
    data.PowerSpec.Awake.thetaBandPower.behavior,data.PowerSpec.Sleep.thetaBandPower.behavior,data.PowerSpec.All.thetaBandPower.behavior);
Theta_PSD01_Table.PSD01 = cat(1,data.PowerSpec.Rest.thetaBandPower.S01,data.PowerSpec.NREM.thetaBandPower.S01,data.PowerSpec.REM.thetaBandPower.S01,...
    data.PowerSpec.Awake.thetaBandPower.S01,data.PowerSpec.Sleep.thetaBandPower.S01,data.PowerSpec.All.thetaBandPower.S01);
Theta_PSD01_FitFormula = 'PSD01 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
Theta_PSD01_Stats = fitglme(Theta_PSD01_Table,Theta_PSD01_FitFormula);
% theta PSD @ 0.01 Hz
Theta_PSD001_tableSize = cat(1,data.PowerSpec.Awake.thetaBandPower.S001,data.PowerSpec.Sleep.thetaBandPower.S001,data.PowerSpec.All.thetaBandPower.S001);
Theta_PSD001_Table = table('Size',[size(Theta_PSD001_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','PSD001'});
Theta_PSD001_Table.Mouse = cat(1,data.PowerSpec.Awake.thetaBandPower.animalID,data.PowerSpec.Sleep.thetaBandPower.animalID,data.PowerSpec.All.thetaBandPower.animalID);
Theta_PSD001_Table.Hemisphere = cat(1,data.PowerSpec.Awake.thetaBandPower.hemisphere,data.PowerSpec.Sleep.thetaBandPower.hemisphere,data.PowerSpec.All.thetaBandPower.hemisphere);
Theta_PSD001_Table.Behavior = cat(1,data.PowerSpec.Awake.thetaBandPower.behavior,data.PowerSpec.Sleep.thetaBandPower.behavior,data.PowerSpec.All.thetaBandPower.behavior);
Theta_PSD001_Table.PSD001 = cat(1,data.PowerSpec.Awake.thetaBandPower.S001,data.PowerSpec.Sleep.thetaBandPower.S001,data.PowerSpec.All.thetaBandPower.S001);
Theta_PSD001_FitFormula = 'PSD001 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
Theta_PSD001_Stats = fitglme(Theta_PSD001_Table,Theta_PSD001_FitFormula);
% alpha PSD @ 0.1 Hz
Alpha_PSD01_tableSize = cat(1,data.PowerSpec.Rest.alphaBandPower.S01,data.PowerSpec.NREM.alphaBandPower.S01,data.PowerSpec.REM.alphaBandPower.S01,...
    data.PowerSpec.Awake.alphaBandPower.S01,data.PowerSpec.Sleep.alphaBandPower.S01,data.PowerSpec.All.alphaBandPower.S01);
Alpha_PSD01_Table = table('Size',[size(Alpha_PSD01_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','PSD01'});
Alpha_PSD01_Table.Mouse = cat(1,data.PowerSpec.Rest.alphaBandPower.animalID,data.PowerSpec.NREM.alphaBandPower.animalID,data.PowerSpec.REM.alphaBandPower.animalID,...
    data.PowerSpec.Awake.alphaBandPower.animalID,data.PowerSpec.Sleep.alphaBandPower.animalID,data.PowerSpec.All.alphaBandPower.animalID);
Alpha_PSD01_Table.Hemisphere = cat(1,data.PowerSpec.Rest.alphaBandPower.hemisphere,data.PowerSpec.NREM.alphaBandPower.hemisphere,data.PowerSpec.REM.alphaBandPower.hemisphere,...
    data.PowerSpec.Awake.alphaBandPower.hemisphere,data.PowerSpec.Sleep.alphaBandPower.hemisphere,data.PowerSpec.All.alphaBandPower.hemisphere);
Alpha_PSD01_Table.Behavior = cat(1,data.PowerSpec.Rest.alphaBandPower.behavior,data.PowerSpec.NREM.alphaBandPower.behavior,data.PowerSpec.REM.alphaBandPower.behavior,...
    data.PowerSpec.Awake.alphaBandPower.behavior,data.PowerSpec.Sleep.alphaBandPower.behavior,data.PowerSpec.All.alphaBandPower.behavior);
Alpha_PSD01_Table.PSD01 = cat(1,data.PowerSpec.Rest.alphaBandPower.S01,data.PowerSpec.NREM.alphaBandPower.S01,data.PowerSpec.REM.alphaBandPower.S01,...
    data.PowerSpec.Awake.alphaBandPower.S01,data.PowerSpec.Sleep.alphaBandPower.S01,data.PowerSpec.All.alphaBandPower.S01);
Alpha_PSD01_FitFormula = 'PSD01 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
Alpha_PSD01_Stats = fitglme(Alpha_PSD01_Table,Alpha_PSD01_FitFormula);
% alpha PSD @ 0.01 Hz
Alpha_PSD001_tableSize = cat(1,data.PowerSpec.Awake.alphaBandPower.S001,data.PowerSpec.Sleep.alphaBandPower.S001,data.PowerSpec.All.alphaBandPower.S001);
Alpha_PSD001_Table = table('Size',[size(Alpha_PSD001_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','PSD001'});
Alpha_PSD001_Table.Mouse = cat(1,data.PowerSpec.Awake.alphaBandPower.animalID,data.PowerSpec.Sleep.alphaBandPower.animalID,data.PowerSpec.All.alphaBandPower.animalID);
Alpha_PSD001_Table.Hemisphere = cat(1,data.PowerSpec.Awake.alphaBandPower.hemisphere,data.PowerSpec.Sleep.alphaBandPower.hemisphere,data.PowerSpec.All.alphaBandPower.hemisphere);
Alpha_PSD001_Table.Behavior = cat(1,data.PowerSpec.Awake.alphaBandPower.behavior,data.PowerSpec.Sleep.alphaBandPower.behavior,data.PowerSpec.All.alphaBandPower.behavior);
Alpha_PSD001_Table.PSD001 = cat(1,data.PowerSpec.Awake.alphaBandPower.S001,data.PowerSpec.Sleep.alphaBandPower.S001,data.PowerSpec.All.alphaBandPower.S001);
Alpha_PSD001_FitFormula = 'PSD001 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
Alpha_PSD001_Stats = fitglme(Alpha_PSD001_Table,Alpha_PSD001_FitFormula);
% beta PSD @ 0.1 Hz
Beta_PSD01_tableSize = cat(1,data.PowerSpec.Rest.betaBandPower.S01,data.PowerSpec.NREM.betaBandPower.S01,data.PowerSpec.REM.betaBandPower.S01,...
    data.PowerSpec.Awake.betaBandPower.S01,data.PowerSpec.Sleep.betaBandPower.S01,data.PowerSpec.All.betaBandPower.S01);
Beta_PSD01_Table = table('Size',[size(Beta_PSD01_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','PSD01'});
Beta_PSD01_Table.Mouse = cat(1,data.PowerSpec.Rest.betaBandPower.animalID,data.PowerSpec.NREM.betaBandPower.animalID,data.PowerSpec.REM.betaBandPower.animalID,...
    data.PowerSpec.Awake.betaBandPower.animalID,data.PowerSpec.Sleep.betaBandPower.animalID,data.PowerSpec.All.betaBandPower.animalID);
Beta_PSD01_Table.Hemisphere = cat(1,data.PowerSpec.Rest.betaBandPower.hemisphere,data.PowerSpec.NREM.betaBandPower.hemisphere,data.PowerSpec.REM.betaBandPower.hemisphere,...
    data.PowerSpec.Awake.betaBandPower.hemisphere,data.PowerSpec.Sleep.betaBandPower.hemisphere,data.PowerSpec.All.betaBandPower.hemisphere);
Beta_PSD01_Table.Behavior = cat(1,data.PowerSpec.Rest.betaBandPower.behavior,data.PowerSpec.NREM.betaBandPower.behavior,data.PowerSpec.REM.betaBandPower.behavior,...
    data.PowerSpec.Awake.betaBandPower.behavior,data.PowerSpec.Sleep.betaBandPower.behavior,data.PowerSpec.All.betaBandPower.behavior);
Beta_PSD01_Table.PSD01 = cat(1,data.PowerSpec.Rest.betaBandPower.S01,data.PowerSpec.NREM.betaBandPower.S01,data.PowerSpec.REM.betaBandPower.S01,...
    data.PowerSpec.Awake.betaBandPower.S01,data.PowerSpec.Sleep.betaBandPower.S01,data.PowerSpec.All.betaBandPower.S01);
Beta_PSD01_FitFormula = 'PSD01 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
Beta_PSD01_Stats = fitglme(Beta_PSD01_Table,Beta_PSD01_FitFormula);
% beta PSD @ 0.01 Hz
Beta_PSD001_tableSize = cat(1,data.PowerSpec.Awake.betaBandPower.S001,data.PowerSpec.Sleep.betaBandPower.S001,data.PowerSpec.All.betaBandPower.S001);
Beta_PSD001_Table = table('Size',[size(Beta_PSD001_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','PSD001'});
Beta_PSD001_Table.Mouse = cat(1,data.PowerSpec.Awake.betaBandPower.animalID,data.PowerSpec.Sleep.betaBandPower.animalID,data.PowerSpec.All.betaBandPower.animalID);
Beta_PSD001_Table.Hemisphere = cat(1,data.PowerSpec.Awake.betaBandPower.hemisphere,data.PowerSpec.Sleep.betaBandPower.hemisphere,data.PowerSpec.All.betaBandPower.hemisphere);
Beta_PSD001_Table.Behavior = cat(1,data.PowerSpec.Awake.betaBandPower.behavior,data.PowerSpec.Sleep.betaBandPower.behavior,data.PowerSpec.All.betaBandPower.behavior);
Beta_PSD001_Table.PSD001 = cat(1,data.PowerSpec.Awake.betaBandPower.S001,data.PowerSpec.Sleep.betaBandPower.S001,data.PowerSpec.All.betaBandPower.S001);
Beta_PSD001_FitFormula = 'PSD001 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
Beta_PSD001_Stats = fitglme(Beta_PSD001_Table,Beta_PSD001_FitFormula);
%% Fig. 7-S3
summaryFigure = figure('Name','Fig7-S3 (a-p');
sgtitle('Figure 7-S3 - Turner et al. 2020')
%% [7-S3a] delta PSD @ 0.1 Hz
ax1 = subplot(4,4,1);
s1 = scatter(ones(1,length(data.PowerSpec.Rest.deltaBandPower.S01))*1,data.PowerSpec.Rest.deltaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.PowerSpec.Rest.deltaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(ones(1,length(data.PowerSpec.NREM.deltaBandPower.S01))*2,data.PowerSpec.NREM.deltaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.PowerSpec.NREM.deltaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(ones(1,length(data.PowerSpec.REM.deltaBandPower.S01))*3,data.PowerSpec.REM.deltaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.PowerSpec.REM.deltaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
s4 = scatter(ones(1,length(data.PowerSpec.Awake.deltaBandPower.S01))*4,data.PowerSpec.Awake.deltaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.PowerSpec.Awake.deltaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
s5 = scatter(ones(1,length(data.PowerSpec.Sleep.deltaBandPower.S01))*5,data.PowerSpec.Sleep.deltaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.PowerSpec.Sleep.deltaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
s6 = scatter(ones(1,length(data.PowerSpec.All.deltaBandPower.S01))*6,data.PowerSpec.All.deltaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.PowerSpec.All.deltaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
title({'[7-S3a] PSD @ 0.1 Hz','Delta-band [1-4 Hz]'})
ylabel('Power (a.u.)')
legend([s1,s2,s3,s4,s5,s6],'Rest','NREM','REM','Alert','Asleep','All')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
axis square
xlim([0,7])
ylim([0.1,100])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [7-S3b] delta PSD @ 0.01 Hz
ax2 = subplot(4,4,2);
scatter(ones(1,length(data.PowerSpec.Awake.deltaBandPower.S001))*1,data.PowerSpec.Awake.deltaBandPower.S001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.PowerSpec.Awake.deltaBandPower.meanS001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Sleep.deltaBandPower.S001))*2,data.PowerSpec.Sleep.deltaBandPower.S001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.PowerSpec.Sleep.deltaBandPower.meanS001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.All.deltaBandPower.S001))*3,data.PowerSpec.All.deltaBandPower.S001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.PowerSpec.All.deltaBandPower.meanS001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[7-S3b] PSD @ 0.01 Hz','Delta-band [1-4 Hz]'})
ylabel('Power (a.u.)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
axis square
xlim([0,4])
ylim([1,1000])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [7-S3c] delta coherence^2 @ 0.1 Hz
ax3 = subplot(4,4,3);
scatter(ones(1,length(data.Coherr.Rest.deltaBandPower.C01))*1,data.Coherr.Rest.deltaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Coherr.Rest.deltaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.Coherr.NREM.deltaBandPower.C01))*2,data.Coherr.NREM.deltaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Coherr.NREM.deltaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.Coherr.REM.deltaBandPower.C01))*3,data.Coherr.REM.deltaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.Coherr.REM.deltaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.Coherr.Awake.deltaBandPower.C01))*4,data.Coherr.Awake.deltaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.Coherr.Awake.deltaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.Coherr.Sleep.deltaBandPower.C01))*5,data.Coherr.Sleep.deltaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.Coherr.Sleep.deltaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(ones(1,length(data.Coherr.All.deltaBandPower.C01))*6,data.Coherr.All.deltaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.Coherr.All.deltaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
title({'[7-S3c] Coherence^2 @ 0.1 Hz','Delta-band [1-4 Hz]'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,7])
ylim([0,1])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [7-S3d] delta coherence^2 @ 0.01 Hz
ax4 = subplot(4,4,4);
scatter(ones(1,length(data.Coherr.Awake.deltaBandPower.C001))*1,data.Coherr.Awake.deltaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Coherr.Awake.deltaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.Coherr.Sleep.deltaBandPower.C001))*2,data.Coherr.Sleep.deltaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Coherr.Sleep.deltaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.Coherr.All.deltaBandPower.C001))*3,data.Coherr.All.deltaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.Coherr.All.deltaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[7-S3d] Coherence^2 @ 0.01 Hz','Delta-band [1-4 Hz]'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,4])
ylim([0,1])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [7-S3e] theta PSD @ 0.1 Hz
ax5 = subplot(4,4,5);
scatter(ones(1,length(data.PowerSpec.Rest.thetaBandPower.S01))*1,data.PowerSpec.Rest.thetaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.PowerSpec.Rest.thetaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.NREM.thetaBandPower.S01))*2,data.PowerSpec.NREM.thetaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.PowerSpec.NREM.thetaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.REM.thetaBandPower.S01))*3,data.PowerSpec.REM.thetaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.PowerSpec.REM.thetaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Awake.thetaBandPower.S01))*4,data.PowerSpec.Awake.thetaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.PowerSpec.Awake.thetaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Sleep.thetaBandPower.S01))*5,data.PowerSpec.Sleep.thetaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.PowerSpec.Sleep.thetaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.All.thetaBandPower.S01))*6,data.PowerSpec.All.thetaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.PowerSpec.All.thetaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
title({'[7-S3e] PSD @ 0.1 Hz','Theta-band [4-10 Hz]'})
ylabel('Power (a.u.)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
axis square
xlim([0,7])
ylim([0.1,1000])
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [S17f] theta PSD @ 0.01 Hz
ax6 = subplot(4,4,6);
scatter(ones(1,length(data.PowerSpec.Awake.thetaBandPower.S001))*1,data.PowerSpec.Awake.thetaBandPower.S001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.PowerSpec.Awake.thetaBandPower.meanS001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Sleep.thetaBandPower.S001))*2,data.PowerSpec.Sleep.thetaBandPower.S001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.PowerSpec.Sleep.thetaBandPower.meanS001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.All.thetaBandPower.S001))*3,data.PowerSpec.All.thetaBandPower.S001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.PowerSpec.All.thetaBandPower.meanS001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[7-S3f] PSD @ 0.01 Hz','Theta-band [4-10 Hz]'})
ylabel('Power (a.u.)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
axis square
xlim([0,4])
ylim([1,1000])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% [S17g] theta coherence^2 @ 0.1 Hz
ax7 = subplot(4,4,7);
scatter(ones(1,length(data.Coherr.Rest.thetaBandPower.C01))*1,data.Coherr.Rest.thetaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Coherr.Rest.thetaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.Coherr.NREM.thetaBandPower.C01))*2,data.Coherr.NREM.thetaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Coherr.NREM.thetaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.Coherr.REM.thetaBandPower.C01))*3,data.Coherr.REM.thetaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.Coherr.REM.thetaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.Coherr.Awake.thetaBandPower.C01))*4,data.Coherr.Awake.thetaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.Coherr.Awake.thetaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.Coherr.Sleep.thetaBandPower.C01))*5,data.Coherr.Sleep.thetaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.Coherr.Sleep.thetaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(ones(1,length(data.Coherr.All.thetaBandPower.C01))*6,data.Coherr.All.thetaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.Coherr.All.thetaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
title({'[7-S3g] Coherence^2 @ 0.1 Hz','Theta-band [4-10 Hz]'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,7])
ylim([0,1])
set(gca,'box','off')
ax7.TickLength = [0.03,0.03];
%% [7-S3h] theta coherence^2 @ 0.01 Hz
ax8 = subplot(4,4,8);
scatter(ones(1,length(data.Coherr.Awake.thetaBandPower.C001))*1,data.Coherr.Awake.thetaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Coherr.Awake.thetaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.Coherr.Sleep.thetaBandPower.C001))*2,data.Coherr.Sleep.thetaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Coherr.Sleep.thetaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.Coherr.All.thetaBandPower.C001))*3,data.Coherr.All.thetaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.Coherr.All.thetaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[7-S3h] Coherence^2 @ 0.01 Hz','Theta-band [4-10 Hz]'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,4])
ylim([0,1])
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];
%% [7-S3i] alpha PSD @ 0.1 Hz
ax9 = subplot(4,4,9);
scatter(ones(1,length(data.PowerSpec.Rest.alphaBandPower.S01))*1,data.PowerSpec.Rest.alphaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.PowerSpec.Rest.alphaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.NREM.alphaBandPower.S01))*2,data.PowerSpec.NREM.alphaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.PowerSpec.NREM.alphaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.REM.alphaBandPower.S01))*3,data.PowerSpec.REM.alphaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.PowerSpec.REM.alphaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Awake.alphaBandPower.S01))*4,data.PowerSpec.Awake.alphaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.PowerSpec.Awake.alphaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Sleep.alphaBandPower.S01))*5,data.PowerSpec.Sleep.alphaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.PowerSpec.Sleep.alphaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.All.alphaBandPower.S01))*6,data.PowerSpec.All.alphaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.PowerSpec.All.alphaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
title({'[7-S3i] PSD @ 0.1 Hz','Alpha-band [10-13 Hz]'})
ylabel('Power (a.u.)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
axis square
xlim([0,7])
ylim([0.1,10000])
set(gca,'box','off')
ax9.TickLength = [0.03,0.03];
%% [7-S3j] alpha PSD @ 0.01 Hz
ax10 = subplot(4,4,10);
scatter(ones(1,length(data.PowerSpec.Awake.alphaBandPower.S001))*1,data.PowerSpec.Awake.alphaBandPower.S001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.PowerSpec.Awake.alphaBandPower.meanS001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Sleep.alphaBandPower.S001))*2,data.PowerSpec.Sleep.alphaBandPower.S001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.PowerSpec.Sleep.alphaBandPower.meanS001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.All.alphaBandPower.S001))*3,data.PowerSpec.All.alphaBandPower.S001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.PowerSpec.All.alphaBandPower.meanS001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[7-S3j] PSD @ 0.01 Hz','Alpha-band [10-13 Hz]'})
ylabel('Power (a.u.)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
axis square
xlim([0,4])
ylim([1,1000])
set(gca,'box','off')
ax10.TickLength = [0.03,0.03];
%% [7-S3k] alpha coherence^2 @ 0.1 Hz
ax11 = subplot(4,4,11);
scatter(ones(1,length(data.Coherr.Rest.alphaBandPower.C01))*1,data.Coherr.Rest.alphaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Coherr.Rest.alphaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.Coherr.NREM.alphaBandPower.C01))*2,data.Coherr.NREM.alphaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Coherr.NREM.alphaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.Coherr.REM.alphaBandPower.C01))*3,data.Coherr.REM.alphaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.Coherr.REM.alphaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.Coherr.Awake.alphaBandPower.C01))*4,data.Coherr.Awake.alphaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.Coherr.Awake.alphaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.Coherr.Sleep.alphaBandPower.C01))*5,data.Coherr.Sleep.alphaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.Coherr.Sleep.alphaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(ones(1,length(data.Coherr.All.alphaBandPower.C01))*6,data.Coherr.All.alphaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.Coherr.All.alphaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
title({'[7-S3k] Coherence^2 @ 0.1 Hz','Alpha-band [10-13 Hz]'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,7])
ylim([0,1])
set(gca,'box','off')
ax11.TickLength = [0.03,0.03];
%% [7-S3l] alpha coherence^2 @ 0.01 Hz
ax12 = subplot(4,4,12);
scatter(ones(1,length(data.Coherr.Awake.alphaBandPower.C001))*1,data.Coherr.Awake.alphaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Coherr.Awake.alphaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.Coherr.Sleep.alphaBandPower.C001))*2,data.Coherr.Sleep.alphaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Coherr.Sleep.alphaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.Coherr.All.alphaBandPower.C001))*3,data.Coherr.All.alphaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.Coherr.All.alphaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[7-S3l] Coherence^2 @ 0.01 Hz','Alpha-band [10-13 Hz]'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,4])
ylim([0,1])
set(gca,'box','off')
ax12.TickLength = [0.03,0.03];
%% [7-S3m] beta PSD @ 0.1 Hz
ax13 = subplot(4,4,13);
scatter(ones(1,length(data.PowerSpec.Rest.betaBandPower.S01))*1,data.PowerSpec.Rest.betaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.PowerSpec.Rest.betaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.NREM.betaBandPower.S01))*2,data.PowerSpec.NREM.betaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.PowerSpec.NREM.betaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.REM.betaBandPower.S01))*3,data.PowerSpec.REM.betaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.PowerSpec.REM.betaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Awake.betaBandPower.S01))*4,data.PowerSpec.Awake.betaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.PowerSpec.Awake.betaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Sleep.betaBandPower.S01))*5,data.PowerSpec.Sleep.betaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.PowerSpec.Sleep.betaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.All.betaBandPower.S01))*6,data.PowerSpec.All.betaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.PowerSpec.All.betaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
title({'[7-S3m] PSD @ 0.1 Hz','Beta-band [13-30 Hz]'})
ylabel('Power (a.u.)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
axis square
xlim([0,7])
ylim([0.1,10000])
set(gca,'box','off')
ax13.TickLength = [0.03,0.03];
%% [7-S3n] beta PSD @ 0.01 Hz
ax14 = subplot(4,4,14);
scatter(ones(1,length(data.PowerSpec.Awake.betaBandPower.S001))*1,data.PowerSpec.Awake.betaBandPower.S001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.PowerSpec.Awake.betaBandPower.meanS001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Sleep.betaBandPower.S001))*2,data.PowerSpec.Sleep.betaBandPower.S001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.PowerSpec.Sleep.betaBandPower.meanS001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.All.betaBandPower.S001))*3,data.PowerSpec.All.betaBandPower.S001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.PowerSpec.All.betaBandPower.meanS001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[7-S3n] PSD @ 0.01 Hz','Beta-band [13-30 Hz]'})
ylabel('Power (a.u.)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
axis square
xlim([0,4])
ylim([1,1000])
set(gca,'box','off')
ax14.TickLength = [0.03,0.03];
%% [7-S3o] delta coherence^2 @ 0.1 Hz
ax15 = subplot(4,4,15);
scatter(ones(1,length(data.Coherr.Rest.betaBandPower.C01))*1,data.Coherr.Rest.betaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Coherr.Rest.betaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.Coherr.NREM.betaBandPower.C01))*2,data.Coherr.NREM.betaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Coherr.NREM.betaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.Coherr.REM.betaBandPower.C01))*3,data.Coherr.REM.betaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.Coherr.REM.betaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.Coherr.Awake.betaBandPower.C01))*4,data.Coherr.Awake.betaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.Coherr.Awake.betaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.Coherr.Sleep.betaBandPower.C01))*5,data.Coherr.Sleep.betaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.Coherr.Sleep.betaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(ones(1,length(data.Coherr.All.betaBandPower.C01))*6,data.Coherr.All.betaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.Coherr.All.betaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
title({'[7-S3o] Coherence^2 @ 0.1 Hz','Beta-band [13-30 Hz]'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,7])
ylim([0,1])
set(gca,'box','off')
ax15.TickLength = [0.03,0.03];
%% [7-S3p] beta coherence^2 @ 0.01 Hz
ax16 = subplot(4,4,16);
scatter(ones(1,length(data.Coherr.Awake.betaBandPower.C001))*1,data.Coherr.Awake.betaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Coherr.Awake.betaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.Coherr.Sleep.betaBandPower.C001))*2,data.Coherr.Sleep.betaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Coherr.Sleep.betaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.Coherr.All.betaBandPower.C001))*3,data.Coherr.All.betaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.Coherr.All.betaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[7-S3p] Coherence^2 @ 0.01 Hz','Beta-band [13-30 Hz]'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,4])
ylim([0,1])
set(gca,'box','off')
ax16.TickLength = [0.03,0.03];
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig7-S3']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig7-S3'])
    %% statistical diary
    diaryFile = [dirpath 'Fig7-S3_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    % delta-band 0.1 Hz PSD statistical diary
    disp('======================================================================================================================')
    disp('[7-S3a] Generalized linear mixed-effects model statistics for delta-band PSD @ 0.1 Hz for Rest, NREM, REM, Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Delta_PSD01_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  Delta 0.1 Hz PSD: ' num2str(round(data.PowerSpec.Rest.deltaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.Rest.deltaBandPower.stdS01,1))]); disp(' ')
    disp(['NREM  Delta 0.1 Hz PSD: ' num2str(round(data.PowerSpec.NREM.deltaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.NREM.deltaBandPower.stdS01,1))]); disp(' ')
    disp(['REM   Delta 0.1 Hz PSD: ' num2str(round(data.PowerSpec.REM.deltaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.REM.deltaBandPower.stdS01,1))]); disp(' ')
    disp(['Awake Delta 0.1 Hz PSD: ' num2str(round(data.PowerSpec.Awake.deltaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.Awake.deltaBandPower.stdS01,1))]); disp(' ')
    disp(['Sleep Delta 0.1 Hz PSD: ' num2str(round(data.PowerSpec.Sleep.deltaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.Sleep.deltaBandPower.stdS01,1))]); disp(' ')
    disp(['All   Delta 0.1 Hz PSD: ' num2str(round(data.PowerSpec.All.deltaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.All.deltaBandPower.stdS01,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % delta-band 0.01 Hz PSD statistical diary
    disp('======================================================================================================================')
    disp('[7-S3b] Generalized linear mixed-effects model statistics for delta-band PSD @ 0.01 Hz for Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Delta_PSD001_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake Delta 0.01 Hz PSD: ' num2str(round(data.PowerSpec.Awake.deltaBandPower.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.Awake.deltaBandPower.stdS001,1))]); disp(' ')
    disp(['Sleep Delta 0.01 Hz PSD: ' num2str(round(data.PowerSpec.Sleep.deltaBandPower.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.Sleep.deltaBandPower.stdS001,1))]); disp(' ')
    disp(['All   Delta 0.01 Hz PSD: ' num2str(round(data.PowerSpec.All.deltaBandPower.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.All.deltaBandPower.stdS001,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % delta-band 0.1 Hz coherence^2 statistical diary
    disp('======================================================================================================================')
    disp('[7-S3c] Generalized linear mixed-effects model statistics for delta-band coherence^2 @ 0.1 Hz for Rest, NREM, REM, Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Delta_Coh01_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  Delta 0.1 Hz Coh2: ' num2str(round(data.Coherr.Rest.deltaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.Rest.deltaBandPower.stdC01,2))]); disp(' ')
    disp(['NREM  Delta 0.1 Hz Coh2: ' num2str(round(data.Coherr.NREM.deltaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.NREM.deltaBandPower.stdC01,2))]); disp(' ')
    disp(['REM   Delta 0.1 Hz Coh2: ' num2str(round(data.Coherr.REM.deltaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.REM.deltaBandPower.stdC01,2))]); disp(' ')
    disp(['Awake Delta 0.1 Hz Coh2: ' num2str(round(data.Coherr.Awake.deltaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.Awake.deltaBandPower.stdC01,2))]); disp(' ')
    disp(['Sleep Delta 0.1 Hz Coh2: ' num2str(round(data.Coherr.Sleep.deltaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.Sleep.deltaBandPower.stdC01,2))]); disp(' ')
    disp(['All   Delta 0.1 Hz Coh2: ' num2str(round(data.Coherr.All.deltaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.All.deltaBandPower.stdC01,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % delta-band 0.01 Hz coherence^2 statistical diary
    disp('======================================================================================================================')
    disp('[7-S3d] Generalized linear mixed-effects model statistics for delta-band coherence^2 @ 0.01 Hz for Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Delta_Coh001_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake Delta 0.01 Hz Coh2: ' num2str(round(data.Coherr.Awake.deltaBandPower.meanC001,2)) ' +/- ' num2str(round(data.Coherr.Awake.deltaBandPower.stdC001,2))]); disp(' ')
    disp(['Sleep Delta 0.01 Hz Coh2: ' num2str(round(data.Coherr.Sleep.deltaBandPower.meanC001,2)) ' +/- ' num2str(round(data.Coherr.Sleep.deltaBandPower.stdC001,2))]); disp(' ')
    disp(['All   Delta 0.01 Hz Coh2: ' num2str(round(data.Coherr.All.deltaBandPower.meanC001,2)) ' +/- ' num2str(round(data.Coherr.All.deltaBandPower.stdC001,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % theta-band 0.1 Hz PSD statistical diary
    disp('======================================================================================================================')
    disp('[7-S3e] Generalized linear mixed-effects model statistics for theta-band PSD @ 0.1 Hz for Rest, NREM, REM, Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Theta_PSD01_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  Theta 0.1 Hz PSD: ' num2str(round(data.PowerSpec.Rest.thetaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.Rest.thetaBandPower.stdS01,1))]); disp(' ')
    disp(['NREM  Theta 0.1 Hz PSD: ' num2str(round(data.PowerSpec.NREM.thetaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.NREM.thetaBandPower.stdS01,1))]); disp(' ')
    disp(['REM   Theta 0.1 Hz PSD: ' num2str(round(data.PowerSpec.REM.thetaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.REM.thetaBandPower.stdS01,1))]); disp(' ')
    disp(['Awake Theta 0.1 Hz PSD: ' num2str(round(data.PowerSpec.Awake.thetaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.Awake.thetaBandPower.stdS01,1))]); disp(' ')
    disp(['Sleep Theta 0.1 Hz PSD: ' num2str(round(data.PowerSpec.Sleep.thetaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.Sleep.thetaBandPower.stdS01,1))]); disp(' ')
    disp(['All   Theta 0.1 Hz PSD: ' num2str(round(data.PowerSpec.All.thetaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.All.thetaBandPower.stdS01,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % theta-band 0.01 Hz PSD statistical diary
    disp('======================================================================================================================')
    disp('[7-S3f] Generalized linear mixed-effects model statistics for theta-band PSD @ 0.01 Hz for Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Theta_PSD001_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake Theta 0.01 Hz PSD: ' num2str(round(data.PowerSpec.Awake.thetaBandPower.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.Awake.thetaBandPower.stdS001,1))]); disp(' ')
    disp(['Sleep Theta 0.01 Hz PSD: ' num2str(round(data.PowerSpec.Sleep.thetaBandPower.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.Sleep.thetaBandPower.stdS001,1))]); disp(' ')
    disp(['All   Theta 0.01 Hz PSD: ' num2str(round(data.PowerSpec.All.thetaBandPower.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.All.thetaBandPower.stdS001,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % theta-band 0.1 Hz coherence^2 statistical diary
    disp('======================================================================================================================')
    disp('[7-S3g] Generalized linear mixed-effects model statistics for theta-band coherence^2 @ 0.1 Hz for Rest, NREM, REM, Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Theta_Coh01_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  Theta 0.1 Hz Coh2: ' num2str(round(data.Coherr.Rest.thetaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.Rest.thetaBandPower.stdC01,2))]); disp(' ')
    disp(['NREM  Theta 0.1 Hz Coh2: ' num2str(round(data.Coherr.NREM.thetaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.NREM.thetaBandPower.stdC01,2))]); disp(' ')
    disp(['REM   Theta 0.1 Hz Coh2: ' num2str(round(data.Coherr.REM.thetaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.REM.thetaBandPower.stdC01,2))]); disp(' ')
    disp(['Awake Theta 0.1 Hz Coh2: ' num2str(round(data.Coherr.Awake.thetaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.Awake.thetaBandPower.stdC01,2))]); disp(' ')
    disp(['Sleep Theta 0.1 Hz Coh2: ' num2str(round(data.Coherr.Sleep.thetaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.Sleep.thetaBandPower.stdC01,2))]); disp(' ')
    disp(['All   Theta 0.1 Hz Coh2: ' num2str(round(data.Coherr.All.thetaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.All.thetaBandPower.stdC01,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % theta-band 0.01 Hz coherence^2 statistical diary
    disp('======================================================================================================================')
    disp('[7-S3h] Generalized linear mixed-effects model statistics for theta-band coherence^2 @ 0.01 Hz for Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Theta_Coh001_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake Theta 0.01 Hz Coh2: ' num2str(round(data.Coherr.Awake.thetaBandPower.meanC001,2)) ' +/- ' num2str(round(data.Coherr.Awake.thetaBandPower.stdC001,2))]); disp(' ')
    disp(['Sleep Theta 0.01 Hz Coh2: ' num2str(round(data.Coherr.Sleep.thetaBandPower.meanC001,2)) ' +/- ' num2str(round(data.Coherr.Sleep.thetaBandPower.stdC001,2))]); disp(' ')
    disp(['All   Theta 0.01 Hz Coh2: ' num2str(round(data.Coherr.All.thetaBandPower.meanC001,2)) ' +/- ' num2str(round(data.Coherr.All.thetaBandPower.stdC001,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % alpha-band 0.1 Hz PSD statistical diary
    disp('======================================================================================================================')
    disp('[7-S3i] Generalized linear mixed-effects model statistics for alpha-band PSD @ 0.1 Hz for Rest, NREM, REM, Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Alpha_PSD01_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  Alpha 0.1 Hz PSD: ' num2str(round(data.PowerSpec.Rest.alphaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.Rest.alphaBandPower.stdS01,1))]); disp(' ')
    disp(['NREM  Alpha 0.1 Hz PSD: ' num2str(round(data.PowerSpec.NREM.alphaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.NREM.alphaBandPower.stdS01,1))]); disp(' ')
    disp(['REM   Alpha 0.1 Hz PSD: ' num2str(round(data.PowerSpec.REM.alphaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.REM.alphaBandPower.stdS01,1))]); disp(' ')
    disp(['Awake Alpha 0.1 Hz PSD: ' num2str(round(data.PowerSpec.Awake.alphaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.Awake.alphaBandPower.stdS01,1))]); disp(' ')
    disp(['Sleep Alpha 0.1 Hz PSD: ' num2str(round(data.PowerSpec.Sleep.alphaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.Sleep.alphaBandPower.stdS01,1))]); disp(' ')
    disp(['All   Alpha 0.1 Hz PSD: ' num2str(round(data.PowerSpec.All.alphaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.All.alphaBandPower.stdS01,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % alpha-band 0.01 Hz PSD statistical diary
    disp('======================================================================================================================')
    disp('[7-S3j] Generalized linear mixed-effects model statistics for alpha-band PSD @ 0.01 Hz for Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Alpha_PSD001_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake Alpha 0.01 Hz PSD: ' num2str(round(data.PowerSpec.Awake.alphaBandPower.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.Awake.alphaBandPower.stdS001,1))]); disp(' ')
    disp(['Sleep Alpha 0.01 Hz PSD: ' num2str(round(data.PowerSpec.Sleep.alphaBandPower.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.Sleep.alphaBandPower.stdS001,1))]); disp(' ')
    disp(['All   Alpha 0.01 Hz PSD: ' num2str(round(data.PowerSpec.All.alphaBandPower.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.All.alphaBandPower.stdS001,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % alpha-band 0.1 Hz coherence^2 statistical diary
    disp('======================================================================================================================')
    disp('[7-S3k] Generalized linear mixed-effects model statistics for alpha-band coherence^2 @ 0.1 Hz for Rest, NREM, REM, Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Alpha_Coh01_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  Alpha 0.1 Hz Coh2: ' num2str(round(data.Coherr.Rest.alphaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.Rest.alphaBandPower.stdC01,2))]); disp(' ')
    disp(['NREM  Alpha 0.1 Hz Coh2: ' num2str(round(data.Coherr.NREM.alphaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.NREM.alphaBandPower.stdC01,2))]); disp(' ')
    disp(['REM   Alpha 0.1 Hz Coh2: ' num2str(round(data.Coherr.REM.alphaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.REM.alphaBandPower.stdC01,2))]); disp(' ')
    disp(['Awake Alpha 0.1 Hz Coh2: ' num2str(round(data.Coherr.Awake.alphaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.Awake.alphaBandPower.stdC01,2))]); disp(' ')
    disp(['Sleep Alpha 0.1 Hz Coh2: ' num2str(round(data.Coherr.Sleep.alphaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.Sleep.alphaBandPower.stdC01,2))]); disp(' ')
    disp(['All   Alpha 0.1 Hz Coh2: ' num2str(round(data.Coherr.All.alphaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.All.alphaBandPower.stdC01,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % alpha-band 0.01 Hz coherence^2 statistical diary
    disp('======================================================================================================================')
    disp('[7-S3l] Generalized linear mixed-effects model statistics for alpha-band coherence^2 @ 0.01 Hz for Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Alpha_Coh001_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake Alpha 0.01 Hz Coh2: ' num2str(round(data.Coherr.Awake.alphaBandPower.meanC001,2)) ' +/- ' num2str(round(data.Coherr.Awake.alphaBandPower.stdC001,2))]); disp(' ')
    disp(['Sleep Alpha 0.01 Hz Coh2: ' num2str(round(data.Coherr.Sleep.alphaBandPower.meanC001,2)) ' +/- ' num2str(round(data.Coherr.Sleep.alphaBandPower.stdC001,2))]); disp(' ')
    disp(['All   Alpha 0.01 Hz Coh2: ' num2str(round(data.Coherr.All.alphaBandPower.meanC001,2)) ' +/- ' num2str(round(data.Coherr.All.alphaBandPower.stdC001,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % beta-band 0.1 Hz PSD statistical diary
    disp('======================================================================================================================')
    disp('[7-S3m] Generalized linear mixed-effects model statistics for beta-band PSD @ 0.1 Hz for Rest, NREM, REM, Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Beta_PSD01_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  Beta 0.1 Hz PSD: ' num2str(round(data.PowerSpec.Rest.betaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.Rest.betaBandPower.stdS01,1))]); disp(' ')
    disp(['NREM  Beta 0.1 Hz PSD: ' num2str(round(data.PowerSpec.NREM.betaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.NREM.betaBandPower.stdS01,1))]); disp(' ')
    disp(['REM   Beta 0.1 Hz PSD: ' num2str(round(data.PowerSpec.REM.betaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.REM.betaBandPower.stdS01,1))]); disp(' ')
    disp(['Awake Beta 0.1 Hz PSD: ' num2str(round(data.PowerSpec.Awake.betaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.Awake.betaBandPower.stdS01,1))]); disp(' ')
    disp(['Sleep Beta 0.1 Hz PSD: ' num2str(round(data.PowerSpec.Sleep.betaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.Sleep.betaBandPower.stdS01,1))]); disp(' ')
    disp(['All   Beta 0.1 Hz PSD: ' num2str(round(data.PowerSpec.All.betaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.All.betaBandPower.stdS01,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % beta-band 0.01 Hz PSD statistical diary
    disp('======================================================================================================================')
    disp('[7-S3n] Generalized linear mixed-effects model statistics for beta-band PSD @ 0.01 Hz for Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Beta_PSD001_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake Beta 0.01 Hz PSD: ' num2str(round(data.PowerSpec.Awake.betaBandPower.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.Awake.betaBandPower.stdS001,1))]); disp(' ')
    disp(['Sleep Beta 0.01 Hz PSD: ' num2str(round(data.PowerSpec.Sleep.betaBandPower.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.Sleep.betaBandPower.stdS001,1))]); disp(' ')
    disp(['All   Beta 0.01 Hz PSD: ' num2str(round(data.PowerSpec.All.betaBandPower.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.All.betaBandPower.stdS001,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % beta-band 0.1 Hz coherence^2 statistical diary
    disp('======================================================================================================================')
    disp('[7-S3o] Generalized linear mixed-effects model statistics for beta-band coherence^2 @ 0.1 Hz for Rest, NREM, REM, Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Beta_Coh01_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  Beta 0.1 Hz Coh2: ' num2str(round(data.Coherr.Rest.betaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.Rest.betaBandPower.stdC01,2))]); disp(' ')
    disp(['NREM  Beta 0.1 Hz Coh2: ' num2str(round(data.Coherr.NREM.betaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.NREM.betaBandPower.stdC01,2))]); disp(' ')
    disp(['REM   Beta 0.1 Hz Coh2: ' num2str(round(data.Coherr.REM.betaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.REM.betaBandPower.stdC01,2))]); disp(' ')
    disp(['Awake Beta 0.1 Hz Coh2: ' num2str(round(data.Coherr.Awake.betaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.Awake.betaBandPower.stdC01,2))]); disp(' ')
    disp(['Sleep Beta 0.1 Hz Coh2: ' num2str(round(data.Coherr.Sleep.betaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.Sleep.betaBandPower.stdC01,2))]); disp(' ')
    disp(['All   Beta 0.1 Hz Coh2: ' num2str(round(data.Coherr.All.betaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.All.betaBandPower.stdC01,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % beta-band 0.01 Hz coherence^2 statistical diary
    disp('======================================================================================================================')
    disp('[7-S3p] Generalized linear mixed-effects model statistics for beta-band coherence^2 @ 0.01 Hz for Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Beta_Coh001_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake Beta 0.01 Hz Coh2: ' num2str(round(data.Coherr.Awake.betaBandPower.meanC001,2)) ' +/- ' num2str(round(data.Coherr.Awake.betaBandPower.stdC001,2))]); disp(' ')
    disp(['Sleep Beta 0.01 Hz Coh2: ' num2str(round(data.Coherr.Sleep.betaBandPower.meanC001,2)) ' +/- ' num2str(round(data.Coherr.Sleep.betaBandPower.stdC001,2))]); disp(' ')
    disp(['All   Beta 0.01 Hz Coh2: ' num2str(round(data.Coherr.All.betaBandPower.meanC001,2)) ' +/- ' num2str(round(data.Coherr.All.betaBandPower.stdC001,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    diary off
    %% organized for supplemental table
    % variable names
    ColumnNames = {'Rest','NREM','REM','Awake','Sleep','All'};
    % delta-band 0.1 Hz PSD power
    for aa = 1:length(ColumnNames)
        Delta_S01_MeanStD{1,aa} = [num2str(round(data.PowerSpec.(ColumnNames{1,aa}).deltaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.(ColumnNames{1,aa}).deltaBandPower.stdS01,1))]; %#ok<*AGROW>
    end
    % delta-band 0.1 Hz PSD p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Rest') == true
            Delta_S01_pVal{1,aa} = {' '};
        else
            Delta_S01_pVal{1,aa} = ['p < ' num2str(Delta_PSD01_Stats.Coefficients.pValue(aa,1))];
        end
    end
    % delta-band 0.01 Hz PSD power
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Awake') == true || strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Delta_S001_MeanStD{1,aa} = [num2str(round(data.PowerSpec.(ColumnNames{1,aa}).deltaBandPower.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.(ColumnNames{1,aa}).deltaBandPower.stdS001,1))];
        else
            Delta_S001_MeanStD{1,aa} = {' '};
        end
    end
    % delta-band 0.01 Hz PSD p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Delta_S001_pVal{1,aa} = ['p < ' num2str(Delta_PSD001_Stats.Coefficients.pValue(aa - 3,1))];
        else
            Delta_S001_pVal{1,aa} = {' '};
        end
    end
    % delta-band 0.1 Hz Coh2 power
    for aa = 1:length(ColumnNames)
        Delta_C01_MeanStD{1,aa} = [num2str(round(data.Coherr.(ColumnNames{1,aa}).deltaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.(ColumnNames{1,aa}).deltaBandPower.stdC01,2))]; %#ok<*AGROW>
    end
    % delta-band 0.1 Hz Coh2 p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Rest') == true
            Delta_C01_pVal{1,aa} = {' '};
        else
            Delta_C01_pVal{1,aa} = ['p < ' num2str(Delta_Coh01_Stats.Coefficients.pValue(aa,1))];
        end
    end
    % delta-band 0.01 Hz Coh2 power
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Awake') == true || strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Delta_C001_MeanStD{1,aa} = [num2str(round(data.Coherr.(ColumnNames{1,aa}).deltaBandPower.meanC001,2)) ' +/- ' num2str(round(data.Coherr.(ColumnNames{1,aa}).deltaBandPower.stdC001,2))];
        else
            Delta_C001_MeanStD{1,aa} = {' '};
        end
    end
    % delta-band 0.01 Hz Coh2 p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Delta_C001_pVal{1,aa} = ['p < ' num2str(Delta_Coh001_Stats.Coefficients.pValue(aa - 3,1))];
        else
            Delta_C001_pVal{1,aa} = {' '};
        end
    end
    % theta-band 0.1 Hz PSD power
    for aa = 1:length(ColumnNames)
        Theta_S01_MeanStD{1,aa} = [num2str(round(data.PowerSpec.(ColumnNames{1,aa}).thetaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.(ColumnNames{1,aa}).thetaBandPower.stdS01,1))]; %#ok<*AGROW>
    end
    % theta-band 0.1 Hz PSD p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Rest') == true
            Theta_S01_pVal{1,aa} = {' '};
        else
            Theta_S01_pVal{1,aa} = ['p < ' num2str(Theta_PSD01_Stats.Coefficients.pValue(aa,1))];
        end
    end
    % theta-band 0.01 Hz PSD power
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Awake') == true || strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Theta_S001_MeanStD{1,aa} = [num2str(round(data.PowerSpec.(ColumnNames{1,aa}).thetaBandPower.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.(ColumnNames{1,aa}).thetaBandPower.stdS001,1))];
        else
            Theta_S001_MeanStD{1,aa} = {' '};
        end
    end
    % theta-band 0.01 Hz PSD p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Theta_S001_pVal{1,aa} = ['p < ' num2str(Theta_PSD001_Stats.Coefficients.pValue(aa - 3,1))];
        else
            Theta_S001_pVal{1,aa} = {' '};
        end
    end
    % theta-band 0.1 Hz Coh2 power
    for aa = 1:length(ColumnNames)
        Theta_C01_MeanStD{1,aa} = [num2str(round(data.Coherr.(ColumnNames{1,aa}).thetaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.(ColumnNames{1,aa}).thetaBandPower.stdC01,2))]; %#ok<*AGROW>
    end
    % theta-band 0.1 Hz Coh2 p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Rest') == true
            Theta_C01_pVal{1,aa} = {' '};
        else
            Theta_C01_pVal{1,aa} = ['p < ' num2str(Theta_Coh01_Stats.Coefficients.pValue(aa,1))];
        end
    end
    % theta-band 0.01 Hz Coh2 power
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Awake') == true || strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Theta_C001_MeanStD{1,aa} = [num2str(round(data.Coherr.(ColumnNames{1,aa}).thetaBandPower.meanC001,2)) ' +/- ' num2str(round(data.Coherr.(ColumnNames{1,aa}).thetaBandPower.stdC001,2))];
        else
            Theta_C001_MeanStD{1,aa} = {' '};
        end
    end
    % theta-band 0.01 Hz Coh2 p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Theta_C001_pVal{1,aa} = ['p < ' num2str(Theta_Coh001_Stats.Coefficients.pValue(aa - 3,1))];
        else
            Theta_C001_pVal{1,aa} = {' '};
        end
    end
    % alpha-band 0.1 Hz PSD power
    for aa = 1:length(ColumnNames)
        Alpha_S01_MeanStD{1,aa} = [num2str(round(data.PowerSpec.(ColumnNames{1,aa}).alphaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.(ColumnNames{1,aa}).alphaBandPower.stdS01,1))]; %#ok<*AGROW>
    end
    % alpha-band 0.1 Hz PSD p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Rest') == true
            Alpha_S01_pVal{1,aa} = {' '};
        else
            Alpha_S01_pVal{1,aa} = ['p < ' num2str(Alpha_PSD01_Stats.Coefficients.pValue(aa,1))];
        end
    end
    % alpha-band 0.01 Hz PSD power
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Awake') == true || strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Alpha_S001_MeanStD{1,aa} = [num2str(round(data.PowerSpec.(ColumnNames{1,aa}).alphaBandPower.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.(ColumnNames{1,aa}).alphaBandPower.stdS001,1))];
        else
            Alpha_S001_MeanStD{1,aa} = {' '};
        end
    end
    % alpha-band 0.01 Hz PSD p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Alpha_S001_pVal{1,aa} = ['p < ' num2str(Alpha_PSD001_Stats.Coefficients.pValue(aa - 3,1))];
        else
            Alpha_S001_pVal{1,aa} = {' '};
        end
    end
    % alpha-band 0.1 Hz Coh2 power
    for aa = 1:length(ColumnNames)
        Alpha_C01_MeanStD{1,aa} = [num2str(round(data.Coherr.(ColumnNames{1,aa}).alphaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.(ColumnNames{1,aa}).alphaBandPower.stdC01,2))]; %#ok<*AGROW>
    end
    % alpha-band 0.1 Hz Coh2 p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Rest') == true
            Alpha_C01_pVal{1,aa} = {' '};
        else
            Alpha_C01_pVal{1,aa} = ['p < ' num2str(Alpha_Coh01_Stats.Coefficients.pValue(aa,1))];
        end
    end
    % alpha-band 0.01 Hz Coh2 power
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Awake') == true || strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Alpha_C001_MeanStD{1,aa} = [num2str(round(data.Coherr.(ColumnNames{1,aa}).alphaBandPower.meanC001,2)) ' +/- ' num2str(round(data.Coherr.(ColumnNames{1,aa}).alphaBandPower.stdC001,2))];
        else
            Alpha_C001_MeanStD{1,aa} = {' '};
        end
    end
    % alpha-band 0.01 Hz Coh2 p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Alpha_C001_pVal{1,aa} = ['p < ' num2str(Alpha_Coh001_Stats.Coefficients.pValue(aa - 3,1))];
        else
            Alpha_C001_pVal{1,aa} = {' '};
        end
    end
    % beta-band 0.1 Hz PSD power
    for aa = 1:length(ColumnNames)
        Beta_S01_MeanStD{1,aa} = [num2str(round(data.PowerSpec.(ColumnNames{1,aa}).betaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.(ColumnNames{1,aa}).betaBandPower.stdS01,1))]; %#ok<*AGROW>
    end
    % beta-band 0.1 Hz PSD p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Rest') == true
            Beta_S01_pVal{1,aa} = {' '};
        else
            Beta_S01_pVal{1,aa} = ['p < ' num2str(Beta_PSD01_Stats.Coefficients.pValue(aa,1))];
        end
    end
    % beta-band 0.01 Hz PSD power
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Awake') == true || strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Beta_S001_MeanStD{1,aa} = [num2str(round(data.PowerSpec.(ColumnNames{1,aa}).betaBandPower.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.(ColumnNames{1,aa}).betaBandPower.stdS001,1))];
        else
            Beta_S001_MeanStD{1,aa} = {' '};
        end
    end
    % beta-band 0.01 Hz PSD p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Beta_S001_pVal{1,aa} = ['p < ' num2str(Beta_PSD001_Stats.Coefficients.pValue(aa - 3,1))];
        else
            Beta_S001_pVal{1,aa} = {' '};
        end
    end
    % beta-band 0.1 Hz Coh2 power
    for aa = 1:length(ColumnNames)
        Beta_C01_MeanStD{1,aa} = [num2str(round(data.Coherr.(ColumnNames{1,aa}).betaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.(ColumnNames{1,aa}).betaBandPower.stdC01,2))]; %#ok<*AGROW>
    end
    % beta-band 0.1 Hz Coh2 p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Rest') == true
            Beta_C01_pVal{1,aa} = {' '};
        else
            Beta_C01_pVal{1,aa} = ['p < ' num2str(Beta_Coh01_Stats.Coefficients.pValue(aa,1))];
        end
    end
    % beta-band 0.01 Hz Coh2 power
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Awake') == true || strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Beta_C001_MeanStD{1,aa} = [num2str(round(data.Coherr.(ColumnNames{1,aa}).betaBandPower.meanC001,2)) ' +/- ' num2str(round(data.Coherr.(ColumnNames{1,aa}).betaBandPower.stdC001,2))];
        else
            Beta_C001_MeanStD{1,aa} = {' '};
        end
    end
    % beta-band 0.01 Hz Coh2 p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Beta_C001_pVal{1,aa} = ['p < ' num2str(Beta_Coh001_Stats.Coefficients.pValue(aa - 3,1))];
        else
            Beta_C001_pVal{1,aa} = {' '};
        end
    end
    %% save table data
    if isfield(AnalysisResults,'PSD') == false
        AnalysisResults.PSD = [];
    end
    if isfield(AnalysisResults.PSD,'deltaBandPower') == false
        AnalysisResults.PSD.columnNames = ColumnNames;
        AnalysisResults.PSD.deltaBandPower.meanStD01 = Delta_S01_MeanStD;
        AnalysisResults.PSD.deltaBandPower.p01 = Delta_S01_pVal;
        AnalysisResults.PSD.thetaBandPower.meanStD01 = Theta_S01_MeanStD;
        AnalysisResults.PSD.thetaBandPower.p01 = Theta_S01_pVal;
        AnalysisResults.PSD.alphaBandPower.meanStD01 = Alpha_S01_MeanStD;
        AnalysisResults.PSD.alphaBandPower.p01 = Alpha_S01_pVal;
        AnalysisResults.PSD.betaBandPower.meanStD01 = Beta_S01_MeanStD;
        AnalysisResults.PSD.betaBandPower.p01 = Beta_S01_pVal;
        AnalysisResults.PSD.deltaBandPower.meanStD001 = Delta_S001_MeanStD;
        AnalysisResults.PSD.deltaBandPower.p001 = Delta_S001_pVal;
        AnalysisResults.PSD.thetaBandPower.meanStD001 = Theta_S001_MeanStD;
        AnalysisResults.PSD.thetaBandPower.p001 = Theta_S001_pVal;
        AnalysisResults.PSD.alphaBandPower.meanStD001 = Alpha_S001_MeanStD;
        AnalysisResults.PSD.alphaBandPower.p001 = Alpha_S001_pVal;
        AnalysisResults.PSD.betaBandPower.meanStD001 = Beta_S001_MeanStD;
        AnalysisResults.PSD.betaBandPower.p001 = Beta_S001_pVal;
        AnalysisResults.Coherr.columnNames = ColumnNames;
        AnalysisResults.Coherr.deltaBandPower.meanStD01 = Delta_C01_MeanStD;
        AnalysisResults.Coherr.deltaBandPower.p01 = Delta_C01_pVal;
        AnalysisResults.Coherr.thetaBandPower.meanStD01 = Theta_C01_MeanStD;
        AnalysisResults.Coherr.thetaBandPower.p01 = Theta_C01_pVal;
        AnalysisResults.Coherr.alphaBandPower.meanStD01 = Alpha_C01_MeanStD;
        AnalysisResults.Coherr.alphaBandPower.p01 = Alpha_C01_pVal;
        AnalysisResults.Coherr.betaBandPower.meanStD01 = Beta_C01_MeanStD;
        AnalysisResults.Coherr.betaBandPower.p01 = Beta_C01_pVal;
        AnalysisResults.Coherr.deltaBandPower.meanStD001 = Delta_C001_MeanStD;
        AnalysisResults.Coherr.deltaBandPower.p001 = Delta_C001_pVal;
        AnalysisResults.Coherr.thetaBandPower.meanStD001 = Theta_C001_MeanStD;
        AnalysisResults.Coherr.thetaBandPower.p001 = Theta_C001_pVal;
        AnalysisResults.Coherr.alphaBandPower.meanStD001 = Alpha_C001_MeanStD;
        AnalysisResults.Coherr.alphaBandPower.p001 = Alpha_C001_pVal;
        AnalysisResults.Coherr.betaBandPower.meanStD001 = Beta_C001_MeanStD;
        AnalysisResults.Coherr.betaBandPower.p001 = Beta_C001_pVal;
        cd(rootFolder)
        save('AnalysisResults.mat','AnalysisResults')
    end
end

end
