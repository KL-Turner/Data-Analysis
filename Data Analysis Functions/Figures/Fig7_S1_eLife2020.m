function [AnalysisResults] = Fig7_S1_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 7-S1 for Turner_Kederasetti_Gheres_Proctor_Costanzo_Drew_eLife2020
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
TwoP_animalIDs = {'T115','T116','T117','T118','T125','T126'};
behavFields = {'Rest','NREM','REM','Awake','Sleep','All'};
behavFields2 = {'Rest','NREM','REM','Awake','All'};
dataTypes = {'CBV_HbT','gammaBandPower'};
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
for ee = 1:length(behavFields)
    behavField = behavFields{1,ee};
    for ff = 1:length(dataTypes)
        dataType = dataTypes{1,ff};
        for gg = 1:size(data.Coherr.(behavField).(dataType).C,2)
            if strcmp(behavField,'Rest') == true
                f_short = data.Coherr.(behavField).(dataType).f(gg,:);
                C = data.Coherr.(behavField).(dataType).C(:,gg);
                f_long = 0:0.01:0.5;
                C_long = interp1(f_short,C,f_long);
                index01 = find(f_long == 0.1);
                data.Coherr.(behavField).(dataType).C01(gg,1) = C_long(index01).^2; %#ok<FNDSB>
            elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
                F = round(data.Coherr.(behavField).(dataType).f(gg,:),2);
                C = data.Coherr.(behavField).(dataType).C(:,gg);
                index01 = find(F == 0.1);
                data.Coherr.(behavField).(dataType).C01(gg,1) = C(index01(1)).^2;
            else
                F = round(data.Coherr.(behavField).(dataType).f(gg,:),3);
                C = data.Coherr.(behavField).(dataType).C(:,gg);
                index01 = find(F == 0.1);
                index001 = find(F == 0.01);
                data.Coherr.(behavField).(dataType).C01(gg,1) = C(index01(1)).^2;
                data.Coherr.(behavField).(dataType).C001(gg,1) = C(index001(1)).^2;
            end
        end
    end
end
% take mean/StD of peak C
for ee = 1:length(behavFields)
    behavField = behavFields{1,ee};
    for ff = 1:length(dataTypes)
        dataType = dataTypes{1,ff};
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
% gamma PSD @ 0.1 Hz
Gamma_Coh01_tableSize = cat(1,data.Coherr.Rest.gammaBandPower.C01,data.Coherr.NREM.gammaBandPower.C01,data.Coherr.REM.gammaBandPower.C01,...
    data.Coherr.Awake.gammaBandPower.C01,data.Coherr.Sleep.gammaBandPower.C01,data.Coherr.All.gammaBandPower.C01);
Gamma_Coh01_Table = table('Size',[size(Gamma_Coh01_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','Coh01'});
Gamma_Coh01_Table.Mouse = cat(1,data.Coherr.Rest.gammaBandPower.animalID,data.Coherr.NREM.gammaBandPower.animalID,data.Coherr.REM.gammaBandPower.animalID,...
    data.Coherr.Awake.gammaBandPower.animalID,data.Coherr.Sleep.gammaBandPower.animalID,data.Coherr.All.gammaBandPower.animalID);
Gamma_Coh01_Table.Behavior = cat(1,data.Coherr.Rest.gammaBandPower.behavior,data.Coherr.NREM.gammaBandPower.behavior,data.Coherr.REM.gammaBandPower.behavior,...
    data.Coherr.Awake.gammaBandPower.behavior,data.Coherr.Sleep.gammaBandPower.behavior,data.Coherr.All.gammaBandPower.behavior);
Gamma_Coh01_Table.Coh01 = cat(1,data.Coherr.Rest.gammaBandPower.C01,data.Coherr.NREM.gammaBandPower.C01,data.Coherr.REM.gammaBandPower.C01,...
    data.Coherr.Awake.gammaBandPower.C01,data.Coherr.Sleep.gammaBandPower.C01,data.Coherr.All.gammaBandPower.C01);
Gamma_Coh01_FitFormula = 'Coh01 ~ 1 + Behavior + (1|Mouse)';
Gamma_Coh01_Stats = fitglme(Gamma_Coh01_Table,Gamma_Coh01_FitFormula);
% gamma PSD @ 0.01 Hz
Gamma_Coh001_tableSize = cat(1,data.Coherr.Awake.gammaBandPower.C001,data.Coherr.Sleep.gammaBandPower.C001,data.Coherr.All.gammaBandPower.C001);
Gamma_Coh001_Table = table('Size',[size(Gamma_Coh001_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','Coh001'});
Gamma_Coh001_Table.Mouse = cat(1,data.Coherr.Awake.gammaBandPower.animalID,data.Coherr.Sleep.gammaBandPower.animalID,data.Coherr.All.gammaBandPower.animalID);
Gamma_Coh001_Table.Behavior = cat(1,data.Coherr.Awake.gammaBandPower.behavior,data.Coherr.Sleep.gammaBandPower.behavior,data.Coherr.All.gammaBandPower.behavior);
Gamma_Coh001_Table.Coh001 = cat(1,data.Coherr.Awake.gammaBandPower.C001,data.Coherr.Sleep.gammaBandPower.C001,data.Coherr.All.gammaBandPower.C001);
Gamma_Coh001_FitFormula = 'Coh001 ~ 1 + Behavior + (1|Mouse)';
Gamma_Coh001_Stats = fitglme(Gamma_Coh001_Table,Gamma_Coh001_FitFormula);
% HbT PSD @ 0.1 Hz
HbT_Coh01_tableSize = cat(1,data.Coherr.Rest.CBV_HbT.C01,data.Coherr.NREM.CBV_HbT.C01,data.Coherr.REM.CBV_HbT.C01,...
    data.Coherr.Awake.CBV_HbT.C01,data.Coherr.Sleep.CBV_HbT.C01,data.Coherr.All.CBV_HbT.C01);
HbT_Coh01_Table = table('Size',[size(HbT_Coh01_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','Coh01'});
HbT_Coh01_Table.Mouse = cat(1,data.Coherr.Rest.CBV_HbT.animalID,data.Coherr.NREM.CBV_HbT.animalID,data.Coherr.REM.CBV_HbT.animalID,...
    data.Coherr.Awake.CBV_HbT.animalID,data.Coherr.Sleep.CBV_HbT.animalID,data.Coherr.All.CBV_HbT.animalID);
HbT_Coh01_Table.Behavior = cat(1,data.Coherr.Rest.CBV_HbT.behavior,data.Coherr.NREM.CBV_HbT.behavior,data.Coherr.REM.CBV_HbT.behavior,...
    data.Coherr.Awake.CBV_HbT.behavior,data.Coherr.Sleep.CBV_HbT.behavior,data.Coherr.All.CBV_HbT.behavior);
HbT_Coh01_Table.Coh01 = cat(1,data.Coherr.Rest.CBV_HbT.C01,data.Coherr.NREM.CBV_HbT.C01,data.Coherr.REM.CBV_HbT.C01,...
    data.Coherr.Awake.CBV_HbT.C01,data.Coherr.Sleep.CBV_HbT.C01,data.Coherr.All.CBV_HbT.C01);
HbT_Coh01_FitFormula = 'Coh01 ~ 1 + Behavior + (1|Mouse)';
HbT_Coh01_Stats = fitglme(HbT_Coh01_Table,HbT_Coh01_FitFormula);
% HbT PSD @ 0.01 Hz
HbT_Coh001_tableSize = cat(1,data.Coherr.Awake.CBV_HbT.C001,data.Coherr.Sleep.CBV_HbT.C001,data.Coherr.All.CBV_HbT.C001);
HbT_Coh001_Table = table('Size',[size(HbT_Coh001_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','Coh001'});
HbT_Coh001_Table.Mouse = cat(1,data.Coherr.Awake.CBV_HbT.animalID,data.Coherr.Sleep.CBV_HbT.animalID,data.Coherr.All.CBV_HbT.animalID);
HbT_Coh001_Table.Behavior = cat(1,data.Coherr.Awake.CBV_HbT.behavior,data.Coherr.Sleep.CBV_HbT.behavior,data.Coherr.All.CBV_HbT.behavior);
HbT_Coh001_Table.Coh001 = cat(1,data.Coherr.Awake.CBV_HbT.C001,data.Coherr.Sleep.CBV_HbT.C001,data.Coherr.All.CBV_HbT.C001);
HbT_Coh001_FitFormula = 'Coh001 ~ 1 + Behavior + (1|Mouse)';
HbT_Coh001_Stats = fitglme(HbT_Coh001_Table,HbT_Coh001_FitFormula);
%% power spectra during different behaviors
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
for ee = 1:length(behavFields)
    behavField = behavFields{1,ee};
    for ff = 1:length(dataTypes)
        dataType = dataTypes{1,ff};
        data.PowerSpec.(behavField).(dataType).cat_S = [];
        data.PowerSpec.(behavField).(dataType).cat_f = [];
        data.PowerSpec.(behavField).(dataType).animalID = {};
        data.PowerSpec.(behavField).(dataType).behavior = {};
        data.PowerSpec.(behavField).(dataType).hemisphere = {};
        for z = 1:length(data.PowerSpec.(behavField).(dataType).normLH)
            data.PowerSpec.(behavField).(dataType).cat_S = cat(2,data.PowerSpec.(behavField).(dataType).cat_S,data.PowerSpec.(behavField).(dataType).normLH{z,1},data.PowerSpec.(behavField).(dataType).normRH{z,1});
            data.PowerSpec.(behavField).(dataType).cat_f = cat(1,data.PowerSpec.(behavField).(dataType).cat_f,data.PowerSpec.(behavField).(dataType).adjLH.f{z,1},data.PowerSpec.(behavField).(dataType).adjRH.f{z,1});
            if isempty(data.PowerSpec.(behavField).(dataType).normLH{z,1}) == false
                data.PowerSpec.(behavField).(dataType).animalID = cat(1,data.PowerSpec.(behavField).(dataType).animalID,animalID,animalID);
                data.PowerSpec.(behavField).(dataType).behavior = cat(1,data.PowerSpec.(behavField).(dataType).behavior,behavField,behavField);
                data.PowerSpec.(behavField).(dataType).hemisphere = cat(1,data.PowerSpec.(behavField).(dataType).hemisphere,'LH','RH');
            end
        end
    end
end
% find 0.1/0.01 Hz peaks in PSD
for ee = 1:length(behavFields)
    behavField = behavFields{1,ee};
    for ff = 1:length(dataTypes)
        dataType = dataTypes{1,ff};
        for gg = 1:size(data.PowerSpec.(behavField).(dataType).cat_S,2)
            if strcmp(behavField,'Rest') == true
                f_short = data.PowerSpec.(behavField).(dataType).cat_f(gg,:);
                S = data.PowerSpec.(behavField).(dataType).cat_S(:,gg);
                f_long = 0:0.01:0.5;
                S_long = interp1(f_short,S,f_long);
                index01 = find(f_long == 0.1);
                data.PowerSpec.(behavField).(dataType).S01(gg,1) = S_long(index01); %#ok<FNDSB>
            elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
                F = round(data.PowerSpec.(behavField).(dataType).cat_f(gg,:),2);
                S = data.PowerSpec.(behavField).(dataType).cat_S(:,gg);
                index01 = find(F == 0.1);
                data.PowerSpec.(behavField).(dataType).S01(gg,1) = S(index01(1));
            else
                F = round(data.PowerSpec.(behavField).(dataType).cat_f(gg,:),3);
                S = data.PowerSpec.(behavField).(dataType).cat_S(:,gg);
                index01 = find(F == 0.1);
                index001 = find(F == 0.01);
                data.PowerSpec.(behavField).(dataType).S01(gg,1) = S(index01(1));
                data.PowerSpec.(behavField).(dataType).S001(gg,1) = S(index001(1));
            end
        end
    end
end
% take mean/StD of peak S
for ee = 1:length(behavFields)
    behavField = behavFields{1,ee};
    for ff = 1:length(dataTypes)
        dataType = dataTypes{1,ff};
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
% gamma PSD @ 0.1 Hz
Gamma_PSD01_tableSize = cat(1,data.PowerSpec.Rest.gammaBandPower.S01,data.PowerSpec.NREM.gammaBandPower.S01,data.PowerSpec.REM.gammaBandPower.S01,...
    data.PowerSpec.Awake.gammaBandPower.S01,data.PowerSpec.Sleep.gammaBandPower.S01,data.PowerSpec.All.gammaBandPower.S01);
Gamma_PSD01_Table = table('Size',[size(Gamma_PSD01_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','PSD01'});
Gamma_PSD01_Table.Mouse = cat(1,data.PowerSpec.Rest.gammaBandPower.animalID,data.PowerSpec.NREM.gammaBandPower.animalID,data.PowerSpec.REM.gammaBandPower.animalID,...
    data.PowerSpec.Awake.gammaBandPower.animalID,data.PowerSpec.Sleep.gammaBandPower.animalID,data.PowerSpec.All.gammaBandPower.animalID);
Gamma_PSD01_Table.Hemisphere = cat(1,data.PowerSpec.Rest.gammaBandPower.hemisphere,data.PowerSpec.NREM.gammaBandPower.hemisphere,data.PowerSpec.REM.gammaBandPower.hemisphere,...
    data.PowerSpec.Awake.gammaBandPower.hemisphere,data.PowerSpec.Sleep.gammaBandPower.hemisphere,data.PowerSpec.All.gammaBandPower.hemisphere);
Gamma_PSD01_Table.Behavior = cat(1,data.PowerSpec.Rest.gammaBandPower.behavior,data.PowerSpec.NREM.gammaBandPower.behavior,data.PowerSpec.REM.gammaBandPower.behavior,...
    data.PowerSpec.Awake.gammaBandPower.behavior,data.PowerSpec.Sleep.gammaBandPower.behavior,data.PowerSpec.All.gammaBandPower.behavior);
Gamma_PSD01_Table.PSD01 = cat(1,data.PowerSpec.Rest.gammaBandPower.S01,data.PowerSpec.NREM.gammaBandPower.S01,data.PowerSpec.REM.gammaBandPower.S01,...
    data.PowerSpec.Awake.gammaBandPower.S01,data.PowerSpec.Sleep.gammaBandPower.S01,data.PowerSpec.All.gammaBandPower.S01);
Gamma_PSD01_FitFormula = 'PSD01 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
Gamma_PSD01_Stats = fitglme(Gamma_PSD01_Table,Gamma_PSD01_FitFormula);
% gamma PSD @ 0.01 Hz
Gamma_PSD001_tableSize = cat(1,data.PowerSpec.Awake.gammaBandPower.S001,data.PowerSpec.Sleep.gammaBandPower.S001,data.PowerSpec.All.gammaBandPower.S001);
Gamma_PSD001_Table = table('Size',[size(Gamma_PSD001_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','PSD001'});
Gamma_PSD001_Table.Mouse = cat(1,data.PowerSpec.Awake.gammaBandPower.animalID,data.PowerSpec.Sleep.gammaBandPower.animalID,data.PowerSpec.All.gammaBandPower.animalID);
Gamma_PSD001_Table.Hemisphere = cat(1,data.PowerSpec.Awake.gammaBandPower.hemisphere,data.PowerSpec.Sleep.gammaBandPower.hemisphere,data.PowerSpec.All.gammaBandPower.hemisphere);
Gamma_PSD001_Table.Behavior = cat(1,data.PowerSpec.Awake.gammaBandPower.behavior,data.PowerSpec.Sleep.gammaBandPower.behavior,data.PowerSpec.All.gammaBandPower.behavior);
Gamma_PSD001_Table.PSD001 = cat(1,data.PowerSpec.Awake.gammaBandPower.S001,data.PowerSpec.Sleep.gammaBandPower.S001,data.PowerSpec.All.gammaBandPower.S001);
Gamma_PSD001_FitFormula = 'PSD001 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
Gamma_PSD001_Stats = fitglme(Gamma_PSD001_Table,Gamma_PSD001_FitFormula);
% HbT PSD @ 0.1 Hz
HbT_PSD01_tableSize = cat(1,data.PowerSpec.Rest.CBV_HbT.S01,data.PowerSpec.NREM.CBV_HbT.S01,data.PowerSpec.REM.CBV_HbT.S01,...
    data.PowerSpec.Awake.CBV_HbT.S01,data.PowerSpec.Sleep.CBV_HbT.S01,data.PowerSpec.All.CBV_HbT.S01);
HbT_PSD01_Table = table('Size',[size(HbT_PSD01_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','PSD01'});
HbT_PSD01_Table.Mouse = cat(1,data.PowerSpec.Rest.CBV_HbT.animalID,data.PowerSpec.NREM.CBV_HbT.animalID,data.PowerSpec.REM.CBV_HbT.animalID,...
    data.PowerSpec.Awake.CBV_HbT.animalID,data.PowerSpec.Sleep.CBV_HbT.animalID,data.PowerSpec.All.CBV_HbT.animalID);
HbT_PSD01_Table.Hemisphere = cat(1,data.PowerSpec.Rest.CBV_HbT.hemisphere,data.PowerSpec.NREM.CBV_HbT.hemisphere,data.PowerSpec.REM.CBV_HbT.hemisphere,...
    data.PowerSpec.Awake.CBV_HbT.hemisphere,data.PowerSpec.Sleep.CBV_HbT.hemisphere,data.PowerSpec.All.CBV_HbT.hemisphere);
HbT_PSD01_Table.Behavior = cat(1,data.PowerSpec.Rest.CBV_HbT.behavior,data.PowerSpec.NREM.CBV_HbT.behavior,data.PowerSpec.REM.CBV_HbT.behavior,...
    data.PowerSpec.Awake.CBV_HbT.behavior,data.PowerSpec.Sleep.CBV_HbT.behavior,data.PowerSpec.All.CBV_HbT.behavior);
HbT_PSD01_Table.PSD01 = cat(1,data.PowerSpec.Rest.CBV_HbT.S01,data.PowerSpec.NREM.CBV_HbT.S01,data.PowerSpec.REM.CBV_HbT.S01,...
    data.PowerSpec.Awake.CBV_HbT.S01,data.PowerSpec.Sleep.CBV_HbT.S01,data.PowerSpec.All.CBV_HbT.S01);
HbT_PSD01_FitFormula = 'PSD01 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
HbT_PSD01_Stats = fitglme(HbT_PSD01_Table,HbT_PSD01_FitFormula);
% HbT PSD @ 0.01 Hz
HbT_PSD001_tableSize = cat(1,data.PowerSpec.Awake.CBV_HbT.S001,data.PowerSpec.Sleep.CBV_HbT.S001,data.PowerSpec.All.CBV_HbT.S001);
HbT_PSD001_Table = table('Size',[size(HbT_PSD001_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','PSD001'});
HbT_PSD001_Table.Mouse = cat(1,data.PowerSpec.Awake.CBV_HbT.animalID,data.PowerSpec.Sleep.CBV_HbT.animalID,data.PowerSpec.All.CBV_HbT.animalID);
HbT_PSD001_Table.Hemisphere = cat(1,data.PowerSpec.Awake.CBV_HbT.hemisphere,data.PowerSpec.Sleep.CBV_HbT.hemisphere,data.PowerSpec.All.CBV_HbT.hemisphere);
HbT_PSD001_Table.Behavior = cat(1,data.PowerSpec.Awake.CBV_HbT.behavior,data.PowerSpec.Sleep.CBV_HbT.behavior,data.PowerSpec.All.CBV_HbT.behavior);
HbT_PSD001_Table.PSD001 = cat(1,data.PowerSpec.Awake.CBV_HbT.S001,data.PowerSpec.Sleep.CBV_HbT.S001,data.PowerSpec.All.CBV_HbT.S001);
HbT_PSD001_FitFormula = 'PSD001 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
HbT_PSD001_Stats = fitglme(HbT_PSD001_Table,HbT_PSD001_FitFormula);
%% vessel power spectra of different behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.VesselPowerSpec = [];
for aa = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,aa};
    for bb = 1:length(behavFields2)
        behavField = behavFields2{1,bb};
        if isfield(AnalysisResults.(animalID).PowerSpectra,behavField) == true
            vesselIDs = fieldnames(AnalysisResults.(animalID).PowerSpectra.(behavField));
            for cc = 1:length(vesselIDs)
                vesselID = vesselIDs{cc,1};
                data.VesselPowerSpec.(animalID).(behavField).(vesselID).S = AnalysisResults.(animalID).PowerSpectra.(behavField).(vesselID).S;
                data.VesselPowerSpec.(animalID).(behavField).(vesselID).f = AnalysisResults.(animalID).PowerSpectra.(behavField).(vesselID).f;
            end
        end
    end
end
% find the peak of the resting PSD for each arteriole
for aa = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,aa};
    vesselIDs = fieldnames(data.VesselPowerSpec.(animalID).Rest);
    for cc = 1:length(vesselIDs)
        vesselID = vesselIDs{cc,1};
        data.VesselPowerSpec.(animalID).baseline.(vesselID) = max(data.VesselPowerSpec.(animalID).Rest.(vesselID).S);
    end
end
% DC-shift each arteriole's PSD with respect to the resting peak
for aa = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,aa};
    for bb = 1:length(behavFields2)
        behavField = behavFields2{1,bb};
        if isfield(data.VesselPowerSpec.(animalID),behavField) == true
            vesselIDs = fieldnames(data.VesselPowerSpec.(animalID).(behavField));
            for cc = 1:length(vesselIDs)
                vesselID = vesselIDs{cc,1};
                data.VesselPowerSpec.(animalID).(behavField).(vesselID).normS = data.VesselPowerSpec.(animalID).(behavField).(vesselID).S*(1/data.VesselPowerSpec.(animalID).baseline.(vesselID));
            end
        end
    end
end
% concatenate the data from all arterioles - removes any empty data
for aa = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,aa};
    for bb = 1:length(behavFields2)
        behavField = behavFields2{1,bb};
        if isfield(data.VesselPowerSpec,behavField) == false
            data.VesselPowerSpec.(behavField).S = [];
            data.VesselPowerSpec.(behavField).f = [];
            data.VesselPowerSpec.(behavField).animalID = {};
            data.VesselPowerSpec.(behavField).behavior = {};
            data.VesselPowerSpec.(behavField).vessel = {};
        end
        if isfield(data.VesselPowerSpec.(animalID),behavField) == true
            vesselIDs = fieldnames(data.VesselPowerSpec.(animalID).(behavField));
            for cc = 1:length(vesselIDs)
                vesselID = vesselIDs{cc,1};
                data.VesselPowerSpec.(behavField).S = cat(2,data.VesselPowerSpec.(behavField).S,data.VesselPowerSpec.(animalID).(behavField).(vesselID).normS);
                data.VesselPowerSpec.(behavField).f = cat(1,data.VesselPowerSpec.(behavField).f,data.VesselPowerSpec.(animalID).(behavField).(vesselID).f);
                data.VesselPowerSpec.(behavField).animalID = cat(1,data.VesselPowerSpec.(behavField).animalID,animalID);
                data.VesselPowerSpec.(behavField).behavior = cat(1,data.VesselPowerSpec.(behavField).behavior,behavField);
                data.VesselPowerSpec.(behavField).vessel = cat(1,data.VesselPowerSpec.(behavField).vessel,vesselID);
            end
        end
    end
end
% find 0.1/0.01 Hz peaks in PSD
for ee = 1:length(behavFields2)
    behavField = behavFields2{1,ee};
    for gg = 1:size(data.VesselPowerSpec.(behavField).S,2)
        if strcmp(behavField,'Rest') == true
            f_short = data.VesselPowerSpec.(behavField).f(gg,:);
            S = data.VesselPowerSpec.(behavField).S(:,gg);
            f_long = 0:0.01:0.5;
            S_long = interp1(f_short,S,f_long);
            index01 = find(f_long == 0.1);
            data.VesselPowerSpec.(behavField).S01(gg,1) = S_long(index01); %#ok<FNDSB>
        elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
            F = round(data.VesselPowerSpec.(behavField).f(gg,:),2);
            S = data.VesselPowerSpec.(behavField).S(:,gg);
            index01 = find(F == 0.1);
            data.VesselPowerSpec.(behavField).S01(gg,1) = S(index01(1));
        else
            F = round(data.VesselPowerSpec.(behavField).f(gg,:),3);
            S = data.VesselPowerSpec.(behavField).S(:,gg);
            index01 = find(F == 0.1);
            index001 = find(F == 0.01);
            data.VesselPowerSpec.(behavField).S01(gg,1) = S(index01(1));
            data.VesselPowerSpec.(behavField).S001(gg,1) = S(index001(1));
        end
    end
end
% take mean/StD of peak S
for ee = 1:length(behavFields2)
    behavField = behavFields2{1,ee};
    if strcmp(behavField,'Rest') == true || strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
        data.VesselPowerSpec.(behavField).meanS01 = mean(data.VesselPowerSpec.(behavField).S01,1);
        data.VesselPowerSpec.(behavField).stdS01 = std(data.VesselPowerSpec.(behavField).S01,0,1);
    else
        data.VesselPowerSpec.(behavField).meanS01 = mean(data.VesselPowerSpec.(behavField).S01,1);
        data.VesselPowerSpec.(behavField).stdS01 = std(data.VesselPowerSpec.(behavField).S01,0,1);
        data.VesselPowerSpec.(behavField).meanS001 = mean(data.VesselPowerSpec.(behavField).S001,1);
        data.VesselPowerSpec.(behavField).stdS001 = std(data.VesselPowerSpec.(behavField).S001,0,1);
    end
end
%% statistics - generalized linear mixed effects model
% vessel PSD @ 0.1 Hz
TwoP_PSD01_tableSize = cat(1,data.VesselPowerSpec.Rest.S01,data.VesselPowerSpec.NREM.S01,data.VesselPowerSpec.REM.S01,...
    data.VesselPowerSpec.Awake.S01,data.VesselPowerSpec.All.S01);
TwoP_PSD01_Table = table('Size',[size(TwoP_PSD01_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','PSD01'});
TwoP_PSD01_Table.Mouse = cat(1,data.VesselPowerSpec.Rest.animalID,data.VesselPowerSpec.NREM.animalID,data.VesselPowerSpec.REM.animalID,...
    data.VesselPowerSpec.Awake.animalID,data.VesselPowerSpec.All.animalID);
TwoP_PSD01_Table.Vessel = cat(1,data.VesselPowerSpec.Rest.vessel,data.VesselPowerSpec.NREM.vessel,data.VesselPowerSpec.REM.vessel,...
    data.VesselPowerSpec.Awake.vessel,data.VesselPowerSpec.All.vessel);
TwoP_PSD01_Table.Behavior = cat(1,data.VesselPowerSpec.Rest.behavior,data.VesselPowerSpec.NREM.behavior,data.VesselPowerSpec.REM.behavior,...
    data.VesselPowerSpec.Awake.behavior,data.VesselPowerSpec.All.behavior);
TwoP_PSD01_Table.PSD01 = cat(1,data.VesselPowerSpec.Rest.S01,data.VesselPowerSpec.NREM.S01,data.VesselPowerSpec.REM.S01,...
    data.VesselPowerSpec.Awake.S01,data.VesselPowerSpec.All.S01);
TwoP_PSD01_FitFormula = 'PSD01 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Vessel)';
TwoP_PSD01_Stats = fitglme(TwoP_PSD01_Table,TwoP_PSD01_FitFormula);
% vessel PSD @ 0.01 Hz
TwoP_PSD001_tableSize = cat(1,data.VesselPowerSpec.Awake.S001,data.VesselPowerSpec.All.S001);
TwoP_PSD001_Table = table('Size',[size(TwoP_PSD001_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Vessel','Behavior','PSD001'});
TwoP_PSD001_Table.Mouse = cat(1,data.VesselPowerSpec.Awake.animalID,data.VesselPowerSpec.All.animalID);
TwoP_PSD001_Table.Vessel = cat(1,data.VesselPowerSpec.Awake.vessel,data.VesselPowerSpec.All.vessel);
TwoP_PSD001_Table.Behavior = cat(1,data.VesselPowerSpec.Awake.behavior,data.VesselPowerSpec.All.behavior);
TwoP_PSD001_Table.PSD001 = cat(1,data.VesselPowerSpec.Awake.S001,data.VesselPowerSpec.All.S001);
TwoP_PSD001_FitFormula = 'PSD001 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Vessel)';
TwoP_PSD001_Stats = fitglme(TwoP_PSD001_Table,TwoP_PSD001_FitFormula);
%% Fig. 7-S1
summaryFigure = figure('Name','Fig7-S1 (a-j)');
sgtitle('Figure 7-S1 - Turner et al. 2020')
%% [7-S1a] gamma PSD @ 0.1 Hz
ax1 = subplot(3,4,1);
s1 = scatter(ones(1,length(data.PowerSpec.Rest.gammaBandPower.S01))*1,data.PowerSpec.Rest.gammaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.PowerSpec.Rest.gammaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(ones(1,length(data.PowerSpec.NREM.gammaBandPower.S01))*2,data.PowerSpec.NREM.gammaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.PowerSpec.NREM.gammaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(ones(1,length(data.PowerSpec.REM.gammaBandPower.S01))*3,data.PowerSpec.REM.gammaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.PowerSpec.REM.gammaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
s4 = scatter(ones(1,length(data.PowerSpec.Awake.gammaBandPower.S01))*4,data.PowerSpec.Awake.gammaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.PowerSpec.Awake.gammaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
s5 = scatter(ones(1,length(data.PowerSpec.Sleep.gammaBandPower.S01))*5,data.PowerSpec.Sleep.gammaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.PowerSpec.Sleep.gammaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
s6 = scatter(ones(1,length(data.PowerSpec.All.gammaBandPower.S01))*6,data.PowerSpec.All.gammaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.PowerSpec.All.gammaBandPower.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
title({'[7-S1a] PSD @ 0.1 Hz','Gamma-band [30-100 Hz]'})
ylabel('Power (a.u.)')
legend([s1,s2,s3,s4,s5,s6],'Rest','NREM','REM','Alert','Asleep','All')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
axis square
xlim([0,7])
ylim([0.1,1000])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [7-S1b] gamma PSD @ 0.01 Hz
ax2 = subplot(3,4,2);
scatter(ones(1,length(data.PowerSpec.Awake.gammaBandPower.S001))*1,data.PowerSpec.Awake.gammaBandPower.S001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.PowerSpec.Awake.gammaBandPower.meanS001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Sleep.gammaBandPower.S001))*2,data.PowerSpec.Sleep.gammaBandPower.S001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.PowerSpec.Sleep.gammaBandPower.meanS001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.All.gammaBandPower.S001))*3,data.PowerSpec.All.gammaBandPower.S001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.PowerSpec.All.gammaBandPower.meanS001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[7-S1b] PSD @ 0.01 Hz','Gamma-band [30-100 Hz]'})
ylabel('Power (a.u.)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
axis square
xlim([0,4])
ylim([0.1,1000])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [7-S1c] gamma coherence^2 @ 0.1 Hz
ax3 = subplot(3,4,3);
scatter(ones(1,length(data.Coherr.Rest.gammaBandPower.C01))*1,data.Coherr.Rest.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Coherr.Rest.gammaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.Coherr.NREM.gammaBandPower.C01))*2,data.Coherr.NREM.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Coherr.NREM.gammaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.Coherr.REM.gammaBandPower.C01))*3,data.Coherr.REM.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.Coherr.REM.gammaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.Coherr.Awake.gammaBandPower.C01))*4,data.Coherr.Awake.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.Coherr.Awake.gammaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.Coherr.Sleep.gammaBandPower.C01))*5,data.Coherr.Sleep.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.Coherr.Sleep.gammaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(ones(1,length(data.Coherr.All.gammaBandPower.C01))*6,data.Coherr.All.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.Coherr.All.gammaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
title({'[7-S1c] Coherence^2 @ 0.1 Hz','Gamma-band [30-100 Hz]'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,7])
ylim([0,1])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [7-S1d] gamma coherence^2 @ 0.01 Hz
ax4 = subplot(3,4,4);
scatter(ones(1,length(data.Coherr.Awake.gammaBandPower.C001))*1,data.Coherr.Awake.gammaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Coherr.Awake.gammaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.Coherr.Sleep.gammaBandPower.C001))*2,data.Coherr.Sleep.gammaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Coherr.Sleep.gammaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.Coherr.All.gammaBandPower.C001))*3,data.Coherr.All.gammaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.Coherr.All.gammaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[7-S1d] Coherence^2 @ 0.01 Hz','Gamma-band [30-100 Hz]'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,4])
ylim([0,1])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [7-S1e] HbT PSD @ 0.1 Hz
ax5 = subplot(3,4,5);
scatter(ones(1,length(data.PowerSpec.Rest.CBV_HbT.S01))*1,data.PowerSpec.Rest.CBV_HbT.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.PowerSpec.Rest.CBV_HbT.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.NREM.CBV_HbT.S01))*2,data.PowerSpec.NREM.CBV_HbT.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.PowerSpec.NREM.CBV_HbT.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.REM.CBV_HbT.S01))*3,data.PowerSpec.REM.CBV_HbT.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.PowerSpec.REM.CBV_HbT.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Awake.CBV_HbT.S01))*4,data.PowerSpec.Awake.CBV_HbT.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.PowerSpec.Awake.CBV_HbT.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Sleep.CBV_HbT.S01))*5,data.PowerSpec.Sleep.CBV_HbT.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.PowerSpec.Sleep.CBV_HbT.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.All.CBV_HbT.S01))*6,data.PowerSpec.All.CBV_HbT.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.PowerSpec.All.CBV_HbT.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
title({'[7-S1e] PSD @ 0.1 Hz','\Delta[HbT] (\muM)'})
ylabel('Power (a.u.)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
axis square
xlim([0,7])
ylim([0.1,100])
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [7-S1f] HbT PSD @ 0.01 Hz
ax6 = subplot(3,4,6);
scatter(ones(1,length(data.PowerSpec.Awake.CBV_HbT.S001))*1,data.PowerSpec.Awake.CBV_HbT.S001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.PowerSpec.Awake.CBV_HbT.meanS001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Sleep.CBV_HbT.S001))*2,data.PowerSpec.Sleep.CBV_HbT.S001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.PowerSpec.Sleep.CBV_HbT.meanS001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.All.CBV_HbT.S001))*3,data.PowerSpec.All.CBV_HbT.S001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.PowerSpec.All.CBV_HbT.meanS001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[7-S1f] PSD @ 0.01 Hz','\Delta[HbT] (\muM)'})
ylabel('Power (a.u.)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
axis square
xlim([0,4])
ylim([1,1000])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% [7-S1g] HbT coherence^2 @ 0.1 Hz
ax7 = subplot(3,4,7);
scatter(ones(1,length(data.Coherr.Rest.CBV_HbT.C01))*1,data.Coherr.Rest.CBV_HbT.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Coherr.Rest.CBV_HbT.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.Coherr.NREM.CBV_HbT.C01))*2,data.Coherr.NREM.CBV_HbT.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Coherr.NREM.CBV_HbT.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.Coherr.REM.CBV_HbT.C01))*3,data.Coherr.REM.CBV_HbT.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.Coherr.REM.CBV_HbT.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.Coherr.Awake.CBV_HbT.C01))*4,data.Coherr.Awake.CBV_HbT.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.Coherr.Awake.CBV_HbT.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.Coherr.Sleep.CBV_HbT.C01))*5,data.Coherr.Sleep.CBV_HbT.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.Coherr.Sleep.CBV_HbT.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(ones(1,length(data.Coherr.All.CBV_HbT.C01))*6,data.Coherr.All.CBV_HbT.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.Coherr.All.CBV_HbT.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
title({'[7-S1g] Coherence^2 @ 0.1 Hz','\Delta[HbT] (\muM)'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,7])
ylim([0,1])
set(gca,'box','off')
ax7.TickLength = [0.03,0.03];
%% [7-S1h] HbT coherence^2 @ 0.01 Hz
ax8 = subplot(3,4,8);
scatter(ones(1,length(data.Coherr.Awake.CBV_HbT.C001))*1,data.Coherr.Awake.CBV_HbT.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Coherr.Awake.CBV_HbT.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.Coherr.Sleep.CBV_HbT.C001))*2,data.Coherr.Sleep.CBV_HbT.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Coherr.Sleep.CBV_HbT.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.Coherr.All.CBV_HbT.C001))*3,data.Coherr.All.CBV_HbT.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.Coherr.All.CBV_HbT.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[7-S1h] Coherence^2 @ 0.01 Hz','\Delta[HbT] (\muM)'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,4])
ylim([0,1])
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];
%% [7-S1i] vessel PSD @ 0.1 Hz
ax9 = subplot(3,4,9);
scatter(ones(1,length(data.VesselPowerSpec.Rest.S01))*1,data.VesselPowerSpec.Rest.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.VesselPowerSpec.Rest.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.VesselPowerSpec.NREM.S01))*2,data.VesselPowerSpec.NREM.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.VesselPowerSpec.NREM.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.VesselPowerSpec.REM.S01))*3,data.VesselPowerSpec.REM.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.VesselPowerSpec.REM.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.VesselPowerSpec.Awake.S01))*4,data.VesselPowerSpec.Awake.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.VesselPowerSpec.Awake.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.VesselPowerSpec.All.S01))*5,data.VesselPowerSpec.All.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.VesselPowerSpec.All.meanS01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
title({'[7-S1i] PSD @ 0.1 Hz','\DeltaD/D (%)'})
ylabel('PSD (a.u.)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
axis square
xlim([0,6])
ylim([0.1,100])
set(gca,'box','off')
ax9.TickLength = [0.03,0.03];
%% [7-S1j] vessel PSD @ 0.01 Hz
ax10 = subplot(3,4,10);
scatter(ones(1,length(data.VesselPowerSpec.Awake.S001))*1,data.VesselPowerSpec.Awake.S001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.VesselPowerSpec.Awake.meanS001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.VesselPowerSpec.All.S001))*2,data.VesselPowerSpec.All.S001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.VesselPowerSpec.All.meanS001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
title({'[7-S1j] PSD @ 0.01 Hz','\DeltaD/D (%)'})
ylabel('PSD (a.u.)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
axis square
xlim([0,3])
ylim([1,1000])
set(gca,'box','off')
ax10.TickLength = [0.03,0.03];
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig7-S1']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig7-S1'])
    %% statistical diary
    diaryFile = [dirpath 'Fig7-S1_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    % gamma-band 0.1 Hz PSD statistical diary
    disp('======================================================================================================================')
    disp('[7-S1a] Generalized linear mixed-effects model statistics for gamma-band PSD @ 0.1 Hz for Rest, NREM, REM, Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Gamma_PSD01_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  Gamma 0.1 Hz PSD: ' num2str(round(data.PowerSpec.Rest.gammaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.Rest.gammaBandPower.stdS01,1))]); disp(' ')
    disp(['NREM  Gamma 0.1 Hz PSD: ' num2str(round(data.PowerSpec.NREM.gammaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.NREM.gammaBandPower.stdS01,1))]); disp(' ')
    disp(['REM   Gamma 0.1 Hz PSD: ' num2str(round(data.PowerSpec.REM.gammaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.REM.gammaBandPower.stdS01,1))]); disp(' ')
    disp(['Awake Gamma 0.1 Hz PSD: ' num2str(round(data.PowerSpec.Awake.gammaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.Awake.gammaBandPower.stdS01,1))]); disp(' ')
    disp(['Sleep Gamma 0.1 Hz PSD: ' num2str(round(data.PowerSpec.Sleep.gammaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.Sleep.gammaBandPower.stdS01,1))]); disp(' ')
    disp(['All   Gamma 0.1 Hz PSD: ' num2str(round(data.PowerSpec.All.gammaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.All.gammaBandPower.stdS01,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % gamma-band 0.01 Hz PSD statistical diary
    disp('======================================================================================================================')
    disp('[7-S1b] Generalized linear mixed-effects model statistics for gamma-band PSD @ 0.01 Hz for Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Gamma_PSD001_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake Gamma 0.01 Hz PSD: ' num2str(round(data.PowerSpec.Awake.gammaBandPower.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.Awake.gammaBandPower.stdS001,1))]); disp(' ')
    disp(['Sleep Gamma 0.01 Hz PSD: ' num2str(round(data.PowerSpec.Sleep.gammaBandPower.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.Sleep.gammaBandPower.stdS001,1))]); disp(' ')
    disp(['All   Gamma 0.01 Hz PSD: ' num2str(round(data.PowerSpec.All.gammaBandPower.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.All.gammaBandPower.stdS001,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % gamma-band 0.1 Hz coherence^2 statistical diary
    disp('======================================================================================================================')
    disp('[7-S1c] Generalized linear mixed-effects model statistics for gamma-band coherence^2 @ 0.1 Hz for Rest, NREM, REM, Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Gamma_Coh01_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  Gamma 0.1 Hz Coh2: ' num2str(round(data.Coherr.Rest.gammaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.Rest.gammaBandPower.stdC01,2))]); disp(' ')
    disp(['NREM  Gamma 0.1 Hz Coh2: ' num2str(round(data.Coherr.NREM.gammaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.NREM.gammaBandPower.stdC01,2))]); disp(' ')
    disp(['REM   Gamma 0.1 Hz Coh2: ' num2str(round(data.Coherr.REM.gammaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.REM.gammaBandPower.stdC01,2))]); disp(' ')
    disp(['Awake Gamma 0.1 Hz Coh2: ' num2str(round(data.Coherr.Awake.gammaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.Awake.gammaBandPower.stdC01,2))]); disp(' ')
    disp(['Sleep Gamma 0.1 Hz Coh2: ' num2str(round(data.Coherr.Sleep.gammaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.Sleep.gammaBandPower.stdC01,2))]); disp(' ')
    disp(['All   Gamma 0.1 Hz Coh2: ' num2str(round(data.Coherr.All.gammaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.All.gammaBandPower.stdC01,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % gamma-band 0.01 Hz coherence^2 statistical diary
    disp('======================================================================================================================')
    disp('[7-S1d] Generalized linear mixed-effects model statistics for gamma-band coherence^2 @ 0.01 Hz for Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Gamma_Coh001_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake Gamma 0.01 Hz Coh2: ' num2str(round(data.Coherr.Awake.gammaBandPower.meanC001,2)) ' +/- ' num2str(round(data.Coherr.Awake.gammaBandPower.stdC001,2))]); disp(' ')
    disp(['Sleep Gamma 0.01 Hz Coh2: ' num2str(round(data.Coherr.Sleep.gammaBandPower.meanC001,2)) ' +/- ' num2str(round(data.Coherr.Sleep.gammaBandPower.stdC001,2))]); disp(' ')
    disp(['All   Gamma 0.01 Hz Coh2: ' num2str(round(data.Coherr.All.gammaBandPower.meanC001,2)) ' +/- ' num2str(round(data.Coherr.All.gammaBandPower.stdC001,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % HbT 0.1 Hz PSD statistical diary
    disp('======================================================================================================================')
    disp('[7-S1e] Generalized linear mixed-effects model statistics for HbT PSD @ 0.1 Hz for Rest, NREM, REM, Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(HbT_PSD01_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  [HbT] 0.1 Hz PSD: ' num2str(round(data.PowerSpec.Rest.CBV_HbT.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.Rest.CBV_HbT.stdS01,1))]); disp(' ')
    disp(['NREM  [HbT] 0.1 Hz PSD: ' num2str(round(data.PowerSpec.NREM.CBV_HbT.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.NREM.CBV_HbT.stdS01,1))]); disp(' ')
    disp(['REM   [HbT] 0.1 Hz PSD: ' num2str(round(data.PowerSpec.REM.CBV_HbT.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.REM.CBV_HbT.stdS01,1))]); disp(' ')
    disp(['Awake [HbT] 0.1 Hz PSD: ' num2str(round(data.PowerSpec.Awake.CBV_HbT.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.Awake.CBV_HbT.stdS01,1))]); disp(' ')
    disp(['Sleep [HbT] 0.1 Hz PSD: ' num2str(round(data.PowerSpec.Sleep.CBV_HbT.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.Sleep.CBV_HbT.stdS01,1))]); disp(' ')
    disp(['All   [HbT] 0.1 Hz PSD: ' num2str(round(data.PowerSpec.All.CBV_HbT.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.All.CBV_HbT.stdS01,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % HbT 0.01 Hz PSD statistical diary
    disp('======================================================================================================================')
    disp('[7-S1f] Generalized linear mixed-effects model statistics for HbT PSD @ 0.01 Hz for Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(HbT_PSD001_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake [HbT] 0.01 Hz PSD: ' num2str(round(data.PowerSpec.Awake.CBV_HbT.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.Awake.CBV_HbT.stdS001,1))]); disp(' ')
    disp(['Sleep [HbT] 0.01 Hz PSD: ' num2str(round(data.PowerSpec.Sleep.CBV_HbT.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.Sleep.CBV_HbT.stdS001,1))]); disp(' ')
    disp(['All   [HbT] 0.01 Hz PSD: ' num2str(round(data.PowerSpec.All.CBV_HbT.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.All.CBV_HbT.stdS001,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % HbT 0.1 Hz coherence^2 statistical diary
    disp('======================================================================================================================')
    disp('[7-S1g] Generalized linear mixed-effects model statistics for HbT coherence^2 @ 0.1 Hz for Rest, NREM, REM, Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(HbT_Coh01_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  [HbT] 0.1 Hz Coh2: ' num2str(round(data.Coherr.Rest.CBV_HbT.meanC01,2)) ' +/- ' num2str(round(data.Coherr.Rest.CBV_HbT.stdC01,2))]); disp(' ')
    disp(['NREM  [HbT] 0.1 Hz Coh2: ' num2str(round(data.Coherr.NREM.CBV_HbT.meanC01,2)) ' +/- ' num2str(round(data.Coherr.NREM.CBV_HbT.stdC01,2))]); disp(' ')
    disp(['REM   [HbT] 0.1 Hz Coh2: ' num2str(round(data.Coherr.REM.CBV_HbT.meanC01,2)) ' +/- ' num2str(round(data.Coherr.REM.CBV_HbT.stdC01,2))]); disp(' ')
    disp(['Awake [HbT] 0.1 Hz Coh2: ' num2str(round(data.Coherr.Awake.CBV_HbT.meanC01,2)) ' +/- ' num2str(round(data.Coherr.Awake.CBV_HbT.stdC01,2))]); disp(' ')
    disp(['Sleep [HbT] 0.1 Hz Coh2: ' num2str(round(data.Coherr.Sleep.CBV_HbT.meanC01,2)) ' +/- ' num2str(round(data.Coherr.Sleep.CBV_HbT.stdC01,2))]); disp(' ')
    disp(['All   [HbT] 0.1 Hz Coh2: ' num2str(round(data.Coherr.All.CBV_HbT.meanC01,2)) ' +/- ' num2str(round(data.Coherr.All.CBV_HbT.stdC01,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % HbT 0.01 Hz coherence^2 statistical diary
    disp('======================================================================================================================')
    disp('[7-S1h] Generalized linear mixed-effects model statistics for HbT coherence^2 @ 0.01 Hz for Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(HbT_Coh001_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake [HbT] 0.01 Hz Coh2: ' num2str(round(data.Coherr.Awake.CBV_HbT.meanC001,2)) ' +/- ' num2str(round(data.Coherr.Awake.CBV_HbT.stdC001,2))]); disp(' ')
    disp(['Sleep [HbT] 0.01 Hz Coh2: ' num2str(round(data.Coherr.Sleep.CBV_HbT.meanC001,2)) ' +/- ' num2str(round(data.Coherr.Sleep.CBV_HbT.stdC001,2))]); disp(' ')
    disp(['All   [HbT] 0.01 Hz Coh2: ' num2str(round(data.Coherr.All.CBV_HbT.meanC001,2)) ' +/- ' num2str(round(data.Coherr.All.CBV_HbT.stdC001,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % TwoP 0.1 Hz PSD statistical diary
    disp('======================================================================================================================')
    disp('[7-S1i] Generalized linear mixed-effects model statistics for D/D PSD @ 0.1 Hz for Rest, NREM, REM, Awake, and All')
    disp('======================================================================================================================')
    disp(TwoP_PSD01_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  [D/D] 0.1 Hz PSD: ' num2str(round(data.VesselPowerSpec.Rest.meanS01,1)) ' +/- ' num2str(round(data.VesselPowerSpec.Rest.stdS01,1))]); disp(' ')
    disp(['NREM  [D/D] 0.1 Hz PSD: ' num2str(round(data.VesselPowerSpec.NREM.meanS01,1)) ' +/- ' num2str(round(data.VesselPowerSpec.NREM.stdS01,1))]); disp(' ')
    disp(['REM   [D/D] 0.1 Hz PSD: ' num2str(round(data.VesselPowerSpec.REM.meanS01,1)) ' +/- ' num2str(round(data.VesselPowerSpec.REM.stdS01,1))]); disp(' ')
    disp(['Awake [D/D] 0.1 Hz PSD: ' num2str(round(data.VesselPowerSpec.Awake.meanS01,1)) ' +/- ' num2str(round(data.VesselPowerSpec.Awake.stdS01,1))]); disp(' ')
    disp(['All   [D/D] 0.1 Hz PSD: ' num2str(round(data.VesselPowerSpec.All.meanS01,1)) ' +/- ' num2str(round(data.VesselPowerSpec.All.stdS01,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % TwoP 0.01 Hz PSD statistical diary
    disp('======================================================================================================================')
    disp('[7-S1j] Generalized linear mixed-effects model statistics for D/D PSD @ 0.01 Hz for Awake and All')
    disp('======================================================================================================================')
    disp(TwoP_PSD001_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake [D/D] 0.01 Hz PSD: ' num2str(round(data.VesselPowerSpec.Awake.meanS001,1)) ' +/- ' num2str(round(data.VesselPowerSpec.Awake.stdS001,1))]); disp(' ')
    disp(['All   [D/D] 0.01 Hz PSD: ' num2str(round(data.VesselPowerSpec.All.meanS001,1)) ' +/- ' num2str(round(data.VesselPowerSpec.All.stdS001,1))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    diary off
    %% organized for supplemental table
    % variable names
    ColumnNames = {'Rest','NREM','REM','Awake','Sleep','All'};
    % gamma-band 0.1 Hz PSD power
    for aa = 1:length(ColumnNames)
        Gamma_S01_MeanStD{1,aa} = [num2str(round(data.PowerSpec.(ColumnNames{1,aa}).gammaBandPower.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.(ColumnNames{1,aa}).gammaBandPower.stdS01,1))]; %#ok<*AGROW>
    end
    % gamma-band 0.1 Hz PSD p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Rest') == true
            Gamma_S01_pVal{1,aa} = {' '};
        else
            Gamma_S01_pVal{1,aa} = ['p < ' num2str(Gamma_PSD01_Stats.Coefficients.pValue(aa,1))];
        end
    end
    % HbT 0.1 Hz PSD power
    for aa = 1:length(ColumnNames)
        HbT_S01_MeanStD{1,aa} = [num2str(round(data.PowerSpec.(ColumnNames{1,aa}).CBV_HbT.meanS01,1)) ' +/- ' num2str(round(data.PowerSpec.(ColumnNames{1,aa}).CBV_HbT.stdS01,1))];
    end
    % HbT 0.1 Hz PSD p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Rest') == true
            HbT_S01_pVal{1,aa} = {' '};
        else
            HbT_S01_pVal{1,aa} = ['p < ' num2str(HbT_PSD01_Stats.Coefficients.pValue(aa,1))];
        end
    end
    % D/D 0.1 Hz PSD power
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Sleep') == false
            TwoP_S01_MeanStD{1,aa} = [num2str(round(data.VesselPowerSpec.(ColumnNames{1,aa}).meanS01,1)) ' +/- ' num2str(round(data.VesselPowerSpec.(ColumnNames{1,aa}).stdS01,1))];
        else
            TwoP_S01_MeanStD{1,aa} = {' '};
        end
    end
    % D/D 0.1 Hz PSD p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Rest') == true
            TwoP_S01_pVal{1,aa} = {' '};
        elseif strcmp(ColumnNames{1,aa},'Sleep') == false && strcmp(ColumnNames{1,aa},'All') == false
            TwoP_S01_pVal{1,aa} = ['p < ' num2str(TwoP_PSD01_Stats.Coefficients.pValue(aa,1))];
        elseif strcmp(ColumnNames{1,aa},'All') == true
            TwoP_S01_pVal{1,aa} = ['p < ' num2str(TwoP_PSD01_Stats.Coefficients.pValue(aa - 1,1))];
        else
            TwoP_S01_pVal{1,aa} = {' '};
        end
    end
    % gamma-band 0.01 Hz PSD power
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Awake') == true || strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Gamma_S001_MeanStD{1,aa} = [num2str(round(data.PowerSpec.(ColumnNames{1,aa}).gammaBandPower.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.(ColumnNames{1,aa}).gammaBandPower.stdS001,1))];
        else
            Gamma_S001_MeanStD{1,aa} = {' '};
        end
    end
    % gamma-band 0.01 Hz PSD p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Gamma_S001_pVal{1,aa} = ['p < ' num2str(Gamma_PSD001_Stats.Coefficients.pValue(aa - 3,1))];
        else
            Gamma_S001_pVal{1,aa} = {' '};
        end
    end
    % HbT 0.01 Hz PSD power
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Awake') == true || strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            HbT_S001_MeanStD{1,aa} = [num2str(round(data.PowerSpec.(ColumnNames{1,aa}).CBV_HbT.meanS001,1)) ' +/- ' num2str(round(data.PowerSpec.(ColumnNames{1,aa}).CBV_HbT.stdS001,1))];
        else
            HbT_S001_MeanStD{1,aa} = {' '};
        end
    end
    % HbT 0.01 Hz PSD p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            HbT_S001_pVal{1,aa} = ['p < ' num2str(HbT_PSD001_Stats.Coefficients.pValue(aa - 3,1))];
        else
            HbT_S001_pVal{1,aa} = {' '};
        end
    end
    % D/D 0.01 Hz PSD power
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Awake') == true || strcmp(ColumnNames{1,aa},'All') == true
            TwoP_S001_MeanStD{1,aa} = [num2str(round(data.VesselPowerSpec.(ColumnNames{1,aa}).meanS001,1)) ' +/- ' num2str(round(data.VesselPowerSpec.(ColumnNames{1,aa}).stdS001,1))];
        else
            TwoP_S001_MeanStD{1,aa} = {' '};
        end
    end
    % D/D 0.01 Hz PSD p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'All') == true
            TwoP_S001_pVal{1,aa} = ['p < ' num2str(TwoP_PSD001_Stats.Coefficients.pValue(aa - 4,1))];
        else
            TwoP_S001_pVal{1,aa} = {' '};
        end
    end
    % gamma-band 0.1 Hz Coh2 power
    for aa = 1:length(ColumnNames)
        Gamma_C01_MeanStD{1,aa} = [num2str(round(data.Coherr.(ColumnNames{1,aa}).gammaBandPower.meanC01,2)) ' +/- ' num2str(round(data.Coherr.(ColumnNames{1,aa}).gammaBandPower.stdC01,2))]; %#ok<*AGROW>
    end
    % gamma-band 0.1 Hz Coh2 p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Rest') == true
            Gamma_C01_pVal{1,aa} = {' '};
        else
            Gamma_C01_pVal{1,aa} = ['p < ' num2str(Gamma_Coh01_Stats.Coefficients.pValue(aa,1))];
        end
    end
    % HbT 0.1 Hz Coh2 power
    for aa = 1:length(ColumnNames)
        HbT_C01_MeanStD{1,aa} = [num2str(round(data.Coherr.(ColumnNames{1,aa}).CBV_HbT.meanC01,2)) ' +/- ' num2str(round(data.Coherr.(ColumnNames{1,aa}).CBV_HbT.stdC01,2))];
    end
    % HbT 0.1 Hz Coh2 p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Rest') == true
            HbT_C01_pVal{1,aa} = {' '};
        else
            HbT_C01_pVal{1,aa} = ['p < ' num2str(HbT_Coh01_Stats.Coefficients.pValue(aa,1))];
        end
    end
    % gamma-band 0.01 Hz Coh2 power
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Awake') == true || strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Gamma_C001_MeanStD{1,aa} = [num2str(round(data.Coherr.(ColumnNames{1,aa}).gammaBandPower.meanC001,2)) ' +/- ' num2str(round(data.Coherr.(ColumnNames{1,aa}).gammaBandPower.stdC001,2))];
        else
            Gamma_C001_MeanStD{1,aa} = {' '};
        end
    end
    % gamma-band 0.01 Hz Coh2 p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Gamma_C001_pVal{1,aa} = ['p < ' num2str(Gamma_Coh001_Stats.Coefficients.pValue(aa - 3,1))];
        else
            Gamma_C001_pVal{1,aa} = {' '};
        end
    end
    % HbT 0.01 Hz Coh2 power
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Awake') == true || strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            HbT_C001_MeanStD{1,aa} = [num2str(round(data.Coherr.(ColumnNames{1,aa}).CBV_HbT.meanC001,2)) ' +/- ' num2str(round(data.Coherr.(ColumnNames{1,aa}).CBV_HbT.stdC001,2))];
        else
            HbT_C001_MeanStD{1,aa} = {' '};
        end
    end
    % HbT 0.01 Hz Coh2 p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            HbT_C001_pVal{1,aa} = ['p < ' num2str(HbT_Coh001_Stats.Coefficients.pValue(aa - 3,1))];
        else
            HbT_C001_pVal{1,aa} = {' '};
        end
    end
    %% save table data
    if isfield(AnalysisResults,'PSD') == false
        AnalysisResults.PSD = [];
    end
    if isfield(AnalysisResults.PSD,'gammaBandPower') == false
        AnalysisResults.PSD.columnNames = ColumnNames;
        AnalysisResults.PSD.gammaBandPower.meanStD01 = Gamma_S01_MeanStD;
        AnalysisResults.PSD.gammaBandPower.p01 = Gamma_S01_pVal;
        AnalysisResults.PSD.CBV_HbT.meanStD01 = HbT_S01_MeanStD;
        AnalysisResults.PSD.CBV_HbT.p01 = HbT_S01_pVal;
        AnalysisResults.PSD.TwoP.meanStD01 = TwoP_S01_MeanStD;
        AnalysisResults.PSD.TwoP.p01 = TwoP_S01_pVal;
        AnalysisResults.PSD.gammaBandPower.meanStD001 = Gamma_S001_MeanStD;
        AnalysisResults.PSD.gammaBandPower.p001 = Gamma_S001_pVal;
        AnalysisResults.PSD.CBV_HbT.meanStD001 = HbT_S001_MeanStD;
        AnalysisResults.PSD.CBV_HbT.p001 = HbT_S001_pVal;
        AnalysisResults.PSD.TwoP.meanStD001 = TwoP_S001_MeanStD;
        AnalysisResults.PSD.TwoP.p001 = TwoP_S001_pVal;
        AnalysisResults.Coherr.columnNames = ColumnNames;
        AnalysisResults.Coherr.gammaBandPower.meanStD01 = Gamma_C01_MeanStD;
        AnalysisResults.Coherr.gammaBandPower.p01 = Gamma_C01_pVal;
        AnalysisResults.Coherr.CBV_HbT.meanStD01 = HbT_C01_MeanStD;
        AnalysisResults.Coherr.CBV_HbT.p01 = HbT_C01_pVal;
        AnalysisResults.Coherr.gammaBandPower.meanStD001 = Gamma_C001_MeanStD;
        AnalysisResults.Coherr.gammaBandPower.p001 = Gamma_C001_pVal;
        AnalysisResults.Coherr.CBV_HbT.meanStD001 = HbT_C001_MeanStD;
        AnalysisResults.Coherr.CBV_HbT.p001 = HbT_C001_pVal;
        cd(rootFolder)
        save('AnalysisResults.mat','AnalysisResults')
    end
end

end
