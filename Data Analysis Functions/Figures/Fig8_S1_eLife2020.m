function [AnalysisResults] = Fig8_S1_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Generate figure panel 8-S1 for Turner_Gheres_Proctor_Drew_eLife2020
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
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
behavFields = {'Rest','NREM','REM','Awake','Sleep','All'};
dataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
%% average coherence during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.NeuralHemoCoherence = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        % create the behavior folder for the first iteration of the loop
        if isfield(data.NeuralHemoCoherence,behavField) == false
            data.NeuralHemoCoherence.(behavField) = [];
        end
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.C) == false
                % create the data type folder for the first iteration of the loop
                if isfield(data.NeuralHemoCoherence.(behavField),dataType) == false
                    data.NeuralHemoCoherence.(behavField).(dataType).C = [];
                    data.NeuralHemoCoherence.(behavField).(dataType).f = [];
                    data.NeuralHemoCoherence.(behavField).(dataType).confC = [];
                    data.NeuralHemoCoherence.(behavField).(dataType).animalID = {};
                    data.NeuralHemoCoherence.(behavField).(dataType).behavior = {};
                    data.NeuralHemoCoherence.(behavField).(dataType).hemisphere = {};
                end
                % concatenate C/f for existing data - exclude any empty sets
                data.NeuralHemoCoherence.(behavField).(dataType).C = cat(2,data.NeuralHemoCoherence.(behavField).(dataType).C,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.C,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjRH.C);
                data.NeuralHemoCoherence.(behavField).(dataType).f = cat(1,data.NeuralHemoCoherence.(behavField).(dataType).f,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.f,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjRH.f);
                data.NeuralHemoCoherence.(behavField).(dataType).confC = cat(1,data.NeuralHemoCoherence.(behavField).(dataType).confC,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.confC,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjRH.confC);
                if isempty(AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.C) == false
                    data.NeuralHemoCoherence.(behavField).(dataType).animalID = cat(1,data.NeuralHemoCoherence.(behavField).(dataType).animalID,animalID,animalID);
                    data.NeuralHemoCoherence.(behavField).(dataType).behavior = cat(1,data.NeuralHemoCoherence.(behavField).(dataType).behavior,behavField,behavField);
                    data.NeuralHemoCoherence.(behavField).(dataType).hemisphere = cat(1,data.NeuralHemoCoherence.(behavField).(dataType).hemisphere,'LH','RH');
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
        for ff = 1:size(data.NeuralHemoCoherence.(behavField).(dataType).C,2)
            if strcmp(behavField,'Rest') == true
                f_short = data.NeuralHemoCoherence.(behavField).(dataType).f(ff,:);
                C = data.NeuralHemoCoherence.(behavField).(dataType).C(:,ff);
                f_long = 0:0.01:0.5;
                C_long = interp1(f_short,C,f_long);
                index01 = find(f_long == 0.1);
                data.NeuralHemoCoherence.(behavField).(dataType).C01(ff,1) = C_long(index01).^2; %#ok<FNDSB>
            elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
                F = round(data.NeuralHemoCoherence.(behavField).(dataType).f(ff,:),2);
                C = data.NeuralHemoCoherence.(behavField).(dataType).C(:,ff);
                index01 = find(F == 0.1);
                data.NeuralHemoCoherence.(behavField).(dataType).C01(ff,1) = C(index01(1)).^2;
            else
                F = round(data.NeuralHemoCoherence.(behavField).(dataType).f(ff,:),3);
                C = data.NeuralHemoCoherence.(behavField).(dataType).C(:,ff);
                index01 = find(F == 0.1);
                index001 = find(F == 0.01);
                data.NeuralHemoCoherence.(behavField).(dataType).C01(ff,1) = C(index01(1)).^2;
                data.NeuralHemoCoherence.(behavField).(dataType).C001(ff,1) = C(index001(1)).^2;
            end
        end
    end
end
% take mean/StD of C/f and determine confC line
for gg = 1:length(behavFields)
    behavField = behavFields{1,gg};
    for hh = 1:length(dataTypes)
        dataType = dataTypes{1,hh};
        data.NeuralHemoCoherence.(behavField).(dataType).meanC = mean(data.NeuralHemoCoherence.(behavField).(dataType).C,2);
        data.NeuralHemoCoherence.(behavField).(dataType).stdC = std(data.NeuralHemoCoherence.(behavField).(dataType).C,0,2);
        data.NeuralHemoCoherence.(behavField).(dataType).meanf = mean(data.NeuralHemoCoherence.(behavField).(dataType).f,1);
        data.NeuralHemoCoherence.(behavField).(dataType).maxConfC = geomean(data.NeuralHemoCoherence.(behavField).(dataType).confC);
        data.NeuralHemoCoherence.(behavField).(dataType).maxConfC_Y = ones(length(data.NeuralHemoCoherence.(behavField).(dataType).meanf),1)*data.NeuralHemoCoherence.(behavField).(dataType).maxConfC;
        if strcmp(behavField,'Rest') == true || strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
            data.NeuralHemoCoherence.(behavField).(dataType).meanC01 = mean(data.NeuralHemoCoherence.(behavField).(dataType).C01,1);
            data.NeuralHemoCoherence.(behavField).(dataType).stdC01 = std(data.NeuralHemoCoherence.(behavField).(dataType).C01,0,1);
        else
            data.NeuralHemoCoherence.(behavField).(dataType).meanC01 = mean(data.NeuralHemoCoherence.(behavField).(dataType).C01,1);
            data.NeuralHemoCoherence.(behavField).(dataType).stdC01 = std(data.NeuralHemoCoherence.(behavField).(dataType).C01,0,1);
            data.NeuralHemoCoherence.(behavField).(dataType).meanC001 = mean(data.NeuralHemoCoherence.(behavField).(dataType).C001,1);
            data.NeuralHemoCoherence.(behavField).(dataType).stdC001 = std(data.NeuralHemoCoherence.(behavField).(dataType).C001,0,1);
        end
    end
end
%% statistics - generalized linear mixed effects model
% delta PSD @ 0.1 Hz
Delta_Coh01_tableSize = cat(1,data.NeuralHemoCoherence.Rest.deltaBandPower.C01,data.NeuralHemoCoherence.NREM.deltaBandPower.C01,data.NeuralHemoCoherence.REM.deltaBandPower.C01,...
    data.NeuralHemoCoherence.Awake.deltaBandPower.C01,data.NeuralHemoCoherence.Sleep.deltaBandPower.C01,data.NeuralHemoCoherence.All.deltaBandPower.C01);
Delta_Coh01_Table = table('Size',[size(Delta_Coh01_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Hemisphere','Behavior','Coh01'});
Delta_Coh01_Table.Mouse = cat(1,data.NeuralHemoCoherence.Rest.deltaBandPower.animalID,data.NeuralHemoCoherence.NREM.deltaBandPower.animalID,data.NeuralHemoCoherence.REM.deltaBandPower.animalID,...
    data.NeuralHemoCoherence.Awake.deltaBandPower.animalID,data.NeuralHemoCoherence.Sleep.deltaBandPower.animalID,data.NeuralHemoCoherence.All.deltaBandPower.animalID);
Delta_Coh01_Table.Behavior = cat(1,data.NeuralHemoCoherence.Rest.deltaBandPower.behavior,data.NeuralHemoCoherence.NREM.deltaBandPower.behavior,data.NeuralHemoCoherence.REM.deltaBandPower.behavior,...
    data.NeuralHemoCoherence.Awake.deltaBandPower.behavior,data.NeuralHemoCoherence.Sleep.deltaBandPower.behavior,data.NeuralHemoCoherence.All.deltaBandPower.behavior);
Delta_Coh01_Table.Hemisphere = cat(1,data.NeuralHemoCoherence.Rest.deltaBandPower.hemisphere,data.NeuralHemoCoherence.NREM.deltaBandPower.hemisphere,data.NeuralHemoCoherence.REM.deltaBandPower.hemisphere,...
    data.NeuralHemoCoherence.Awake.deltaBandPower.hemisphere,data.NeuralHemoCoherence.Sleep.deltaBandPower.hemisphere,data.NeuralHemoCoherence.All.deltaBandPower.hemisphere);
Delta_Coh01_Table.Coh01 = cat(1,data.NeuralHemoCoherence.Rest.deltaBandPower.C01,data.NeuralHemoCoherence.NREM.deltaBandPower.C01,data.NeuralHemoCoherence.REM.deltaBandPower.C01,...
    data.NeuralHemoCoherence.Awake.deltaBandPower.C01,data.NeuralHemoCoherence.Sleep.deltaBandPower.C01,data.NeuralHemoCoherence.All.deltaBandPower.C01);
Delta_Coh01_FitFormula = 'Coh01 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
Delta_Coh01_Stats = fitglme(Delta_Coh01_Table,Delta_Coh01_FitFormula);
% delta PSD @ 0.01 Hz
Delta_Coh001_tableSize = cat(1,data.NeuralHemoCoherence.Awake.deltaBandPower.C001,data.NeuralHemoCoherence.Sleep.deltaBandPower.C001,data.NeuralHemoCoherence.All.deltaBandPower.C001);
Delta_Coh001_Table = table('Size',[size(Delta_Coh001_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Hemisphere','Behavior','Coh001'});
Delta_Coh001_Table.Mouse = cat(1,data.NeuralHemoCoherence.Awake.deltaBandPower.animalID,data.NeuralHemoCoherence.Sleep.deltaBandPower.animalID,data.NeuralHemoCoherence.All.deltaBandPower.animalID);
Delta_Coh001_Table.Behavior = cat(1,data.NeuralHemoCoherence.Awake.deltaBandPower.behavior,data.NeuralHemoCoherence.Sleep.deltaBandPower.behavior,data.NeuralHemoCoherence.All.deltaBandPower.behavior);
Delta_Coh001_Table.Hemisphere = cat(1,data.NeuralHemoCoherence.Awake.deltaBandPower.hemisphere,data.NeuralHemoCoherence.Sleep.deltaBandPower.hemisphere,data.NeuralHemoCoherence.All.deltaBandPower.hemisphere);
Delta_Coh001_Table.Coh001 = cat(1,data.NeuralHemoCoherence.Awake.deltaBandPower.C001,data.NeuralHemoCoherence.Sleep.deltaBandPower.C001,data.NeuralHemoCoherence.All.deltaBandPower.C001);
Delta_Coh001_FitFormula = 'Coh001 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
Delta_Coh001_Stats = fitglme(Delta_Coh001_Table,Delta_Coh001_FitFormula);
% theta PSD @ 0.1 Hz
Theta_Coh01_tableSize = cat(1,data.NeuralHemoCoherence.Rest.thetaBandPower.C01,data.NeuralHemoCoherence.NREM.thetaBandPower.C01,data.NeuralHemoCoherence.REM.thetaBandPower.C01,...
    data.NeuralHemoCoherence.Awake.thetaBandPower.C01,data.NeuralHemoCoherence.Sleep.thetaBandPower.C01,data.NeuralHemoCoherence.All.thetaBandPower.C01);
Theta_Coh01_Table = table('Size',[size(Theta_Coh01_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Hemisphere','Behavior','Coh01'});
Theta_Coh01_Table.Mouse = cat(1,data.NeuralHemoCoherence.Rest.thetaBandPower.animalID,data.NeuralHemoCoherence.NREM.thetaBandPower.animalID,data.NeuralHemoCoherence.REM.thetaBandPower.animalID,...
    data.NeuralHemoCoherence.Awake.thetaBandPower.animalID,data.NeuralHemoCoherence.Sleep.thetaBandPower.animalID,data.NeuralHemoCoherence.All.thetaBandPower.animalID);
Theta_Coh01_Table.Behavior = cat(1,data.NeuralHemoCoherence.Rest.thetaBandPower.behavior,data.NeuralHemoCoherence.NREM.thetaBandPower.behavior,data.NeuralHemoCoherence.REM.thetaBandPower.behavior,...
    data.NeuralHemoCoherence.Awake.thetaBandPower.behavior,data.NeuralHemoCoherence.Sleep.thetaBandPower.behavior,data.NeuralHemoCoherence.All.thetaBandPower.behavior);
Theta_Coh01_Table.Hemisphere = cat(1,data.NeuralHemoCoherence.Rest.thetaBandPower.hemisphere,data.NeuralHemoCoherence.NREM.thetaBandPower.hemisphere,data.NeuralHemoCoherence.REM.thetaBandPower.hemisphere,...
    data.NeuralHemoCoherence.Awake.thetaBandPower.hemisphere,data.NeuralHemoCoherence.Sleep.thetaBandPower.hemisphere,data.NeuralHemoCoherence.All.thetaBandPower.hemisphere);
Theta_Coh01_Table.Coh01 = cat(1,data.NeuralHemoCoherence.Rest.thetaBandPower.C01,data.NeuralHemoCoherence.NREM.thetaBandPower.C01,data.NeuralHemoCoherence.REM.thetaBandPower.C01,...
    data.NeuralHemoCoherence.Awake.thetaBandPower.C01,data.NeuralHemoCoherence.Sleep.thetaBandPower.C01,data.NeuralHemoCoherence.All.thetaBandPower.C01);
Theta_Coh01_FitFormula = 'Coh01 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
Theta_Coh01_Stats = fitglme(Theta_Coh01_Table,Theta_Coh01_FitFormula);
% theta PSD @ 0.01 Hz
Theta_Coh001_tableSize = cat(1,data.NeuralHemoCoherence.Awake.thetaBandPower.C001,data.NeuralHemoCoherence.Sleep.thetaBandPower.C001,data.NeuralHemoCoherence.All.thetaBandPower.C001);
Theta_Coh001_Table = table('Size',[size(Theta_Coh001_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Hemisphere','Behavior','Coh001'});
Theta_Coh001_Table.Mouse = cat(1,data.NeuralHemoCoherence.Awake.thetaBandPower.animalID,data.NeuralHemoCoherence.Sleep.thetaBandPower.animalID,data.NeuralHemoCoherence.All.thetaBandPower.animalID);
Theta_Coh001_Table.Behavior = cat(1,data.NeuralHemoCoherence.Awake.thetaBandPower.behavior,data.NeuralHemoCoherence.Sleep.thetaBandPower.behavior,data.NeuralHemoCoherence.All.thetaBandPower.behavior);
Theta_Coh001_Table.Hemisphere = cat(1,data.NeuralHemoCoherence.Awake.thetaBandPower.hemisphere,data.NeuralHemoCoherence.Sleep.thetaBandPower.hemisphere,data.NeuralHemoCoherence.All.thetaBandPower.hemisphere);
Theta_Coh001_Table.Coh001 = cat(1,data.NeuralHemoCoherence.Awake.thetaBandPower.C001,data.NeuralHemoCoherence.Sleep.thetaBandPower.C001,data.NeuralHemoCoherence.All.thetaBandPower.C001);
Theta_Coh001_FitFormula = 'Coh001 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
Theta_Coh001_Stats = fitglme(Theta_Coh001_Table,Theta_Coh001_FitFormula);
% alpha PSD @ 0.1 Hz
Alpha_Coh01_tableSize = cat(1,data.NeuralHemoCoherence.Rest.alphaBandPower.C01,data.NeuralHemoCoherence.NREM.alphaBandPower.C01,data.NeuralHemoCoherence.REM.alphaBandPower.C01,...
    data.NeuralHemoCoherence.Awake.alphaBandPower.C01,data.NeuralHemoCoherence.Sleep.alphaBandPower.C01,data.NeuralHemoCoherence.All.alphaBandPower.C01);
Alpha_Coh01_Table = table('Size',[size(Alpha_Coh01_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Hemisphere','Behavior','Coh01'});
Alpha_Coh01_Table.Mouse = cat(1,data.NeuralHemoCoherence.Rest.alphaBandPower.animalID,data.NeuralHemoCoherence.NREM.alphaBandPower.animalID,data.NeuralHemoCoherence.REM.alphaBandPower.animalID,...
    data.NeuralHemoCoherence.Awake.alphaBandPower.animalID,data.NeuralHemoCoherence.Sleep.alphaBandPower.animalID,data.NeuralHemoCoherence.All.alphaBandPower.animalID);
Alpha_Coh01_Table.Behavior = cat(1,data.NeuralHemoCoherence.Rest.alphaBandPower.behavior,data.NeuralHemoCoherence.NREM.alphaBandPower.behavior,data.NeuralHemoCoherence.REM.alphaBandPower.behavior,...
    data.NeuralHemoCoherence.Awake.alphaBandPower.behavior,data.NeuralHemoCoherence.Sleep.alphaBandPower.behavior,data.NeuralHemoCoherence.All.alphaBandPower.behavior);
Alpha_Coh01_Table.Hemisphere = cat(1,data.NeuralHemoCoherence.Rest.alphaBandPower.hemisphere,data.NeuralHemoCoherence.NREM.alphaBandPower.hemisphere,data.NeuralHemoCoherence.REM.alphaBandPower.hemisphere,...
    data.NeuralHemoCoherence.Awake.alphaBandPower.hemisphere,data.NeuralHemoCoherence.Sleep.alphaBandPower.hemisphere,data.NeuralHemoCoherence.All.alphaBandPower.hemisphere);
Alpha_Coh01_Table.Coh01 = cat(1,data.NeuralHemoCoherence.Rest.alphaBandPower.C01,data.NeuralHemoCoherence.NREM.alphaBandPower.C01,data.NeuralHemoCoherence.REM.alphaBandPower.C01,...
    data.NeuralHemoCoherence.Awake.alphaBandPower.C01,data.NeuralHemoCoherence.Sleep.alphaBandPower.C01,data.NeuralHemoCoherence.All.alphaBandPower.C01);
Alpha_Coh01_FitFormula = 'Coh01 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
Alpha_Coh01_Stats = fitglme(Alpha_Coh01_Table,Alpha_Coh01_FitFormula);
% alpha PSD @ 0.01 Hz
Alpha_Coh001_tableSize = cat(1,data.NeuralHemoCoherence.Awake.alphaBandPower.C001,data.NeuralHemoCoherence.Sleep.alphaBandPower.C001,data.NeuralHemoCoherence.All.alphaBandPower.C001);
Alpha_Coh001_Table = table('Size',[size(Alpha_Coh001_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Hemisphere','Behavior','Coh001'});
Alpha_Coh001_Table.Mouse = cat(1,data.NeuralHemoCoherence.Awake.alphaBandPower.animalID,data.NeuralHemoCoherence.Sleep.alphaBandPower.animalID,data.NeuralHemoCoherence.All.alphaBandPower.animalID);
Alpha_Coh001_Table.Behavior = cat(1,data.NeuralHemoCoherence.Awake.alphaBandPower.behavior,data.NeuralHemoCoherence.Sleep.alphaBandPower.behavior,data.NeuralHemoCoherence.All.alphaBandPower.behavior);
Alpha_Coh001_Table.Hemisphere = cat(1,data.NeuralHemoCoherence.Awake.alphaBandPower.hemisphere,data.NeuralHemoCoherence.Sleep.alphaBandPower.hemisphere,data.NeuralHemoCoherence.All.alphaBandPower.hemisphere);
Alpha_Coh001_Table.Coh001 = cat(1,data.NeuralHemoCoherence.Awake.alphaBandPower.C001,data.NeuralHemoCoherence.Sleep.alphaBandPower.C001,data.NeuralHemoCoherence.All.alphaBandPower.C001);
Alpha_Coh001_FitFormula = 'Coh001 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
Alpha_Coh001_Stats = fitglme(Alpha_Coh001_Table,Alpha_Coh001_FitFormula);
% beta PSD @ 0.1 Hz
Beta_Coh01_tableSize = cat(1,data.NeuralHemoCoherence.Rest.betaBandPower.C01,data.NeuralHemoCoherence.NREM.betaBandPower.C01,data.NeuralHemoCoherence.REM.betaBandPower.C01,...
    data.NeuralHemoCoherence.Awake.betaBandPower.C01,data.NeuralHemoCoherence.Sleep.betaBandPower.C01,data.NeuralHemoCoherence.All.betaBandPower.C01);
Beta_Coh01_Table = table('Size',[size(Beta_Coh01_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Hemisphere','Behavior','Coh01'});
Beta_Coh01_Table.Mouse = cat(1,data.NeuralHemoCoherence.Rest.betaBandPower.animalID,data.NeuralHemoCoherence.NREM.betaBandPower.animalID,data.NeuralHemoCoherence.REM.betaBandPower.animalID,...
    data.NeuralHemoCoherence.Awake.betaBandPower.animalID,data.NeuralHemoCoherence.Sleep.betaBandPower.animalID,data.NeuralHemoCoherence.All.betaBandPower.animalID);
Beta_Coh01_Table.Behavior = cat(1,data.NeuralHemoCoherence.Rest.betaBandPower.behavior,data.NeuralHemoCoherence.NREM.betaBandPower.behavior,data.NeuralHemoCoherence.REM.betaBandPower.behavior,...
    data.NeuralHemoCoherence.Awake.betaBandPower.behavior,data.NeuralHemoCoherence.Sleep.betaBandPower.behavior,data.NeuralHemoCoherence.All.betaBandPower.behavior);
Beta_Coh01_Table.Hemisphere = cat(1,data.NeuralHemoCoherence.Rest.betaBandPower.hemisphere,data.NeuralHemoCoherence.NREM.betaBandPower.hemisphere,data.NeuralHemoCoherence.REM.betaBandPower.hemisphere,...
    data.NeuralHemoCoherence.Awake.betaBandPower.hemisphere,data.NeuralHemoCoherence.Sleep.betaBandPower.hemisphere,data.NeuralHemoCoherence.All.betaBandPower.hemisphere);
Beta_Coh01_Table.Coh01 = cat(1,data.NeuralHemoCoherence.Rest.betaBandPower.C01,data.NeuralHemoCoherence.NREM.betaBandPower.C01,data.NeuralHemoCoherence.REM.betaBandPower.C01,...
    data.NeuralHemoCoherence.Awake.betaBandPower.C01,data.NeuralHemoCoherence.Sleep.betaBandPower.C01,data.NeuralHemoCoherence.All.betaBandPower.C01);
Beta_Coh01_FitFormula = 'Coh01 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
Beta_Coh01_Stats = fitglme(Beta_Coh01_Table,Beta_Coh01_FitFormula);
% beta PSD @ 0.01 Hz
Beta_Coh001_tableSize = cat(1,data.NeuralHemoCoherence.Awake.betaBandPower.C001,data.NeuralHemoCoherence.Sleep.betaBandPower.C001,data.NeuralHemoCoherence.All.betaBandPower.C001);
Beta_Coh001_Table = table('Size',[size(Beta_Coh001_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Hemisphere','Behavior','Coh001'});
Beta_Coh001_Table.Mouse = cat(1,data.NeuralHemoCoherence.Awake.betaBandPower.animalID,data.NeuralHemoCoherence.Sleep.betaBandPower.animalID,data.NeuralHemoCoherence.All.betaBandPower.animalID);
Beta_Coh001_Table.Behavior = cat(1,data.NeuralHemoCoherence.Awake.betaBandPower.behavior,data.NeuralHemoCoherence.Sleep.betaBandPower.behavior,data.NeuralHemoCoherence.All.betaBandPower.behavior);
Beta_Coh001_Table.Hemisphere = cat(1,data.NeuralHemoCoherence.Awake.betaBandPower.hemisphere,data.NeuralHemoCoherence.Sleep.betaBandPower.hemisphere,data.NeuralHemoCoherence.All.betaBandPower.hemisphere);
Beta_Coh001_Table.Coh001 = cat(1,data.NeuralHemoCoherence.Awake.betaBandPower.C001,data.NeuralHemoCoherence.Sleep.betaBandPower.C001,data.NeuralHemoCoherence.All.betaBandPower.C001);
Beta_Coh001_FitFormula = 'Coh001 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
Beta_Coh001_Stats = fitglme(Beta_Coh001_Table,Beta_Coh001_FitFormula);
% gamma PSD @ 0.1 Hz
Gamma_Coh01_tableSize = cat(1,data.NeuralHemoCoherence.Rest.gammaBandPower.C01,data.NeuralHemoCoherence.NREM.gammaBandPower.C01,data.NeuralHemoCoherence.REM.gammaBandPower.C01,...
    data.NeuralHemoCoherence.Awake.gammaBandPower.C01,data.NeuralHemoCoherence.Sleep.gammaBandPower.C01,data.NeuralHemoCoherence.All.gammaBandPower.C01);
Gamma_Coh01_Table = table('Size',[size(Gamma_Coh01_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Hemisphere','Behavior','Coh01'});
Gamma_Coh01_Table.Mouse = cat(1,data.NeuralHemoCoherence.Rest.gammaBandPower.animalID,data.NeuralHemoCoherence.NREM.gammaBandPower.animalID,data.NeuralHemoCoherence.REM.gammaBandPower.animalID,...
    data.NeuralHemoCoherence.Awake.gammaBandPower.animalID,data.NeuralHemoCoherence.Sleep.gammaBandPower.animalID,data.NeuralHemoCoherence.All.gammaBandPower.animalID);
Gamma_Coh01_Table.Behavior = cat(1,data.NeuralHemoCoherence.Rest.gammaBandPower.behavior,data.NeuralHemoCoherence.NREM.gammaBandPower.behavior,data.NeuralHemoCoherence.REM.gammaBandPower.behavior,...
    data.NeuralHemoCoherence.Awake.gammaBandPower.behavior,data.NeuralHemoCoherence.Sleep.gammaBandPower.behavior,data.NeuralHemoCoherence.All.gammaBandPower.behavior);
Gamma_Coh01_Table.Hemisphere = cat(1,data.NeuralHemoCoherence.Rest.gammaBandPower.hemisphere,data.NeuralHemoCoherence.NREM.gammaBandPower.hemisphere,data.NeuralHemoCoherence.REM.gammaBandPower.hemisphere,...
    data.NeuralHemoCoherence.Awake.gammaBandPower.hemisphere,data.NeuralHemoCoherence.Sleep.gammaBandPower.hemisphere,data.NeuralHemoCoherence.All.gammaBandPower.hemisphere);
Gamma_Coh01_Table.Coh01 = cat(1,data.NeuralHemoCoherence.Rest.gammaBandPower.C01,data.NeuralHemoCoherence.NREM.gammaBandPower.C01,data.NeuralHemoCoherence.REM.gammaBandPower.C01,...
    data.NeuralHemoCoherence.Awake.gammaBandPower.C01,data.NeuralHemoCoherence.Sleep.gammaBandPower.C01,data.NeuralHemoCoherence.All.gammaBandPower.C01);
Gamma_Coh01_FitFormula = 'Coh01 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
Gamma_Coh01_Stats = fitglme(Gamma_Coh01_Table,Gamma_Coh01_FitFormula);
% gamma PSD @ 0.01 Hz
Gamma_Coh001_tableSize = cat(1,data.NeuralHemoCoherence.Awake.gammaBandPower.C001,data.NeuralHemoCoherence.Sleep.gammaBandPower.C001,data.NeuralHemoCoherence.All.gammaBandPower.C001);
Gamma_Coh001_Table = table('Size',[size(Gamma_Coh001_tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Hemisphere','Behavior','Coh001'});
Gamma_Coh001_Table.Mouse = cat(1,data.NeuralHemoCoherence.Awake.gammaBandPower.animalID,data.NeuralHemoCoherence.Sleep.gammaBandPower.animalID,data.NeuralHemoCoherence.All.gammaBandPower.animalID);
Gamma_Coh001_Table.Behavior = cat(1,data.NeuralHemoCoherence.Awake.gammaBandPower.behavior,data.NeuralHemoCoherence.Sleep.gammaBandPower.behavior,data.NeuralHemoCoherence.All.gammaBandPower.behavior);
Gamma_Coh001_Table.Hemisphere = cat(1,data.NeuralHemoCoherence.Awake.gammaBandPower.hemisphere,data.NeuralHemoCoherence.Sleep.gammaBandPower.hemisphere,data.NeuralHemoCoherence.All.gammaBandPower.hemisphere);
Gamma_Coh001_Table.Coh001 = cat(1,data.NeuralHemoCoherence.Awake.gammaBandPower.C001,data.NeuralHemoCoherence.Sleep.gammaBandPower.C001,data.NeuralHemoCoherence.All.gammaBandPower.C001);
Gamma_Coh001_FitFormula = 'Coh001 ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
Gamma_Coh001_Stats = fitglme(Gamma_Coh001_Table,Gamma_Coh001_FitFormula);
%% Fig. 8-S1
summaryFigure = figure('Name','Fig8-S1 (a-o)');
sgtitle('Figure 8-S1 - Turner et al. 2020')
%% [8-S1a] coherence between delta-band power and HbT during different arousal-states
ax1 = subplot(5,3,1);
s1 = semilogx(data.NeuralHemoCoherence.Rest.deltaBandPower.meanf,data.NeuralHemoCoherence.Rest.deltaBandPower.meanC.^2,'color',colorRest,'LineWidth',2);
hold on
rectangle('Position',[0.005,0.05,0.1 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
s2 = semilogx(data.NeuralHemoCoherence.NREM.deltaBandPower.meanf,data.NeuralHemoCoherence.NREM.deltaBandPower.meanC.^2,'color',colorNREM,'LineWidth',2);
rectangle('Position',[0.005,0.05,1/30 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
s3 = semilogx(data.NeuralHemoCoherence.REM.deltaBandPower.meanf,data.NeuralHemoCoherence.REM.deltaBandPower.meanC.^2,'color',colorREM,'LineWidth',2);
rectangle('Position',[0.005,0.05,1/60 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
s4 = semilogx(data.NeuralHemoCoherence.Awake.deltaBandPower.meanf,data.NeuralHemoCoherence.Awake.deltaBandPower.meanC.^2,'color',colorAlert,'LineWidth',2);
s5 = semilogx(data.NeuralHemoCoherence.Sleep.deltaBandPower.meanf,data.NeuralHemoCoherence.Sleep.deltaBandPower.meanC.^2,'color',colorAsleep,'LineWidth',2);
s6 = semilogx(data.NeuralHemoCoherence.All.deltaBandPower.meanf,data.NeuralHemoCoherence.All.deltaBandPower.meanC.^2,'color',colorAll,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[8-S1a] Neural-hemo coherence^2','Delta-[HbT]'})
legend([s1,s2,s3,s4,s5,s6],'Rest','NREM','REM','Alert','Asleep','All','Location','SouthEast')
axis square
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [8-S1b] coherence between delta-band power and HbT during different arousal-states @ 0.1 Hz
ax2 = subplot(5,3,2);
scatter(ones(1,length(data.NeuralHemoCoherence.Rest.deltaBandPower.C01))*1,data.NeuralHemoCoherence.Rest.deltaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.NeuralHemoCoherence.Rest.deltaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.NREM.deltaBandPower.C01))*2,data.NeuralHemoCoherence.NREM.deltaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.NeuralHemoCoherence.NREM.deltaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.REM.deltaBandPower.C01))*3,data.NeuralHemoCoherence.REM.deltaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NeuralHemoCoherence.REM.deltaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.Awake.deltaBandPower.C01))*4,data.NeuralHemoCoherence.Awake.deltaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.NeuralHemoCoherence.Awake.deltaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.Sleep.deltaBandPower.C01))*5,data.NeuralHemoCoherence.Sleep.deltaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.NeuralHemoCoherence.Sleep.deltaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.All.deltaBandPower.C01))*6,data.NeuralHemoCoherence.All.deltaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.NeuralHemoCoherence.All.deltaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
title({'[8-S1b] Coherence^2 @ 0.1 Hz','Delta-[HbT]'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,7])
ylim([0,1])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [8-S1c] coherence between delta-band power and HbT during different arousal-states @ 0.01 Hz
ax3 = subplot(5,3,3);
scatter(ones(1,length(data.NeuralHemoCoherence.Awake.deltaBandPower.C001))*1,data.NeuralHemoCoherence.Awake.deltaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.NeuralHemoCoherence.Awake.deltaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.Sleep.deltaBandPower.C001))*2,data.NeuralHemoCoherence.Sleep.deltaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.NeuralHemoCoherence.Sleep.deltaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.All.deltaBandPower.C001))*3,data.NeuralHemoCoherence.All.deltaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NeuralHemoCoherence.All.deltaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[8-S1c] Coherence^2 @ 0.01 Hz','Delta-[HbT]'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,4])
ylim([0,1])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [8-S1d] coherence between theta-band power and HbT during different arousal-states
ax4 = subplot(5,3,4);
semilogx(data.NeuralHemoCoherence.Rest.thetaBandPower.meanf,data.NeuralHemoCoherence.Rest.thetaBandPower.meanC.^2,'color',colorRest,'LineWidth',2);
hold on
rectangle('Position',[0.005,0.05,0.1 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.NeuralHemoCoherence.NREM.thetaBandPower.meanf,data.NeuralHemoCoherence.NREM.thetaBandPower.meanC.^2,'color',colorNREM,'LineWidth',2);
rectangle('Position',[0.005,0.05,1/30 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.NeuralHemoCoherence.REM.thetaBandPower.meanf,data.NeuralHemoCoherence.REM.thetaBandPower.meanC.^2,'color',colorREM,'LineWidth',2);
rectangle('Position',[0.005,0.05,1/60 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.NeuralHemoCoherence.Awake.thetaBandPower.meanf,data.NeuralHemoCoherence.Awake.thetaBandPower.meanC.^2,'color',colorAlert,'LineWidth',2);
semilogx(data.NeuralHemoCoherence.Sleep.thetaBandPower.meanf,data.NeuralHemoCoherence.Sleep.thetaBandPower.meanC.^2,'color',colorAsleep,'LineWidth',2);
semilogx(data.NeuralHemoCoherence.All.thetaBandPower.meanf,data.NeuralHemoCoherence.All.thetaBandPower.meanC.^2,'color',colorAll,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[8-S1d] Neural-hemo coherence^2','Theta-[HbT]'})
axis square
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [8-S1e] coherence between theta-band power and HbT during different arousal-states @ 0.1 Hz
ax5 = subplot(5,3,5);
scatter(ones(1,length(data.NeuralHemoCoherence.Rest.thetaBandPower.C01))*1,data.NeuralHemoCoherence.Rest.thetaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.NeuralHemoCoherence.Rest.thetaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.NREM.thetaBandPower.C01))*2,data.NeuralHemoCoherence.NREM.thetaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.NeuralHemoCoherence.NREM.thetaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.REM.thetaBandPower.C01))*3,data.NeuralHemoCoherence.REM.thetaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NeuralHemoCoherence.REM.thetaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.Awake.thetaBandPower.C01))*4,data.NeuralHemoCoherence.Awake.thetaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.NeuralHemoCoherence.Awake.thetaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.Sleep.thetaBandPower.C01))*5,data.NeuralHemoCoherence.Sleep.thetaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.NeuralHemoCoherence.Sleep.thetaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.All.thetaBandPower.C01))*6,data.NeuralHemoCoherence.All.thetaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.NeuralHemoCoherence.All.thetaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
title({'[8-S1e] Coherence^2 @ 0.1 Hz','Theta-[HbT]'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,7])
ylim([0,1])
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [8-S1f] coherence between theta-band power and HbT during different arousal-states @ 0.01 Hz
ax6 = subplot(5,3,6);
scatter(ones(1,length(data.NeuralHemoCoherence.Awake.thetaBandPower.C001))*1,data.NeuralHemoCoherence.Awake.thetaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.NeuralHemoCoherence.Awake.thetaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.Sleep.thetaBandPower.C001))*2,data.NeuralHemoCoherence.Sleep.thetaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.NeuralHemoCoherence.Sleep.thetaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.All.thetaBandPower.C001))*3,data.NeuralHemoCoherence.All.thetaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NeuralHemoCoherence.All.thetaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[8-S1f] Coherence^2 @ 0.01 Hz','Theta-[HbT]'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,4])
ylim([0,1])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% [8-S1g] coherence between alpha-band power and HbT during different arousal-states
ax7 = subplot(5,3,7);
semilogx(data.NeuralHemoCoherence.Rest.alphaBandPower.meanf,data.NeuralHemoCoherence.Rest.alphaBandPower.meanC.^2,'color',colorRest,'LineWidth',2);
hold on
rectangle('Position',[0.005,0.05,0.1 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.NeuralHemoCoherence.NREM.alphaBandPower.meanf,data.NeuralHemoCoherence.NREM.alphaBandPower.meanC.^2,'color',colorNREM,'LineWidth',2);
rectangle('Position',[0.005,0.05,1/30 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.NeuralHemoCoherence.REM.alphaBandPower.meanf,data.NeuralHemoCoherence.REM.alphaBandPower.meanC.^2,'color',colorREM,'LineWidth',2);
rectangle('Position',[0.005,0.05,1/60 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.NeuralHemoCoherence.Awake.alphaBandPower.meanf,data.NeuralHemoCoherence.Awake.alphaBandPower.meanC.^2,'color',colorAlert,'LineWidth',2);
semilogx(data.NeuralHemoCoherence.Sleep.alphaBandPower.meanf,data.NeuralHemoCoherence.Sleep.alphaBandPower.meanC.^2,'color',colorAsleep,'LineWidth',2);
semilogx(data.NeuralHemoCoherence.All.alphaBandPower.meanf,data.NeuralHemoCoherence.All.alphaBandPower.meanC.^2,'color',colorAll,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[8-S1g] Neural-hemo coherence^2','Alpha-[HbT]'})
axis square
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
ax7.TickLength = [0.03,0.03];
%% [8-S1h] coherence between alpha-band power and HbT during different arousal-states @ 0.1 Hz
ax8 = subplot(5,3,8);
scatter(ones(1,length(data.NeuralHemoCoherence.Rest.alphaBandPower.C01))*1,data.NeuralHemoCoherence.Rest.alphaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.NeuralHemoCoherence.Rest.alphaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.NREM.alphaBandPower.C01))*2,data.NeuralHemoCoherence.NREM.alphaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.NeuralHemoCoherence.NREM.alphaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.REM.alphaBandPower.C01))*3,data.NeuralHemoCoherence.REM.alphaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NeuralHemoCoherence.REM.alphaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.Awake.alphaBandPower.C01))*4,data.NeuralHemoCoherence.Awake.alphaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.NeuralHemoCoherence.Awake.alphaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.Sleep.alphaBandPower.C01))*5,data.NeuralHemoCoherence.Sleep.alphaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.NeuralHemoCoherence.Sleep.alphaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.All.alphaBandPower.C01))*6,data.NeuralHemoCoherence.All.alphaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.NeuralHemoCoherence.All.alphaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
title({'[8-S1h] Coherence^2 @ 0.1 Hz','Alpha-[HbT]'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,7])
ylim([0,1])
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];
%% [8-S1i] coherence between alpha-band power and HbT during different arousal-states @ 0.01 Hz
ax9 = subplot(5,3,9);
scatter(ones(1,length(data.NeuralHemoCoherence.Awake.alphaBandPower.C001))*1,data.NeuralHemoCoherence.Awake.alphaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.NeuralHemoCoherence.Awake.alphaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.Sleep.alphaBandPower.C001))*2,data.NeuralHemoCoherence.Sleep.alphaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.NeuralHemoCoherence.Sleep.alphaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.All.alphaBandPower.C001))*3,data.NeuralHemoCoherence.All.alphaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NeuralHemoCoherence.All.alphaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[8-S1i] Coherence^2 @ 0.01 Hz','Alpha-[HbT]'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,4])
ylim([0,1])
set(gca,'box','off')
ax9.TickLength = [0.03,0.03];
%% [8-S1j] coherence between beta-band power and HbT during different arousal-states
ax10 = subplot(5,3,10);
semilogx(data.NeuralHemoCoherence.Rest.betaBandPower.meanf,data.NeuralHemoCoherence.Rest.betaBandPower.meanC.^2,'color',colorRest,'LineWidth',2);
hold on
rectangle('Position',[0.005,0.05,0.1 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.NeuralHemoCoherence.NREM.betaBandPower.meanf,data.NeuralHemoCoherence.NREM.betaBandPower.meanC.^2,'color',colorNREM,'LineWidth',2);
rectangle('Position',[0.005,0.05,1/30 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.NeuralHemoCoherence.REM.betaBandPower.meanf,data.NeuralHemoCoherence.REM.betaBandPower.meanC.^2,'color',colorREM,'LineWidth',2);
rectangle('Position',[0.005,0.05,1/60 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.NeuralHemoCoherence.Awake.betaBandPower.meanf,data.NeuralHemoCoherence.Awake.betaBandPower.meanC.^2,'color',colorAlert,'LineWidth',2);
semilogx(data.NeuralHemoCoherence.Sleep.betaBandPower.meanf,data.NeuralHemoCoherence.Sleep.betaBandPower.meanC.^2,'color',colorAsleep,'LineWidth',2);
semilogx(data.NeuralHemoCoherence.All.betaBandPower.meanf,data.NeuralHemoCoherence.All.betaBandPower.meanC.^2,'color',colorAll,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[8-S1j] Neural-hemo coherence^2','Beta-[HbT]'})
axis square
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
ax10.TickLength = [0.03,0.03];
%% [8-S1k] coherence between beta-band power and HbT during different arousal-states @ 0.1 Hz
ax11 = subplot(5,3,11);
scatter(ones(1,length(data.NeuralHemoCoherence.Rest.betaBandPower.C01))*1,data.NeuralHemoCoherence.Rest.betaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.NeuralHemoCoherence.Rest.betaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.NREM.betaBandPower.C01))*2,data.NeuralHemoCoherence.NREM.betaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.NeuralHemoCoherence.NREM.betaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.REM.betaBandPower.C01))*3,data.NeuralHemoCoherence.REM.betaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NeuralHemoCoherence.REM.betaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.Awake.betaBandPower.C01))*4,data.NeuralHemoCoherence.Awake.betaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.NeuralHemoCoherence.Awake.betaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.Sleep.betaBandPower.C01))*5,data.NeuralHemoCoherence.Sleep.betaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.NeuralHemoCoherence.Sleep.betaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.All.betaBandPower.C01))*6,data.NeuralHemoCoherence.All.betaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.NeuralHemoCoherence.All.betaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
title({'[8-S1k] Coherence^2 @ 0.1 Hz','Beta-[HbT]'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,7])
ylim([0,1])
set(gca,'box','off')
ax11.TickLength = [0.03,0.03];
%% [8-S1l] coherence between beta-band power and HbT during different arousal-states @ 0.01 Hz
ax12 = subplot(5,3,12);
scatter(ones(1,length(data.NeuralHemoCoherence.Awake.betaBandPower.C001))*1,data.NeuralHemoCoherence.Awake.betaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.NeuralHemoCoherence.Awake.betaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.Sleep.betaBandPower.C001))*2,data.NeuralHemoCoherence.Sleep.betaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.NeuralHemoCoherence.Sleep.betaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.All.betaBandPower.C001))*3,data.NeuralHemoCoherence.All.betaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NeuralHemoCoherence.All.betaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[8-S1l] Coherence^2 @ 0.01 Hz','Beta-[HbT]'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,4])
ylim([0,1])
set(gca,'box','off')
ax12.TickLength = [0.03,0.03];
%% [8-S1m] coherence between gamma-band power and HbT during different arousal-states
ax13 = subplot(5,3,13);
semilogx(data.NeuralHemoCoherence.Rest.gammaBandPower.meanf,data.NeuralHemoCoherence.Rest.gammaBandPower.meanC.^2,'color',colorRest,'LineWidth',2);
hold on
rectangle('Position',[0.005,0.05,0.1 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.NeuralHemoCoherence.NREM.gammaBandPower.meanf,data.NeuralHemoCoherence.NREM.gammaBandPower.meanC.^2,'color',colorNREM,'LineWidth',2);
rectangle('Position',[0.005,0.05,1/30 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.NeuralHemoCoherence.REM.gammaBandPower.meanf,data.NeuralHemoCoherence.REM.gammaBandPower.meanC.^2,'color',colorREM,'LineWidth',2);
rectangle('Position',[0.005,0.05,1/60 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.NeuralHemoCoherence.Awake.gammaBandPower.meanf,data.NeuralHemoCoherence.Awake.gammaBandPower.meanC.^2,'color',colorAlert,'LineWidth',2);
semilogx(data.NeuralHemoCoherence.Sleep.gammaBandPower.meanf,data.NeuralHemoCoherence.Sleep.gammaBandPower.meanC.^2,'color',colorAsleep,'LineWidth',2);
semilogx(data.NeuralHemoCoherence.All.gammaBandPower.meanf,data.NeuralHemoCoherence.All.gammaBandPower.meanC.^2,'color',colorAll,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[8-S1m] Neural-hemo coherence^2','Gamma-[HbT]'})
axis square
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
ax13.TickLength = [0.03,0.03];
%% [8-S1n] coherence between gamma-band power and HbT during different arousal-states @ 0.1 Hz
ax14 = subplot(5,3,14);
scatter(ones(1,length(data.NeuralHemoCoherence.Rest.gammaBandPower.C01))*1,data.NeuralHemoCoherence.Rest.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.NeuralHemoCoherence.Rest.gammaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.NREM.gammaBandPower.C01))*2,data.NeuralHemoCoherence.NREM.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.NeuralHemoCoherence.NREM.gammaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.REM.gammaBandPower.C01))*3,data.NeuralHemoCoherence.REM.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NeuralHemoCoherence.REM.gammaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.Awake.gammaBandPower.C01))*4,data.NeuralHemoCoherence.Awake.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.NeuralHemoCoherence.Awake.gammaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.Sleep.gammaBandPower.C01))*5,data.NeuralHemoCoherence.Sleep.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.NeuralHemoCoherence.Sleep.gammaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.All.gammaBandPower.C01))*6,data.NeuralHemoCoherence.All.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.NeuralHemoCoherence.All.gammaBandPower.meanC01,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
title({'[8-S1n] Coherence^2 @ 0.1 Hz','Gamma-[HbT]'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,7])
ylim([0,1])
set(gca,'box','off')
ax14.TickLength = [0.03,0.03];
%% [8-S1o] coherence between gamma-band power and HbT during different arousal-states @ 0.01 Hz
ax15 = subplot(5,3,15);
scatter(ones(1,length(data.NeuralHemoCoherence.Awake.gammaBandPower.C001))*1,data.NeuralHemoCoherence.Awake.gammaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.NeuralHemoCoherence.Awake.gammaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.Sleep.gammaBandPower.C001))*2,data.NeuralHemoCoherence.Sleep.gammaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.NeuralHemoCoherence.Sleep.gammaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.NeuralHemoCoherence.All.gammaBandPower.C001))*3,data.NeuralHemoCoherence.All.gammaBandPower.C001,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NeuralHemoCoherence.All.gammaBandPower.meanC001,0,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[8-S1o] Coherence^2 @ 0.01 Hz','Gamma-[HbT]'})
ylabel('Coherence^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,4])
ylim([0,1])
set(gca,'box','off')
ax15.TickLength = [0.03,0.03];
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig8-S1']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig8-S1'])
    %% statistical diary
    diaryFile = [dirpath 'Fig8-S1_Statistics.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on
    % delta-HbT 0.1 Hz coherence statistical diary
    disp('======================================================================================================================')
    disp('[8-S1b] Generalized linear mixed-effects model statistics for delta-HbT coherence^2 @ 0.1 Hz for Rest, NREM, REM, Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Delta_Coh01_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  Delta-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Rest.deltaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Rest.deltaBandPower.stdC01,2))]); disp(' ')
    disp(['NREM  Delta-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.NREM.deltaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.NREM.deltaBandPower.stdC01,2))]); disp(' ')
    disp(['REM   Delta-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.REM.deltaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.REM.deltaBandPower.stdC01,2))]); disp(' ')
    disp(['Awake Delta-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Awake.deltaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Awake.deltaBandPower.stdC01,2))]); disp(' ')
    disp(['Sleep Delta-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Sleep.deltaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Sleep.deltaBandPower.stdC01,2))]); disp(' ')
    disp(['All   Delta-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.All.deltaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.All.deltaBandPower.stdC01,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % delta-HbT 0.01 Hz coherence statistical diary
    disp('======================================================================================================================')
    disp('[8-S1c] Generalized linear mixed-effects model statistics for delta-HbT coherence^2 @ 0.01 Hz for Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Delta_Coh001_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake Delta-[HbT] 0.01 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Awake.deltaBandPower.meanC001,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Awake.deltaBandPower.stdC001,2))]); disp(' ')
    disp(['Sleep Delta-[HbT] 0.01 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Sleep.deltaBandPower.meanC001,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Sleep.deltaBandPower.stdC001,2))]); disp(' ')
    disp(['All   Delta-[HbT] 0.01 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.All.deltaBandPower.meanC001,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.All.deltaBandPower.stdC001,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % Theta-[HbT] 0.1 Hz coherence statistical diary
    disp('======================================================================================================================')
    disp('[8-S1e] Generalized linear mixed-effects model statistics for Theta-[HbT] coherence^2 @ 0.1 Hz for Rest, NREM, REM, Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Theta_Coh01_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  Theta-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Rest.thetaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Rest.thetaBandPower.stdC01,2))]); disp(' ')
    disp(['NREM  Theta-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.NREM.thetaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.NREM.thetaBandPower.stdC01,2))]); disp(' ')
    disp(['REM   Theta-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.REM.thetaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.REM.thetaBandPower.stdC01,2))]); disp(' ')
    disp(['Awake Theta-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Awake.thetaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Awake.thetaBandPower.stdC01,2))]); disp(' ')
    disp(['Sleep Theta-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Sleep.thetaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Sleep.thetaBandPower.stdC01,2))]); disp(' ')
    disp(['All   Theta-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.All.thetaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.All.thetaBandPower.stdC01,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % Theta-[HbT] 0.01 Hz coherence statistical diary
    disp('======================================================================================================================')
    disp('[8-S1f] Generalized linear mixed-effects model statistics for Theta-[HbT] coherence^2 @ 0.01 Hz for Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Theta_Coh001_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake Theta-[HbT] 0.01 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Awake.thetaBandPower.meanC001,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Awake.thetaBandPower.stdC001,2))]); disp(' ')
    disp(['Sleep Theta-[HbT] 0.01 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Sleep.thetaBandPower.meanC001,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Sleep.thetaBandPower.stdC001,2))]); disp(' ')
    disp(['All   Theta-[HbT] 0.01 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.All.thetaBandPower.meanC001,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.All.thetaBandPower.stdC001,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % Alpha-[HbT] 0.1 Hz coherence statistical diary
    disp('======================================================================================================================')
    disp('[8-S1h] Generalized linear mixed-effects model statistics for Alpha-[HbT] coherence^2 @ 0.1 Hz for Rest, NREM, REM, Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Alpha_Coh01_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  Alpha-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Rest.alphaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Rest.alphaBandPower.stdC01,2))]); disp(' ')
    disp(['NREM  Alpha-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.NREM.alphaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.NREM.alphaBandPower.stdC01,2))]); disp(' ')
    disp(['REM   Alpha-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.REM.alphaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.REM.alphaBandPower.stdC01,2))]); disp(' ')
    disp(['Awake Alpha-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Awake.alphaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Awake.alphaBandPower.stdC01,2))]); disp(' ')
    disp(['Sleep Alpha-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Sleep.alphaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Sleep.alphaBandPower.stdC01,2))]); disp(' ')
    disp(['All   Alpha-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.All.alphaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.All.alphaBandPower.stdC01,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % Alpha-[HbT] 0.01 Hz coherence statistical diary
    disp('======================================================================================================================')
    disp('[8-S1i] Generalized linear mixed-effects model statistics for Alpha-[HbT] coherence^2 @ 0.01 Hz for Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Alpha_Coh001_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake Alpha-[HbT] 0.01 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Awake.alphaBandPower.meanC001,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Awake.alphaBandPower.stdC001,2))]); disp(' ')
    disp(['Sleep Alpha-[HbT] 0.01 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Sleep.alphaBandPower.meanC001,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Sleep.alphaBandPower.stdC001,2))]); disp(' ')
    disp(['All   Alpha-[HbT] 0.01 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.All.alphaBandPower.meanC001,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.All.alphaBandPower.stdC001,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % Beta-[HbT] 0.1 Hz coherence statistical diary
    disp('======================================================================================================================')
    disp('[8-S1k] Generalized linear mixed-effects model statistics for Beta-[HbT] coherence^2 @ 0.1 Hz for Rest, NREM, REM, Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Beta_Coh01_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  Beta-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Rest.betaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Rest.betaBandPower.stdC01,2))]); disp(' ')
    disp(['NREM  Beta-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.NREM.betaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.NREM.betaBandPower.stdC01,2))]); disp(' ')
    disp(['REM   Beta-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.REM.betaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.REM.betaBandPower.stdC01,2))]); disp(' ')
    disp(['Awake Beta-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Awake.betaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Awake.betaBandPower.stdC01,2))]); disp(' ')
    disp(['Sleep Beta-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Sleep.betaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Sleep.betaBandPower.stdC01,2))]); disp(' ')
    disp(['All   Beta-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.All.betaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.All.betaBandPower.stdC01,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % Beta-[HbT] 0.01 Hz coherence statistical diary
    disp('======================================================================================================================')
    disp('[8-S1l] Generalized linear mixed-effects model statistics for Beta-[HbT] coherence^2 @ 0.01 Hz for Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Beta_Coh001_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake Beta-[HbT] 0.01 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Awake.betaBandPower.meanC001,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Awake.betaBandPower.stdC001,2))]); disp(' ')
    disp(['Sleep Beta-[HbT] 0.01 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Sleep.betaBandPower.meanC001,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Sleep.betaBandPower.stdC001,2))]); disp(' ')
    disp(['All   Beta-[HbT] 0.01 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.All.betaBandPower.meanC001,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.All.betaBandPower.stdC001,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % Gamma-[HbT] 0.1 Hz coherence statistical diary
    disp('======================================================================================================================')
    disp('[8-S1n] Generalized linear mixed-effects model statistics for Gamma-[HbT] coherence^2 @ 0.1 Hz for Rest, NREM, REM, Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Gamma_Coh01_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Rest  Gamma-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Rest.gammaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Rest.gammaBandPower.stdC01,2))]); disp(' ')
    disp(['NREM  Gamma-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.NREM.gammaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.NREM.gammaBandPower.stdC01,2))]); disp(' ')
    disp(['REM   Gamma-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.REM.gammaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.REM.gammaBandPower.stdC01,2))]); disp(' ')
    disp(['Awake Gamma-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Awake.gammaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Awake.gammaBandPower.stdC01,2))]); disp(' ')
    disp(['Sleep Gamma-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Sleep.gammaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Sleep.gammaBandPower.stdC01,2))]); disp(' ')
    disp(['All   Gamma-[HbT] 0.1 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.All.gammaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.All.gammaBandPower.stdC01,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    % Gamma-[HbT] 0.01 Hz coherence statistical diary
    disp('======================================================================================================================')
    disp('[8-S1o] Generalized linear mixed-effects model statistics for Gamma-[HbT] coherence^2 @ 0.01 Hz for Awake, Sleep, and All')
    disp('======================================================================================================================')
    disp(Gamma_Coh001_Stats)
    disp('----------------------------------------------------------------------------------------------------------------------')
    disp(['Awake Gamma-[HbT] 0.01 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Awake.gammaBandPower.meanC001,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Awake.gammaBandPower.stdC001,2))]); disp(' ')
    disp(['Sleep Gamma-[HbT] 0.01 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.Sleep.gammaBandPower.meanC001,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.Sleep.gammaBandPower.stdC001,2))]); disp(' ')
    disp(['All   Gamma-[HbT] 0.01 Hz Coh2: ' num2str(round(data.NeuralHemoCoherence.All.gammaBandPower.meanC001,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.All.gammaBandPower.stdC001,2))]); disp(' ')
    disp('----------------------------------------------------------------------------------------------------------------------')
    diary off
    %% organized for supplemental table
    % variable names
    ColumnNames = {'Rest','NREM','REM','Awake','Sleep','All'};
    % delta-band 0.1 Hz Coh power
    for aa = 1:length(ColumnNames)
        Delta_C01_MeanStD{1,aa} = [num2str(round(data.NeuralHemoCoherence.(ColumnNames{1,aa}).deltaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.(ColumnNames{1,aa}).deltaBandPower.stdC01,2))]; %#ok<*AGROW>
    end
    % delta-band 0.1 Hz Coh p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Rest') == true
            Delta_C01_pVal{1,aa} = {' '};
        else
            Delta_C01_pVal{1,aa} = ['p < ' num2str(Delta_Coh01_Stats.Coefficients.pValue(aa,1))];
        end
    end
    % delta-band 0.01 Hz Coh power
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Awake') == true || strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Delta_C001_MeanStD{1,aa} = [num2str(round(data.NeuralHemoCoherence.(ColumnNames{1,aa}).deltaBandPower.meanC001,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.(ColumnNames{1,aa}).deltaBandPower.stdC001,2))];
        else
            Delta_C001_MeanStD{1,aa} = {' '};
        end
    end
    % delta-band 0.01 Hz Coh p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Delta_C001_pVal{1,aa} = ['p < ' num2str(Delta_Coh001_Stats.Coefficients.pValue(aa - 3,1))];
        else
            Delta_C001_pVal{1,aa} = {' '};
        end
    end
    % theta-band 0.1 Hz Coh power
    for aa = 1:length(ColumnNames)
        Theta_C01_MeanStD{1,aa} = [num2str(round(data.NeuralHemoCoherence.(ColumnNames{1,aa}).thetaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.(ColumnNames{1,aa}).thetaBandPower.stdC01,2))]; %#ok<*AGROW>
    end
    % theta-band 0.1 Hz Coh p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Rest') == true
            Theta_C01_pVal{1,aa} = {' '};
        else
            Theta_C01_pVal{1,aa} = ['p < ' num2str(Theta_Coh01_Stats.Coefficients.pValue(aa,1))];
        end
    end
    % theta-band 0.01 Hz Coh power
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Awake') == true || strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Theta_C001_MeanStD{1,aa} = [num2str(round(data.NeuralHemoCoherence.(ColumnNames{1,aa}).thetaBandPower.meanC001,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.(ColumnNames{1,aa}).thetaBandPower.stdC001,2))];
        else
            Theta_C001_MeanStD{1,aa} = {' '};
        end
    end
    % theta-band 0.01 Hz Coh p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Theta_C001_pVal{1,aa} = ['p < ' num2str(Theta_Coh001_Stats.Coefficients.pValue(aa - 3,1))];
        else
            Theta_C001_pVal{1,aa} = {' '};
        end
    end
    % alpha-band 0.1 Hz Coh power
    for aa = 1:length(ColumnNames)
        Alpha_C01_MeanStD{1,aa} = [num2str(round(data.NeuralHemoCoherence.(ColumnNames{1,aa}).alphaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.(ColumnNames{1,aa}).alphaBandPower.stdC01,2))]; %#ok<*AGROW>
    end
    % alpha-band 0.1 Hz Coh p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Rest') == true
            Alpha_C01_pVal{1,aa} = {' '};
        else
            Alpha_C01_pVal{1,aa} = ['p < ' num2str(Alpha_Coh01_Stats.Coefficients.pValue(aa,1))];
        end
    end
    % alpha-band 0.01 Hz Coh power
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Awake') == true || strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Alpha_C001_MeanStD{1,aa} = [num2str(round(data.NeuralHemoCoherence.(ColumnNames{1,aa}).alphaBandPower.meanC001,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.(ColumnNames{1,aa}).alphaBandPower.stdC001,2))];
        else
            Alpha_C001_MeanStD{1,aa} = {' '};
        end
    end
    % alpha-band 0.01 Hz Coh p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Alpha_C001_pVal{1,aa} = ['p < ' num2str(Alpha_Coh001_Stats.Coefficients.pValue(aa - 3,1))];
        else
            Alpha_C001_pVal{1,aa} = {' '};
        end
    end
    % beta-band 0.1 Hz Coh power
    for aa = 1:length(ColumnNames)
        Beta_C01_MeanStD{1,aa} = [num2str(round(data.NeuralHemoCoherence.(ColumnNames{1,aa}).betaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.(ColumnNames{1,aa}).betaBandPower.stdC01,2))]; %#ok<*AGROW>
    end
    % beta-band 0.1 Hz Coh p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Rest') == true
            Beta_C01_pVal{1,aa} = {' '};
        else
            Beta_C01_pVal{1,aa} = ['p < ' num2str(Beta_Coh01_Stats.Coefficients.pValue(aa,1))];
        end
    end
    % beta-band 0.01 Hz Coh power
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Awake') == true || strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Beta_C001_MeanStD{1,aa} = [num2str(round(data.NeuralHemoCoherence.(ColumnNames{1,aa}).betaBandPower.meanC001,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.(ColumnNames{1,aa}).betaBandPower.stdC001,2))];
        else
            Beta_C001_MeanStD{1,aa} = {' '};
        end
    end
    % beta-band 0.01 Hz Coh p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Beta_C001_pVal{1,aa} = ['p < ' num2str(Beta_Coh001_Stats.Coefficients.pValue(aa - 3,1))];
        else
            Beta_C001_pVal{1,aa} = {' '};
        end
    end
    % gamma-band 0.1 Hz Coh power
    for aa = 1:length(ColumnNames)
        Gamma_C01_MeanStD{1,aa} = [num2str(round(data.NeuralHemoCoherence.(ColumnNames{1,aa}).gammaBandPower.meanC01,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.(ColumnNames{1,aa}).gammaBandPower.stdC01,2))]; %#ok<*AGROW>
    end
    % gamma-band 0.1 Hz Coh p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Rest') == true
            Gamma_C01_pVal{1,aa} = {' '};
        else
            Gamma_C01_pVal{1,aa} = ['p < ' num2str(Gamma_Coh01_Stats.Coefficients.pValue(aa,1))];
        end
    end
    % gamma-band 0.01 Hz Coh power
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Awake') == true || strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Gamma_C001_MeanStD{1,aa} = [num2str(round(data.NeuralHemoCoherence.(ColumnNames{1,aa}).gammaBandPower.meanC001,2)) ' +/- ' num2str(round(data.NeuralHemoCoherence.(ColumnNames{1,aa}).gammaBandPower.stdC001,2))];
        else
            Gamma_C001_MeanStD{1,aa} = {' '};
        end
    end
    % gamma-band 0.01 Hz Coh p-values
    for aa = 1:length(ColumnNames)
        if strcmp(ColumnNames{1,aa},'Sleep') == true || strcmp(ColumnNames{1,aa},'All') == true
            Gamma_C001_pVal{1,aa} = ['p < ' num2str(Gamma_Coh001_Stats.Coefficients.pValue(aa - 3,1))];
        else
            Gamma_C001_pVal{1,aa} = {' '};
        end
    end
    %% save table data
    if isfield(AnalysisResults,'NeuralHemoCoherence') == false
        AnalysisResults.NeuralHemoCoherence = [];
    end
    if isfield(AnalysisResults.NeuralHemoCoherence,'deltaBandPower') == false
        AnalysisResults.NeuralHemoCoherence.columnNames = ColumnNames;
        AnalysisResults.NeuralHemoCoherence.deltaBandPower.meanStD01 = Delta_C01_MeanStD;
        AnalysisResults.NeuralHemoCoherence.deltaBandPower.p01 = Delta_C01_pVal;
        AnalysisResults.NeuralHemoCoherence.thetaBandPower.meanStD01 = Theta_C01_MeanStD;
        AnalysisResults.NeuralHemoCoherence.thetaBandPower.p01 = Theta_C01_pVal;
        AnalysisResults.NeuralHemoCoherence.alphaBandPower.meanStD01 = Alpha_C01_MeanStD;
        AnalysisResults.NeuralHemoCoherence.alphaBandPower.p01 = Alpha_C01_pVal;
        AnalysisResults.NeuralHemoCoherence.betaBandPower.meanStD01 = Beta_C01_MeanStD;
        AnalysisResults.NeuralHemoCoherence.betaBandPower.p01 = Beta_C01_pVal;
        AnalysisResults.NeuralHemoCoherence.gammaBandPower.meanStD01 = Gamma_C01_MeanStD;
        AnalysisResults.NeuralHemoCoherence.gammaBandPower.p01 = Gamma_C01_pVal;
        AnalysisResults.NeuralHemoCoherence.deltaBandPower.meanStD001 = Delta_C001_MeanStD;
        AnalysisResults.NeuralHemoCoherence.deltaBandPower.p001 = Delta_C001_pVal;
        AnalysisResults.NeuralHemoCoherence.thetaBandPower.meanStD001 = Theta_C001_MeanStD;
        AnalysisResults.NeuralHemoCoherence.thetaBandPower.p001 = Theta_C001_pVal;
        AnalysisResults.NeuralHemoCoherence.alphaBandPower.meanStD001 = Alpha_C001_MeanStD;
        AnalysisResults.NeuralHemoCoherence.alphaBandPower.p001 = Alpha_C001_pVal;
        AnalysisResults.NeuralHemoCoherence.betaBandPower.meanStD001 = Beta_C001_MeanStD;
        AnalysisResults.NeuralHemoCoherence.betaBandPower.p001 = Beta_C001_pVal;
        AnalysisResults.NeuralHemoCoherence.gammaBandPower.meanStD001 = Gamma_C001_MeanStD;
        AnalysisResults.NeuralHemoCoherence.gammaBandPower.p001 = Gamma_C001_pVal;
        cd(rootFolder)
        save('AnalysisResults.mat','AnalysisResults')
    end
end

end
