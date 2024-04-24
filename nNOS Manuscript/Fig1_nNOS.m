function [] = Fig1_nNOS(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
path = [rootFolder delim 'Results_Turner'];
cd(path)
% conversion from circular ROI to cubic mm
radius = 0.5; % 1 mm diameter circle counting ROI;
circArea = pi*radius^2;
squareRatio = 1/circArea;

%% NADPH diaphorase counts
% setup and pull data from excel sheet
msExcelFile = 'Diaphorase_Counts.xlsx';
[~,~,allNADPHdata] = xlsread(msExcelFile); %#ok<XLSRD>
groups = {'Naive','SSP_SAP','Blank_SAP'};
% pre-allocate for concatenation
for aa = 1:length(groups)
    group = groups{1,aa};
    NADPHdata.(group).AnimalID = {};
    NADPHdata.(group).Sex = [];
    NADPHdata.(group).Group = {};
    NADPHdata.(group).LH = [];
    NADPHdata.(group).RH = [];
    NADPHdata.(group).hemLH = {};
    NADPHdata.(group).hemRH = {};
end
% concatenate data for each group/hemishpere
for aa = 2:size(allNADPHdata,1)
    group = allNADPHdata{aa,3};
    NADPHdata.(group).AnimalID = cat(1,NADPHdata.(group).AnimalID,allNADPHdata{aa,1});
    NADPHdata.(group).Sex = cat(1,NADPHdata.(group).Sex,allNADPHdata{aa,2});
    NADPHdata.(group).Group = cat(1,NADPHdata.(group).Group,allNADPHdata{aa,3});
    NADPHdata.(group).LH = cat(1,NADPHdata.(group).LH,allNADPHdata{aa,4}*squareRatio);
    NADPHdata.(group).RH = cat(1,NADPHdata.(group).RH,allNADPHdata{aa,5}*squareRatio);
    NADPHdata.(group).hemLH = cat(1,NADPHdata.(group).hemLH,'LH');
    NADPHdata.(group).hemRH = cat(1,NADPHdata.(group).hemRH,'RH');
end
% mean/std of each hemisphere
for aa = 1:length(groups)
    group = groups{1,aa};
    NADPHdata.(group).LH_Mean = mean(NADPHdata.(group).LH,1);
    NADPHdata.(group).LH_StD = std(NADPHdata.(group).LH,0,1);
    NADPHdata.(group).RH_Mean = mean(NADPHdata.(group).RH,1);
    NADPHdata.(group).RH_StD = std(NADPHdata.(group).RH,0,1);
end
% statistics - generalized linear mixed effects model
for aa = 1:length(groups)
    group = groups{1,aa};
    NADPHstats.(group).tableSize = cat(1,NADPHdata.(group).LH,NADPHdata.(group).RH);
    NADPHstats.(group).Table = table('Size',[size(NADPHstats.(group).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Hemisphere','Count'});
    NADPHstats.(group).Table.Mouse = cat(1,NADPHdata.(group).AnimalID,NADPHdata.(group).AnimalID);
    NADPHstats.(group).Table.Hemisphere = cat(1,NADPHdata.(group).hemLH,NADPHdata.(group).hemRH);
    NADPHstats.(group).Table.Count = cat(1,NADPHdata.(group).LH,NADPHdata.(group).RH);
    NADPHstats.(group).FitFormula = 'Count ~ 1 + Hemisphere + (1|Mouse)';
    NADPHstats.(group).Stats = fitglme(NADPHstats.(group).Table,NADPHstats.(group).FitFormula);
end
% Naive vs blank RH
NADPHstats.NaiveBlank.tableSize = cat(1,NADPHdata.Naive.RH,NADPHdata.Blank_SAP.RH);
NADPHstats.NaiveBlank.Table = table('Size',[size(NADPHstats.NaiveBlank.tableSize,1),2],'VariableTypes',{'string','double'},'VariableNames',{'Group','Count'});
NADPHstats.NaiveBlank.Table.Group = cat(1,NADPHdata.Naive.Group,NADPHdata.Blank_SAP.Group);
NADPHstats.NaiveBlank.Table.Count = cat(1,NADPHdata.Naive.RH,NADPHdata.Blank_SAP.RH);
NADPHstats.NaiveBlank.FitFormula = 'Count ~ 1 + Group';
NADPHstats.NaiveBlank.Stats = fitglme(NADPHstats.NaiveBlank.Table,NADPHstats.NaiveBlank.FitFormula);
% Blank vs SSP RH
NADPHstats.BlankSSP.tableSize = cat(1,NADPHdata.Blank_SAP.RH,NADPHdata.SSP_SAP.RH);
NADPHstats.BlankSSP.Table = table('Size',[size(NADPHstats.BlankSSP.tableSize,1),2],'VariableTypes',{'string','double'},'VariableNames',{'Group','Count'});
NADPHstats.BlankSSP.Table.Group = cat(1,NADPHdata.Blank_SAP.Group,NADPHdata.SSP_SAP.Group);
NADPHstats.BlankSSP.Table.Count = cat(1,NADPHdata.Blank_SAP.RH,NADPHdata.SSP_SAP.RH);
NADPHstats.BlankSSP.FitFormula = 'Count ~ 1 + Group';
NADPHstats.BlankSSP.Stats = fitglme(NADPHstats.BlankSSP.Table,NADPHstats.BlankSSP.FitFormula);

%% DAPI IHC counts
path = [rootFolder delim 'Results_Turner'];
cd(path)
% setup and pull data from excel sheet
msExcelFile = 'DAPI_Counts.xlsx';
[~,~,allDAPIdata] = xlsread(msExcelFile); %#ok<XLSRD>
groups = {'Naive','SSP_SAP','Blank_SAP'};
% pre-allocate for concatenation
for aa = 1:length(groups)
    group = groups{1,aa};
    DAPIdata.(group).AnimalID = {};
    DAPIdata.(group).Group = {};
    DAPIdata.(group).Count = [];
end
% conversion from circular ROI to cubic mm
% concatenate data for each group/hemishpere
for aa = 2:size(allDAPIdata,1)
    group = allDAPIdata{aa,2};
    DAPIdata.(group).AnimalID = cat(1,DAPIdata.(group).AnimalID,allDAPIdata{aa,1});
    DAPIdata.(group).Group = cat(1,DAPIdata.(group).Group,allDAPIdata{aa,2});
    DAPIdata.(group).Count = cat(1,DAPIdata.(group).Count,allDAPIdata{aa,6}*(1/(.650*.650)));
end
% mean/std of each hemisphere
for aa = 1:length(groups)
    group = groups{1,aa};
    DAPIdata.(group).Mean = mean(DAPIdata.(group).Count,1);
    DAPIdata.(group).StD = std(DAPIdata.(group).Count,0,1);
end
% Naive vs blank RH
DAPIstats.NaiveBlank.tableSize = cat(1,DAPIdata.Naive.Count,DAPIdata.Blank_SAP.Count);
DAPIstats.NaiveBlank.Table = table('Size',[size(DAPIstats.NaiveBlank.tableSize,1),2],'VariableTypes',{'string','double'},'VariableNames',{'Group','Count'});
DAPIstats.NaiveBlank.Table.Group = cat(1,DAPIdata.Naive.Group,DAPIdata.Blank_SAP.Group);
DAPIstats.NaiveBlank.Table.Count = cat(1,DAPIdata.Naive.Count,DAPIdata.Blank_SAP.Count);
DAPIstats.NaiveBlank.FitFormula = 'Count ~ 1 + Group';
DAPIstats.NaiveBlank.Stats = fitglme(DAPIstats.NaiveBlank.Table,DAPIstats.NaiveBlank.FitFormula);
% Blank vs SSP RH
DAPIstats.BlankSSP.tableSize = cat(1,DAPIdata.Blank_SAP.Count,DAPIdata.SSP_SAP.Count);
DAPIstats.BlankSSP.Table = table('Size',[size(DAPIstats.BlankSSP.tableSize,1),2],'VariableTypes',{'string','double'},'VariableNames',{'Group','Count'});
DAPIstats.BlankSSP.Table.Group = cat(1,DAPIdata.Blank_SAP.Group,DAPIdata.SSP_SAP.Group);
DAPIstats.BlankSSP.Table.Count = cat(1,DAPIdata.Blank_SAP.Count,DAPIdata.SSP_SAP.Count);
DAPIstats.BlankSSP.FitFormula = 'Count ~ 1 + Group';
DAPIstats.BlankSSP.Stats = fitglme(DAPIstats.BlankSSP.Table,DAPIstats.BlankSSP.FitFormula);

%% IBA1 IHC counts
% setup and pull data from excel sheet
msExcelFile = 'IBA1_Counts.xlsx';
[~,~,allIBA1data] = xlsread(msExcelFile); %#ok<XLSRD>
groups = {'SSP_SAP','Blank_SAP'};
% pre-allocate for concatenation
for aa = 1:length(groups)
    group = groups{1,aa};
    IBA1data.(group).AnimalID = {};
    IBA1data.(group).Sex = [];
    IBA1data.(group).Group = {};
    IBA1data.(group).LH = [];
    IBA1data.(group).RH = [];
    IBA1data.(group).hemLH = {};
    IBA1data.(group).hemRH = {};
end
% concatenate data for each group/hemishpere
for aa = 2:size(allIBA1data,1)
    group = allIBA1data{aa,3};
    IBA1data.(group).AnimalID = cat(1,IBA1data.(group).AnimalID,allIBA1data{aa,1});
    IBA1data.(group).Sex = cat(1,IBA1data.(group).Sex,allIBA1data{aa,2});
    IBA1data.(group).Group = cat(1,IBA1data.(group).Group,allIBA1data{aa,3});
    IBA1data.(group).LH = cat(1,IBA1data.(group).LH,allIBA1data{aa,4}*squareRatio);
    IBA1data.(group).RH = cat(1,IBA1data.(group).RH,allIBA1data{aa,5}*squareRatio);
    IBA1data.(group).hemLH = cat(1,IBA1data.(group).hemLH,'LH');
    IBA1data.(group).hemRH = cat(1,IBA1data.(group).hemRH,'RH');
end
% mean/std of each hemisphere
for aa = 1:length(groups)
    group = groups{1,aa};
    IBA1data.(group).LH_Mean = mean(IBA1data.(group).LH,1);
    IBA1data.(group).LH_StD = std(IBA1data.(group).LH,0,1);
    IBA1data.(group).RH_Mean = mean(IBA1data.(group).RH,1);
    IBA1data.(group).RH_StD = std(IBA1data.(group).RH,0,1);
end
% statistics - generalized linear mixed effects model
for aa = 1:length(groups)
    group = groups{1,aa};
    IBA1stats.(group).tableSize = cat(1,IBA1data.(group).LH,IBA1data.(group).RH);
    IBA1stats.(group).Table = table('Size',[size(IBA1stats.(group).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Hemisphere','Count'});
    IBA1stats.(group).Table.Mouse = cat(1,IBA1data.(group).AnimalID,IBA1data.(group).AnimalID);
    IBA1stats.(group).Table.Hemisphere = cat(1,IBA1data.(group).hemLH,IBA1data.(group).hemRH);
    IBA1stats.(group).Table.Count = cat(1,IBA1data.(group).LH,IBA1data.(group).RH);
    IBA1stats.(group).FitFormula = 'Count ~ 1 + Hemisphere + (1|Mouse)';
    IBA1stats.(group).Stats = fitglme(IBA1stats.(group).Table,IBA1stats.(group).FitFormula);
end
% Blank vs SSP RH
IBA1stats.BlankSSP.tableSize = cat(1,IBA1data.Blank_SAP.RH,IBA1data.SSP_SAP.RH);
IBA1stats.BlankSSP.Table = table('Size',[size(IBA1stats.BlankSSP.tableSize,1),2],'VariableTypes',{'string','double'},'VariableNames',{'Group','Count'});
IBA1stats.BlankSSP.Table.Group = cat(1,IBA1data.Blank_SAP.Group,IBA1data.SSP_SAP.Group);
IBA1stats.BlankSSP.Table.Count = cat(1,IBA1data.Blank_SAP.RH,IBA1data.SSP_SAP.RH);
IBA1stats.BlankSSP.FitFormula = 'Count ~ 1 + Group';
IBA1stats.BlankSSP.Stats = fitglme(IBA1stats.BlankSSP.Table,IBA1stats.BlankSSP.FitFormula);

%% nNOS IHC counts
% setup and pull data from excel sheet
msExcelFile = 'nNOS_Counts.xlsx';
[~,~,allnNOSdata] = xlsread(msExcelFile); %#ok<XLSRD>
groups = {'SSP_SAP','Blank_SAP'};
% pre-allocate for concatenation
for aa = 1:length(groups)
    group = groups{1,aa};
    nNOSdata.(group).AnimalID = {};
    nNOSdata.(group).Sex = [];
    nNOSdata.(group).Group = {};
    nNOSdata.(group).LH = [];
    nNOSdata.(group).RH = [];
    nNOSdata.(group).hemLH = {};
    nNOSdata.(group).hemRH = {};
end
% concatenate data for each group/hemishpere
for aa = 2:size(allnNOSdata,1)
    group = allnNOSdata{aa,3};
    nNOSdata.(group).AnimalID = cat(1,nNOSdata.(group).AnimalID,allnNOSdata{aa,1});
    nNOSdata.(group).Sex = cat(1,nNOSdata.(group).Sex,allnNOSdata{aa,2});
    nNOSdata.(group).Group = cat(1,nNOSdata.(group).Group,allnNOSdata{aa,3});
    nNOSdata.(group).LH = cat(1,nNOSdata.(group).LH,allnNOSdata{aa,4}*squareRatio);
    nNOSdata.(group).RH = cat(1,nNOSdata.(group).RH,allnNOSdata{aa,5}*squareRatio);
    nNOSdata.(group).hemLH = cat(1,nNOSdata.(group).hemLH,'LH');
    nNOSdata.(group).hemRH = cat(1,nNOSdata.(group).hemRH,'RH');
end
% mean/std of each hemisphere
for aa = 1:length(groups)
    group = groups{1,aa};
    nNOSdata.(group).LH_Mean = mean(nNOSdata.(group).LH,1);
    nNOSdata.(group).LH_StD = std(nNOSdata.(group).LH,0,1);
    nNOSdata.(group).RH_Mean = mean(nNOSdata.(group).RH,1);
    nNOSdata.(group).RH_StD = std(nNOSdata.(group).RH,0,1);
end
% statistics - generalized linear mixed effects model
for aa = 1:length(groups)
    group = groups{1,aa};
    nNOSstats.(group).tableSize = cat(1,nNOSdata.(group).LH,nNOSdata.(group).RH);
    nNOSstats.(group).Table = table('Size',[size(nNOSstats.(group).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Hemisphere','Count'});
    nNOSstats.(group).Table.Mouse = cat(1,nNOSdata.(group).AnimalID,nNOSdata.(group).AnimalID);
    nNOSstats.(group).Table.Hemisphere = cat(1,nNOSdata.(group).hemLH,nNOSdata.(group).hemRH);
    nNOSstats.(group).Table.Count = cat(1,nNOSdata.(group).LH,nNOSdata.(group).RH);
    nNOSstats.(group).FitFormula = 'Count ~ 1 + Hemisphere + (1|Mouse)';
    nNOSstats.(group).Stats = fitglme(nNOSstats.(group).Table,nNOSstats.(group).FitFormula);
end
% Blank vs SSP RH
nNOSstats.BlankSSP.tableSize = cat(1,nNOSdata.Blank_SAP.RH,nNOSdata.SSP_SAP.RH);
nNOSstats.BlankSSP.Table = table('Size',[size(nNOSstats.BlankSSP.tableSize,1),2],'VariableTypes',{'string','double'},'VariableNames',{'Group','Count'});
nNOSstats.BlankSSP.Table.Group = cat(1,nNOSdata.Blank_SAP.Group,nNOSdata.SSP_SAP.Group);
nNOSstats.BlankSSP.Table.Count = cat(1,nNOSdata.Blank_SAP.RH,nNOSdata.SSP_SAP.RH);
nNOSstats.BlankSSP.FitFormula = 'Count ~ 1 + Group';
nNOSstats.BlankSSP.Stats = fitglme(nNOSstats.BlankSSP.Table,nNOSstats.BlankSSP.FitFormula);

%% DAPI IHC fluorescence
% setup and pull data from excel sheet
msExcelFile = 'DAPI_Fluorescence.xlsx';
[~,~,allDAPIFluordata] = xlsread(msExcelFile); %#ok<XLSRD>
groups = {'SSP_SAP','Blank_SAP'};
% pre-allocate for concatenation
for aa = 1:length(groups)
    group = groups{1,aa};
    DAPIFluordata.(group).AnimalID = {};
    DAPIFluordata.(group).Sex = [];
    DAPIFluordata.(group).Group = {};
    DAPIFluordata.(group).LH = [];
    DAPIFluordata.(group).RH = [];
    DAPIFluordata.(group).hemLH = {};
    DAPIFluordata.(group).hemRH = {};
end
% concatenate data for each group/hemishpere
for aa = 2:size(allDAPIFluordata,1)
    group = allDAPIFluordata{aa,3};
    DAPIFluordata.(group).AnimalID = cat(1,DAPIFluordata.(group).AnimalID,allDAPIFluordata{aa,1});
    DAPIFluordata.(group).Sex = cat(1,DAPIFluordata.(group).Sex,allDAPIFluordata{aa,2});
    DAPIFluordata.(group).Group = cat(1,DAPIFluordata.(group).Group,allDAPIFluordata{aa,3});
    DAPIFluordata.(group).LH = cat(1,DAPIFluordata.(group).LH,allDAPIFluordata{aa,4}*squareRatio);
    DAPIFluordata.(group).RH = cat(1,DAPIFluordata.(group).RH,allDAPIFluordata{aa,5}*squareRatio);
    DAPIFluordata.(group).hemLH = cat(1,DAPIFluordata.(group).hemLH,'LH');
    DAPIFluordata.(group).hemRH = cat(1,DAPIFluordata.(group).hemRH,'RH');
end
% mean/std of each hemisphere
for aa = 1:length(groups)
    group = groups{1,aa};
    DAPIFluordata.(group).LH_Mean = mean(DAPIFluordata.(group).LH,1);
    DAPIFluordata.(group).LH_StD = std(DAPIFluordata.(group).LH,0,1);
    DAPIFluordata.(group).RH_Mean = mean(DAPIFluordata.(group).RH,1);
    DAPIFluordata.(group).RH_StD = std(DAPIFluordata.(group).RH,0,1);
end
% statistics - generalized linear mixed effects model
for aa = 1:length(groups)
    group = groups{1,aa};
    DAPIstats.(group).tableSize = cat(1,DAPIFluordata.(group).LH,DAPIFluordata.(group).RH);
    DAPIstats.(group).Table = table('Size',[size(DAPIstats.(group).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Hemisphere','Count'});
    DAPIstats.(group).Table.Mouse = cat(1,DAPIFluordata.(group).AnimalID,DAPIFluordata.(group).AnimalID);
    DAPIstats.(group).Table.Hemisphere = cat(1,DAPIFluordata.(group).hemLH,DAPIFluordata.(group).hemRH);
    DAPIstats.(group).Table.Count = cat(1,DAPIFluordata.(group).LH,DAPIFluordata.(group).RH);
    DAPIstats.(group).FitFormula = 'Count ~ 1 + Hemisphere + (1|Mouse)';
    DAPIstats.(group).Stats = fitglme(DAPIstats.(group).Table,DAPIstats.(group).FitFormula);
end
% Blank vs SSP RH
DAPIstats.BlankSSP.tableSize = cat(1,DAPIFluordata.Blank_SAP.RH,DAPIFluordata.SSP_SAP.RH);
DAPIstats.BlankSSP.Table = table('Size',[size(DAPIstats.BlankSSP.tableSize,1),2],'VariableTypes',{'string','double'},'VariableNames',{'Group','Count'});
DAPIstats.BlankSSP.Table.Group = cat(1,DAPIFluordata.Blank_SAP.Group,DAPIFluordata.SSP_SAP.Group);
DAPIstats.BlankSSP.Table.Count = cat(1,DAPIFluordata.Blank_SAP.RH,DAPIFluordata.SSP_SAP.RH);
DAPIstats.BlankSSP.FitFormula = 'Count ~ 1 + Group';
DAPIstats.BlankSSP.Stats = fitglme(DAPIstats.BlankSSP.Table,DAPIstats.BlankSSP.FitFormula);

%% GFAP IHC fluorescence
% setup and pull data from excel sheet
msExcelFile = 'GFAP_Fluorescence.xlsx';
[~,~,allGFAPdata] = xlsread(msExcelFile); %#ok<XLSRD>
groups = {'SSP_SAP','Blank_SAP'};
% pre-allocate for concatenation
for aa = 1:length(groups)
    group = groups{1,aa};
    GFAPdata.(group).AnimalID = {};
    GFAPdata.(group).Sex = [];
    GFAPdata.(group).Group = {};
    GFAPdata.(group).LH = [];
    GFAPdata.(group).RH = [];
    GFAPdata.(group).hemLH = {};
    GFAPdata.(group).hemRH = {};
end
% concatenate data for each group/hemishpere
for aa = 2:size(allGFAPdata,1)
    group = allGFAPdata{aa,3};
    GFAPdata.(group).AnimalID = cat(1,GFAPdata.(group).AnimalID,allGFAPdata{aa,1});
    GFAPdata.(group).Sex = cat(1,GFAPdata.(group).Sex,allGFAPdata{aa,2});
    GFAPdata.(group).Group = cat(1,GFAPdata.(group).Group,allGFAPdata{aa,3});
    GFAPdata.(group).LH = cat(1,GFAPdata.(group).LH,allGFAPdata{aa,4}*squareRatio);
    GFAPdata.(group).RH = cat(1,GFAPdata.(group).RH,allGFAPdata{aa,5}*squareRatio);
    GFAPdata.(group).hemLH = cat(1,GFAPdata.(group).hemLH,'LH');
    GFAPdata.(group).hemRH = cat(1,GFAPdata.(group).hemRH,'RH');
end
% mean/std of each hemisphere
for aa = 1:length(groups)
    group = groups{1,aa};
    GFAPdata.(group).LH_Mean = mean(GFAPdata.(group).LH,1);
    GFAPdata.(group).LH_StD = std(GFAPdata.(group).LH,0,1);
    GFAPdata.(group).RH_Mean = mean(GFAPdata.(group).RH,1);
    GFAPdata.(group).RH_StD = std(GFAPdata.(group).RH,0,1);
end
% statistics - generalized linear mixed effects model
for aa = 1:length(groups)
    group = groups{1,aa};
    GFAPstats.(group).tableSize = cat(1,GFAPdata.(group).LH,GFAPdata.(group).RH);
    GFAPstats.(group).Table = table('Size',[size(GFAPstats.(group).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Hemisphere','Count'});
    GFAPstats.(group).Table.Mouse = cat(1,GFAPdata.(group).AnimalID,GFAPdata.(group).AnimalID);
    GFAPstats.(group).Table.Hemisphere = cat(1,GFAPdata.(group).hemLH,GFAPdata.(group).hemRH);
    GFAPstats.(group).Table.Count = cat(1,GFAPdata.(group).LH,GFAPdata.(group).RH);
    GFAPstats.(group).FitFormula = 'Count ~ 1 + Hemisphere + (1|Mouse)';
    GFAPstats.(group).Stats = fitglme(GFAPstats.(group).Table,GFAPstats.(group).FitFormula);
end
% Blank vs SSP RH
GFAPstats.BlankSSP.tableSize = cat(1,GFAPdata.Blank_SAP.RH,GFAPdata.SSP_SAP.RH);
GFAPstats.BlankSSP.Table = table('Size',[size(GFAPstats.BlankSSP.tableSize,1),2],'VariableTypes',{'string','double'},'VariableNames',{'Group','Count'});
GFAPstats.BlankSSP.Table.Group = cat(1,GFAPdata.Blank_SAP.Group,GFAPdata.SSP_SAP.Group);
GFAPstats.BlankSSP.Table.Count = cat(1,GFAPdata.Blank_SAP.RH,GFAPdata.SSP_SAP.RH);
GFAPstats.BlankSSP.FitFormula = 'Count ~ 1 + Group';
GFAPstats.BlankSSP.Stats = fitglme(GFAPstats.BlankSSP.Table,GFAPstats.BlankSSP.FitFormula);

%% NeuN IHC fluorescence
% setup and pull data from excel sheet
msExcelFile = 'NeuN_Fluorescence.xlsx';
[~,~,allNeuNdata] = xlsread(msExcelFile); %#ok<XLSRD>
groups = {'SSP_SAP','Blank_SAP'};
% pre-allocate for concatenation
for aa = 1:length(groups)
    group = groups{1,aa};
    NeuNdata.(group).AnimalID = {};
    NeuNdata.(group).Sex = [];
    NeuNdata.(group).Group = {};
    NeuNdata.(group).LH = [];
    NeuNdata.(group).RH = [];
    NeuNdata.(group).hemLH = {};
    NeuNdata.(group).hemRH = {};
end
% concatenate data for each group/hemishpere
for aa = 2:size(allNeuNdata,1)
    group = allNeuNdata{aa,3};
    NeuNdata.(group).AnimalID = cat(1,NeuNdata.(group).AnimalID,allNeuNdata{aa,1});
    NeuNdata.(group).Sex = cat(1,NeuNdata.(group).Sex,allNeuNdata{aa,2});
    NeuNdata.(group).Group = cat(1,NeuNdata.(group).Group,allNeuNdata{aa,3});
    NeuNdata.(group).LH = cat(1,NeuNdata.(group).LH,allNeuNdata{aa,4}*squareRatio);
    NeuNdata.(group).RH = cat(1,NeuNdata.(group).RH,allNeuNdata{aa,5}*squareRatio);
    NeuNdata.(group).hemLH = cat(1,NeuNdata.(group).hemLH,'LH');
    NeuNdata.(group).hemRH = cat(1,NeuNdata.(group).hemRH,'RH');
end
% mean/std of each hemisphere
for aa = 1:length(groups)
    group = groups{1,aa};
    NeuNdata.(group).LH_Mean = mean(NeuNdata.(group).LH,1);
    NeuNdata.(group).LH_StD = std(NeuNdata.(group).LH,0,1);
    NeuNdata.(group).RH_Mean = mean(NeuNdata.(group).RH,1);
    NeuNdata.(group).RH_StD = std(NeuNdata.(group).RH,0,1);
end
% statistics - generalized linear mixed effects model
for aa = 1:length(groups)
    group = groups{1,aa};
    NeuNstats.(group).tableSize = cat(1,NeuNdata.(group).LH,NeuNdata.(group).RH);
    NeuNstats.(group).Table = table('Size',[size(NeuNstats.(group).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Hemisphere','Count'});
    NeuNstats.(group).Table.Mouse = cat(1,NeuNdata.(group).AnimalID,NeuNdata.(group).AnimalID);
    NeuNstats.(group).Table.Hemisphere = cat(1,NeuNdata.(group).hemLH,NeuNdata.(group).hemRH);
    NeuNstats.(group).Table.Count = cat(1,NeuNdata.(group).LH,NeuNdata.(group).RH);
    NeuNstats.(group).FitFormula = 'Count ~ 1 + Hemisphere + (1|Mouse)';
    NeuNstats.(group).Stats = fitglme(NeuNstats.(group).Table,NeuNstats.(group).FitFormula);
end
% Blank vs SSP RH
NeuNstats.BlankSSP.tableSize = cat(1,NeuNdata.Blank_SAP.RH,NeuNdata.SSP_SAP.RH);
NeuNstats.BlankSSP.Table = table('Size',[size(NeuNstats.BlankSSP.tableSize,1),2],'VariableTypes',{'string','double'},'VariableNames',{'Group','Count'});
NeuNstats.BlankSSP.Table.Group = cat(1,NeuNdata.Blank_SAP.Group,NeuNdata.SSP_SAP.Group);
NeuNstats.BlankSSP.Table.Count = cat(1,NeuNdata.Blank_SAP.RH,NeuNdata.SSP_SAP.RH);
NeuNstats.BlankSSP.FitFormula = 'Count ~ 1 + Group';
NeuNstats.BlankSSP.Stats = fitglme(NeuNstats.BlankSSP.Table,NeuNstats.BlankSSP.FitFormula);

%% Figure 1
Fig1 = figure('Name','Figure 1','units','normalized','outerposition',[0 0 1 1]);

% NADPH diaphorase count
subplot(2,4,1)
xInds = ones(1,length(NADPHdata.Naive.RH));
s1 = scatter(xInds*1,NADPHdata.Naive.RH,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','off','jitterAmount',0.25);
hold on
e1 = errorbar(1,NADPHdata.Naive.RH_Mean,NADPHdata.Naive.RH_StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(NADPHdata.Blank_SAP.RH));
s2 = scatter(xInds*2,NADPHdata.Blank_SAP.RH,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(2,NADPHdata.Blank_SAP.RH_Mean,NADPHdata.Blank_SAP.RH_StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(NADPHdata.SSP_SAP.RH));
s3 = scatter(xInds*3,NADPHdata.SSP_SAP.RH,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(3,NADPHdata.SSP_SAP.RH_Mean,NADPHdata.SSP_SAP.RH_StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('nNOS Cells/mm^2 ')
legend([s1,s2,s3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'box','off')
set(gca,'xtick',[])
axis square
xlim([0,4]);
ylim([0,30])

% DAPI IHC count
subplot(2,4,2)
xInds = ones(1,length(DAPIdata.Naive.Count));
scatter(xInds*1,DAPIdata.Naive.Count,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','off','jitterAmount',0.25);
hold on
e1 = errorbar(1,DAPIdata.Naive.Mean,DAPIdata.Naive.StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(DAPIdata.Blank_SAP.Count));
scatter(xInds*2,DAPIdata.Blank_SAP.Count,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(2,DAPIdata.Blank_SAP.Mean,DAPIdata.Blank_SAP.StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
xInds = ones(1,length(DAPIdata.SSP_SAP.Count));
scatter(xInds*3,DAPIdata.SSP_SAP.Count,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e3 = errorbar(3,DAPIdata.SSP_SAP.Mean,DAPIdata.SSP_SAP.StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
ylabel('DAPI labeled Cells/mm^2')
set(gca,'box','off')
set(gca,'xtick',[])
axis square
xlim([0,4]);
ylim([0,3000])

% nNOS IHC count
subplot(2,4,3)
xInds = ones(1,length(nNOSdata.Blank_SAP.RH));
s1 = scatter(xInds*1,nNOSdata.Blank_SAP.RH,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e1 = errorbar(1,nNOSdata.Blank_SAP.RH_Mean,nNOSdata.Blank_SAP.RH_StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(nNOSdata.SSP_SAP.RH));
s2 = scatter(xInds*2,nNOSdata.SSP_SAP.RH,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(2,nNOSdata.SSP_SAP.RH_Mean,nNOSdata.SSP_SAP.RH_StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
ylabel('nNOS Cells/mm^2 ')
legend([s1,s2],'Blank-SAP','SSP-SAP')
set(gca,'box','off')
set(gca,'xtick',[])
axis square
xlim([0,3]);
ylim([0,15])

% IBA1 count
subplot(2,4,4)
xInds = ones(1,length(IBA1data.Blank_SAP.RH));
s1 = scatter(xInds*1,IBA1data.Blank_SAP.RH,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e1 = errorbar(1,IBA1data.Blank_SAP.RH_Mean,IBA1data.Blank_SAP.RH_StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(IBA1data.SSP_SAP.RH));
s2 = scatter(xInds*2,IBA1data.SSP_SAP.RH,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(2,IBA1data.SSP_SAP.RH_Mean,IBA1data.SSP_SAP.RH_StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
ylabel('IBA1 Cells/mm^2 ')
legend([s1,s2],'Blank-SAP','SSP-SAP')
set(gca,'box','off')
set(gca,'xtick',[])
axis square
xlim([0,3]);
ylim([0,500])

% GFAP fluorescence
subplot(2,4,5)
xInds = ones(1,length(GFAPdata.Blank_SAP.RH));
s1 = scatter(xInds*1,GFAPdata.Blank_SAP.RH,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e1 = errorbar(1,GFAPdata.Blank_SAP.RH_Mean,GFAPdata.Blank_SAP.RH_StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(GFAPdata.SSP_SAP.RH));
s2 = scatter(xInds*2,GFAPdata.SSP_SAP.RH,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(2,GFAPdata.SSP_SAP.RH_Mean,GFAPdata.SSP_SAP.RH_StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
ylabel('GFAP Fluor/mm^2 (a.u.)')
legend([s1,s2],'Blank-SAP','SSP-SAP')
set(gca,'box','off')
set(gca,'xtick',[])
axis square
xlim([0,3]);
ylim([0,2000])

% DAPI fluorescence
subplot(2,4,6)
xInds = ones(1,length(DAPIFluordata.Blank_SAP.RH));
s1 = scatter(xInds*1,DAPIFluordata.Blank_SAP.RH,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e1 = errorbar(1,DAPIFluordata.Blank_SAP.RH_Mean,DAPIFluordata.Blank_SAP.RH_StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(DAPIFluordata.SSP_SAP.RH));
s2 = scatter(xInds*2,DAPIFluordata.SSP_SAP.RH,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(2,DAPIFluordata.SSP_SAP.RH_Mean,DAPIFluordata.SSP_SAP.RH_StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
ylabel('DAPI Fluor/mm^2 (a.u.)')
legend([s1,s2],'Blank-SAP','SSP-SAP')
set(gca,'box','off')
set(gca,'xtick',[])
axis square
xlim([0,3]);
ylim([0,10000])

% NeuN fluorescence
subplot(2,4,7)
xInds = ones(1,length(NeuNdata.Blank_SAP.RH));
s1 = scatter(xInds*1,NeuNdata.Blank_SAP.RH,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off','jitterAmount',0.25);
hold on
e1 = errorbar(1,NeuNdata.Blank_SAP.RH_Mean,NeuNdata.Blank_SAP.RH_StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
xInds = ones(1,length(NeuNdata.SSP_SAP.RH));
s2 = scatter(xInds*2,NeuNdata.SSP_SAP.RH,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off','jitterAmount',0.25);
e2 = errorbar(2,NeuNdata.SSP_SAP.RH_Mean,NeuNdata.SSP_SAP.RH_StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
ylabel('NeuN Fluor/mm^2 (a.u.)')
legend([s1,s2],'Blank-SAP','SSP-SAP')
set(gca,'box','off')
set(gca,'xtick',[])
axis square
xlim([0,3]);
ylim([0,5000])

%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'MATLAB Figure Panels' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(Fig1,[dirpath 'Fig1']);
    set(Fig1,'PaperPositionMode','auto');
    print('-vector','-dpdf','-fillpage',[dirpath 'Fig1'])
    diaryFile = [dirpath 'Fig1_Readout.txt'];
    if exist(diaryFile,'file') == 2
        delete(diaryFile)
    end
    diary(diaryFile)
    diary on

    % NADPH diaphorase quantification
    disp('NADPH diaphorase counts (RH) 9 mice per group, mean +/- std'); disp(' ')
    disp(['Naive: ' num2str(NADPHdata.Naive.RH_Mean) ' +/- ' num2str(NADPHdata.Naive.RH_StD)]); disp(' ')
    disp(['Blank: ' num2str(NADPHdata.Blank_SAP.RH_Mean) ' +/- ' num2str(NADPHdata.Blank_SAP.RH_StD)]); disp(' ')
    disp(['SSP: ' num2str(NADPHdata.SSP_SAP.RH_Mean) ' +/- ' num2str(NADPHdata.SSP_SAP.RH_StD)]); disp(' ')
    disp('======================================================================================================================')
    disp('GLME statistics for R/R Naive vs. Blank cell counts')
    disp('======================================================================================================================')
    disp(NADPHstats.NaiveBlank.Stats)
    disp('======================================================================================================================')
    disp('GLME statistics for R/R Blank vs. SSP cell counts')
    disp('======================================================================================================================')
    disp(NADPHstats.BlankSSP.Stats)
    
    % DAPI quantification
    disp('DAPI counts (RH) 5-8 mice per group, mean +/- std'); disp(' ')
    disp(['Naive: ' num2str(DAPIdata.Naive.Mean) ' +/- ' num2str(DAPIdata.Naive.StD)]); disp(' ')
    disp(['Blank: ' num2str(DAPIdata.Blank_SAP.Mean) ' +/- ' num2str(DAPIdata.Blank_SAP.StD)]); disp(' ')
    disp(['SSP: ' num2str(DAPIdata.SSP_SAP.Mean) ' +/- ' num2str(DAPIdata.SSP_SAP.StD)]); disp(' ')
    disp('======================================================================================================================')
    disp('GLME statistics for R/R Naive vs. Blank cell counts')
    disp('======================================================================================================================')
    disp(DAPIstats.NaiveBlank.Stats)
    disp('======================================================================================================================')
    disp('GLME statistics for R/R Blank vs. SSP cell counts')
    disp('======================================================================================================================')
    disp(DAPIstats.BlankSSP.Stats)
    
    diary off
end