function [] = Fig1_nNOS(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
path = [rootFolder delim 'Results_Turner'];
cd(path)
% setup and pull data from excel sheet
msExcelFile = 'DiaphoraseCellCounts.xlsx';
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
% conversion from circular ROI to cubic mm
% height = 70/1000; % 70 micron section -> mm
radius = 0.5; % 1 mm diameter circle counting ROI;
circArea = pi*radius^2;
% sliceVolume = pi*radius^2*height;
% cubicRatio = 1/sliceVolume;
cubicRatio = 1/circArea;
% concatenate data for each group/hemishpere
for aa = 2:size(allNADPHdata,1)
    group = allNADPHdata{aa,3};
    NADPHdata.(group).AnimalID = cat(1,NADPHdata.(group).AnimalID,allNADPHdata{aa,1});
    NADPHdata.(group).Sex = cat(1,NADPHdata.(group).Sex,allNADPHdata{aa,2});
    NADPHdata.(group).Group = cat(1,NADPHdata.(group).Group,allNADPHdata{aa,3});
    NADPHdata.(group).LH = cat(1,NADPHdata.(group).LH,allNADPHdata{aa,4}*cubicRatio);
    NADPHdata.(group).RH = cat(1,NADPHdata.(group).RH,allNADPHdata{aa,5}*cubicRatio);
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

%% DAPI quantification
path = [rootFolder delim 'Results_Turner'];
cd(path)
% setup and pull data from excel sheet
msExcelFile = 'DAPICellCounts.xlsx';
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
squareRatio = 1/(.650*.650);
% concatenate data for each group/hemishpere
for aa = 2:size(allDAPIdata,1)
    group = allDAPIdata{aa,2};
    DAPIdata.(group).AnimalID = cat(1,DAPIdata.(group).AnimalID,allDAPIdata{aa,1});
    DAPIdata.(group).Group = cat(1,DAPIdata.(group).Group,allDAPIdata{aa,2});
    DAPIdata.(group).Count = cat(1,DAPIdata.(group).Count,allDAPIdata{aa,6}*squareRatio);
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

%% Figure 1
Fig1 = figure('Name','Figure 1','units','normalized','outerposition',[0 0 1 1]);

% NADPH diaphorase quantification
subplot(1,2,1)
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

% DAPI quantification
subplot(1,2,2)
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