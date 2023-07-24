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
[~,~,alldata] = xlsread(msExcelFile); %#ok<XLSRD>
groups = {'Naive','SSP_SAP','Blank_SAP'};
% pre-allocate for concatenation
for aa = 1:length(groups)
    group = groups{1,aa};
    data.(group).AnimalID = {};
    data.(group).Sex = [];
    data.(group).Group = {};
    data.(group).LH = [];
    data.(group).RH = [];
    data.(group).hemLH = {};
    data.(group).hemRH = {};
end
% conversion from circular ROI to cubic mm
height = 70/1000; % 70 micron section -> mm
radius = 0.5; % 1 mm diameter circle counting ROI;
sliceVolume = pi*radius^2*height;
cubicRatio = 1/sliceVolume;
% concatenate data for each group/hemishpere
for aa = 2:size(alldata,1)
    group = alldata{aa,3};
    data.(group).AnimalID = cat(1,data.(group).AnimalID,alldata{aa,1});
    data.(group).Sex = cat(1,data.(group).Sex,alldata{aa,2});
    data.(group).Group = cat(1,data.(group).Group,alldata{aa,3});
    data.(group).LH = cat(1,data.(group).LH,alldata{aa,4}*cubicRatio);
    data.(group).RH = cat(1,data.(group).RH,alldata{aa,5}*cubicRatio);
    data.(group).hemLH = cat(1,data.(group).hemLH,'LH');
    data.(group).hemRH = cat(1,data.(group).hemRH,'RH');
end
% mean/std of each hemisphere
for aa = 1:length(groups)
    group = groups{1,aa};
    data.(group).LH_Mean = mean(data.(group).LH,1);
    data.(group).LH_StD = std(data.(group).LH,0,1);
    data.(group).RH_Mean = mean(data.(group).RH,1);
    data.(group).RH_StD = std(data.(group).RH,0,1);
end
% statistics - generalized linear mixed effects model
for aa = 1:length(groups)
    group = groups{1,aa};
    Stats.(group).tableSize = cat(1,data.(group).LH,data.(group).RH);
    Stats.(group).Table = table('Size',[size(Stats.(group).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Hemisphere','Count'});
    Stats.(group).Table.Mouse = cat(1,data.(group).AnimalID,data.(group).AnimalID);
    Stats.(group).Table.Hemisphere = cat(1,data.(group).hemLH,data.(group).hemRH);
    Stats.(group).Table.Count = cat(1,data.(group).LH,data.(group).RH);
    Stats.(group).FitFormula = 'Count ~ 1 + Hemisphere + (1|Mouse)';
    Stats.(group).Stats = fitglme(Stats.(group).Table,Stats.(group).FitFormula);
end
% Naive vs blank RH
Stats.NaiveBlank.tableSize = cat(1,data.Naive.RH,data.Blank_SAP.RH);
Stats.NaiveBlank.Table = table('Size',[size(Stats.Blank_SAP.tableSize,1),2],'VariableTypes',{'string','double'},'VariableNames',{'Group','Count'});
Stats.NaiveBlank.Table.Group = cat(1,data.Naive.Group,data.Blank_SAP.Group);
Stats.NaiveBlank.Table.Count = cat(1,data.Naive.RH,data.Blank_SAP.RH);
Stats.NaiveBlank.FitFormula = 'Count ~ 1 + Group';
Stats.NaiveBlank.Stats = fitglme(Stats.NaiveBlank.Table,Stats.NaiveBlank.FitFormula);
% Blank vs SSP RH
Stats.BlankSSP.tableSize = cat(1,data.Blank_SAP.RH,data.SSP_SAP.RH);
Stats.BlankSSP.Table = table('Size',[size(Stats.Blank_SAP.tableSize,1),2],'VariableTypes',{'string','double'},'VariableNames',{'Group','Count'});
Stats.BlankSSP.Table.Group = cat(1,data.Blank_SAP.Group,data.SSP_SAP.Group);
Stats.BlankSSP.Table.Count = cat(1,data.Blank_SAP.RH,data.SSP_SAP.RH);
Stats.BlankSSP.FitFormula = 'Count ~ 1 + Group';
Stats.BlankSSP.Stats = fitglme(Stats.BlankSSP.Table,Stats.BlankSSP.FitFormula);
% indeces for scatter plot
C57_LH_inds = ones(length(data.Naive.LH),1)*1;
C57_RH_inds = ones(length(data.Naive.RH),1)*2;
Blank_LH_inds = ones(length(data.Blank_SAP.LH),1)*3;
Blank_RH_inds = ones(length(data.Blank_SAP.RH),1)*4;
SSP_LH_inds = ones(length(data.SSP_SAP.LH),1)*5;
SSP_RH_inds = ones(length(data.SSP_SAP.RH),1)*6;
% Figure 1
Fig1 = figure('Name','Figure 1','units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
b1 = bar(1,data.Naive.LH_Mean,'FaceColor',colors('sapphire'));
hold on
bar(2,data.Naive.RH_Mean,'FaceColor',colors('sapphire'))
for aa = 1:length(data.Naive.LH)
    x = [C57_LH_inds(aa,1),C57_RH_inds(aa,1)];
    y = [data.Naive.LH(aa,1),data.Naive.RH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','off', 'jitterAmount',0.25)
end
% Blank-SAP - plot each data point, connect L/R hemispheres
b2 = bar(3,data.Blank_SAP.LH_Mean,'FaceColor',colors('north texas green'));
bar(4,data.Blank_SAP.LH_Mean,'FaceColor',colors('north texas green'))
for aa = 1:length(data.Blank_SAP.LH)
    x = [Blank_LH_inds(aa,1),Blank_RH_inds(aa,1)];
    y = [data.Blank_SAP.LH(aa,1),data.Blank_SAP.RH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off', 'jitterAmount',0.25)
end
% SSP-SAP - plot each data point, connect L/R hemispheres
b3 = bar(5,data.SSP_SAP.LH_Mean,'FaceColor',colors('electric purple'));
bar(6,data.SSP_SAP.RH_Mean,'FaceColor',colors('electric purple'))
for aa = 1:length(data.SSP_SAP.LH)
    x = [SSP_LH_inds(aa,1),SSP_RH_inds(aa,1)];
    y = [data.SSP_SAP.LH(aa,1),data.SSP_SAP.RH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off', 'jitterAmount',0.25)
end
% figure characteristics
ylabel('nNOS Cells/mm^3 ')
legend([b1,b2,b3],'Naive','Blank-SAP','SSP-SAP')
set(gca,'xtick',[1,2,3,4,5,6])
set(gca,'xticklabel',{'LH','RH','LH','RH (Rx)','LH','RH (Rx)'})
xtickangle(45)
axis square
xlim([0,7])
set(gca,'box','off')
subplot(1,2,2)
text(0.5,0.5,'DAPI quantification');
axis off
% save figure(s)
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
    disp('======================================================================================================================')
    disp('GLME statistics for L/R Naive cell counts')
    disp('======================================================================================================================')
    disp(Stats.Naive.Stats)
    disp('======================================================================================================================')
    disp('GLME statistics for L/R Blank-SAP cell counts')
    disp('======================================================================================================================')
    disp(Stats.Blank_SAP.Stats)
    disp('======================================================================================================================')
    disp('GLME statistics for L/R SSP-SAP cell counts')
    disp('======================================================================================================================')
    disp(Stats.SSP_SAP.Stats)
    disp('======================================================================================================================')
    disp('GLME statistics for R/R Naive vs. Blank cell counts')
    disp('======================================================================================================================')
    disp(Stats.NaiveBlank.Stats)
    disp('======================================================================================================================')
    disp('GLME statistics for R/R Blank vs. SSP cell counts')
    disp('======================================================================================================================')
    disp(Stats.BlankSSP.Stats)
    diary off
end