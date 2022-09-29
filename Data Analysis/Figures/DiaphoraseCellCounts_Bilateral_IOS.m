function [] = DiaphoraseCellCounts_Bilateral_IOS(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

%% setup and pull data from excel sheet
msExcelFile = 'DiaphoraseCellCounts_Bilateral_IOS.xlsx';
[~,~,alldata] = xlsread(msExcelFile);
treatments = {'Naive','SSP_SAP','Blank_SAP'};
% pre-allocate for concatenation
for aa = 1:length(treatments)
    treatment = treatments{1,aa};
    data.(treatment).LH = [];
    data.(treatment).RH = [];
    data.(treatment).animalID = {};
    data.(treatment).treatment = {};
    data.(treatment).hemLH = {};
    data.(treatment).hemRH = {};
end
% conversion from circular ROI to cubic mm
height = 70/1000; % 70 micron section -> mm
radius = 0.5; % 1 mm diameter circle counting ROI;
sliceVolume = pi*radius^2*height;
cubicRatio = 1/sliceVolume;
% concatenate data for each treatment/hemishpere
for bb = 2:size(alldata,1)
    treatment = alldata{bb,2};
    data.(treatment).LH = cat(1,data.(treatment).LH,alldata{bb,3}*cubicRatio);
    data.(treatment).RH = cat(1,data.(treatment).RH,alldata{bb,4}*cubicRatio);
    data.(treatment).animalID = cat(1,data.(treatment).animalID,alldata{bb,1});
    data.(treatment).treatment = cat(1,data.(treatment).treatment,treatment);
    data.(treatment).hemLH = cat(1,data.(treatment).hemLH,'LH');
    data.(treatment).hemRH = cat(1,data.(treatment).hemRH,'RH');
end
% mean/std of each hemisphere
for cc = 1:length(treatments)
    treatment = treatments{1,cc};
    data.(treatment).LH_Mean = mean(data.(treatment).LH,1);
    data.(treatment).LH_StD = std(data.(treatment).LH,0,1);
    data.(treatment).RH_Mean = mean(data.(treatment).RH,1);
    data.(treatment).RH_StD = std(data.(treatment).RH,0,1);
end
%% statistics - generalized linear mixed effects model
for qq = 1:length(treatments)
    treatment = treatments{1,qq};
    Stats.(treatment).tableSize = cat(1,data.(treatment).LH,data.(treatment).RH);
    Stats.(treatment).Table = table('Size',[size(Stats.(treatment).tableSize,1),3],'VariableTypes',{'string','string','double'},'VariableNames',{'Mouse','Hemisphere','Count'});
    Stats.(treatment).Table.Mouse = cat(1,data.(treatment).animalID,data.(treatment).animalID);
    Stats.(treatment).Table.Hemisphere = cat(1,data.(treatment).hemLH,data.(treatment).hemRH);
    Stats.(treatment).Table.Count = cat(1,data.(treatment).LH,data.(treatment).RH);
    Stats.(treatment).FitFormula = 'Count ~ 1 + Hemisphere + (1|Mouse)';
    Stats.(treatment).Stats = fitglme(Stats.(treatment).Table,Stats.(treatment).FitFormula);
end
Stats.Comp.tableSize = cat(1,data.Blank_SAP.LH,data.SSP_SAP.LH);
Stats.Comp.Table = table('Size',[size(Stats.Blank_SAP.tableSize,1),2],'VariableTypes',{'string','double'},'VariableNames',{'Treatment','Count'});
Stats.Comp.Table.Treatment = cat(1,data.Blank_SAP.treatment,data.SSP_SAP.treatment);
Stats.Comp.Table.Count = cat(1,data.Blank_SAP.LH,data.SSP_SAP.LH);
Stats.Comp.FitFormula = 'Count ~ 1 + Treatment';
Stats.Comp.Stats = fitglme(Stats.Comp.Table,Stats.Comp.FitFormula);
% indeces for scatter plot
C57_LH_inds = ones(length(data.Naive.LH),1)*1;
C57_RH_inds = ones(length(data.Naive.RH),1)*2;
SSP_LH_inds = ones(length(data.SSP_SAP.LH),1)*3;
SSP_RH_inds = ones(length(data.SSP_SAP.RH),1)*4;
Blank_LH_inds = ones(length(data.Blank_SAP.LH),1)*5;
Blank_RH_inds = ones(length(data.Blank_SAP.RH),1)*6;
%% cell counting figure
summaryFigure = figure;
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
% SSP-SAP - plot each data point, connect L/R hemispheres
b2 = bar(3,data.SSP_SAP.LH_Mean,'FaceColor',colors('electric purple'));
bar(4,data.SSP_SAP.RH_Mean,'FaceColor',colors('electric purple'))
for aa = 1:length(data.SSP_SAP.LH)
    x = [SSP_LH_inds(aa,1),SSP_RH_inds(aa,1)];
    y = [data.SSP_SAP.LH(aa,1),data.SSP_SAP.RH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off', 'jitterAmount',0.25)
end
% Blank-SAP - plot each data point, connect L/R hemispheres
b3 = bar(5,data.Blank_SAP.LH_Mean,'FaceColor',colors('north texas green'));
bar(6,data.Blank_SAP.LH_Mean,'FaceColor',colors('north texas green'))
for aa = 1:length(data.Blank_SAP.LH)
    x = [Blank_LH_inds(aa,1),Blank_RH_inds(aa,1)];
    y = [data.Blank_SAP.LH(aa,1),data.Blank_SAP.RH(aa,1)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off', 'jitterAmount',0.25)
end
%% figure characteristics
ylabel('Cell density (per cubic mm of somatosensory cortical tissue)')
legend([b1,b2,b3],'Naive','SSP-SAP','Blank-SAP')
set(gca,'xtick',[1,2,3,4,5,6])
set(gca,'xticklabel',{'LH','RH','LH','RH (Rx)','LH','RH (Rx)'})
xtickangle(45)
axis square
xlim([0,7])
set(gca,'box','off')
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Cell Counts - Bilateral IOS' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'NADPH_DiaphoraseCellCounts_Bilateral_IOS']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'NADPH_DiaphoraseCellCounts_Bilateral_IOS'])
    %% statistical diary
    diaryFile = [dirpath 'NADPH_DiaphoraseCellCounts_Statistics_Bilateral_IOS.txt'];
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
    disp('GLME statistics for R/R Blank vs. SSP cell counts')
    disp('======================================================================================================================')
    disp(Stats.Comp.Stats)
end

end
