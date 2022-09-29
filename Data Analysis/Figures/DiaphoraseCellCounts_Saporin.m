function [] = DiaphoraseCellCounts_Saporin(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

%% setup and pull data from excel sheet
msExcelFile = 'Diaphorase_Cell_Counts.xlsx';
[~,~,alldata] = xlsread(msExcelFile);
treatments = {'C57BL6J','SSP_SAP','Blank_SAP'};
% pre-allocate for concatenation
for aa = 1:length(treatments)
    treatment = treatments{1,aa};
    data.(treatment).LH = [];
    data.(treatment).RH = [];
end
% concatenate data for each treatment/hemishpere 
for bb = 2:size(alldata,1)
    treatment = alldata{bb,2};
    data.(treatment).LH = cat(2,data.(treatment).LH,alldata{bb,3});
    data.(treatment).RH = cat(2,data.(treatment).RH,alldata{bb,4});
end
% mean/std of each hemisphere
for cc = 1:length(treatments)
    treatment = treatments{1,cc};
    data.(treatment).LH_Mean = mean(data.(treatment).LH,2);
    data.(treatment).LH_StD = std(data.(treatment).LH,0,2);
    data.(treatment).RH_Mean = mean(data.(treatment).RH,2);
    data.(treatment).RH_StD = std(data.(treatment).RH,0,2);
end
% indeces for scatter plot
C57_LH_inds = ones(1,length(data.C57BL6J.LH))*1;
C57_RH_inds = ones(1,length(data.C57BL6J.RH))*2;
SSP_LH_inds = ones(1,length(data.SSP_SAP.LH))*3;
SSP_RH_inds = ones(1,length(data.SSP_SAP.RH))*4;
Blank_LH_inds = ones(1,length(data.Blank_SAP.LH))*5;
Blank_RH_inds = ones(1,length(data.Blank_SAP.RH))*6;
%% cell counting figure
summaryFigure = figure;
b1 = bar(1,data.C57BL6J.LH_Mean,'FaceColor',colors('sapphire'));
hold on
bar(2,data.C57BL6J.RH_Mean,'FaceColor',colors('sapphire'))
for aa = 1:length(data.C57BL6J.LH)
    x = [C57_LH_inds(1,aa),C57_RH_inds(1,aa)];
    y = [data.C57BL6J.LH(1,aa),data.C57BL6J.RH(1,aa)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','k','MarkerFaceColor',colors('sapphire'),'jitter','off', 'jitterAmount',0.25)
end
% SSP-SAP - plot each data point, connect L/R hemispheres
b2 = bar(3,data.SSP_SAP.LH_Mean,'FaceColor',colors('electric purple'));
bar(4,data.SSP_SAP.RH_Mean,'FaceColor',colors('electric purple'))
for aa = 1:length(data.SSP_SAP.LH)
    x = [SSP_LH_inds(1,aa),SSP_RH_inds(1,aa)];
    y = [data.SSP_SAP.LH(1,aa),data.SSP_SAP.RH(1,aa)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','k','MarkerFaceColor',colors('electric purple'),'jitter','off', 'jitterAmount',0.25)
end
% Blank-SAP - plot each data point, connect L/R hemispheres
b3 = bar(5,data.Blank_SAP.LH_Mean,'FaceColor',colors('north texas green'));
bar(6,data.Blank_SAP.LH_Mean,'FaceColor',colors('north texas green'))
for aa = 1:length(data.Blank_SAP.LH)
    x = [Blank_LH_inds(1,aa),Blank_RH_inds(1,aa)];
    y = [data.Blank_SAP.LH(1,aa),data.Blank_SAP.RH(1,aa)];
    plot(x,y,'-k')
    scatter(x(1),y(1),150,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off', 'jitterAmount',0.25)
    scatter(x(2),y(2),150,'MarkerEdgeColor','k','MarkerFaceColor',colors('north texas green'),'jitter','off', 'jitterAmount',0.25)
end
%% figure characteristics
sgtitle({'Mean NADPH diaphorase-stained cells per 70 \muM section';'within 1 mm diameter cortical ROI'})
ylabel('Cell Density (0.78 mm^2 of cortical tissue)')
legend([b1,b2,b3],'C57BL/6J','SSP-SAP','Blank-SAP')
set(gca,'xtick',[1,2,3,4,5,6])
set(gca,'xticklabel',{'LH','RH','LH','RH (Rx)','LH','RH (Rx)'})
xtickangle(45)
axis square
xlim([0,7])
set(gca,'box','off')
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'NADPH_Diaphorase_Cell_Counts']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'NADPH_Diaphorase_Cell_Counts'])
end

end
