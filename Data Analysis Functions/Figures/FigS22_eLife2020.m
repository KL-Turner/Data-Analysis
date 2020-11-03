function [AnalysisResults] = FigS22_eLife2020(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate figure panel S22 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% set-up and process data
if isfield(AnalysisResults,'PixelDrift') == true
    x = AnalysisResults.PixelDrift.x;
    filtCatCement_cementData = AnalysisResults.PixelDrift.filtCatCement_cementData;
    Cement_modelFit_Y = AnalysisResults.PixelDrift.Cement_modelFit_Y;
    catLH_CBVdata = AnalysisResults.PixelDrift.catLH_CBVdata;
    catRH_CBVdata = AnalysisResults.PixelDrift.catRH_CBVdata;
    Cement_modelFit_flip = AnalysisResults.PixelDrift.Cement_modelFit_flip;
    LH_adjCatC_CBVdata = AnalysisResults.PixelDrift.LH_adjCatC_CBVdata;
    RH_adjCatC_CBVdata = AnalysisResults.PixelDrift.RH_adjCatC_CBVdata;
else
    dataDir = [rootFolder '\Summary Figures and Structures\Pixel Drift Correction\'];
    cd(dataDir)
    % character list of ProcData files
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
    catLH_CBVdata = [];
    catRH_CBVdata = [];
    catCement_cementData = [];
    % load the processed CBV/cement data from each file and concat it into one array
    for aa = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(aa,:);
        load(procDataFileID,'-mat')
        samplingRate = ProcData.notes.CBVCamSamplingRate;
        LH_CBVdata = ProcData.data.CBV.LH;
        RH_CBVdata = ProcData.data.CBV.RH;
        Cement_cementData = ProcData.data.CBV.Cement;
        catLH_CBVdata = horzcat(catLH_CBVdata,LH_CBVdata); %#ok<*AGROW>
        catRH_CBVdata = horzcat(catRH_CBVdata,RH_CBVdata);
        catCement_cementData = horzcat(catCement_cementData,Cement_cementData);
    end
    % establish whether a slow exponential trend exists for the data
    [B,A] = butter(3,0.01/(samplingRate/2),'low');
    filtCatCement_cementData = filtfilt(B,A,catCement_cementData);
    x = ((1:length(filtCatCement_cementData))/samplingRate)';
    % create a weight vector for the trend
    Cement_weightVec = ones(1,length(x));
    Cement_secondHalfMean = mean(filtCatCement_cementData(floor(length(filtCatCement_cementData/2)):end));
    for bb = 1:length(Cement_weightVec)
        if filtCatCement_cementData(bb) > Cement_secondHalfMean
            Cement_weightVec(bb) = 10;
        end
    end
    % compare weighted models
    Cement_modelFit = fit(x,filtCatCement_cementData','exp2','Weight',Cement_weightVec);
    Cement_modelFit_Y = Cement_modelFit(x);
    Cement_modelFit_norm = (Cement_modelFit_Y - min(Cement_modelFit_Y))./min(Cement_modelFit_Y);
    Cement_modelFit_flip = 1 - Cement_modelFit_norm;
    % apply exponential correction to original data
    LH_adjCatC_CBVdata = catLH_CBVdata.*Cement_modelFit_flip';
    RH_adjCatC_CBVdata = catRH_CBVdata.*Cement_modelFit_flip';
    % update analysis structure
    AnalysisResults.PixelDrift.x = x;
    AnalysisResults.PixelDrift.filtCatCement_cementData = filtCatCement_cementData;
    AnalysisResults.PixelDrift.Cement_modelFit_Y = Cement_modelFit_Y;
    AnalysisResults.PixelDrift.catLH_CBVdata = catLH_CBVdata;
    AnalysisResults.PixelDrift.catRH_CBVdata = catRH_CBVdata;
    AnalysisResults.PixelDrift.Cement_modelFit_flip = Cement_modelFit_flip;
    AnalysisResults.PixelDrift.LH_adjCatC_CBVdata = LH_adjCatC_CBVdata;
    AnalysisResults.PixelDrift.RH_adjCatC_CBVdata = RH_adjCatC_CBVdata;
    % save results
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end
%% Fig. S22
summaryFigure = figure('Name','FigS22 (a-f)');
sgtitle('Figure S22 - Turner et al. 2020')
%% [S22a] data and exponential fit for cement ROI
ax1 = subplot(2,3,1);
p4a = plot(x,filtCatCement_cementData,'color',colors_eLife2020('deep carrot orange'),'LineWidth',1);
hold on
p4b = plot(x,Cement_modelFit_Y,'color',colors_eLife2020('electric purple'),'LineWidth',1);
title('[S22a] Cement ROI pixel drift')
xlabel('Time (s)')
ylabel('Pixel intensity (12-bit)')
legend([p4a,p4b],'Cement ROI','Cement ROI fit')
axis tight
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [S22b] original left hemisphere
ax2 = subplot(2,3,2);
plot(x,catLH_CBVdata,'color','r','LineWidth',0.5)
title({'[S22b] Left hemisphere','original data'})
xlabel('Time (s)')
ylabel('Pixel intensity (12-bit)')
axis tight
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [S22c] original right hemisphere
ax3 = subplot(2,3,3);
plot(x,catRH_CBVdata,'color','c','LineWidth',0.5)
title({'[S22c] Right hemisphere','original data'})
xlabel('Time (s)')
ylabel('Pixel intensity (12-bit)')
axis tight
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [S22d] correction profile
ax4 = subplot(2,3,4);
plot(x,Cement_modelFit_flip,'color',colors_eLife2020('electric purple'),'LineWidth',1)
title('[S22d] Correction profile')
xlabel('Time (s)')
ylabel('Correction profile (%)')
axis tight
axis square
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [S22e] left hemisphere correction
ax5 = subplot(2,3,5);
plot(x,catLH_CBVdata,'color','r','LineWidth',0.5)
hold on
p8 = plot(x,LH_adjCatC_CBVdata,'color',colors_eLife2020('electric purple'),'LineWidth',0.5);
title({'[S22e] Left hemisphere','original vs. corrected data'})
xlabel('Time (s)')
ylabel('Pixel intensity (12-bit)')
legend(p8,'Corrected left ROI')
axis tight
axis square
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [S22f] right hemisphere correction
ax6 = subplot(2,3,6);
plot(x,catRH_CBVdata,'color','c','LineWidth',0.5)
hold on
p9 = plot(x,RH_adjCatC_CBVdata,'color',colors_eLife2020('electric purple'),'LineWidth',0.5);
title({'[S22f] Right hemisphere','original vs. corrected data'})
xlabel('Time (s)')
ylabel('Pixel intensity (12-bit)')
legend(p9,'Corrected right ROI')
axis tight
axis square
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'FigS22']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'FigS22'])
end

end
