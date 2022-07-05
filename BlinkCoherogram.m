function [] = BlinkCoherogram(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Generate figure panel 8 for Turner_Gheres_Proctor_Drew_eLife2020
%________________________________________________________________________________________________________________________

%% variables for loops
resultsStruct = 'Results_BlinkCoherogram';
load(resultsStruct);
animalIDs = fieldnames(Results_BlinkCoherogram);
behavFields = {'Awake','Asleep','All'};
dataTypes = {'HbT','gamma','left','right'};
%% take data from each animal corresponding to the CBV-gamma relationship
for aa = 1:length(animalIDs)
    animalID = animalIDs{aa,1};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        data.(dataType).dummyCheck = 1;
        for cc = 1:length(behavFields)
            behavField = behavFields{1,cc};
            if isfield(data.(dataType),behavField) == false
                data.(dataType).(behavField).C = [];
                data.(dataType).(behavField).f = [];
                data.(dataType).(behavField).t = [];
                data.(dataType).(behavField).leadC = [];
                data.(dataType).(behavField).lagC = [];
                data.(dataType).(behavField).leadf = [];
                data.(dataType).(behavField).lagf = [];
            end
            C = Results_BlinkCoherogram.(animalID).(dataType).(behavField).C;
            meanC = mean(C(:,1:40*10),2);
            matC = meanC.*ones(size(C));
            msC = (C - matC);
            data.(dataType).(behavField).C = cat(3,data.(dataType).(behavField).C,msC);
            data.(dataType).(behavField).t = cat(1,data.(dataType).(behavField).t,Results_BlinkCoherogram.(animalID).(dataType).(behavField).t);
            data.(dataType).(behavField).f = cat(1,data.(dataType).(behavField).f,Results_BlinkCoherogram.(animalID).(dataType).(behavField).f);

            data.(dataType).(behavField).leadC = cat(2,data.(dataType).(behavField).leadC,Results_BlinkCoherogram.(animalID).(dataType).(behavField).leadC);
            data.(dataType).(behavField).lagC = cat(2,data.(dataType).(behavField).lagC,Results_BlinkCoherogram.(animalID).(dataType).(behavField).lagC);
            data.(dataType).(behavField).leadf = cat(1,data.(dataType).(behavField).leadf,Results_BlinkCoherogram.(animalID).(dataType).(behavField).leadf);
            data.(dataType).(behavField).lagf = cat(1,data.(dataType).(behavField).lagf,Results_BlinkCoherogram.(animalID).(dataType).(behavField).lagf');
        end
    end
end
for aa = 1:length(behavFields)
    behavField = behavFields{1,aa};
    data.gammaHbT.(behavField).C = cat(3,data.left.(behavField).C,data.right.(behavField).C);
    data.gammaHbT.(behavField).t = cat(1,data.left.(behavField).t,data.right.(behavField).t);
    data.gammaHbT.(behavField).f = cat(1,data.left.(behavField).f,data.right.(behavField).f);
    data.gammaHbT.(behavField).leadC = cat(1,data.left.(behavField).leadC,data.right.(behavField).leadC);
    data.gammaHbT.(behavField).lagC = cat(1,data.left.(behavField).lagC,data.right.(behavField).lagC);
    data.gammaHbT.(behavField).leadf = cat(2,data.left.(behavField).leadf,data.right.(behavField).leadf);
    data.gammaHbT.(behavField).lagf = cat(2,data.left.(behavField).lagf,data.right.(behavField).lagf);
end
%% mean
dataTypes = {'HbT','gamma','gammaHbT'};
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        data.(dataType).(behavField).meanC = mean(data.(dataType).(behavField).C,3);
        data.(dataType).(behavField).meanT = mean(data.(dataType).(behavField).t,1);
        data.(dataType).(behavField).meanF = mean(data.(dataType).(behavField).f,1);  
        data.(dataType).(behavField).meanLeadC = mean(data.(dataType).(behavField).leadC,2);
        data.(dataType).(behavField).meanLagC = mean(data.(dataType).(behavField).lagC,2);
        data.(dataType).(behavField).meanLeadF = mean(data.(dataType).(behavField).leadf,1);
        data.(dataType).(behavField).meanLagF = mean(data.(dataType).(behavField).lagf,1);
        
        data.(dataType).(behavField).stdLeadC = std(data.(dataType).(behavField).leadC,0,2)./sqrt(size(data.(dataType).(behavField).leadC,2));
        data.(dataType).(behavField).stdLagC = std(data.(dataType).(behavField).lagC,0,2)./sqrt(size(data.(dataType).(behavField).lagC,2));
    end
end
%%
summaryFigure = figure;
sgtitle('HbT Coherogram')
subplot(1,3,1)
Semilog_ImageSC(data.HbT.All.meanT,data.HbT.All.meanF,data.HbT.All.meanC,'y')
axis xy
axis square
title('All blinks')
subplot(1,3,2)
Semilog_ImageSC(data.HbT.Awake.meanT,data.HbT.Awake.meanF,data.HbT.Awake.meanC,'y')
axis xy
axis square
title('Awake blinks')
subplot(1,3,3)
Semilog_ImageSC(data.HbT.Asleep.meanT,data.HbT.Asleep.meanF,data.HbT.Asleep.meanC,'y')
axis xy
axis square
title('Asleep blinks')
%%
summaryFigure = figure;
sgtitle('Gammaband Coherogram')
subplot(1,3,1)
Semilog_ImageSC(data.gamma.All.meanT,data.gamma.All.meanF,data.gamma.All.meanC,'y')
axis xy
axis square
title('All blinks')
subplot(1,3,2)
Semilog_ImageSC(data.gamma.Awake.meanT,data.gamma.Awake.meanF,data.gamma.Awake.meanC,'y')
axis xy
axis square
title('Awake blinks')
subplot(1,3,3)
Semilog_ImageSC(data.gamma.Asleep.meanT,data.gamma.Asleep.meanF,data.gamma.Asleep.meanC,'y')
axis xy
axis square
title('Asleep blinks')
%%
summaryFigure = figure;
sgtitle('Gamma-HbT Coherogram')
subplot(1,3,1)
Semilog_ImageSC(data.gammaHbT.All.meanT,data.gammaHbT.All.meanF,data.gammaHbT.All.meanC,'y')
axis xy
axis square
title('All blinks')
subplot(1,3,2)
Semilog_ImageSC(data.gammaHbT.Awake.meanT,data.gammaHbT.Awake.meanF,data.gammaHbT.Awake.meanC,'y')
axis xy
axis square
title('Awake blinks')
subplot(1,3,3)
Semilog_ImageSC(data.gammaHbT.Asleep.meanT,data.gammaHbT.Asleep.meanF,data.gammaHbT.Asleep.meanC,'y')
axis xy
axis square
title('Asleep blinks')
%%
figure
s1 = semilogx(data.gamma.Awake.meanLeadF,data.gamma.Awake.meanLeadC,'r','LineWidth',2);
hold on
semilogx(data.gamma.Awake.meanLeadF,data.gamma.Awake.meanLeadC + data.gamma.Awake.stdLeadC,'r','LineWidth',0.5)
semilogx(data.gamma.Awake.meanLeadF,data.gamma.Awake.meanLeadC - data.gamma.Awake.stdLeadC,'r','LineWidth',0.5)

s2 = semilogx(data.gamma.Awake.meanLeadF,data.gamma.Awake.meanLagC,'b','LineWidth',2);
semilogx(data.gamma.Awake.meanLeadF,data.gamma.Awake.meanLagC + data.gamma.Awake.stdLagC,'b','LineWidth',0.5)
semilogx(data.gamma.Awake.meanLeadF,data.gamma.Awake.meanLagC - data.gamma.Awake.stdLagC,'b','LineWidth',0.5)
axis tight
title('gamma-band coherence around blinking, n = 22 mice')
ylabel('Coherence')
xlabel('Freq (Hz)')
legend([s1,s2],'leading blink +/- SEM','lagging blink +/- SEM')
axis square
%%
figure
s1 = semilogx(data.HbT.Awake.meanLeadF,data.HbT.Awake.meanLeadC,'r','LineWidth',2);
hold on
semilogx(data.HbT.Awake.meanLeadF,data.HbT.Awake.meanLeadC + data.HbT.Awake.stdLeadC,'r','LineWidth',0.5)
semilogx(data.HbT.Awake.meanLeadF,data.HbT.Awake.meanLeadC - data.HbT.Awake.stdLeadC,'r','LineWidth',0.5)

s2 = semilogx(data.HbT.Awake.meanLeadF,data.HbT.Awake.meanLagC,'b','LineWidth',2);
semilogx(data.HbT.Awake.meanLeadF,data.HbT.Awake.meanLagC + data.HbT.Awake.stdLagC,'b','LineWidth',0.5)
semilogx(data.HbT.Awake.meanLeadF,data.HbT.Awake.meanLagC - data.HbT.Awake.stdLagC,'b','LineWidth',0.5)
axis tight
title('HbT coherence around blinking, n = 22 mice')
ylabel('Coherence')
xlabel('Freq (Hz)')
legend([s1,s2],'leading blink +/- SEM','lagging blink +/- SEM')
axis square
%%
figure
s1 = semilogx(data.gammaHbT.Awake.meanLeadF,data.gammaHbT.Awake.meanLeadC,'r','LineWidth',2);
hold on
semilogx(data.gammaHbT.Awake.meanLeadF,data.gammaHbT.Awake.meanLeadC + data.gammaHbT.Awake.stdLeadC,'r','LineWidth',0.5)
semilogx(data.gammaHbT.Awake.meanLeadF,data.gammaHbT.Awake.meanLeadC - data.gammaHbT.Awake.stdLeadC,'r','LineWidth',0.5)

s2 = semilogx(data.gammaHbT.Awake.meanLeadF,data.gammaHbT.Awake.meanLagC,'b','LineWidth',2);
semilogx(data.gammaHbT.Awake.meanLeadF,data.gammaHbT.Awake.meanLagC + data.gammaHbT.Awake.stdLagC,'b','LineWidth',0.5)
semilogx(data.gammaHbT.Awake.meanLeadF,data.gammaHbT.Awake.meanLagC - data.gammaHbT.Awake.stdLagC,'b','LineWidth',0.5)
axis tight
title('gammaHbT coherence around blinking, n = 22 mice')
ylabel('Coherence')
xlabel('Freq (Hz)')
legend([s1,s2],'leading blink +/- SEM','lagging blink +/- SEM')
axis square
%% save figure(s)
% if saveFigs == true
%     dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'Pupil-Gamma Relationship' delim];
%     if ~exist(dirpath,'dir')
%         mkdir(dirpath);
%     end
%     imwrite(h1Img,[dirpath 'Fig8d_Awake.png'])
%     imwrite(h2Img,[dirpath 'Fig8d_NREM.png'])
%     imwrite(h3Img,[dirpath 'Fig8d_REM.png'])
%     savefig(summaryFigure,[dirpath 'PupilGamma_Relationship']);
%     set(summaryFigure,'PaperPositionMode','auto');
%     print('-painters','-dpdf','-fillpage',[dirpath 'PupilGamma_Relationship'])
% end

end
