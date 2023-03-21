function [] = StimEvoked_Pulse_Figures(rootFolder,saveFigs,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

cd([rootFolder delim 'Results_Turner'])
resultsStruct = 'Results_Evoked_Pulse';
load(resultsStruct);
cd(rootFolder)
% loop variables
groups = {'SSP_SAP','Blank_SAP'};
solenoids = {'LPadSol','RPadSol','AudSol'};
dataTypes = {'HbT','timeVector'};
% extract the analysis results
for aa = 1:length(groups)
    group = groups{1,aa};
    animalIDs = fieldnames(Results_Evoked_Pulse.(group));
    for bb = 1:length(animalIDs)
        animalID = animalIDs{bb,1};
        for cc = 1:length(solenoids)
            solenoid = solenoids{1,cc};
            data.(group).(solenoid).dummCheck = 1;
            for dd = 1:length(dataTypes)
                dataType = dataTypes{1,dd};
                if isfield(data.(group).(solenoid),dataType) == false
                    data.(group).(solenoid).(dataType) = [];
                end
            end
            data.(group).(solenoid).HbT = cat(1,data.(group).(solenoid).HbT,Results_Evoked_Pulse.(group).(animalID).Stim.adjBarrels.(solenoid).HbT);
            data.(group).(solenoid).timeVector = cat(1,data.(group).(solenoid).timeVector,Results_Evoked_Pulse.(group).(animalID).Stim.adjBarrels.(solenoid).timeVector);
        end
    end
end
% pair stimulation types with hemisphere
for aa = 1:length(groups)
    group = groups{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        % contra
        data.(group).contra.(dataType) = data.(group).LPadSol.(dataType);
        data.(group).contra.(['mean_' dataType]) = mean(data.(group).LPadSol.(dataType),1);
        data.(group).contra.(['stdErr_' dataType]) = std(data.(group).LPadSol.(dataType),1)./sqrt(size(data.(group).LPadSol.(dataType),1));
        % ipsi
        data.(group).ipsi.(dataType) = data.(group).RPadSol.(dataType);
        data.(group).ipsi.(['mean_' dataType]) = mean(data.(group).RPadSol.(dataType),1);
        data.(group).ipsi.(['stdErr_' dataType]) = std(data.(group).RPadSol.(dataType),1)./sqrt(size(data.(group).RPadSol.(dataType),1));
        % aud
        data.(group).aud.(dataType) = data.(group).AudSol.(dataType);
        data.(group).aud.(['mean_' dataType]) = mean(data.(group).AudSol.(dataType),1);
        data.(group).aud.(['stdErr_' dataType]) = std(data.(group).AudSol.(dataType),1)./sqrt(size(data.(group).AudSol.(dataType),1));
    end
end
% figure
figure;
x0 =  xline(0,'k');
hold on
x1 = xline(5,'r');
% Blank-SAP
p1 = plot(data.Blank_SAP.contra.mean_timeVector,data.Blank_SAP.contra.mean_HbT,'color',colors('north texas green'),'LineWidth',2);
plot(data.Blank_SAP.contra.mean_timeVector,data.Blank_SAP.contra.mean_HbT + data.Blank_SAP.contra.stdErr_HbT,'color',colors('north texas green'),'LineWidth',0.5)
plot(data.Blank_SAP.contra.mean_timeVector,data.Blank_SAP.contra.mean_HbT - data.Blank_SAP.contra.stdErr_HbT,'color',colors('north texas green'),'LineWidth',0.5)
% SSP-SAP
p2 = plot(data.SSP_SAP.contra.mean_timeVector,data.SSP_SAP.contra.mean_HbT,'color',colors('electric purple'),'LineWidth',2);
plot(data.SSP_SAP.contra.mean_timeVector,data.SSP_SAP.contra.mean_HbT + data.SSP_SAP.contra.stdErr_HbT,'color',colors('electric purple'),'LineWidth',0.5)
plot(data.SSP_SAP.contra.mean_timeVector,data.SSP_SAP.contra.mean_HbT - data.SSP_SAP.contra.stdErr_HbT,'color',colors('electric purple'),'LineWidth',0.5)
title('Pulse contralateral whisker stimlation')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-stimulus time (s)')
legend([p1,p2,x0,x1],'Blank-SAP','SSP-SAP','stimOn','stimOff')
set(gca,'box','off')
xlim([-2,10])
axis square
