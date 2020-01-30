%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%________________________________________________________________________________________________________________________

whiskAnimalIDs = {'T72', 'T73', 'T74', 'T75', 'T76'};
x = 1;
for a = 1:length(whiskAnimalIDs)
    animal = whiskAnimalIDs{a};
    cd(['I:\' animal '\Combined Imaging\']);
    load([animal '_ComparisonData.mat']);
    
    for b = 1:length(ComparisonData.tblVals.vesselIDs)
        animalID{x,1} = animal;
        vesselID{x,1} = ComparisonData.tblVals.vesselIDs{b,1};
        baselineDiam{x,1} = ComparisonData.tblVals.baselines{b,1};
        minutesPerVessel{x,1} = ComparisonData.tblVals.timePerVessel{b,1};
        eventsPerCond1{x,1} = ComparisonData.tblVals.C1events{b,1};
        eventsPerCond2{x,1} = ComparisonData.tblVals.C2events{b,1};
        eventsPerCond3{x,1} = ComparisonData.tblVals.C3events{b,1};
        x = x + 1;
    end 
end

T = table(animalID, vesselID, baselineDiam, minutesPerVessel, eventsPerCond1, eventsPerCond2, eventsPerCond3,...
    'VariableNames', {'Animal_ID', 'Vessel_ID', 'Baseline_diameter_um', 'Total_minutes_per_Vessel', 'Total_05_2_sec_events', 'Total_2_5_sec_events', 'Total_5_10_sec_events'});
figure('Name', 'Individual vessel imaging information', 'NumberTitle', 'off')
u = uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
TString = evalc('disp(T)');
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
FixedWidth = get(0,'FixedWidthFontName');
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex','FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);
set(u,'ColumnWidth',{75})
