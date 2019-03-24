function DetermineVesselStatistics_SlowOscReview2019(ComparisonData)
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

animalIDs = fields(ComparisonData);
x = 1;
for a = 1:length(animalIDs)
    animalID = animalIDs{a};  
    for b = 1:length(ComparisonData.(animalID).tblVals.vesselIDs)
        animal{x,1} = animalID;
        vesselID{x,1} = ComparisonData.(animalID).tblVals.vesselIDs{b,1};
        baselineDiam{x,1} = ComparisonData.(animalID).tblVals.baselines{b,1};
        minutesPerVessel{x,1} = ComparisonData.(animalID).tblVals.timePerVessel{b,1};
        x = x + 1;
    end
end

T = table(animal, vesselID, baselineDiam, minutesPerVessel, 'VariableNames', {'Animal_ID', 'Vessel_ID', ...
    'Baseline_diameter_um', 'Total_minutes_per_Vessel'});
figure('Name', 'Individual vessel imaging information', 'NumberTitle', 'off')
u = uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
set(u,'ColumnWidth',{125})
pause(1)
end
