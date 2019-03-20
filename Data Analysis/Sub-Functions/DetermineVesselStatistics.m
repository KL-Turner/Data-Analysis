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
    cd(['/Volumes/Blackberry/' animal]);
    load([animal '_ComparisonData.mat']);
    
    for b = 1:length(ComparisonData.tblVals.vesselIDs)
        animalID{x,1} = animal;
        vesselID{x,1} = ComparisonData.tblVals.vesselIDs{b,1};
        baselineDiam{x,1} = ComparisonData.tblVals.baselines{b,1};
        minutesPerVessel{x,1} = ComparisonData.tblVals.timePerVessel{b,1};
        minutesPerCond1{x,1} = ComparisonData.tblVals.C1events{b,1};
        minutesPerCond2{x,1} = ComparisonData.tblVals.C2events{b,1};
        minutesPerCond3{x,1} = ComparisonData.tblVals.C3events{b,1};
        x = x + 1;
    end 
end

figure;


