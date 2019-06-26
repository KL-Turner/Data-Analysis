function [SpectrogramData] = NormalizeSpectrograms(animal, dataType, RestingBaselines, SpectrogramData)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: //
%________________________________________________________________________________________________________________________
%
%   Inputs: //
%
%   Outputs: //
%________________________________________________________________________________________________________________________

uniqueFileIDs = SpectrogramData.(dataType).FileIDs;

for ii = 1:length(uniqueFileIDs)
    fileID = uniqueFileIDs{ii, :};
    date = ConvertDate(fileID);
    disp(['Normalizing ' (dataType) ' spectrogram ' num2str(ii) ' of ' num2str(length(uniqueFileIDs)) '...']); disp(' ')
    baseLine1 = RestingBaselines.Spectrograms.(dataType).OneSec.(date);
    baseLine5 = RestingBaselines.Spectrograms.(dataType).FiveSec.(date);

    S1 = SpectrogramData.(dataType).OneSec.S{ii};
    S5 = SpectrogramData.(dataType).FiveSec.S{ii};
    
    hold_matrix1 = baseLine1.*ones(size(S1));
    hold_matrix5 = baseLine5.*ones(size(S5));
    
    S1_Norm = (S1 - hold_matrix1) ./ hold_matrix1;
    S5_Norm = (S5 - hold_matrix5) ./ hold_matrix5;

    SpectrogramData.(dataType).OneSec.S_Norm{ii, 1} = S1_Norm;
    SpectrogramData.(dataType).FiveSec.S_Norm{ii, 1} = S5_Norm;
end

disp('Saving...'); disp(' ')
save([animal '_SpectrogramData.mat'], '-v7.3', 'SpectrogramData');

end
