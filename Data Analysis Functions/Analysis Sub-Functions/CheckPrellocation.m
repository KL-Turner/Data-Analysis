function [dataStruct] = CheckPrellocation(dataStruct,variables)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Check if folders already exist for pre-allocated variable
%________________________________________________________________________________________________________________________

for aa = 1:length(variables)
    if isfield(dataStruct,variables{1,aa}) == false
        dataStruct.(variables{1,aa}) = [];
    end
end

end