function [] = FigS1_nNOS(rootFolder,saveFigs,delim)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
path = [rootFolder delim 'Results_Turner'];
cd(path)
% Naive NADPH Nissl Histology
Naive_NADPH_Nissl = imread('Naive_NADPH_Nissl.tif');
% Blank NADPH Nissl Histology
Blank_NADPH_Nissl = imread('Blank_NADPH_Nissl.tif');
% SSP NADPH Nissl Histology
SSP_NADPH_Nissl = imread('SSP_NADPH_Nissl.tif');
% Figure S1
FigS1 = figure('Name','Figure S1');
subplot(2,2,1)
imshow(Naive_NADPH_Nissl)
axis image
title('a')
subplot(2,2,2)
imshow(Blank_NADPH_Nissl)
axis image
title('b')
subplot(2,2,3)
imshow(SSP_NADPH_Nissl)
axis image
title('c')
subplot(2,2,4)
text(0.5,0.5,'Cell counting example');
axis off
title('d')
% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'Figure Panels' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(FigS1,[dirpath 'FigS1']);
end