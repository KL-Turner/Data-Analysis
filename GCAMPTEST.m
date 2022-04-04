clear; clc;
curDir = cd;
animalIDs = {'T253','T254','T255'};
% animalIDs = {'T224'};
blueChannel = [];
greenChannel = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    cd([animalID '/Combined Imaging/']);
%     cd([animalID '/Bilateral Imaging/']);
    procDataFileIDs = ls('*_ProcData.mat');
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        disp(['Pulling data from ' animalID ' ProcData file number ' num2str(bb) '/' num2str(size(procDataFileIDs,1))]); disp(' ')
        load(procDataFileID)
        blueChannel = cat(2,blueChannel,rescale(ProcData.data.GCaMP7s.adjLH,0,1),rescale(ProcData.data.GCaMP7s.adjRH,0,1));
        greenChannel = cat(2,greenChannel,rescale(ProcData.data.CBV.adjLH,0,1),rescale(ProcData.data.CBV.adjRH,0,1));
    end
    cd(curDir)
end
figure;
h1 = histogram2(greenChannel,blueChannel);
h1Vals = h1.Values;
s = pcolor(h1Vals');
s.FaceColor = 'interp';
set(s,'EdgeColor','none');
figure
h1 = histogram2(greenChannel,blueChannel);
x = sum(h1.BinCounts,1);
y = repmat(x,100,1);
z = h1.BinCounts./y;
figure; imagesc(z); axis xy; caxis([0,0.1]);
[f,gof] = fit(greenChannel',blueChannel','poly1');
figure; plot(f,greenChannel,blueChannel)
%% current GCaMP correction
% [~,fileDate,~] = GetFileInfo_IOS(procDataFileID);
% strDay = ConvertDate_IOS(fileDate);
% LH_scale = RestingBaselines.manualSelection.CBV.adjLH.(strDay).mean/RestingBaselines.manualSelection.GCaMP7s.adjLH.(strDay).mean;
% RH_scale = RestingBaselines.manualSelection.CBV.adjRH.(strDay).mean/RestingBaselines.manualSelection.GCaMP7s.adjRH.(strDay).mean;
% ProcData.data.GCaMP7s.corLH = (ProcData.data.GCaMP7s.adjLH./ProcData.data.CBV.adjLH)*LH_scale;
% ProcData.data.GCaMP7s.corRH = (ProcData.data.GCaMP7s.adjRH./ProcData.data.CBV.adjRH)*RH_scale;
