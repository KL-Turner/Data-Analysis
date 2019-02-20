function [] = CheckElectrodeQuality(rawDataFile, procDataFile)
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

load(rawDataFile)
load(procDataFile)

[animal, ~, fileDate, ~] = GetFileInfo(rawDataFile);
strDay = ConvertDate(fileDate);

electrodeCheck = figure;
ax1 = subplot(4,1,1);
plot((1:length(RawData.Data.Neural_LH)) / RawData.Notes.analogSamplingRate, RawData.Data.Neural_LH, 'k')
hold on;
plot((1:length(RawData.Data.Neural_RH)) / RawData.Notes.analogSamplingRate, RawData.Data.Neural_RH, 'b')
title({[animal ' ' strDay ' Electrode Check']; 'Raw Neural Data'})
legend('LH', 'RH')
ylabel('Volts')
set(gca, 'Ticklength', [0 0])
axis tight;
xlim([0 300])

ax2 = subplot(4,1,2);
plot((1:length(ProcData.Data.DeltaBand_Power.LH)) / ProcData.Notes.deltaBandSamplingRate, ProcData.Data.DeltaBand_Power.LH, 'k')
hold on;
plot((1:length(ProcData.Data.DeltaBand_Power.RH)) / ProcData.Notes.deltaBandSamplingRate, ProcData.Data.DeltaBand_Power.RH, 'b')
title('Delta Band Power')
ylabel('Power')
set(gca, 'Ticklength', [0 0])
axis tight;
xlim([0 300])

ax3 = subplot(4,1,3);
plot((1:length(ProcData.Data.GammaBand_Power.LH)) / ProcData.Notes.gammaBandSamplingRate, ProcData.Data.GammaBand_Power.LH, 'k')
hold on;
plot((1:length(ProcData.Data.GammaBand_Power.RH)) / ProcData.Notes.gammaBandSamplingRate, ProcData.Data.GammaBand_Power.RH, 'b')
title('Gamma Band Power')
ylabel('Power')
set(gca, 'Ticklength', [0 0])
axis tight;
xlim([0 300])

ax4 = subplot(4,1,4);
plot((1:length(ProcData.Data.MUA_Power.LH)) / ProcData.Notes.multiUnitSamplingRate, ProcData.Data.MUA_Power.LH, 'k')
hold on;
plot((1:length(ProcData.Data.MUA_Power.RH)) / ProcData.Notes.multiUnitSamplingRate, ProcData.Data.MUA_Power.RH, 'b')
title('MUA Power')
ylabel('Power')
xlabel('Time (sec)')
set(gca, 'Ticklength', [0 0])
axis tight;
xlim([0 300])

linkaxes([ax1 ax2 ax3 ax4], 'x')

[pathstr, ~, ~] = fileparts(cd);
dirpath = [pathstr '/Figures/Electrode Checks/'];

if ~exist(dirpath, 'dir') 
    mkdir(dirpath); 
end

savefig(electrodeCheck, [dirpath animal '_' strDay '_ElectrodeCheck']);

end