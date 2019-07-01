function [] = startup()

cd 'C:\Users\klt8\Documents\'
addpath(genpath('C:\Users\klt8\Documents\IOS-Data-Analysis\'))
addpath(genpath('C:\Users\klt8\Documents\Two-Photon-Data-Analysis\'))
addpath(genpath('C:\Users\klt8\Documents\Gheres-Manuscript\'))
addpath(genpath('C:\Users\klt8\Documents\TurnerFigs-SlowOscReview2019\'))
addpath(genpath('C:\Users\klt8\Documents\Condensor-Microphone\'))
addpath(genpath('C:\Users\klt8\Documents\Stream-RealSenseD435-Cam\'))

id = 'signal:filtfilt:ParseSOS';
warning('off', id)
clear id 

set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))

dbstop if error

end