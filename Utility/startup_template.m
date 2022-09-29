function [] = startup_template()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: runs automatially upon opening of matlab IDE. keep in Documents/MATLAB folder
%________________________________________________________________________________________________________________________

% go to code directory - generate and add file paths of desired projects
cd 'C:\Users\klt8\Documents\Github Repositories\'
% addpath(genpath('C:\Users\klt8\Documents\Github Repositories\IOS-Data-Analysis\'))
% addpath(genpath('C:\Users\klt8\Documents\Github Repositories\Two-Photon-Data-Analysis\'))
% addpath(genpath('C:\Users\klt8\Documents\Github Repositories\Stream-RealSenseD435-Cam\'))
% addpath(genpath('C:\Users\klt8\Documents\Winder_Echagarruga_Zhang_Drew_2017_Code-master\'))
% addpath(genpath('C:\Users\klt8\Documents\Github Repositories\Drew_Mateo_Turner_Yu_Kleinfeld_Neuron2020\'))
addpath(genpath('C:\Users\klt8\Documents\Github Repositories\Turner_Gheres_Proctor_Drew_eLife2020\'))
% addpath(genpath('C:\Users\klt8\Documents\Github Repositories\Norwood_Gharpure_Turner_Ferrer-Pistone_VanderWal_Drew_Manuscript2020\'))
% addpath(genpath('C:\Users\klt8\Documents\Github Repositories\Gheres_Turner_Drew_Manuscript2020\'))
% addpath(genpath('C:\Users\klt8\Documents\Github Repositories\Kederasetti_Turner_Drew_Costanzo_ManuscriptTBD\'))
addpath(genpath('C:\Users\klt8\Documents\Github Repositories\Turner_Drew_ManuscriptTBD\'))
addpath(genpath('C:\Users\klt8\Documents\Github Repositories\Ultrasound-Stimulation'))
% remove useless warnings
id = 'signal:filtfilt:ParseSOS';
warning('off',id)
clear id 
% add old style for figure toolbar because the new version sucks
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))
% stop and enter debug upon errors
dbstop if error
% email settings
server = 'smtp.gmail.com';
% Apply prefs and props
setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Server',server);
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.port','587');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
%|
%|            the cake is a lie.
%|            the cake is a lie.
%|            the cake is a lie.
%|            the cake is a lie.
%|            the cake is a lie.
%|
%|
%|              ,:/+/-                          
%|	            /A/              .,-=;//;-      
%|         .:/= ;IS/,    ,=/+EISALIETH:         
%|        -KEISALIETHECAKEISALIETH:.    -/HE    
%|   .,LIET ECAKEISA -ETHECAKE-     -+ECAKEI    
%|    .,ISAL;      +EISAL/,     =EISALIE;-      
%|  TH-  :EISALIETHECA.    .:ETHECAK:            
%|  CAKE,   +ETHEC/-.  ,;IETHECA,          -    
%|  LIETHE=,,---,.-ETHECAKEI:          ,ALIE    
%|  SALIETHECAKEISALIETH/.         :HECALE-                 
%|  HECAKEISALIETHECA,         ;CAKEIS=         
%|  THECAKEISALIETHEC.    .=IETHEC=             
%|  EISALIETHECAKEISA..;THECAKEI=          .:+  
%|  IETHECAKEISALIETHECAKE=           =+SALI    
%|  AKEISALIETHECAKEISA/.         =+LIET=       
%|  =+LIETHECAKEISALIET,      ,/LIE+:,          
%|    .;ISALIETHECAKEI=   ,/HEC+:;              
%|       .=+AKEISALIET+/+AKE+=.                 
%|           ,:/CAKEISAL/.                      
%|                ,.:=-.                                                                                                            
%|
end
