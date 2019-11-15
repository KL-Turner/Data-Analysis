function [DataStruct,FiltArray,Criteria] = SelectBehavioralEvents2(DataStruct,Behavior,hem)
%   [DataStruct,FiltArray,Criteria] = SelectBehavioralEvents(DataStruct,Behavior)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Sets definitions of various categories of behavior and
%   selects the data that match each category
%   
%_______________________________________________________________
%   PARAMETERS:             
%                   DataStruct - [struct] structure containing data and
%                   behavioral information for 'Behavior'.
%
%                   Behavior - [string] behavioral category desired               
%_______________________________________________________________
%   RETURN:                     
%                   DataStruct - [struct] structure containing data and 
%                       behavioral information about selected category.
%
%                   FiltArray - [array] logical array used to filter the
%                   events that satisfy the conditions in "Criteria"
%
%                   Criteria - [Struc] contains the conditions used to
%                   filter the behavioral data.
%_______________________________________________________________
BehaviorBufferTime = 4; % Seconds, 4 seconds allows CBV signals to stabalize
switch Behavior
    case 'Rest'
        Criteria.Fieldname = {'PuffDistance','duration'};
        Criteria.Comparison = {'gt','gt'};
        Criteria.Value = {5,6+BehaviorBufferTime};
        Criteria.Min_Duration = 5;
        Criteria.Min_PuffDist = 5;
    case 'EndofWhisk'
        Criteria.Fieldname = {'WhiskDurs','duration','PuffDistance'};
        Criteria.Comparison = {'gt','gt','gt'};
        Criteria.Value = {1,2,5};
        DataStruct = DataStruct.rest;
    case 'VW'
        Criteria.Fieldname = {'duration','duration','restTime','PuffDistance'};
        Criteria.Comparison = {'lt','gt','gt','gt'};
        Criteria.Value = {2,0.5,BehaviorBufferTime,5};
        DataStruct = DataStruct.whisk;
    case 'Str'
        Criteria.Fieldname = {'duration','restTime','PuffDistance'};
        Criteria.Comparison = {'gt','gt','gt'};
        Criteria.Value = {2,BehaviorBufferTime,5};
        DataStruct = DataStruct.whisk;
    case 'Contra'
        Criteria.Fieldname = {'Name','whiskScore_Pre', 'movementScore_Pre','PuffDistance'};
        Criteria.Comparison = {'equal','lt','lt','gt'};
        Criteria.Value = {'Contra',0.01, 0.01, 5};
        DataStruct = DataStruct.stim;
    case 'Ipsi'
        Criteria.Fieldname = {'Name','whiskScore_Pre', 'movementScore_Pre','PuffDistance'};
        Criteria.Comparison = {'equal','lt','lt','gt'};
        Criteria.Value = {'Ipsi',0.01, 0.01, 5};
        DataStruct = DataStruct.stim;
    case 'Control'
        Criteria.Fieldname = {'Name','whiskScore_Pre', 'movementScore_Pre','PuffDistance'};
        Criteria.Comparison = {'equal','lt','lt','gt'};
        Criteria.Value = {'Control',0.01, 0.01, 5};
        DataStruct = DataStruct.stim;
end

FiltArray = FilterEvents2(DataStruct,Criteria,hem,Behavior);
