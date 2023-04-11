function [RestData] = NormRestDataStruct_IOS(animalID,RestData,RestingBaselines,baselineType)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataTypes = fieldnames(RestData);
for dT = 1:length(dataTypes)
    dataType = char(dataTypes(dT));
    hemisphereDataTypes = fieldnames(RestData.(dataType));
    for hDT = 1:length(hemisphereDataTypes)
        hemDataType = char(hemisphereDataTypes(hDT));
        normData = [];
        if any(strcmp(dataType,{'HbT','HbO','HbR'})) == true % don't normalize, but keep field name for convenience
            RestData.(dataType).(hemDataType).NormData = RestData.(dataType).(hemDataType).data;
        else
            [uniqueDays,~,~] = GetUniqueDays_IOS(RestData.(dataType).(hemDataType).fileDates);
            for uD = 1:length(uniqueDays)
                date = uniqueDays{uD};
                strDay = ConvertDate_IOS(date);
                [~,dayInds] = GetDayInds_IOS(RestData.(dataType).(hemDataType).fileDates,date);
                disp(['Normalizing ' (hemDataType) ' ' (dataType) ' for ' (strDay) '...']); disp(' ')
                % calculate the baseline differently depending on data type
                dayData = RestData.(dataType).(hemDataType).data(dayInds);
                normDayData = cell(size(dayData));
                try
                    dayBaseline = RestingBaselines.(baselineType).(dataType).(hemDataType).(strDay).mean;
                catch
                    dayBaseline = 0;
                end
                for dD = 1:size(dayData,1)
                    cellBase = dayBaseline*ones(1,size(dayData{dD},2));
                    normDayData{dD} = dayData{dD}./cellBase - 1;
                end
                normData(dayInds) = normDayData;
            end
            RestData.(dataType).(hemDataType).NormData = normData';
        end
    end
end
save([animalID '_RestData.mat'],'RestData','-v7.3')