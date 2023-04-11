function [EventData] = NormEventDataStruct_IOS(animalID,EventData,RestingBaselines,baselineType)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
dataTypes = fieldnames(EventData);
for dT = 1:length(dataTypes)
    dataType = char(dataTypes(dT));
    hemisphereDataTypes = fieldnames(EventData.(dataType));
    for hDT = 1:length(hemisphereDataTypes)
        hemDataType = char(hemisphereDataTypes(hDT));
        behaviorFields = fieldnames(EventData.(dataType).(hemDataType));
        normData = [];
        for bF = 1:length(behaviorFields)
            behavField = char(behaviorFields(bF));
            if isempty(EventData.(dataType).(hemDataType).(behavField).data) == false
                if any(strcmp(dataType,{'HbT','HbO','HbR'})) == true % don't normalize, but keep field name for convenience
                    EventData.(dataType).(hemDataType).(behavField).NormData = EventData.(dataType).(hemDataType).(behavField).data;
                else
                    [uniqueDays,~,~] = GetUniqueDays_IOS(EventData.(dataType).(hemDataType).(behavField).fileDates);
                    for uD = 1:length(uniqueDays)
                        date = uniqueDays{uD};
                        strDay = ConvertDate_IOS(date);
                        [~,dayInds] = GetDayInds_IOS(EventData.(dataType).(hemDataType).(behavField).fileDates,date);
                        disp(['Normalizing ' (hemDataType) ' ' (dataType) ' ' (behavField) ' for ' (strDay) '...']); disp(' ')
                        % calculate the baseline differently depending on data type
                        try
                            dayBaseline = RestingBaselines.(baselineType).(dataType).(hemDataType).(strDay).mean;
                        catch
                            dayBaseline = 0;
                        end
                        % pre-allocate array and use for permutation
                        normDayData = EventData.(dataType).(hemDataType).(behavField).data(dayInds,:,:);
                        % permute norm_session_data to handle both matrix and array (squeeze
                        % causes a matrix dimension error if not permuted)
                        dayData = permute(normDayData,unique([2,1,ndims(normDayData)],'stable'));
                        for dD = 1:size(dayData,2)
                            normDayData(dD,:,:) = squeeze(dayData(:,dD,:))./(ones(size(dayData,1),1)*dayBaseline) - 1;
                        end
                        normData(dayInds,:,:) = normDayData;
                    end
                    EventData.(dataType).(hemDataType).(behavField).NormData = normData';
                end
            end
        end
    end
end
save([animalID '_EventData.mat'],'EventData','-v7.3')