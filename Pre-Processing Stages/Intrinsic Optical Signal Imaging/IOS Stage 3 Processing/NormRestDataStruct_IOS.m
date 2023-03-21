function [RestData] = NormRestDataStruct_IOS(RestData,RestingBaselines,baselineType)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%
% Purpose: Normalizes data structures based on resting values
%________________________________________________________________________________________________________________________

animalID = RestingBaselines.setDuration.baselineFileInfo.animalID;
dataTypes = fieldnames(RestData);
for dT = 1:length(dataTypes)
    dataType = char(dataTypes(dT));
    if any(strcmp(dataType,{'HbT','HbO','HbR'})) == true
        hemisphereDataTypes = fieldnames(RestData.(dataType));
        if strcmp(dataType,'LH_HbT') == false && strcmp(dataType,'RH_HbT') == false
            for hDT = 1:length(hemisphereDataTypes)
                hemDataType = char(hemisphereDataTypes(hDT));
                if isfield(RestData.(dataType).(hemDataType),'whisk')
                    behaviorFields = fieldnames(RestData.(dataType).(hemDataType));
                    for bF = 1:length(behaviorFields)
                        behField = char(behaviorFields(bF));
                        if ~isempty(RestData.(dataType).(hemDataType).(behField).data)
                            NormData = RestData.(dataType).(hemDataType).(behField).data;
                            [uniqueDays,~,~] = GetUniqueDays_IOS(RestData.(dataType).(hemDataType).(behField).fileDates);
                            for uD = 1:length(uniqueDays)
                                date = uniqueDays{uD};
                                strDay = ConvertDate_IOS(date);
                                [~,dayInds] = GetDayInds_IOS(RestData.(dataType).(hemDataType).(behField).fileDates,date);
                                disp(['Normalizing ' (hemDataType) ' ' (dataType) ' ' (behField) ' for ' (strDay) '...']); disp(' ')
                                % calculate the baseline differently depending on data type
                                if iscell(RestData.(dataType).(hemDataType).(behField).data)
                                    dayData = RestData.(dataType).(hemDataType).(behField).data(dayInds);
                                    normDayData = cell(size(dayData));
                                    dayBaseline = RestingBaselines.(baselineType).(dataType).(hemDataType).(strDay);
                                    for dD = 1:size(dayData,1)
                                        cellBase = dayBaseline*ones(1,size(dayData{dD},2));
                                        normDayData{dD} = dayData{dD}./cellBase - 1;
                                    end
                                    NormData(dayInds) = normDayData;
                                else
                                    try
                                        dayBaseline = RestingBaselines.(baselineType).(dataType).(hemDataType).(strDay).mean;
                                    catch
                                        if strcmp(hemDataType,'LH_gammaBandPower') == true
                                            dayBaseline = RestingBaselines.(baselineType).cortical_LH.gammaBandPower.(strDay).mean;
                                        elseif strcmp(hemDataType,'RH_gammaBandPower') == true
                                            dayBaseline = RestingBaselines.(baselineType).cortical_RH.gammaBandPower.(strDay).mean;
                                        else
                                            dayBaseline = NaN;
                                        end
                                    end
                                    % pre-allocate array and use for permutation
                                    normDayData = RestData.(dataType).(hemDataType).(behField).data(dayInds,:,:);
                                    % permute norm_session_data to handle both matrix and array (squeeze
                                    % causes a matrix dimension error if not permuted)
                                    dayData = permute(normDayData,unique([2,1,ndims(normDayData)],'stable'));
                                    for dD = 1:size(dayData,2)
                                        normDayData(dD,:,:) = squeeze(dayData(:,dD,:))./(ones(size(dayData,1),1)*dayBaseline) - 1;
                                    end
                                    NormData(dayInds,:,:) = normDayData;
                                end
                                RestData.(dataType).(hemDataType).(behField).NormData = NormData;
                            end
                        end
                    end
                else
                    NormData = RestData.(dataType).(hemDataType).data;
                    [uniqueDays,~,~] = GetUniqueDays_IOS(RestData.(dataType).(hemDataType).fileDates);
                    for uD = 1:length(uniqueDays)
                        date = uniqueDays{uD};
                        strDay = ConvertDate_IOS(date);
                        [~,dayInds] = GetDayInds_IOS(RestData.(dataType).(hemDataType).fileDates,date);
                        disp(['Normalizing ' (hemDataType) ' ' (dataType) ' for ' (strDay) '...']); disp(' ')
                        % calculate the baseline differently depending on data type
                        if iscell(RestData.(dataType).(hemDataType).data)
                            dayData = RestData.(dataType).(hemDataType).data(dayInds);
                            normDayData = cell(size(dayData));
                            try
                                dayBaseline = RestingBaselines.(baselineType).(dataType).(hemDataType).(strDay).mean;
                            catch
                                if strcmp(hemDataType,'LH_gammaBandPower') == true
                                    dayBaseline = RestingBaselines.(baselineType).cortical_LH.gammaBandPower.(strDay).mean;
                                elseif strcmp(hemDataType,'RH_gammaBandPower') == true
                                    dayBaseline = RestingBaselines.(baselineType).cortical_RH.gammaBandPower.(strDay).mean;
                                else
                                    dayBaseline = NaN;
                                end
                            end
                            for dD = 1:size(dayData,1)
                                cellBase = dayBaseline*ones(1,size(dayData{dD},2));
                                normDayData{dD} = dayData{dD}./cellBase - 1;
                            end
                            NormData(dayInds) = normDayData;
                        else
                            dayBaseline = RestingBaselines.(dataType).(hemDataType).(strDay);
                            % pre-allocate array and use for permutation
                            normDayData = RestData.(dataType).(hemDataType).data(dayInds,:,:);
                            % permute norm_session_data to handle both matrix and array (squeeze
                            % causes a matrix dimension error if not permuted)
                            dayData = permute(normDayData,unique([2,1,ndims(normDayData)],'stable'));

                            for dD = 1:size(dayData,2)
                                normDayData(dD,:,:) = squeeze(dayData(:,dD,:))./(ones(size(dayData,1),1)*dayBaseline) - 1;
                            end
                            NormData(dayInds,:,:) = normDayData;
                        end
                        RestData.(dataType).(hemDataType).NormData = NormData;
                    end
                end
            end
        end
        RestData.(dataType).(hemDataType).NormData = NormData;
    end
end
save([animalID '_RestData.mat'],'RestData','-v7.3')

end
