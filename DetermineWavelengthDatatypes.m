function [dataTypes] = DetermineWavelengthDatatypes(imagingWavelengths,iteration)

if strcmp(imagingWavelengths,'Red, Green, & Blue') == true
    if iteration == 1
        dataTypes = {'green','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    else
        dataTypes = {'green','blue','red','HbT','HbO','HbR','GCaMP','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    end
elseif strcmp(imagingWavelengths,'Red, Lime, & Blue') == true
    if iteration == 1
        dataTypes = {'lime','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    else
        dataTypes = {'lime','blue','red','HbT','HbO','HbR','GCaMP','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    end
elseif strcmp(imagingWavelengths,'Green & Blue') == true
    if iteration == 1
        dataTypes = {'green','blue','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    else
        dataTypes = {'green','blue','HbT','GCaMP','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    end
elseif strcmp(imagingWavelengths,'Lime & Blue') == true
    if iteration == 1
        dataTypes = {'lime','blue','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    else
        dataTypes = {'lime','blue','HbT','GCaMP','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    end
elseif strcmp(imagingWavelengths,'Green') == true
    if iteration == 1
        dataTypes = {'green','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    else
        dataTypes = {'green','HbT','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    end
elseif strcmp(imagingWavelengths,'Lime') == true
    if iteration == 1
        dataTypes = {'lime','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    else
        dataTypes = {'lime','HbT','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    end
elseif strcmp(imagingWavelengths,'Blue') == true
    if iteration == 1
        dataTypes = {'blue','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    else
        dataTypes = {'blue','HbT','cortical_LH','cortical_RH','hippocampus','EMG','pupil'};
    end
end