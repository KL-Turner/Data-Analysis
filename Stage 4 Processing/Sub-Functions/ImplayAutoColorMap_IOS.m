function [] = ImplayAutoColorMap(image, min, max)

handle = implay(image);

handle.Visual.ColorMap.UserRange = 1; 
handle.Visual.ColorMap.UserRangeMin = min; 
handle.Visual.ColorMap.UserRangeMax = max;

end