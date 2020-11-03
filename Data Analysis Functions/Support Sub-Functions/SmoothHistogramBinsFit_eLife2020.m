function [xCurve,yCurve] = SmoothHistogramBinsFit_eLife2020(histData,bins,fitFunc)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Fits a spline to histogram to create a curve
%   Take from https://www.mathworks.com/help/curvefit/examples/smoothing-a-histogram.html
%________________________________________________________________________________________________________________________

curveFig = figure;
h = histfit(histData,bins,fitFunc);
yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', yt/numel(histData))
curve = h(2);
xCurve = curve.XData;
yCurve = curve.YData/numel(histData);
close(curveFig)

end
