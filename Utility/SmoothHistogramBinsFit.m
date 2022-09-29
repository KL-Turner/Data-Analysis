function [xCurve,yCurve] = SmoothHistogramBinsFit(histData,bins,fitFunc)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Fits a spline to histogram to create a curve
%   Take from https://www.mathworks.com/help/curvefit/examples/smoothing-a-histogram.html
%________________________________________________________________________________________________________________________

% h = 2*iqr(histData)*(length(histData)^(-1/3));
% bins = round((max(histData) - min(histData))/h);
curveFig = figure;
h = histfit(histData,bins,fitFunc);
axis tight
yt = get(gca,'YTick');
set(gca,'YTick',yt,'YTickLabel',yt/numel(histData))
curve = h(2);
xCurve = curve.XData;
yCurve = curve.YData/numel(histData);
close(curveFig)

end
