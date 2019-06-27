function [] = k()
% ESC 555 Homework 5
% Kevin Turner
% 6 Apr. 2019

% Problem Three
data = exp(-abs(randn(1000, 1)));
functionIndex = {@mean, @median, @std};
functionIDs = {'mean', 'median', 'std'};
nboot = 1000;

for a = 1:length(functionIndex)
    [output.(functionIDs{a}).data, output.(functionIDs{a}).value] = BootStrapESC555(nboot, functionIndex{a}, data);
    output.(functionIDs{a}).CI = prctile(output.(functionIDs{a}).data,[2.5 97.5]);
end

figure
for a = 1:length(functionIndex)
    p1 = plot(a,output.(functionIDs{a}).value,'k');
    set(p1, 'LineWidth', 2, 'MarkerSize', 10, 'Marker', 'x')
    hold on
    e1 = errorbar(a, output.(functionIDs{a}).value, output.(functionIDs{a}).value - output.(functionIDs{a}).CI(1),'b');
    set(e1, 'LineStyle', 'none'); 
    eline = get(e1, 'Children');
    set(eline,  'LineWidth', 2, 'Marker', 'x', 'MarkerSize', 4)
end

axis square
xlim([0.5 3.5])
legend([p1 e1], {'value', '95% CI'})
set(gca, 'XTick', [1 2 3], 'XTickLabel' ,functionIDs)
title('bootstrapped mean, median, and std')
ylabel('A.U.')

end

function [samples, meanVal] = BootStrapESC555(nboot, funcID, data)

samples = zeros(1,length(nboot));
for a = 1:nboot
    sampleIndex = randi([1,length(data)],1,length(data));
    sampledData = data(sampleIndex);
    sampledResult = funcID(sampledData);
    samples(a) = sampledResult;
end
meanVal = mean(samples);

end
