function [ComparisonData] = AnalyzeAwakeHRF(targetMinutes, dataType, procDataFiles, RestingBaselines, ComparisonData)
%________________________________________________________________________________________________________________________
% Written by Davis Haselden & Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: 
%________________________________________________________________________________________________________________________
%
%   Inputs: 
%
%   Outputs: //
%________________________________________________________________________________________________________________________
fileNumber = targetMinutes / 5;
clear global testData
global testData
for fN = 1:fileNumber
    % Handle file number, load each successive file
    procDataFile = procDataFiles(fN, :);
    load(procDataFile);
    [animal, ~, fileDate, ~] = GetFileInfo(procDataFile);
    strDay = ConvertDate(fileDate);
    
    %% CBV data - normalize by rest, lowpass filer, detrend
    CBV = ProcData.Data.CBV.(dataType);
    normCBV = (CBV - RestingBaselines.CBV.(dataType).(strDay)) ./ (RestingBaselines.CBV.(dataType).(strDay));
    [B, A] = butter(4, 2 / (30 / 2), 'low');   % 2 Hz lowpass
    testData.filtCBV(fN, :) = detrend(filtfilt(B, A, normCBV), 'constant');
    
    %% Gamma band power data - normalize by rest, lowpass filer, detrend
    gammaPower = ProcData.Data.GammaBand_Power.LH;
    normGamma = (gammaPower - RestingBaselines.GammaBand_Power.(dataType).(strDay)) ./ (RestingBaselines.GammaBand_Power.(dataType).(strDay));
    testData.filtGamma(fN, :) = detrend(filtfilt(B, A, normGamma), 'constant');
end

%% Gamma function and minimizing error using fminsearch 
testData.dt = 1/30;   % 30 Hz sampling rate
% Initialize option parameters for fminsearch
options = optimset('Display','none',...     % change iter-> none to display no output
    'FunValCheck','off',...  % check objective values
    'MaxFunEvals', 2000,...  % max number of function evaluations allowed
    'MaxIter', 2000,...      % max number of iteration allowed
    'TolFun',1e-8,...        % termination tolerance on the function value
    'TolX',1e-8,...          % termination tolerance on x
    'UseParallel','always'); % always use parallel computation

fit_fun = @fit_gammafunction_error;
x0 = [4, .75, -1e-3];   % Initial conditions from Winder et al. 2017
fit_gammafunction_error(x0);   %[A0, B0, amplitude]);
[best_fit_params, fval] = fminsearch(fit_fun, x0, options);

% report optimized values
testData.best_x = best_fit_params;
testData.the_best_kernel = best_fit_params(3)*gampdf(testData.kernel_t, testData.best_x(1), testData.best_x(2));
testData.fval = fval;

testData.halfmax = min(testData.the_best_kernel) / 2;
testData.FWHM = testData.dt*length(find(testData.the_best_kernel > testData.halfmax));
[~, peak_point] = max(abs(testData.the_best_kernel)); % find the peak of the HRF
testData.kernel_peak = peak_point*testData.dt; %time of the peak
dataStruct_out = testData;

ComparisonData.HRFs.(dataType).dataStruct_out = testData;
save([animal '_ComparisonData.mat'], 'ComparisonData');
end

function [unexplained_var] = fit_gammafunction_error(X)
% this function returns the unexplained variance
global testData

kernel_maxt=10;
kernel_t=0:testData.dt:kernel_maxt;
testData.kernel_t = kernel_t;

for ii = 1:size(testData.filtGamma,1)
    middle = length(kernel_t):(size(testData.filtGamma, 2) - length(kernel_t));
    the_kernel = X(3)*gampdf(kernel_t,X(1),X(2));
    CBV_predict_temp = conv(testData.filtGamma(ii,:),(the_kernel));
    testData.CBV_predict(ii,:) = CBV_predict_temp(1:size(testData.filtGamma, 2));
    % r^2 on on test data
    unexplained_var_all(ii) = var(testData.CBV_predict(ii, middle) - testData.filtCBV(ii, middle)) / var(testData.filtCBV(ii, middle));
end

unexplained_var = sum(unexplained_var_all);
testData.kernel_t = kernel_t;
testData.the_kernel = the_kernel;
end
