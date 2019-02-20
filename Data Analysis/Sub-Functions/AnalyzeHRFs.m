function [ComparisonData] = AnalyzeHRFs(procDataFiles, dataType, RestingBaselines, ComparisonData) 
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________
%
%   Inputs:
%
%   Outputs:
%________________________________________________________________________________________________________________________

global datastruct
fs = 30;

for ii = 1:size(procDataFiles, 1)
    disp(['Analyzing procDataFile ' num2str(ii) ' of ' num2str(size(procDataFiles, 1)) '...']); disp(' ') 
    procDataFile = procDataFiles(ii, :);
    load(procDataFile);
    [animal, ~, fileDate, fileID] = GetFileInfo(procDataFile);
    strDay = ConvertDate(fileDate);
    
    gammaPower = ProcData.Data.GammaBand_Power.(dataType);
    CBV = ProcData.Data.CBV.(dataType);
    
    [B, A] = butter(4, 2/(fs/2), 'low');
    datastruct.procCBV = filtfilt(B, A, (CBV - RestingBaselines.CBV.(dataType).(strDay)) ./ RestingBaselines.CBV.(dataType).(strDay));
    datastruct.procGamma = filtfilt(B, A, (gammaPower - RestingBaselines.GammaBand_Power.(dataType).(strDay)) ./ RestingBaselines.GammaBand_Power.(dataType).(strDay));
    
    % Initialize option parameters for fminsearch from Davis
    options = optimset('Display','none',...     % change iter-> none to display no output
        'FunValCheck','off',...  % check objective values
        'MaxFunEvals', 1000,...  % max number of function evaluations allowed
        'MaxIter', 1000,...      % max number of iteration allowed
        'TolFun',1e-8,...        % termination tolerance on the function value
        'TolX',1e-8,...          % termination tolerance on x
        'UseParallel','always'); % always use parallel computation
    
    % options = optimset('MaxFunEvals',2e3,'MaxIter',2e3,'TolFun',1e-7,'TolX',1e-7); % Aaron's options from Winder et al. 2017
    
    fit_fun = @fit_gammafunction_error;
    % x0 = [6, .28, .2];   % Patrick's plausible initial conditions from Mateo Data crunch
    x0 = [1, 0.75, -1e-3]; % Aaron's initial conditions from Winder et al. 2017
    [best_fit_params, fval] = fminsearch(fit_fun, x0, options);
    the_kernel = best_fit_params(3)*gampdf(datastruct.kernel_t, best_fit_params(1), best_fit_params(2));
    r_squared = 1 - fval;
    CBV_prediction = conv(datastruct.procGamma, the_kernel, 'full');
    r = corrcoef(datastruct.procCBV, CBV_prediction(1:length(datastruct.procCBV)));
    r = r(2,1);
    
    ComparisonData.HRFs.fileIDs{ii, 1} = fileID;
    ComparisonData.HRFs.Individual_kernels.(dataType)(ii, :) = the_kernel;
    ComparisonData.HRFs.Individual_r_squared.(dataType)(ii, 1) = r_squared;
    ComparisonData.HRFs.Individual_r.(dataType)(ii, 1) = r;
    ComparisonData.HRFs.time = datastruct.kernel_t;
end

save([animal '_ComparisonData.mat'], 'ComparisonData');

end

function [unexplained_var] = fit_gammafunction_error(x0)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
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

% this function returns the unexplained variance
global datastruct
kernel_maxt = 10;   % Seconds
datastruct.kernel_t = 0:1/30:kernel_maxt;

middle = length(datastruct.kernel_t):(length(datastruct.procGamma) - length(datastruct.kernel_t));
kernel = x0(3)*gampdf(datastruct.kernel_t, x0(1), x0(2));
CBV_predict_temp = conv(datastruct.procGamma, (kernel));
CBV_predict = CBV_predict_temp(1:length(datastruct.procGamma));

% r^2 on on test data
unexplained_var = var(detrend(CBV_predict(middle), 'constant') - detrend(datastruct.procCBV(middle), 'constant'))/var(detrend(datastruct.procCBV(middle), 'constant'));
end
