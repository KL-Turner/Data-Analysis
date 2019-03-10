function [fList] = GetFuncDependencies(functionName)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Outputs all dependencies of a given function, including subfunction dependencies. It also verifies this by
%            looping through each subfunction and further checking that it doesn't have any dependencies.
%________________________________________________________________________________________________________________________
%
%   Inputs: functionName (string) of the function to be tested. i.e. GetFuncDependencies('Example.m')
%
%   Outputs: Table showing the results.
%
%   Last Revised: March 9th, 2019
%________________________________________________________________________________________________________________________

% Detect OS and set delimeter to correct value. This typically isn't necessary in newer versions of Matlab.
if isunix
    delimiter = '/';
elseif ispc
    delimiter = '\';
else
    disp('Platform not currently supported');
    return
end

% Attempt to use Matlab's codetool feature to detect function dependencies.
try
    disp('Analyzing function dependencies...');
    [fList, ~] = matlab.codetools.requiredFilesAndProducts(functionName);
catch
    % Catch the instance where the filename does not exist or was spelled incorrectly.
    disp(['Matlab function ' functionName ' does not appear to exist or to be included in the current filepath(s)']);
    return
end

% Remove the first result, as it corresponds to the function itself.
if size(fList, 2) == 1
    fList = [];
else
    fList = fList(1, 2:end);
end

% If the function has dependencies and was not empty, continue.
searchForDep = true;
while searchForDep == true
    if ~isempty(fList)
        pause(1)
        % Cycle through each identified function dependency, and check if it itself has dependencies.
        for x = 1:size(fList, 2)
            indDepFuncPath = fList{1, x};
            funcDelimiters = strfind(indDepFuncPath, delimiter);
            indDepFuncName = string(strip(indDepFuncPath(funcDelimiters(end):end), delimiter));
            [depSubFunc, ~] = matlab.codetools.requiredFilesAndProducts(indDepFuncName);
            % If the dependent function doesn't itself have any dependencies, coninue to the next.
            % If it does, go through the dependent function's dependencies and see if they are already listed.
            if size(depSubFunc, 2) == 1
                continue
            else
                funcExists = false;
                depSubFunc = depSubFunc(1, 2:end);
                for y = 1:size(depSubFunc, 2)
                    while funcExists == false
                        for z = 1:size(fList, 2)
                            if strcmp(depSubFunc{1, y}, fList{1, z})
                                funcExists = true;
                                break
                            end
                        end
                        if funcExists == false
                            % If a new dependency is found, add it to the list of functions and restart the loop.
                            fList = horzcat(fList, depSubFunc{y, 1}); %#ok<AGROW>
                            x = 1; %#ok<FXSET>
                        end
                        break
                    end
                end
            end
        end
        % This is only reached once all functions have been accounted for and have no further dependencies.
        searchForDep = false;
    else
        % The function has no identified dependencies.
        pause(2)
        disp(['Matlab function ' functionName ' is self-contained and has no user-defined dependencies'])
        searchForDep = false;
    end
end

% Find the unique functions listed.
uniqueFuncPaths = unique(fList);
allDepFuncNames = cell(size(uniqueFuncPaths, 2), 1);
allDepFuncPaths = cell(size(uniqueFuncPaths, 2), 1);
for x = 1:size(uniqueFuncPaths, 2)
    allDepFuncPaths{x, 1} = uniqueFuncPaths{1, x};
    funcDelimiters = strfind(allDepFuncPaths{x, 1}, delimiter);
    allDepFuncNames{x, 1} = char(strip(allDepFuncPaths{x, 1}(funcDelimiters(end):end), delimiter));
end

tableVals = sortrows(horzcat(allDepFuncNames, allDepFuncPaths));
fileNames = tableVals(:, 1);
filePaths = tableVals(:, 2);
T = table(fileNames, filePaths, 'VariableNames', {'File_names', 'Full_file_path'});
figure('Name', ['Function dependencies for ' functionName], 'NumberTitle', 'off')
u = uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
TString = evalc('disp(T)');
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
FixedWidth = get(0,'FixedWidthFontName');
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex','FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);
set(u,'ColumnWidth',{200})

end
