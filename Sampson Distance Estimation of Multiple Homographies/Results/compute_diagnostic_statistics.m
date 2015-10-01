%   Function: compute_diagnostic_statistics
%
%   Generates summary statistics of various diagnostics that are
%   useful in ascertaining the performance and stability of estimation
%   methods.
%
%   Parameters:
%
%      listOfDiagnosticForEachScene   - a data structure containing squared
%                                       diagnostic information such as
%                                       the residuals for each iteration
%                                       of an estimation method, the 
%                                       running time, the total number
%                                       of iterations etc.
%
%
%
%   Returns:
%
%   Summary statistics computer from raw diagnostic information 
%
%   See Also:
%
%
%
%  Zygmunt L. Szpak (c) 2014
%  Last modified 16/09/2014
function [meanIter, medianIter, varIter,stdIter, ...
    meanTiming, medianTiming,varTiming,stdTiming,...
    startResiduals, endResiduals] = ...
                compute_diagnostic_statistics(listOfDiagnosticForEachScene)

nScene = length(listOfDiagnosticForEachScene);
iterations = zeros(1,nScene);
timings = zeros(1,nScene);
startResiduals = zeros(1,nScene);
endResiduals = zeros(1,nScene);
for iScene = 1:nScene
    diagnostic = listOfDiagnosticForEachScene{iScene};
    iterations(iScene) = diagnostic.iterations;
    runningTime(iScene) = diagnostic.runningTime;
    residualHistory = diagnostic.residualHistory;
    startResiduals(iScene) = residualHistory(1);
    endResiduals(iScene) = residualHistory(end);
end

meanIter = mean(iterations);
medianIter = median(iterations);
varIter = var(iterations);
stdIter = std(iterations);
meanTiming = mean(runningTime);
medianTiming = median(runningTime);
varTiming = var(runningTime);
stdTiming = std(runningTime);

end

