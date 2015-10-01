%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to demonstrate what is often (incorrectly) assumed
% to be the gold standard multiple homography estimation method---separate
% bundle adjustment without consistency constraints---discussed in the
% paper:
%
% Z. L. Szpak, W. Chojnacki, A. Eriksson, and A. van den Hengel. 
% Sampson distance based joint estimation of multiple homographies with
% uncalibrated cameras. 
% Comput. Vis. Image Underst., 125:200-213, 2014. 
% http://dx.doi.org/10.1016/j.cviu.2014.04.008
%
% The script generates random synthetic scenes containing the designated 
% number of planes. It then generates corresponding points between two
% views of the planes, and adds user-specified noise to the 
% correspondences. Finally, it estimates the homographies without enforcing
% full consistency constraints. 
% 
% The performance of the method is reported using the symmetric
% transfer error, or the gold-standard reprojection error. Various
% other diagnostic information is also recorded. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randseed = 11;
randnseed = 19;
rand( 'seed', randseed );
randn( 'seed',randnseed );

% number of homographies (must be >= 1)
numOfH = 4 ;

% number of desired data points (expect to get half of this on average)
nPoints = 50;

% number of random scenes
numberOfScenes = 10;



% generate multiple random planar scenes
listOfRandomScenes = generate_multiple_random_planarscenes_new(...
                            numberOfScenes,numOfH,nPoints);

% pick one of the generated scenes at random
sceneData = listOfRandomScenes{randi(numberOfScenes)};
% and visualise it
visualise_groundtruth_scene_and_plot(sceneData)

% set a noise level for measurements with zero mean Gaussian noise and
% standard deviation sigma
sigma = 1 ;

% add noise to each of the scenes
listOfRandomScenesWithNoise = ...
                            add_noise_to_scenes(listOfRandomScenes,sigma);

listOfEstimatedHomographiesForEachScene = cell(1,numberOfScenes);
listOfDiagnosticForEachScene = cell(1,numberOfScenes);
for i = 1:numberOfScenes
    sceneData = listOfRandomScenesWithNoise{i};
    listOfInitialH_DLT = compute_dlt_estimates(sceneData); 
    % refine the DLT homographies using bundle adjustment that does
    % not enforce consistency constraints
    [listOfEstimatedH, diagnostic] = ...
                compute_ba_separate_estimates(...
                                     sceneData,listOfInitialH_DLT);
    listOfEstimatedHomographiesForEachScene{i} =  listOfEstimatedH; 
    listOfDiagnosticForEachScene{i} =  diagnostic;     
end

% you can compute the symmetric transfer error (faster)
listOfErrors = ...
                compute_symmetric_transfer_error(listOfRandomScenes,...
                                  listOfEstimatedHomographiesForEachScene);
 
 % or you can compute the gold-standard reprojection error (slower)                             
 %listOfErrors = ...
 %               compute_reprojection_error(listOfRandomScenes,...
 %                                listOfEstimatedHomographiesForEachScene);                             

 
% finally we compute the root-mean-square errors
 [meanRmsError, rmsErrorForEachPlane] =  ...
                        compute_mean_root_mean_square_error(listOfErrors);
 
% gather summaries of diagnostics
[meanIter, medianIter, varIter,stdIter, ...
    meanTiming, medianTiming,varTiming,stdTiming,...
    startResiduals, endResiduals] = ...
               compute_diagnostic_statistics(listOfDiagnosticForEachScene);
                    
% write output of experiment to data structure
basep.listOfErrors = listOfErrors;
basep.meanRmsError = meanRmsError;
basep.rmsErrorForEachPlane = rmsErrorForEachPlane;
basep.sigma = sigma;
basep.numberOfScenes = numberOfScenes;
basep.nPoints = nPoints;
basep.nHomographies = numOfH;
basep.randseed = randseed;
basep.randnseed = randnseed;

basep.meanIter = meanIter;
basep.medianIter = medianIter;
basep.varIter = varIter;
basep.stdIter = stdIter;
basep.meanTiming = meanTiming;
basep.medianTiming = medianTiming;
basep.stdTiming = stdTiming;
basep.startResiduals = startResiduals;
basep.endResiduals = endResiduals;

basep.listOfDiagnosticForEachScene = listOfDiagnosticForEachScene;

save experiment1-basep basep

fprintf('\n \n Root-Mean-Square Errors for each Homography\n')
basep.rmsErrorForEachPlane
fprintf('Mean Root Mean Square Error \n')
basep.meanRmsError
                   
