%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to demonstrate the new gold standard multiple 
% homography estimation method---fully consistent bundle adjustment---
% described in the paper:
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
% correspondences. Finally, it estimates the homographies and enforces
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

% number of homographies (must be >= 2)
numOfH = 4 ;

% This variable determines whether we factorise the initial homographies
% into latent variables using a procedure proposed by P. Chen et al. or 
% a procedure proposed by Chojnacki et al. The method proposed by Chen et
% al. requires a minimum of 3 homographies, whereas the method proposed
% by Chojnacki et al requires a minimum of two homographies. 
% if numOfH < 3 then useChojnackiInitialisation *must* be set to 1
useChojnackiInitialisation = 1; 

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
listOfRandomScenesWithNoise = add_noise_to_scenes(listOfRandomScenes,sigma);

listOfEstimatedHomographiesForEachScene = cell(1,numberOfScenes);
listOfDiagnosticForEachScene = cell(1,numberOfScenes);
for i = 1:numberOfScenes
    sceneData = listOfRandomScenesWithNoise{i};
    listOfInitialH_DLT = compute_dlt_estimates(sceneData); 
    % enforce consistency constraints on homographies that were estimated 
    % using DLT. Use fully consistent bundle adjustment to refine
    % the homographies. This method should produce the best result.
     [listOfEstimatedH, diagnostic] = ...
                  compute_ba_estimates(sceneData,listOfInitialH_DLT,...
                                              useChojnackiInitialisation);
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
 %                                 listOfEstimatedHomographiesForEachScene);                             

 
% finally we compute the root-mean-square errors
 [meanRmsError, rmsErrorForEachPlane] =  ...
                        compute_mean_root_mean_square_error(listOfErrors);
 
% gather summaries of diagnostics
[meanIter, medianIter, varIter,stdIter, ...
    meanTiming, medianTiming,varTiming,stdTiming,...
    startResiduals, endResiduals] = ...
              compute_diagnostic_statistics(listOfDiagnosticForEachScene);
            

% write output of experiment to data structure
bajoint.listOfErrors = listOfErrors;
bajoint.meanRmsError = meanRmsError;
bajoint.rmsErrorForEachPlane = rmsErrorForEachPlane;
bajoint.sigma = sigma;
bajoint.numberOfScenes = numberOfScenes;
bajoint.nPoints = nPoints;
bajoint.nHomographies = numOfH;
bajoint.randseed = randseed;
bajoint.randnseed = randnseed;

bajoint.meanIter = meanIter;
bajoint.medianIter = medianIter;
bajoint.varIter = varIter;
bajoint.stdIter = stdIter;
bajoint.meanTiming = meanTiming;
bajoint.medianTiming = medianTiming;
bajoint.stdTiming = stdTiming;
bajoint.startResiduals = startResiduals;
bajoint.endResiduals = endResiduals;

bajoint.listOfDiagnosticForEachScene = listOfDiagnosticForEachScene;

save experiment1-bajoint bajoint

fprintf('\n \n Root-Mean-Square Errors for each Homography\n')
bajoint.rmsErrorForEachPlane
fprintf('Mean Root Mean Square Error \n')
bajoint.meanRmsError

 
  