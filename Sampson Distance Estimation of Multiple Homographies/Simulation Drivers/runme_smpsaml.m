%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to demonstrate the new gold standard multiple 
% homography estimation method---fully consistent sampson distance---
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

% this should be set to 1 if we want the  compute_aml_estimates to compute
% initial homography using the FNS scheme as well as covariance matrices
% associated with the FNS estimate. Otherwise, the method will use
% DLT to estimate the homographies together with DLT covariance matrices.
useFNSCovariances = 0;

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
    % enforce consistency constraints on homographies that were estimated 
    % using DLT (if useFNSCovariances = 0), or FNS
    % ( if useFNSCovariance == 1).  Enfroce full consistency using sampson
    % distance on data points to refine the homographies. This method
    % should produce the second best result. However, with the current
    % optimisation method (Levenberg-Marquardt) this algorithm takes a
    % long time to converge.
    % It is a matter of future work to resolve this issue...
    % If your noise level is not too large, and speed is critical, then you
    % should use the method in the "runme_aml" script. Alternatively, if
    % precision is of utmost importance and noise is high  you can use 
    % fully consistent bundle adjustment method in the "runme_bajoint"
    % script. 
    [listOfEstimatedH, diagnostic] = ...
                compute_robustsampsonaml_estimates(...
                                     sceneData,...
                                        useChojnackiInitialisation,...
                                                  useFNSCovariances,Inf);
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
smpsaml.listOfErrors = listOfErrors;
smpsaml.meanRmsError = meanRmsError;
smpsaml.rmsErrorForEachPlane = rmsErrorForEachPlane;
smpsaml.sigma = sigma;
smpsaml.numberOfScenes = numberOfScenes;
smpsaml.nPoints = nPoints;
smpsaml.nHomographies = numOfH;
smpsaml.randseed = randseed;
smpsaml.randnseed = randnseed;

smpsaml.meanIter = meanIter;
smpsaml.medianIter = medianIter;
smpsaml.varIter = varIter;
smpsaml.stdIter = stdIter;
smpsaml.meanTiming = meanTiming;
smpsaml.medianTiming = medianTiming;
smpsaml.stdTiming = stdTiming;
smpsaml.startResiduals = startResiduals;
smpsaml.endResiduals = endResiduals;

smpsaml.listOfDiagnosticForEachScene = listOfDiagnosticForEachScene;

save experiment1-smpsaml smpsaml

fprintf('\n \n Root-Mean-Square Errors for each Homography\n')
smpsaml.rmsErrorForEachPlane
fprintf('Mean Root Mean Square Error \n')
smpsaml.meanRmsError


 
  