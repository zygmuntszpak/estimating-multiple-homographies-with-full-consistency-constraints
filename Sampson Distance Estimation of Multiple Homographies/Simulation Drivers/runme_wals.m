%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to demonstrate a multiple
% homography estimation method---weighted alternating least squares---that
% enforces certain inter-homography constraint originally proposed in
%
% P. Chen and D. Suter. 
% Rank Constraints for Homographies over Two Views: Revisiting the Rank 
% Four Constraint
% Int J Comput Vis (2009) 81: 205–225
% DOI 10.1007/s11263-008-0167-z
%
% and evaluated in
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
% some consistency constraints. 
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

% number of homographies (must be >= 3)
numOfH = 4 ;

% this should be set to 1 if we were to pass initial homographies into the 
% compute_aml_estimates method that were computed using the FNS estimation
% method. Since we are using DLT estimates, we will want to utilise
% homography covariance matrices based on the DLT estimate (hence we set 
% the variable to zero). 
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
    listOfInitialH_DLT = compute_dlt_estimates(sceneData); 
    % Use Chen's estimation procedure with Chen initialisation and 
    % DLT seed and DLT covariances
    [listOfEstimatedH, diagnostic] = ... = ...
                compute_chen_estimates(sceneData,listOfInitialH_DLT,0);
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
wals.listOfErrors = listOfErrors;
wals.meanRmsError = meanRmsError;
wals.rmsErrorForEachPlane = rmsErrorForEachPlane;
wals.sigma = sigma;
wals.numberOfScenes = numberOfScenes;
wals.nPoints = nPoints;
wals.nHomographies = numOfH;
wals.randseed = randseed;
wals.randnseed = randnseed;

wals.meanIter = meanIter;
wals.medianIter = medianIter;
wals.varIter = varIter;
wals.stdIter = stdIter;
wals.meanTiming = meanTiming;
wals.medianTiming = medianTiming;
wals.stdTiming = stdTiming;
wals.startResiduals = startResiduals;
wals.endResiduals = endResiduals;

wals.listOfDiagnosticForEachScene = listOfDiagnosticForEachScene;

save experiment1-wals wals


fprintf('\n \n Root-Mean-Square Errors for each Homography\n')
wals.rmsErrorForEachPlane
fprintf('Mean Root Mean Square Error \n')
wals.meanRmsError
 
  