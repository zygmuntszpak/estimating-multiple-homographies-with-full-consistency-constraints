%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to demonstrate the most popular single homography
% estimation method---the direct linear transform---discussed in many 
% papers, including:
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
% consistency constraints.
% 
% The performance of the method is reported using the symmetric
% transfer error, or the gold-standard reprojection error. 
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
for i = 1:numberOfScenes
    sceneData = listOfRandomScenesWithNoise{i};
    listOfInitialH_DLT = compute_dlt_estimates(sceneData); 
    listOfEstimatedHomographiesForEachScene{i} =  listOfInitialH_DLT ;  
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
 

% write output of experiment to data structure
dlt.listOfErrors = listOfErrors;
dlt.meanRmsError = meanRmsError;
dlt.rmsErrorForEachPlane = rmsErrorForEachPlane;
dlt.sigma = sigma;
dlt.numberOfScenes = numberOfScenes;
dlt.nPoints = nPoints;
dlt.nHomographies = numOfH;
dlt.randseed = randseed;
dlt.randnseed = randnseed;

% We do not compute diagnostics such as running times etc., because this
% method is for all intents and purposes "instantaneous" and non-iterative


save experiment1-dlt dlt

fprintf('\n \n Root-Mean-Square Errors for each Homography\n')
dlt.rmsErrorForEachPlane
fprintf('Mean Root Mean Square Error \n')
dlt.meanRmsError

 
  