%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to demonstrate the multiple homography estimation
% method described in the paper:
%
% Chojnacki, W.; Szpak, Z.L.; Brooks, M.J.; van den Hengel, A.; ,
% "Multiple Homography Estimation with Full Consistency Constraints,"
% Digital Image Computing: Techniques and Applications (DICTA), 2010
% International Conference on , vol., no., pp.480-485, 1-3 Dec. 2010
% doi: 10.1109/DICTA.2010.87
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
    % using DLT (if useFNSCovariances = 0), 
    % or FNS ( if useFNSCovariance == 1).    	
    [listOfEstimatedH, diagnostic] = ...
                compute_aml_estimates(sceneData,...
                                     useChojnackiInitialisation,...
                                     useFNSCovariances);
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
covaml.listOfErrors = listOfErrors;
covaml.meanRmsError = meanRmsError;
covaml.rmsErrorForEachPlane = rmsErrorForEachPlane;
covaml.sigma = sigma;
covaml.numberOfScenes = numberOfScenes;
covaml.nPoints = nPoints;
covaml.nHomographies = numOfH;
covaml.randseed = randseed;
covaml.randnseed = randnseed;

covaml.meanIter = meanIter;
covaml.medianIter = medianIter;
covaml.varIter = varIter;
covaml.stdIter = stdIter;
covaml.meanTiming = meanTiming;
covaml.medianTiming = medianTiming;
covaml.stdTiming = stdTiming;
covaml.startResiduals = startResiduals;
covaml.endResiduals = endResiduals;

covaml.listOfDiagnosticForEachScene = listOfDiagnosticForEachScene;

save experiment1-covaml covaml

fprintf('\n \n Root-Mean-Square Errors for each Homography\n')
covaml.rmsErrorForEachPlane
fprintf('Mean Root Mean Square Error \n')
covaml.meanRmsError
                    
      