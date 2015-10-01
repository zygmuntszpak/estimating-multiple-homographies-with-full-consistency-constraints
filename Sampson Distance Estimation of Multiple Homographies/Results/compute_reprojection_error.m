%   Function: compute_reprojection_error
%
%   Computes the reprojection error (with euclidean distance metric) 
%   for each scene and each homography. The reprojection error is computed
%   for each estimation method by fixing the homography but allowing the 
%   corresponding points to vary.  This provides a gold standard measure of
%   how good a particular homography estimation method is.
%
%   Parameters:
%
%      listOfRandomScenes           - cell array of ground truth sceneData
%   								  structs
%
%      listOfEstimatedHomographiesForEachScene  - cell array of estimated
%                                                  homographies
%
%
%   Returns:
%
%    A cell array of symmetric transfer errors for each scene, 
%    and each homography.
%
%   See Also:
%
%  visualise_symmetric_transfer_errors_for_fixed_noise
%
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 15/5/2012
function listOfErrorsForEachSceneAndHomography  = ...
    compute_reprojection_error(...
    listOfRandomScenes,...
    listOfEstimatedHomographiesForEachScene)

numberOfScenes = length(listOfEstimatedHomographiesForEachScene);
numOfH = length(listOfEstimatedHomographiesForEachScene{1}) ;

% note the structure of the cell array, the first parameter indexes the
% scene, the second determines which plane/homography we are referring to
listOfErrorsForEachSceneAndHomography = cell(numberOfScenes,numOfH);

options = optimset('FinDiffType ','central','Diagnostics','on',...
    'Display','iter','Algorithm',{'levenberg-marquardt',10},...
    'Jacobian','off','DerivativeCheck','off');


for i=1:numberOfScenes
    % listOfRandomScenes containts the ground truth corresponding points
    sceneData = listOfRandomScenes{i};
    
    for k=1:numOfH
        x = sceneData(k).ptsInViewOne;
        xp = sceneData(k).ptsInViewTwo;
        % gold-standard method will optimise eta to find best
        % correspondences
        eta = [x(1,:)' ; x(2,:)'];
        % Ground Truth Homography (useful for debugging)
        %H = sceneData(k).homographies;
        listOfH = listOfEstimatedHomographiesForEachScene{i};
        H = listOfH{k};
        % cost_function_reprojection_error(eta,sceneData(k), H);
        f = @(bestEta)cost_function_reprojection_error(bestEta,...
                                                        sceneData(k),H);
        [etaFinal,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                            lsqnonlin(f,eta,[],[],options);
   
        listOfErrorsForEachSceneAndHomography{i,k}.errorsForEachPlane = ...
                                                    residual.^2; %resnorm;
    end
end
end


function errors = cost_function_reprojection_error(eta,sceneData, H)       
        numberOfPoints = length(sceneData.ptsInViewOne);
        x = reshape(eta,numberOfPoints,2)';
   
        n = numberOfPoints;
        % Generate the corresponding points xp in image 2
        xp = H * [ x; ones( 1, n ) ];
        xp = xp ./ repmat( xp(3,:), 3, 1 );
        xp = xp(1:2,:);
         
        errorsViewOne = sceneData.ptsInViewOne-x;
        errorsViewTwo = sceneData.ptsInViewTwo-xp; 
        combinedErrors = [errorsViewOne; errorsViewTwo];
        
        errors = [ combinedErrors(:) ];
end

