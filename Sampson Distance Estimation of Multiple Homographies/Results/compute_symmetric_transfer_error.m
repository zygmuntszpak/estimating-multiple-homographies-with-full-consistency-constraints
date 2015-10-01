%   Function: compute_symmetric_transfer
%
%   Computes the symmetric transfer error for each scene and each 
%   homography. The symmetric transfer error is  computed by applying
%   the estimated homography the to ground truth corresponding points.
%
%   Parameters:
%
%      listOfRandomScenes  		- cell array of ground truth sceneData
%   							  structs
%
%      listOfEstimatedHomographiesForEachScene  - cell array of estimated
%                                                 homographies
%
%
%   Returns:
%
%    A cell array of symmetric transfer errors for each scene, and each 
%    homography.
%
%   See Also:
%
%  visualise_symmetric_transfer_errors_for_fixed_noise
%
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 15/5/2012

function listOfErrorsForEachSceneAndHomography  = ...
                compute_symmetric_transfer_error(listOfRandomScenes,...
                                   listOfEstimatedHomographiesForEachScene)
numberOfScenes = length(listOfEstimatedHomographiesForEachScene);
numOfH = length(listOfEstimatedHomographiesForEachScene{1}) ;

% note the structure of the cell array, the first parameter indexes the
% scene, the second determines which plane/homography we are referring to
listOfErrorsForEachSceneAndHomography = cell(numberOfScenes,numOfH);

for i=1:numberOfScenes
    % listOfRandomScenes containts the ground truth corresponding points
    sceneData = listOfRandomScenes{i};
    for k=1:numOfH
        x = sceneData(k).ptsInViewOne;
        xp = sceneData(k).ptsInViewTwo;
        % Ground Truth Homography
        % useful for debugging... if we plug in the ground truth ...
        % homography we should get
        % zero symmetric transfer error.
        % H = sceneData(k).homographies;
        
        listOfH = listOfEstimatedHomographiesForEachScene{i};
        H = listOfH{k};
        n = length(x);
        % convert to homogenous coordinates
        x = [ x; ones( 1, n ) ];
        xp = [ xp; ones( 1, n ) ];
        xpp = H * x; xpp = xpp ./ repmat( xpp(3,:), 3, 1 );
        xx = H \ xp; xx = xx ./ repmat( xx(3,:), 3, 1 );
        
        transfer_error_one=  x(1:2,:)-xx(1:2,:);
        transfer_error_one = transfer_error_one(:).^2;
        transfer_error_two = xp(1:2,:)-xpp(1:2,:);
        transfer_error_two = transfer_error_two(:).^2;
        
        listOfErrorsForEachSceneAndHomography{i,k}.errorsForEachPlane = ...
                                    transfer_error_one + transfer_error_two;

    end
        
end

end

