%   Function: add_outliers_to_scenes
%
%   Given a cell array of random scenes and ground truth corresponding points, the function
%   adds randomly generated and uniformely distributed false corresponding points (i.e. outliers) to
%   every scene in the cell array. The number of points that are added is specified as
%   as a fraction of the number of true correspondences (i.e. inliers).
%
%   Parameters:
%
%      listOfRandomScenes - A cell array consisting of sceneData structs.
%                           The sceneData struct contain ground truth corresponding
%							points that may or may not have been perturbed by Gaussian noise.
%
%     fractionOfOutliers  - the number of outliers that should be added to
%                           to the corresponding points specified as a
%                           fraction of the number of existing inlier
%                           corresponding points
%							
%
%     bounds             -  a 2D array containing the minium and maximum
%                           values for the x and y coordinates of the
%                           outliers (we assume that the image is square)
%
%   Returns: 
%
%      A cell array containing all planes and scenes with noisy corresponding points 
%	   between two views and with additional outliers (false correspondences).
%
%   See Also:
%
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 15/5/2012 
function listOfRandomScenesWithOutliers = add_outliers_to_scenes( listOfRandomScenes,fractionOfOutliers,bounds)
    numberOfScenes = length(listOfRandomScenes);
    listOfRandomScenesWithOutliers = cell(1,numberOfScenes);
    for i = 1:numberOfScenes
        sceneData = listOfRandomScenes{i};
        % we copy the sceneData structure
        sceneDataWithOutlier = sceneData;
        numOfH = length(sceneData);
            for j = 1:numOfH
                x = sceneData(j).ptsInViewOne;
                xp = sceneData(j).ptsInViewTwo;
                n = length(x);
                nOutliers = floor(fractionOfOutliers*n);
                xOutlier = randi([bounds(1),bounds(2)],[2,nOutliers]);
                xpOutlier = randi([bounds(1),bounds(2)],[2,nOutliers]);
                % and overwrite parts of the sceneData strucure with new data           
                sceneDataWithOutlier(j).ptsInViewOne = [x, xOutlier];
                sceneDataWithOutlier(j).ptsInViewTwo = [xp, xpOutlier];
            end
        listOfRandomScenesWithOutliers{i} = sceneDataWithOutlier;
    end
end

