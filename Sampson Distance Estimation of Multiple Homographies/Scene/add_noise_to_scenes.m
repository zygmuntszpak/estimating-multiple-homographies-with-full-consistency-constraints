%   Function: add_noise_to_scenes
%
%   Given a cell array of random scenes and ground truth corresponding
%   points, the function applies zero-mean Gaussian noise with specified
%   standard deviation to the corresponding  points of every scene in the
%   cell array.
%
%   Parameters:
%
%      listOfRandomScenes - A cell array consisting of sceneData structs.
%                           The sceneData structs contain ground truth
%                           corresponding points.
%
%       sigma             - the standard deviation of Gaussian noise that
%                           is to be addedd to the ground truth
%                           corresponding points. 
%
%   Returns: 
%
%      A cell array containing all planes and scenes with noisy
%     corresponding points  between two views.
%
%   See Also:
%
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 15/5/2012 
function listOfRandomScenesWithNoise = ...
                             add_noise_to_scenes(listOfRandomScenes,sigma)
    numberOfScenes = length(listOfRandomScenes);
    listOfRandomScenesWithNoise = cell(1,numberOfScenes);
    for i = 1:numberOfScenes
        sceneData = listOfRandomScenes{i};
        sceneDataWithNoise = sceneData;
        numOfH = length(sceneData);
            for j = 1:numOfH
                x = sceneData(j).ptsInViewOne;
                xp = sceneData(j).ptsInViewTwo;
                n = length(x);
                % add zero mean noise with standard deviation sigma
                x = [  x + sigma*randn( 2, n )];
                xp = [  xp(1:2,:) + sigma*randn( 2, n )];
                sceneDataWithNoise(j).ptsInViewOne = x;
                sceneDataWithNoise(j).ptsInViewTwo = xp;
            end
        listOfRandomScenesWithNoise{i} = sceneDataWithNoise;
    end
end

