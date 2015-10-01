%   Function: generate_multiple_random_planarscenes
%
%   Generates multiple random planar scenes each consisting of a 
%   specified number or planes. A random scene consists of homographies 
%    together with corresponding points.
%
%   Parameters:
%
%      numberOfScenes     - desired number of random scenes
%
%      numH               - number of 3D planes, and hence the number of
%							homographies to generate
%				
%
%
%   Returns: 
%
%     A cell array of sceneData structs that contains the homographies
%     and corresponding  points from view one to view two.
%
%   See Also: 
%
%    generate_groundtruth_scene
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 16/9/2014 
function listOfRandomScenes = generate_multiple_random_planarscenes_new(...
                numberOfScenes,numOfH,nPoints)
listOfRandomScenes = cell(1,numberOfScenes);
for i = 1:numberOfScenes
	sceneData = generate_groundtruth_scene_new(...
                                numOfH,nPoints);
	listOfRandomScenes{i} = sceneData;
end

end