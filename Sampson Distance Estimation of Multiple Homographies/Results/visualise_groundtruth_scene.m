%   Function: visualise_groundtruth_scene
%
%   Plots lines connecting the corresponding points between two views
%   for all homographies in the scene. The points in the first view and
%   the points in the second view are overlayed on the same image, and
%   a line joins the points.  
%
%   Parameters:
%
%      sceneData 	- a struct containing the ground truth homography and
%					  corresponding points between two views
% 									 

%	
%
%   Returns: 
% 
%
%   See Also:
%
%  generate_groundtruth_scene
%  
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 15/5/2012 

function [ output_args ] = visualise_groundtruth_scene(sceneData)

%  sceneData is a struct that contains
%  list of homographies for the scene
%  list of points on planes in the first view
%  list of points on planes in the second view
numOfH = length(sceneData);

differentColors = {{'g'},{'r'},{'m'},{'b'}}
 figure

for i=1:numOfH
    H = sceneData(i).homographies;
    x = sceneData(i).ptsInViewOne;
    xp = sceneData(i).ptsInViewTwo;
    
    
    %figure
    %scatter(x(1,:),x(2,:),'r')
   
    scatter(xp(1,:),xp(2,:),'black')
    hold on
 
    color = rand(3,1);
    for i = 1:length(x)
        plot( [x(1,i),xp(1,i)], [x(2,i),xp(2,i)],'-','Color',color);
    end
    
%     axis([0 imageWidth 0 imageHeight])
 axis([-250 250 -250 250])
    
    
end

end

