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

function visualise_groundtruth_scene_and_plot(sceneData)

%  sceneData is a struct that contains
%  list of homographies for the scene
%  list of points on planes in the first view
%  list of points on planes in the second view
numOfH = length(sceneData);

differentMarkers = {'o','o','o','o'};
figure

% need to add more colors if you wish to draw more than four 
% homographies
colors = [166, 206, 227;
           31, 120, 180;
          178, 223 ,138;
           51,160,44]./256;
       
[nColors, ~] = size(colors);

for i=1:numOfH

    x = sceneData(i).ptsInViewOne;
    xp = sceneData(i).ptsInViewTwo;    
    
    if (i <= nColors)
        col = colors(i,:);
    else     
        col = hsv2rgb([(randi(255)/256) 1 1]);
    end
    
    h = scatter(xp(1,:),xp(2,:),10,col,'filled',...
                                            'Marker',differentMarkers{i});
    hold on     
    for n = 1:length(x)
        plot( [x(1,n),xp(1,n)], [x(2,n),xp(2,n)],'-','Color',col,...
                                                      'LineWidth', 1.5 );
    end

 % we are assuming that the image dimensions are 500 x 500, ranging from
 % -250 to 250
 axis([-250 250 -250 250]) 
 hChildren = get(h, 'Children');
 set(hChildren, 'Markersize', 7)        
end

% Annotations
hTitle  = title ('Flow Vectors for Points');

hXLabel = xlabel('x'                     );
hYLabel = ylabel('y' );

% Adjust font and graphics options
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde');
% set([hLegend, gca]             , ...
%     'FontSize'   , 12           );
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 12          );
set( hTitle                    , ...
    'FontSize'   , 14          , ...
    'FontWeight' , 'bold'      );


set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XGrid'       , 'on'      , ...
  'XColor'      , [.0 .0 .0], ...
  'YColor'      , [.0 .0 .0], ...  
  'LineWidth'   , 1.2         );

set(gcf,'color','w');

% print2eps flowVectorsFourHomographiesOverlap.eps
% eps2pdf flowVectorsFourHomographiesOverlap.eps flowVectorsFourHomographiesOverlap.pdf

%print2eps flowVectorsFourHomographiesSpread.eps
%eps2pdf flowVectorsFourHomographiesSpread.eps flowVectorsFourHomographiesSpread.pdf

end

