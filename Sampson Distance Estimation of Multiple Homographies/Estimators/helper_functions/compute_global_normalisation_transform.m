%   Function: compute_global_normalisation_transform
%
% Given a collection of corresponding points between two views
% associated with one or more planes, this function takes all the data
% points and transforms them all jointly into a new coordinate system
% such that the points lie inside a unit box. If a set of initial 
% homographies are passed into the function, then the homographies
% are also transformed into the new coordinate system.
%
%
%
%   Parameters:
%
%      sceneData 	   			   - a struct containing noisy 
%                                    corresponding points between two views
%
%      listOfInitialH  			   - initial homography estimates
%
%										
%
%   Returns: sceneData with points in a new normalised coordinate system
%           
%            listOfInitialH transformed into a new coordinate system
%            if the user passed it as a parameter, otherwise an
%            empty cell array is returned
%           
% 
%            T the transformation that maps all data points in the first
%            image into a unit box.
%
%            Tp the transformation that maps all data points in the second
%            image into a unit box.
%
%   See Also:
%
% 
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 15/5/2012 
function [sceneData, listOfInitialH, T, Tp] = ...
          compute_global_normalisation_transform(sceneData, listOfInitialH)
   
    numOfH = length(sceneData);
    
    
    % collect all the pts in view one into a matrix
    % collect all the pts in view two into a matrix
    x = sceneData(1).ptsInViewOne;
    xp = sceneData(1).ptsInViewTwo;
    for i = 2:numOfH
      x = horzcat(x,sceneData(i).ptsInViewOne);
      xp = horzcat(xp,sceneData(i).ptsInViewTwo);
    end
    
    nn = length(x);
    x = [ x; ones( 1, nn ) ];
    xp = [ xp; ones( 1, nn ) ];
    T = normalise_data_transform( x );
    Tp = normalise_data_transform( xp );
    
    % convert the points and homographys in each plane into the new
    % joint coordinate system    
    for i = 1:numOfH
      
        % Check to see if the user passed in a list of homographies 
        if (exist('listOfInitialH', 'var'))
            % Transform the homographies into the new global coordinate
            % system
           listOfInitialH{i} = Tp * listOfInitialH{i} / T;
           
           % the formula below would be used to transform from the
           % normalised back to the original coordinate system
           %listOfInitialH{i} = Tp \ listOfInitialH{i} * T;
        end
        
        x = sceneData(i).ptsInViewOne;
        xp = sceneData(i).ptsInViewTwo;
        n = length(x);
        x = [ x; ones( 1, n ) ];
        xp = [ xp; ones( 1, n ) ];
        
        xx =  T * x;
        xpp = Tp * xp;
        
        xx = xx ./ repmat( xx(3,:), 3, 1 );
        xpp = xpp ./ repmat( xpp(3,:), 3, 1 );
        sceneData(i).ptsInViewOne = xx(1:2,:);
        sceneData(i).ptsInViewTwo = xpp(1:2,:);       
    end
    
    % if the user did not pass in a list of initial homographies
    % then we just pass back an empty cell array
    if (~exist('listOfInitialH', 'var'))
        listOfInitialH = cell(1,1);
    end

end

