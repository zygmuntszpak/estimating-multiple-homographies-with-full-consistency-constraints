%   Function: compute_rank_estimates
%
%   Given a series of corresponding points and 6 or more planes, this 
%   function computes DLT homographies and upgrades them  so that they 
%   satisfy certain rank-4 constraints that arise in two-view geometry
%   due to the rigidityof the scene. 
%
%
%   Parameters:
%
%      sceneData 	   			   - a struct containing noisy 
%                                    corresponding points between two views
%
%      listOfInitialH  			   - initial homography estimates, either 
%                                    DLT or FNS based
%
% 
%										
%
%   Returns: A cell array of rank-4 consistent homographies.
% 
%
%   See Also:
%
%  compute_aml_estimates
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 15/5/2012 

function  [listOfEstimatedH] = compute_rank_estimates( ...
    sceneData, listOfInitialH)


% transfer the pts in all views and the list of initial homographies to
% a common coordinate system for numerical stability
[~, listOfInitialH, T, Tp] = ...
                       normalise_joint_transform(sceneData,listOfInitialH);


numOfH = length(sceneData);
listOfEstimatedH = cell(1,numOfH);

H = zeros(9,length(listOfInitialH));
for i = 1:length(listOfInitialH)
    H(:,i) = vec(listOfInitialH{i})/norm(vec(listOfInitialH{i}));
    H(:,i) = H(:,i);
end

% enforce the rank-4 constraint
[U, D, V] = svd(H);
total  = min(length(listOfInitialH),9);
for i =  5:total
    D(i,i) = 0;
end
H = U*D*V';
for i = 1:length(listOfInitialH)
    H_norm(:,i) = H(:,i)/sign(H(end,i));
    H_norm(:,i) = H_norm(:,i) / norm(H_norm(:,i));
    Htilde = reshape(H_norm(:,i),3,3);    
    listOfEstimatedH{i}=  Tp \ Htilde  * T;    
end


end


function [sceneData, listOfInitialH, T, Tp] = ...
                       normalise_joint_transform(sceneData, listOfInitialH)
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
    T = normalize_transform( x );
    Tp = normalize_transform( xp );
    
    % convert the points and homographys in each plane into the new
    % joint coordinate system    
    for i = 1:numOfH
        %listOfInitialH{i} = Tp \ listOfInitialH{i} * T;
        listOfInitialH{i} = Tp * listOfInitialH{i} / T;
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

end

function T = normalize_transform( x )

  % Transform taking x's centroid to the origin

  Ttrans = [ 1 0 -mean( x(1,:) ) ; 0 1 -mean( x(2,:) ) ; 0 0 1 ];

  % Calculate appropriate scaling factor

  x = Ttrans * x;
  lengths = sqrt( sum( x(1:2,:).^2 ));
  s = sqrt(2) / mean(lengths);

  % Transform scaling x to an average length of sqrt(2)

  Tscale = [ s 0 0 ; 0 s 0 ; 0 0 1 ];

  % Compose the transforms

  T = Tscale * Ttrans;
end

