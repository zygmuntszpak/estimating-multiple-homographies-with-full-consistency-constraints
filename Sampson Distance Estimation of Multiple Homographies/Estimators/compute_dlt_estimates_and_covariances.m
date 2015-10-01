%   Function: compute_dlt_estimates_and_covariances
%
%   Computes homography estimates using the Direct Linear Transform (DLT)
%   method and additionally compute covariance matrices associated with
%   the estimates. Note that this method assumes that the data points have
%   been globally normalised using a Hartley data transformation. Hence,
%   the estimated homographies and covariance matrices are only valid in
%   normalised space.
%   If you subsequently wish to 'de-normalise' the homographies,
%   then you will have to 'de-normalise' the covariance matrices too.  
%
%   Parameters:
%
%      sceneData 	- a struct containing noisy corresponding points 
%                     between two views
%									
%
%   Returns: A cell array of estimated homographies and a cell array 
%           of associated covariance matrices both in a globaly
%           *normalised* coordinate system.
% 
%
%   See Also:
%
%  compute_dlt_covariances
%  compute_dlt_estimates
%
%  Acknowledgement: I recycled code from Matt Dailey for the DLT estimates.
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 15/5/2012 

% Assumes that the points have been globally normalised
function  [listOfEstimatedH, listOfCovariances ]  =  ...
                compute_dlt_estimates_and_covariances(sceneData)

numOfH = length(sceneData);
listOfEstimatedH = cell(1,numOfH);

% Compute DLT Estimates from globablly normalised data
for i =1:numOfH
    
    x = sceneData(i).ptsInViewOne;
    xp = sceneData(i).ptsInViewTwo;
    n = length(x);
    
    % DLT requires homogenous coordinates
    x = [ x; ones( 1, n ) ];
    xp = [ xp; ones( 1, n ) ];    
    
    H = dlt(  x,  xp );
    H = H / H(3,3);
    listOfEstimatedH{i} = H;
end

% Compute DLT Covariances from globally normalised data
listOfCovariances = cell(1,numOfH);
for k = 1:numOfH
    
    x = sceneData(k).ptsInViewOne;
    xp = sceneData(k).ptsInViewTwo;
    
    n = length(x);
    x = [ x; ones( 1, n ) ];
    xp = [ xp; ones( 1, n ) ];
    
    H = listOfEstimatedH{k};
    theta = reshape(H,1,9)';
    theta = theta / norm(theta,'fro');
    
    Ai = zeros(9,9);
    n = length(x);
        
   
    for i = 1:n
        
        pone = [x(1,i) x(2,i) 1]';
        ptwo = [xp(1,i) xp(2,i) 1]';
        anti = [0 -ptwo(3) ptwo(2); ptwo(3) 0 -ptwo(1); -ptwo(2) ptwo(1) 0]  ;
        Ui = kron(-pone,anti) ;
        Ai = Ai + (Ui *Ui');
        
    end    
    
    Nx = (1 / (theta'*theta)) * Ai;
    
    % Truncated Moore-Penrose Psuedo Inverse
    [U, S, V] = svd(Nx);
    for j = 1:8
        S(j,j) = 1 / S(j,j);
    end
    S(9,9) = 0;
    hCovariance = V*S*U';
    listOfCovariances{k} = hCovariance;
end

end


function H = dlt( x, xp )

  n = size( x, 2 );
  if n < 4
    error( 'DLT requires at least 4 points' );
  end;
  if ( size( x, 1 ) ~= 3 | size( xp, 1 ) ~= 3 )
    error( 'DLT requres homogeneous coordinates' );
  end;

  A = [];

  for i = 1:n

    xip = xp( 1, i );
    yip = xp( 2, i );
    wip = xp( 3, i );

    xi = x( :, i );

    Ai = [ 0, 0, 0,    -wip * xi',   yip * xi' ;
           wip * xi',     0, 0, 0,  -xip * xi' ];

    A = [ A ; Ai ];
  end;

  [~,~,V] = svd( A ,0);

  % In Octave, the SVD is sorted with decreasing singular values
  % so we want the last column of V

  H = reshape( V(:,9), 3, 3 )';
  H = H / H(3,3);

end


