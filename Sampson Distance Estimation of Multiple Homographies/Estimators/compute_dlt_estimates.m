%   Function: compute_dlt_estimates
%
%   Computes homography estimates using the Direct Linear Transform (DLT) 
%   method. 
%
%   Parameters:
%
%      sceneData 	- a struct containing noisy corresponding points 
%                     between two views
%									
%
%   Returns: A cell array of estimated homographies
% 
%
%   See Also:
%
%  compute_dlt_covariances
%  compute_fns_estimates
%
%  Acknowledgement: I recycled code from Matt Dailey for most of these
% 					functions.
%
%  Last modified 15/5/2012 
function  listOfEstimatedH  = compute_dlt_estimates( sceneData )

 
numOfH = length(sceneData);
listOfEstimatedH = cell(1,numOfH);

for i = 1:numOfH
    
    x = sceneData(i).ptsInViewOne;
    xp = sceneData(i).ptsInViewTwo;
    n = length(x);
    
    % DLT requires homogenous coordinates
    x = [ x; ones( 1, n ) ];
    xp = [ xp; ones( 1, n ) ];
    
    T = normalise_data_transform( x );
    Tp = normalise_data_transform( xp );
    
    Htilde = dlt( T * x, Tp * xp );
    H = inv( Tp ) * Htilde * T;
    H = H / H(3,3);
    listOfEstimatedH{i} = H;
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



