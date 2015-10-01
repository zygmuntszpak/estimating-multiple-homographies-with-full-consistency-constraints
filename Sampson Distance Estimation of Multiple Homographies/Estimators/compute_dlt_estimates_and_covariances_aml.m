%   Function: compute_dlt_estimates_and_covariances_aml
%
%   Computes homography estimates using the Direct Linear Transform (DLT) 
%   method and additionally compute covariance matrices associated with the 
%   estimates. Note that this method assumes that the data points have been
%   globally normalised using a Hartley data transformation. Hence, the 
%   estimated homographies and covariance matrices are only valid in
%   normalised space. If you subsequently wish to 'de-normalise' the 
%   homographies, then you will have to 'de-normalise' the covariance
%   matrices too.  
%   The covariance matrices produced by this method are more accurate 
%   than the covariance matrices produce by
%   the compute_dlt_estimates_and_covariances method.
%
%   Parameters:
%
%      sceneData 	- a struct containing noisy corresponding points 
%                     between two views
%									
%
%   Returns: A cell array of estimated homographies and a cell array of 
%            associated covariance matrices both in a normalised coordinate
%            system.
% 
%
%   See Also:
%
%  compute_dlt_covariances
%  compute_dlt_estimates
%  compute_dlt_estimates_and_covariances
%
%  Acknowledgement: I recycled code from Matt Dailey for the DLT estimates.
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 15/5/2012 

% Assumes that the points have been globally normalised
function  [listOfEstimatedH, listOfCovariances ]  = ...
             compute_dlt_estimates_and_covariances_aml( sceneData)

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
    Ei = zeros(9,9);
    Di = zeros(9,9);
    n = length(x);
    
    
    for i = 1:n

        pone = [x(1,i) x(2,i) 1]';
        ptwo = [xp(1,i) xp(2,i) 1]';
        anti = ...
            [0 -ptwo(3) ptwo(2); ptwo(3) 0 -ptwo(1); -ptwo(2) ptwo(1) 0] ;
        Ui = kron(-pone,anti) ;
        u = x(1,i);
        v = x(2,i);
        u2 = xp(1,i);
        v2 = xp(2,i);
        
        dxUi =  [ 0, -1, v2, 0, 0, 0, 0, 0, 0, 1, 0, -u2, 0, 0, 0, 0, 0,...
                                           0, -v2, u2, 0, 0, 0, 0, 0, 0, 0;
                  0, 0, 0, 0, -1, v2, 0, 0, 0, 0, 0, 0, 1, 0, -u2, 0, 0,...
                                           0, 0, 0, 0, -v2, u2, 0, 0, 0, 0; 
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -u, 0, 0, -v, 0, 0,...
                                             -1, 0, u, 0, 0, v, 0, 0, 1, 0; 
                 0, 0, u, 0, 0, v, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
                                            -u, 0, 0, -v, 0, 0, -1, 0, 0]';     

          
        covariancei = eye(4,4);
        % 27 x 27 matrix
        Bi = dxUi * covariancei * dxUi';                   
        Sigmai = kron(eye(3,3),theta') * Bi * kron(eye(3,3),theta);         
        Ai = Ai + (Ui *Ui');   
        Di = Di + (Ui *Sigmai*Ui');        

    end
      
    
    % Add an extra normalisation step to avoid numerical problems    
    Nx = (1 / (theta'*theta)) * Ai;    
    Di = (1 / (theta'*theta)^2) * Di; 
    Nx = (1/n) * Nx;
    Di = (1/n)*Di;
    % Truncated Moore-Penrose Psuedo Inverse
    [U, S, V] = svd(Nx);
    for j = 1:8
        S(j,j) = 1 / S(j,j);
    end
    S(9,9) = 0;
    hCovariance = V*S*U';    
    Mx = (1/n)*hCovariance*Di*hCovariance;
        
    hCovariance = Mx;
    
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

  [U,D,V] = svd( A ,0);

  % In Octave, the SVD is sorted with decreasing singular values
  % so we want the last column of V

  H = reshape( V(:,9), 3, 3 )';
  H = H / H(3,3);

end


