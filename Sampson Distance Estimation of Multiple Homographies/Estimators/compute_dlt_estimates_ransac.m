%   Function: compute_dlt_estimates_ransac
%
%   Computes homography estimates using the Direct Linear Transform (DLT)
%   method and random sampling and consensus. 
%
%   Note that this function is currently *one big hack*--I have hard
%   coded the fact that there are TWO planar structures that we seek
%   and I use sequential ransac to extract the two structures.
%
%   This function should be rewritten so that you can specify the
%   number of strutures that you seek, and it should take as parameters
%   the requisite ransac thresholds etc. 
%
%   Parameters:
%
%      sceneData 	- a struct containing noisy corresponding 
%                     points between two views with potential 
%                     outliers. 
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
%                   I also used Peter Kovesi's *fantastic* functions
%                   for ransac: 
%
%  Last modified 15/5/2012 
function  [listOfEstimatedH, sceneData]  = ...
                            compute_dlt_estimates_ransac( sceneData )

 
numOfH = length(sceneData);
listOfEstimatedH = cell(1,numOfH);

sceneDataClone = sceneData;

nPointsStructureOne = length(sceneData(1).ptsInViewOne);
nPointsStructureTwo = length(sceneData(2).ptsInViewOne);


    
z = sceneData(1).ptsInViewOne;
zp = sceneData(1).ptsInViewTwo;
    
z = [z, sceneData(2).ptsInViewOne];
zp = [zp, sceneData(2).ptsInViewTwo];   

% Convert to homogeneous coordinates    
z = [ z; ones( 1, length(z) ) ];
zp = [ zp; ones( 1, length(z) ) ];
    
   
T = normalise_data_transform( z );
Tp = normalise_data_transform( zp );
    
t = .005;  % Distance threshold for deciding outliers
[H1, inliers] = ransacfithomography(T*z, Tp*zp, t);
H1 = H1 / H1(3,3);
H1 = inv( Tp ) * H1 * T;
H1 = H1 / H1(3,3);
inliers1view1 = z(:,inliers);
inliers1view2 = zp(:,inliers);
medianInlierIndex1 = median(inliers);
    

% remove the inliers corresponding to the first structure
z(:,inliers) = [];
zp(:,inliers) = []; 
    
sceneDataClone(1).ptsInViewOne = inliers1view1; 
sceneDataClone(1).ptsInViewTwo =  inliers1view2;

% ransac on the remaining points
t = .005;  % Distance threshold for deciding outliers
[H2, inliers] = ransacfithomography(T*z, Tp*zp, t);
H2 = inv( Tp ) * H2 * T;
H2 = H2 / H2(3,3);
        
inliers2view1 = z(:,inliers);
inliers2view2 = zp(:,inliers);
    
sceneDataClone(2).ptsInViewOne = inliers2view1; 
sceneDataClone(2).ptsInViewTwo =  inliers2view2;
    
medianInlierIndex2 = median(inliers);
    
% Since later on we will compare the performance of the homographies
% on ground truth data points, we need to be sure that the first
% structure and inliers that ransac found actually correspond to 
% to data points in our sceneData structure that originated from the
% first plane. 
    if (medianInlierIndex1 < nPointsStructureOne)
        listOfEstimatedH{1} = H1;
        listOfEstimatedH{2} = H2;
        sceneData(1).ptsInViewOne = sceneDataClone(1).ptsInViewOne(1:2,:);
        sceneData(1).ptsInViewTwo = sceneDataClone(1).ptsInViewTwo(1:2,:);
        sceneData(2).ptsInViewOne = sceneDataClone(2).ptsInViewOne(1:2,:);
        sceneData(2).ptsInViewTwo = sceneDataClone(2).ptsInViewTwo(1:2,:);
    else
        listOfEstimatedH{1} = H2;
        listOfEstimatedH{2} = H1;
        sceneData(2).ptsInViewOne = sceneDataClone(1).ptsInViewOne(1:2,:);
        sceneData(2).ptsInViewTwo = sceneDataClone(1).ptsInViewTwo(1:2,:);
        sceneData(1).ptsInViewOne = sceneDataClone(2).ptsInViewOne(1:2,:);
        sceneData(1).ptsInViewTwo = sceneDataClone(2).ptsInViewTwo(1:2,:);
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



