
%   Function: compute_dlt_covariances
%
%   Computes a first order approximation of the covariance matrix 
%   associated with a homography that was estimated with the DLT method. 
%   The resulting covariance estimate is valid only for the homography
%   in normalised space.  If you decide to de-normalise the homography
%   using transformation matrices to go back  to your original data space,
%   then the resulting homography covariance matrices will also need to be 
%   transformed by transformation matrices.  We operate in normalised space
%   since the compute_aml_estimates operates in normalised space and uses
%   these covariances. 
%
%   Parameters:
%
%      sceneData 	   - a struct containing noisy corresponding points 
%                        between two views
%
%      listOfT  	   - Cell Array of transformation matrices that 
%                        normalised data points in the first view. 
%                        There is one transformation matrix per 
%						 plane/homography.
%
%      listOfTp  	   - Cell Array of transformation matrices that 
%                        normalised data points in the second view. 
%                        There is one transformation matrix per 
%						 plane/homography.
%
%      jointT  		   - A transformation matrix that normalised data 
%                        points in the first view. This is a global
%                        transformation matrix in that all points,
%                        regardless of what plane/homography they
%						 belong to, are normalised. 
%
%      jointTp 		   - A transformation matrix that normalised data
%                        points in the first view. This is a global 
%                        transformation matrix in that all points, 
%                        regardless of what plane/homography they
%						 belong to, are normalised. 
%
%      listOfH		   - The DLT homography estimates for which we want
%                        to compute  the covarianc matrices.
% 
%										
%
%   Returns: A cell array of covariance matrices
% 
%
%   See Also:
%
%  compute_dlt_estimates
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 15/5/2012 
function listOfCovariances = ...
                    compute_dlt_covariances(sceneData,listOfT,...
                                         listOfTp,jointT,jointTp, listOfH)

numOfH = length(listOfH);
listOfCovariances = cell(1,numOfH);

for k = 1:numOfH
    
    x = sceneData(k).ptsInViewOne;
    xp = sceneData(k).ptsInViewTwo;
    
    n = length(x);
    x = [ x; ones( 1, n ) ];
    xp = [ xp; ones( 1, n ) ];
    
    H = listOfH{k};
    theta = reshape(H,1,9)';
    theta = theta / norm(theta,'fro');
    
    Ai = zeros(9,9);
    n = length(x);
    
    SiOne =  listOfT{k} / jointT;
    SiTwo =  listOfTp{k} / jointTp;
    
    Zi = det(SiTwo)*inv(SiTwo)';
    
    for i = 1:n
        
        pone = [x(1,i) x(2,i) 1]';
        ptwo = [xp(1,i) xp(2,i) 1]';
        anti = [0 -ptwo(3) ptwo(2); ptwo(3) 0 -ptwo(1); ...
                                                    -ptwo(2) ptwo(1) 0]  ;
        Ui = kron(-pone,anti) ;
        Ai = Ai + (Ui *(Zi'*Zi)*Ui');
        
    end
    
    Ci = kron((inv(SiOne)*inv(SiOne)'), (SiTwo'*SiTwo) );
    
    Nx = ((norm(theta,2)^2) / (theta'*Ci*theta)) * Ai;
    
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

