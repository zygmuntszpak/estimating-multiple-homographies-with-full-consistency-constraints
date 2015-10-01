%   Function: compute_fns_estimates
%
%   Computes Approximate Maximum Likelihood homography estimates using 
%   the Fundamental Numerical Scheme (FNS) method.  Since it is an 
%   iterative method, it requires and initial seed value which is
%   usually the DLT estimate. 
%
%   Parameters:
%
%      sceneData 			- a struct containing noisy corresponding 
%                             points between two views
%
%      numberOfIterations 	- a maximum number of iterations for the 
%                            iterative estimation scheme
%
%      listOfInitialH    	- a cell array containing initial homography
%                             estimates
%									
%
%   Returns: A cell array of estimated homographies
% 
%
%   See Also:
%
%  compute_fns_covariances
%  compute_dlt_estimates
%
%
%  Last modified 15/5/2012 
function listOfEstimatedH  = ...
    compute_fns_estimates( sceneData, numberOfIterations, listOfInitialH )

numOfH = length(sceneData);
listOfEstimatedH = cell(1,numOfH);
listOfCovariancesOfH = cell(1,numOfH);
for i=1:numOfH
    x = sceneData(i).ptsInViewOne;
    xp = sceneData(i).ptsInViewTwo;
    n = length(x);
    
    % Need homogenous coordinates for normalise_transform function
    x = [ x; ones( 1, n ) ];
    xp = [ xp; ones( 1, n ) ];
    T = normalise_data_transform( x );
    Tp = normalise_data_transform( xp );
         
    Htilde = ...
        fns( T * x, Tp * xp, numberOfIterations,Tp*listOfInitialH{i} /(T));
        
    
    H =  Tp \ Htilde * T;
    H = H / H(3,3);
    listOfEstimatedH{i} = H;    
end

end

function  H = fns( x, xp , numberOfIterations,initialH)
     % Uses notation of the paper:
     % 'A multi-objective parameter estimation for image mosaicing,
     %  T. Scoleri, W. Chojnacki and M. J. Brooks '            
     
     %theta is the estimate of the homography represented as a vector
     theta = reshape(initialH',1,9)';
     
     oldTheta = zeros(1,9)';
     n = length(x);
     
     
     for k = 1:numberOfIterations
         Xthetai = zeros(9,9);
         % covariance of FNS estimate
         Mi = zeros(9,9);
         for i=1:n
             % the x-coordinate in the first image
             u = x(1,i);
             % the y-coordinate in the first image
             v = x(2,i);
             % the x-coordinate in the second image
             u2 = xp(1,i);
             % the y-coordinate in the second image
             v2 = xp(2,i);
             
             % Constructr Data Carrier matrix
             ux1 = [ 0    0    0  -u  -v -1  u*v2   v*v2   v2 ]';
             ux2 = [ u    v    1   0   0  0  -u*u2 -v*u2  -u2 ]';
             ux3 = [-u*v2 -v*v2 -v2 u*u2 v*u2 u2  0    0      0 ]';
             % 9 x 3 matrix
             Ui = [ux1 ux2 ux3];
             
             % 27 x 4 matrix
             dxUi = [0 0 0 -1 0 0 v2 0 0 1 0 0 0 0 0 -u2 0 ...
                                                0 -v2 0  0  u2 0  0  0 0 0;
                 0 0 0 0 -1 0 0 v2 0 0 1 0 0 0 0  0 -u2 ...
                                                0  0 -v2 0  0  u2 0  0 0 0;
                 0 0 0 0  0 0 0 0  0 0 0 0 0 0 0 -u -v -1 ...
                                                   0  0  0  u  v  1  0 0 0;
                 0 0 0 0  0 0 u v  1 0 0 0 0 0 0  0  0  0 ...
                           -u -v -1  0  0  0  0 0 0]' ;%note the transpose
             % 4 x 4 matrix (assume identity for covariance of data points)
             covariancei = eye(4,4);
             % 27 x 27 matrix
             Bi = dxUi * covariancei * dxUi';
             
             Sigmai = kron(eye(3,3),theta') * Bi * kron(eye(3,3),theta);
             
             % Find the inverse of Sigmai
             [U, S, V] = svd(Sigmai);
             S(1,1) = 1 / S(1,1);
             S(2,2) = 1 / S(2,2);
             S(3,3) = 0;
             inverseSigmai = U*S*V';
             
             etai = inverseSigmai * Ui' * theta    ;
             Xthetai = Xthetai + (Ui *inverseSigmai*Ui' - ...
                          kron(etai',eye(9,9)) * Bi * kron(etai,eye(9,9)));

         end
         
         % Find the smallest eigenvalue and choose its corresponding
         % eigenvector as the estimate of theta
         [eigenvectors eigenvalues] =  eig(Xthetai);
         [val ind] = min(abs(diag(eigenvalues)));
         
         theta = eigenvectors(:,ind);
         theta = theta / norm(theta);
         %difference = theta - oldTheta
         oldTheta = theta;
                           
     end
     
     % return a 9x9 matrix 
     H = reshape(theta,3,3)';     
    
          
end

function H = fns_alternative( x, xp , numberOfIterations,initialH)
 %theta is the estimate of the homography represented as a vector
     %theta = rand(9,1);    
     
     theta = reshape(initialH',1,9)';
     
     oldTheta = zeros(1,9)';
     n = length(x);
     
     for k = 1:numberOfIterations
         Xthetai = zeros(9,9);
         
         for i=1:n
             % the x-coordinate in the first image
             u = x(1,i);
             % the y-coordinate in the first image
             v = x(2,i);
             % the x-coordinate in the second image
             u2 = xp(1,i);
             % the y-coordinate in the second image
             v2 = xp(2,i);
             
             % Constructr Data Carrier matrix
             ux1 = [ 0    0    0  -u  -v -1  u*v2   v*v2   v2 ]';
             ux2 = [ u    v    1   0   0  0  -u*u2 -v*u2  -u2 ]';
             ux3 = [-u*v2 -v*v2 -v2 u*u2 v*u2 u2  0 0 0]';
             % 9 x 3 matrix
             Ui = [ux1 ux2 ux3];
             
             % 27 x 4 matrix
             dxUi = [ 0, 1, -v2, 0, 0, 0, 0, 0, 0, -1, 0, u2, 0, 0,...
                                  0, 0, 0, 0, v2, -u2, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 1, -v2, 0, 0, 0, 0, 0, 0, -1, 0, u2, 0, ...
                                  0, 0, 0, 0, 0, v2, -u2, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, u, 0, 0, v, 0, ...
                                  0, 1, 0, -u, 0, 0, -v, 0, 0, -1, 0;
                 0, 0, -u, 0, 0, -v, 0, 0, -1, 0, 0, 0, 0, 0,0, 0, 0,...
                                  0, u, 0, 0, v, 0, 0, 1, 0, 0]';
             % 4 x 4 matrix (assume identity for covariance of data points)
             covariancei = eye(4,4);
             % 27 x 27 matrix
             Bi = dxUi * covariancei * dxUi';
             
             Sigmai = kron(theta',eye(3,3)) * Bi * kron(theta,eye(3,3));
             
             % Find the inverse of Sigmai
             [U, S, V] = svd(Sigmai);
             S(1,1) = 1 / S(1,1);
             S(2,2) = 1 / S(2,2);
             S(3,3) = 0;
             inverseSigmai = U*S*V';
             
             etai = inverseSigmai * Ui' * theta    ;
             Xthetai = Xthetai + (Ui *inverseSigmai*Ui' - ...
                         kron(eye(9,9),etai') * Bi * kron(eye(9,9),etai));
         end
         
    
         % Find the smallest eigenvalue and choose its corresponding
         % eigenvector as the estimate of theta
         [eigenvectors eigenvalues] =  eig(Xthetai);
         [val, ind] = min(abs(diag(eigenvalues)));
         
         theta = eigenvectors(:,ind);
         theta = theta / norm(theta);
         fprintf('\n iteration %f \n',k)
             
         oldTheta = theta;
     end
      %fprintf('\n ')
     % return a 9x9 matrix 
     H = reshape(theta,3,3)'; 
end



