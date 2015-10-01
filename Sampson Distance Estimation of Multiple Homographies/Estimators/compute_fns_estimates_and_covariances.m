
%   Function: compute_fns_estimates_and_covariances
%
%   Computes approximate maximum likelihood homography estimates using the 
%   Fundamental Numerical Scheme (FNS) method  and additionally computes
%   covariance matrices associated with the estimates. Note that this 
%   method assumes that the data points have been globally normalised
%   using a Hartley data transformation. Hence, the estimated homographies
%   and covariance matrices are only valid in normalised space.
%   If you subsequently wish to 'de-normalise' the homographies, then you
%   will have to 'de-normalise' the covariance matrices too.  
%   The covariance matrices produced by this method are more accurate 
%   than the covariance matrices produce by the 
%   compute_dlt_estimates_and_covariances method.
%   Since it is an iterative method, it requires and initial seed 
%   value which is usually the DLT estimate. 
%
%   Parameters:
%
%      sceneData 	- a struct containing noisy corresponding points
%                     between two views
%
%      numberOfIterations 	- a maximum number of iterations for the
%                             iterative estimation scheme
%
%      listOfInitialH    	- a cell array containing initial 
%                             homography estimates
%									
%
%   Returns: A cell array of estimated homographies and a cell array
%            of associated covariance matrices.
% 
%
%   See Also:
%
%  compute_fns_covariances
%  compute_fns_estimates
%
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 15/5/2012 

function [ listOfEstimatedH, listOfCovariances ]  = ...
                    compute_fns_estimates_and_covariances( sceneData,...
                                         numberOfIterations,listOfInitialH)

%   Assumes that the data has already been globablly normalised

numOfH = length(sceneData);
listOfEstimatedH = cell(1,numOfH);
listOfCovariances = cell(1,numOfH);

for i=1:numOfH
    x = sceneData(i).ptsInViewOne;
    xp = sceneData(i).ptsInViewTwo;
    n = length(x);
    
    % Need homogenous coordinates for normalize_transform function
    x = [ x; ones( 1, n ) ];
    xp = [ xp; ones( 1, n ) ];
    
    H = fns(  x, xp, numberOfIterations,listOfInitialH{i});

    H = H / H(3,3);
    listOfEstimatedH{i} = H;
    
end

for k = 1:numOfH
    x = sceneData(k).ptsInViewOne;
    xp = sceneData(k).ptsInViewTwo;
    
    % Need homogenous coordinates for normalize_transform function
    n = length(x);
    x = [ x; ones( 1, n ) ];
    xp = [ xp; ones( 1, n ) ];
   
    H = listOfEstimatedH{k};
    theta = reshape(H',1,9)';
    theta = theta / norm(theta,'fro');
    1;
    
    n = length(x);
    

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
        1;
        
        % 27 x 4 matrix
        dxUi = [0 0 0 -1 0 0 v2 0 0 1 0 0 0 0 0 -u2 0  0 -v2 0 ...
                                                        0  u2 0  0  0 0 0;
            0 0 0 0 -1 0 0 v2 0 0 1 0 0 0 0  0 -u2 0  0 -v2 0 ...
                                                           0  u2 0  0 0 0;
            0 0 0 0  0 0 0 0  0 0 0 0 0 0 0 -u -v -1  0  0 ...
                                                        0  u  v  1  0 0 0;
            0 0 0 0  0 0 u v  1 0 0 0 0 0 0  0  0  0 -u -v -1 ...
                                      0  0  0  0 0 0]' ;%note the transpose
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
        
        % Used for the computation of the covariance of the FNS
        % estimate
        Mi = Mi +  Ui* inverseSigmai * Ui';   

    end
        
    % compute the covariance of the homography
    thetaNormi = norm(theta,'fro')^2;
    Mi = thetaNormi * Mi;
    % Truncated Moore-Penrose Psuedo Inverse
    [U, S, V] = svd(Mi);
    for j = 1:8
        S(j,j) = 1 / S(j,j);
    end
    S(9,9) = 0;
    homographyCovariance = V*S*U';
    % because theta = reshape(H',1,9)' and not reshape(H,1,9)' we need to
	% use a commutation matrix in the last step. It is possible to derive
    % everything for the case of theta = reshape(H,1,9)' as well, but
    % this function follows the notation of the paper
    % Implementation uses notation of the paper:
    % 'A multi-objective parameter estimation for image mosaicing, ...
    % T. Scoleri, W. Chojnacki and M. J. Brooks '
    hCovariance = commutation(3,3) * homographyCovariance * ...
                                                        commutation(3,3);
    listOfCovariances{k} = hCovariance;
end

end

function  H = fns( x, xp , numberOfIterations,initialH)
% Uses notation of the paper:
% 'A multi-objective parameter estimation for image mosaicing, T. Scoleri,
% W. Chojnacki and M. J. Brooks '

%theta is the estimate of the homography represented as a vector
%theta = rand(9,1);
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
        dxUi = [0 0 0 -1 0 0 v2 0 0 1 0 0 0 0 0 -u2 0  0 -v2 ...
                                                     0  0  u2 0  0  0 0 0;
            0 0 0 0 -1 0 0 v2 0 0 1 0 0 0 0  0 -u2 0  0 -v2 ...
                                                        0  0  u2 0  0 0 0;
            0 0 0 0  0 0 0 0  0 0 0 0 0 0 0 -u -v -1  0  0  0 ...
                                                           u  v  1  0 0 0;
            0 0 0 0  0 0 u v  1 0 0 0 0 0 0  0  0  0 -u -v -1 ...
                                      0  0  0  0 0 0]' ;%note the transpose
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
    [eigenvectors, eigenvalues] =  eig(Xthetai);
    [val, ind] = min(abs(diag(eigenvalues)));
    
    theta = eigenvectors(:,ind);
    theta = theta / norm(theta);
    %difference = theta - oldTheta
    oldTheta = theta;
    
end

% return a 9x9 matrix
H = reshape(theta,3,3)';

end

