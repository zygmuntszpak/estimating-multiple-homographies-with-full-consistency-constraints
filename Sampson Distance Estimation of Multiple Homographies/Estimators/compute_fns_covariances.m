%   Function: compute_fns_covariances
%
%    Computes a first order approximation of the covariance matrix 
%    associated with a homography that was estimated with the FNS method. 
%    It is assumed that the data points have been globally normalised, 
%    and hence that the homography is in normalised space. However, 
%    this is striclty speaking not necessary. If everything is in
%    normalised space, then the resulting covariance estimate is valid only 
%    for the homography in normalised space. If you decide to de-normalise
%    the homography using transformation matrices to go back to your 
%    original data space, then the resulting homography covariance matrices 
%    will also need to be transformed by transformation matrices.  
%    We operate in normalised space since the compute_aml_estimates
%    operates in normalised space and uses these covariances. 
%
%   Parameters:
%
%      sceneData  - a struct containing noisy corresponding points between
%                   two views
%
%
%      listOfH	  -  The FNS homography estimates for which we want to 
%                    compute the covarianc matrices.
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
function listOfCovariances = compute_fns_covariances(sceneData, listOfH)

% Implementation uses notation of the paper:
% 'A multi-objective parameter estimation for image mosaicing, ...
% T. Scoleri, W. Chojnacki and M. J. Brooks '

numOfH = length(listOfH);
listOfCovariances = cell(1,numOfH);

for k = 1:numOfH
    x = sceneData(k).ptsInViewOne;
    xp = sceneData(k).ptsInViewTwo;
    
    % Need homogenous coordinates for normalize_transform function
    n = length(x);
    x = [ x; ones( 1, n ) ];
    xp = [ xp; ones( 1, n ) ];
   
    H = listOfH{k};
    theta = reshape(H',1,9)';
    theta = theta / norm(theta,'fro');
    1;
    
    n = length(x);
    
    Mi = zeros(9,9);
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
        dxUi = [0 0 0 -1 0 0 v2 0 0 1 0 0 0 0 0 -u2 0  0 -v2 0  0  u2 0  0  0 0 0;
            0 0 0 0 -1 0 0 v2 0 0 1 0 0 0 0  0 -u2 0  0 -v2 0  0  u2 0  0 0 0;
            0 0 0 0  0 0 0 0  0 0 0 0 0 0 0 -u -v -1  0  0  0  u  v  1  0 0 0;
            0 0 0 0  0 0 u v  1 0 0 0 0 0 0  0  0  0 -u -v -1  0  0  0  0 0 0]' ;%note the transpose
        % 4 x 4 matrix (assume identity for covariance of data points)
        covariancei = eye(4,4);
        % 27 x 27 matrix
        Bi = dxUi * covariancei * dxUi';
       
        
        Sigmai = kron(eye(3,3),theta') * Bi * kron(eye(3,3),theta);
        
        % Find the inverse of Sigmai
        [U S V] = svd(Sigmai);
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
    hCovariance = commutation(3,3) * homographyCovariance * commutation(3,3);
    listOfCovariances{k} = hCovariance;
	

end

end

