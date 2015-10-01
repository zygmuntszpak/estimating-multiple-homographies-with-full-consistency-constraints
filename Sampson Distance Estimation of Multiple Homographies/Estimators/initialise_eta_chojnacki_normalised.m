
%   Function: initialise_eta_chojnacki_normalised
%
%   Takes a collection of homographies and computes a vector of latent
%   variables encoding new fully consistent homographies that are 'close' 
%   to the original homographies. This is a direct method utilising 
%   eigenvalue decomposition and singular value decomposition. 
%   With noiseless data it gives the a perfect decomposition. 
%   With noisy data, the resulting fully consistent homographies
%   are not very accurate, but are good enough to serve as initial seeds 
%   for joint bundle adjustment, the compute_aml_estimates method or
%   the compute_robustsampsonaml_estimated method which refine the
%   homographies.  

%   With noisy data, the accuracy of the decomposition depends on the
%   covariance associated with each homography. This is because the SVD 
%   operations use the 'Frobernius' norm, which assumes that every component
%   of a homography matrix is estimated with the same level of uncertainty,
%   and that each homography matrix is also estimated with the same level
%   of uncertainty. 
%   This is incorrect. Ideally, the SVD would take into account the
%   covariance of each homography matrix. By transfering all the
%   homographies into a globally normalised co-ordinate system (Hartley
%   type normalisation) the covariance associated with each homography
%   matrix becomes more 'identity' like (very roughly speaking)  so that
%   difference in the level of uncertaint beween homographies is somewhat
%   reduced. The global normalisation is a blunt instrument and the
%   covariancea certainly aren't identity, however,
%   estimates are definitly better when the SVD is performed in the
%   normalised space. Hence, this function will first convert the
%   homographies to the globally normalised space. The vector of latent
%   variables that the function returns  will then be in *normalised* space.
%   However, the  the list of valid initial homography matrices that
%   the function return is in the *original* (unnormalised) space. 
% 
%
%   Parameters:
%
%      sceneData            - A sceneData struct that contains the 
%                             corresponding points from view one to view 
%                             two
%
%      listOfInitialH    	- a cell array containing initial homography 
%                             estimates
%									
%
%   Returns: A vector of latent variables and corresponding fully
%            consistent homographies. The latent variables are
%            strucured as 
%
%            eta = [a,  b, v_1,w_1,... ,v_n, w_n]
%            
%            where a is length-9, b is length-3, v_i is length-3 and w_i
%            is a scalar.
%
%            Note that this is slightly different to the notation
%            used in the paper
%
%            Z. L. Szpak, W. Chojnacki, A. Eriksson, and A. van den Hengel. 
%            Sampson distance based joint estimation of multiple 
%            homographies with uncalibrated cameras. 
%            Comput. Vis. Image Underst., 125:200-213, 2014. 
%            http://dx.doi.org/10.1016/j.cviu.2014.04.008
%
%            which structures eta as:
%
%
%            eta = [a,  b, v_1, ... ,v_n, w_1,...,w_n ].
%
%            The latent variables are associated with a global
%            normalised space.
%
%            On the other hand, the fully consistent homographies
%            that are returned are in the original *unnormalised*
%            coordinate system. 
% 
%
%   See Also:
%
%  initialise_eta_chojnacki
%
%
%  Last modified 15/5/2012 
function [eta, listOfValidInitialH] = ...
              initialise_eta_chojnacki_normalised(sceneData,listOfInitialH)

numOfH = length(listOfInitialH);

[~, listOfInitialH, T, Tp] = ...
          compute_global_normalisation_transform(sceneData,listOfInitialH);

% stores 'corrected' homographies that satisfy all the constraints
% this is a "shotgun" approach at enforcing the constraints so if
% the initial homographies are poor, even though the upgraded homographies
% will be geometrically valid, they may not fit the data well
listOfValidInitialH = cell(1,numOfH);

X_1 = listOfInitialH{1};


% There have to be at least two planes in the scene
for k = 2:numOfH
    X_k = listOfInitialH{k};
    
    eigenvalues = eig(X_1,X_k);
    %eigenvaluesAlternate = eig(X_k\X_1);
    1;
    
   
    res1 = eigenvalues(1);
    res2 = eigenvalues(2);
    res3 = eigenvalues(3);
    
    
    d1 = sum((res1-res2).^2).^0.5;
    d2 = sum((res1-res3).^2).^0.5;
    d3 = sum((res2-res3).^2).^0.5;
    
    distances = [d1 d2 d3];
    min(distances);
    
    ind = find(distances == min(distances));
    
    if (ind == 3)
        combinedMatrices =  [(eigenvalues(2) * listOfInitialH{k} - X_1)...
                             (eigenvalues(3)* listOfInitialH{k} - X_1) ]  ;
        selectedEigenvalues{k} = ...
                          (real(eigenvalues(2))+ real(eigenvalues(3))) / 2;
    elseif (ind == 2)
        combinedMatrices = [(eigenvalues(1) * listOfInitialH{k} - X_1) ...
                              (eigenvalues(3)* listOfInitialH{k} - X_1) ] ;
        selectedEigenvalues{k} = ...
                          (real(eigenvalues(1))+ real(eigenvalues(3))) / 2;
    elseif (ind == 1)
        combinedMatrices =  [(eigenvalues(1) * listOfInitialH{k} - X_1) ...
                              (eigenvalues(2)* listOfInitialH{k} - X_1) ] ;
        selectedEigenvalues{k} = ...
                          (real(eigenvalues(1))+ real(eigenvalues(2))) / 2;
    end
    
    if k == 2
        JUXTAPOSE = combinedMatrices;
    else
        JUXTAPOSE = horzcat(JUXTAPOSE,combinedMatrices);
    end
end

% We will now construct a decomposition w*A + b*vn' which will 
% ensure that the initial homographies satisfy all constraints
% and are valid

% take for b_hat the left singular vector obtained by
% juxtapositioning the matrices:
[U,S,V] = svd(JUXTAPOSE);
            
% take the first column of u
b_hat = real(U(:,1));


% We have several possible choices for computing A_hat and vn_1
% one option is to choose:
%A_hat = X_1;            
%vn_1 = [0 ; 0 ; 0];

% another option is to choose:
Pb = eye(3) - (b_hat*b_hat')/(norm(b_hat,2)^2) ;
A_hat = Pb * X_1;
vn_1 = (X_1'*b_hat)/(norm(b_hat,2)^2);

% either way, we shall set v1 to 1
w1 = 1;

% vectorise A_hat
a = reshape(A_hat,9,1);

% vectorise A_hat
a = reshape(A_hat,9,1);

% build up eta, our parameter vector       
eta = a';
eta =  horzcat(eta,b_hat');
eta = horzcat(eta,vn_1');
eta = horzcat(eta,w1);

listOfValidInitialH{1} = (w1*A_hat + b_hat * vn_1');



% Now compute the remaining homographies and continue to build
% up eta
for k = 2:numOfH
    mu = selectedEigenvalues{k};
    X_k = listOfInitialH{k};
    vn_k = vn_1 +  (mu * X_k - X_1)' * b_hat / (norm(b_hat,2)^2);
    wk = 1 ;    
    listOfValidInitialH{k} = (wk*A_hat + b_hat * vn_k');
        
    % build up our eta parameter vector
    eta = horzcat(eta,vn_k');
    eta = horzcat(eta,wk);
end

% denormalise
for i = 1:length(listOfValidInitialH)
        listOfValidInitialH{i} =  Tp \ listOfValidInitialH{i} * T;
end

eta = eta';

end



function [sceneData, listOfInitialH, T ,Tp] = ...
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