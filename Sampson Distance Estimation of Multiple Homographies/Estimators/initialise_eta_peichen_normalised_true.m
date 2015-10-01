%   Function: initialise_eta_peichen_normalised_true
%
%   Takes a collection of homographies and computes a vector of latent 
%   variables encoding new  fully consistent homographies that are 'close' 
%   to the original homographies. This is a direct method utilising
%   singular value decomposition.  With noiseless data it gives a perfect
%   decomposition. With noisy data, the resulting fully consistent
%   homographies are not very accurate, but are good enough to serve as 
%   initial seeds for joint bundle adjustment, the compute_aml_estimates
%   method or the compute_robustsampsonaml_estimated method which refine 
%   the homographies. 
%
%   This function requires at least three homographies. For the case of
%   two or more homographies, use the initialise_eta_chojnacki_normalised
%   function instead.
%
%   This function is based on work presented in 
%
%   Pei Chen, David Suter
%   Rank Constraints for Homographies over Two Views: Revisiting the
%   Rank Four Constraint
%   International Journal of Computer Vision
%   February 2009, Volume 81, Issue 2, pp 205-225
%   10.1007/s11263-008-0167-z
%
%   although a succint algorithm detailing all of the initialisation
%   steps is not presented in that paper.  
% 
%
%   Parameters:
%
%      listOfInitialH  - a cell array containing initial homography
%                        estimates
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
%            eta = [a,  b, v_1, ... ,v_n, w_1,...,w_n ]
% 
%            The latent variables are associated with a global
%            normalised space.
%
%            On the other hand, the fully consistent homographies
%            that are returned are in the original *unnormalised*
%            coordinate system. 
%
%   See Also:
%
%  initialise_eta_chojnacki_normalised
%
%
%  Last modified 15/5/2012 
function [eta, listOfValidInitialH] = ...
           initialise_eta_peichen_normalised_true(sceneData,listOfInitialH)

numOfH = length(listOfInitialH);

[~, listOfInitialH, T, Tp] = ...
                      normalise_joint_transform(sceneData,listOfInitialH);



e1 = [1,0,0]';
e2 = [0,1,0]';
e3 = [0,0,1]';

E1 = [e2*e1' + e3*e2', e3*e3' - e1*e1', -e1*e2' - e2*e3'];
E2 = [e1*e3' - e2*e2' + e3*e1'];

D = [1 0 0 0 0 0;
     0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 1 0 0;
     0 1 0 0 0 0;
     0 0 0 0 0 1;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     0 0 1 0 0 0];
 
 counter = 1 ;
 listOfM = cell(1,numOfH);
 for i=1:numOfH
  for j =i+1:numOfH
      X_i = listOfInitialH{i};
      X_j = listOfInitialH{j};
      C = -E1*(kron(eye(3),E2))*kron(X_i',X_j');
      listOfM{counter} = C*D;
      counter = counter + 1;
  end
 end
 M  = cell2mat(listOfM');
[U, Sing, V] = svd(M);
n1 = V(:,end-2);
n2 = V(:,end-1);
n3 = V(:,end);

 
 S = [create_matrix_N(n1) create_matrix_N(n2) create_matrix_N(n3)];
 
 [U, Sing, V] = svd(S);
 b = U(:,end);

 % Now compute the remaining homographies and continue to build
% up eta
X_All = [];
for k = 1:numOfH    
    X_All = horzcat(X_All,vec(listOfInitialH{k}));
end

Y = (eye(9) -  kron(eye(3),b*b')/norm(b)^2)*X_All;
[U, D, V] = svd(Y);
A = reshape(U(:,1),3,3);
a = reshape(A,9,1);
ws = D(1,1)*V(:,1);

vs = pinv(kron(eye(3),b))*(X_All-a*ws');

 % vectorise A_hat
a = reshape(A,9,1);

% build up eta, our parameter vector       
eta = a';
eta = horzcat(eta,b');
%eta = horzcat(eta,vs(:,1)');
%eta = horzcat(eta,ws(1));

% Now compute the remaining homographies and continue to build
% up eta
for k = 1:numOfH    
    X_k = listOfInitialH{k};   
    wk = ws(k);
    vn_k = vs(:,k);     
    listOfValidInitialH{k} = (wk*A + b * vn_k');
    % build up our eta parameter vector
    eta = horzcat(eta,vn_k');
    eta = horzcat(eta,wk);
end

% denormalise
for i = 1:length(listOfValidInitialH)
        listOfValidInitialH{i} =  Tp \ listOfValidInitialH{i} * T;
end

eta = eta';

 
 1;

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
