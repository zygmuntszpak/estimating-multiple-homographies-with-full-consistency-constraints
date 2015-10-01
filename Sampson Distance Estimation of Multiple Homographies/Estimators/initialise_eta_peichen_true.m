%   Function: initialise_eta_peichen_true
%
%   Takes a collection of homographies and computes a vector of latent 
%   variables encoding new  fully consistent homographies that are 'close' 
%   to the original homographies. This is a direct method utilising
%   eigenvalue decomposition and singular value decomposition. 
%   With noiseless data it gives the a perfect decomposition.
%   With noisy data, the resulting fully consistent homographies
%   are not very accurate, but are good enough to serve as initial seeds 
%   for joint bundle adjustment, the compute_aml_estimates method or
%   the compute_robustsampsonaml_estimated method which refine the
%   homographies. 
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
%
%           This function also returns the latent variables 
%           in separate groups, in particular as a matrix A,
%           vector b, matrix vs (a concatednation of the individual v_i)
%           and a vector ws. This structuring of latent variables is 
%           required for the WALS algorithm. 
%
%   See Also:
%
%  initialise_eta_chojnacki_normalised
%
%
%  Last modified 15/5/2012 

function [eta, listOfValidInitialH, A, b, vs, ws] = ...
                                initialise_eta_peichen_true(listOfInitialH)
% the method requires at least 3 homographies to produce an initialisation
numOfH = length(listOfInitialH);



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

1;
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

eta = eta'; 
1;
end



